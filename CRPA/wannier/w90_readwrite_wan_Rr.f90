!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE w90_readwrite_wan_Rr(iRs, wan_Rr, io, tmp_dir_wfc)
  !----------------------------------------------------------------------------
  !! Read and write Wannier functions in real space.
  !! If io = +1, write wan_Rr to   file tmp_dir_wfc/prefix.wan_Rr.iR1.iR2.iR3.
  !! If io = -1, read  wan_Rr from file tmp_dir_wfc/prefix.wan_Rr.iR1.iR2.iR3.
  !----------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE mp_pools,         ONLY : me_pool, root_pool
  USE io_files,         ONLY : prefix, diropn
  USE scatter_mod,      ONLY : gather_grid, scatter_grid
  USE fft_base,         ONLY : dffts
  USE noncollin_module, ONLY : npol
  USE w90_interface,    ONLY : num_wann
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: iRs(3)
  !! Real space lattic vector indices
  COMPLEX(DP), INTENT(inout) :: wan_Rr(dffts%nnr, npol, num_wann)
  !! Wavefunction in the Wannier gauge in real space
  INTEGER, INTENT(in) :: io
  !! +1: Write. -1: Read.
  CHARACTER(LEN=256), INTENT(in) :: tmp_dir_wfc
  !! Directory for wavefunctions
  !! If io == +1 (writing), tmp_dir_wfc can be any directory (e.g. outdir/_ph0/).
  !! If io == -1 (reading), tmp_dir_wfc must be the directory where the formatted
  !! wavefunctions are stored.
  !
  INTEGER :: iw
  !! Wannier function index
  INTEGER :: ipol
  !! Spin index
  INTEGER :: lrwfc
  !! Record length
  INTEGER :: iun
  !! IO unit
  INTEGER :: ios
  !! IO status
  INTEGER*8 :: unf_recl
  !! double precision to prevent integer overflow
  INTEGER :: direct_io_factor
  !! Direct IO factor
  REAL(DP) :: dummy
  !! Dummy variable for direct_io_factor
  CHARACTER(LEN=256) :: filename
  !! File name
  COMPLEX(DP), ALLOCATABLE :: wan_Rr_gathered(:, :)
  !! Gathered wavefunctions
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  CALL start_clock("w90_readwrite")
  !
  iun = find_free_unit()
  !
  ALLOCATE(wan_Rr_gathered(dffts%nr1x*dffts%nr2x*dffts%nr3x, npol))
  !
  lrwfc = 2 * dffts%nr1x*dffts%nr2x*dffts%nr3x * npol
  INQUIRE (IOLENGTH=direct_io_factor) dummy
  unf_recl = direct_io_factor * INT(lrwfc, kind=kind(unf_recl))
  !
  WRITE(filename, '(A,I0,".",I0,".",I0)') TRIM(tmp_dir_wfc) // TRIM(prefix) &
      & // '.wan_Rr.', iRs(1), iRs(2), iRs(3)
  !
  IF (me_pool == root_pool) THEN
    ! FIXME: Using NEWUNIT gives a negative unit...
    OPEN(UNIT = iun, FILE = TRIM(ADJUSTL(filename)), IOSTAT = ios, FORM = 'unformatted', &
          STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
    IF (ios /= 0) CALL errore('w90_readwrite_wan_Rr', 'error opening '//TRIM(filename), 1)
  ENDIF
  !
  IF (io == +1) THEN
    !
    ! Write wan_Rr wavefunction to file
    !
    DO iw = 1, num_wann
      DO ipol = 1, npol
        CALL gather_grid(dffts, wan_Rr(:, ipol, iw), wan_Rr_gathered(:, ipol))
      ENDDO
      !
      IF (me_pool == root_pool) THEN
        CALL davcio(wan_Rr_gathered, lrwfc, iun, iw, +1)
      ENDIF
    ENDDO
    !
  ELSEIF (io == -1) THEN
    !
    ! Read wan_Rr wavefunction from file
    !
    wan_Rr = (0.0_DP, 0.0_DP)
    !
    DO iw = 1, num_wann
      IF (me_pool == root_pool) THEN
        CALL davcio(wan_Rr_gathered, lrwfc, iun, iw, -1)
      ENDIF
      !
      DO ipol = 1, npol
        CALL scatter_grid(dffts, wan_Rr_gathered(:, ipol), wan_Rr(:, ipol, iw))
      ENDDO
    ENDDO
    !
  ELSE
    CALL errore('w90_readwrite_wan_Rr', 'invalid io', 1)
  ENDIF
  !
  IF (me_pool == root_pool) THEN
    CLOSE(UNIT=iun, STATUS='keep')
  ENDIF
  !
  CALL stop_clock("w90_readwrite")
  !
END SUBROUTINE w90_readwrite_wan_Rr
