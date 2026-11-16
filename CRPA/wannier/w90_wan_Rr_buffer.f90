!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE w90_wan_Rr_buffer
  ! ----------------------------------------------------------------------------
  !! Module for reading and writing real-space Wannier functions in a distributed buffer.
  !-----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER :: iuwan_r
  !! Unit number for distributed Wannier functions buffer
  INTEGER :: lrwan_r
  !! Record length for distributed Wannier functions buffer (dffts%nnr * npol * num_wann)
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_wan_Rr_save_buffer(tmp_dir_collected, tmp_dir_buffer)
    !----------------------------------------------------------------------------
    !! Read real-space Wannier functions from the collected file and write them to a buffer.
    !! This is done so that during the calculation, each processor can access only the
    !! part of the Wannier functions it needs, without having to read the entire
    !! collected file every time.
    !!
    !! Read two types of files:
    !!   - List of real-space lattice vectors : prefix.wan_R.dat
    !!   - Wannier functions in real space    : prefix.wan_Rr.iR1.iR2.iR3
    !!
    !! The files are written in w90_write_real_space.
    !----------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE mp,               ONLY : mp_bcast
    USE mp_images,        ONLY : intra_image_comm
    USE io_global,        ONLY : stdout, ionode, ionode_id
    USE io_files,         ONLY : prefix
    USE buffers,          ONLY : open_buffer, save_buffer
    USE fft_base,         ONLY : dffts
    USE control_flags,    ONLY : io_level
    USE noncollin_module, ONLY : npol
    USE w90_interface,    ONLY : num_wann
    USE w90_real_space,   ONLY : nR_ws, iRlist_ws
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(in) :: tmp_dir_collected
    !! Directory where the real-space Wannier functions are stored in collected format.
    CHARACTER(LEN=256), INTENT(in) :: tmp_dir_buffer
    !! Directory where the real-space Wannier functions are stored in distributed format.
    !
    CHARACTER(LEN=256) :: filename
    !! File name for the Wigner-Seitz grid points
    LOGICAL :: exst_mem
    !! logical variable to check file exists
    LOGICAL :: exst
    !! logical variable to check file exists in memory
    INTEGER :: iR
    !! Real space lattice vector index
    INTEGER :: iun
    !! Unit number for writing the Wigner-Seitz grid points
    INTEGER :: ios
    !! IO status
    COMPLEX(DP), ALLOCATABLE :: wan_Rr(:, :, :)
    !! Wavefunction in the Wannier gauge in real space
    !
    CALL start_clock("w90_buf_save")
    !
    ! Read the list of real-space lattice vectors
    !
    filename = TRIM(tmp_dir_collected) // TRIM(prefix) // '.wan_R.dat'
    !
    IF (ionode) THEN
      OPEN(NEWUNIT=iun, FILE=filename, STATUS='old', FORM='formatted', IOSTAT=ios)
      IF (ios /= 0) CALL errore('w90_write_wan_Rr_from_collected_to_buffer', &
      'opening file ' // TRIM(filename), ABS(ios))
      !
      READ(iun, '(I8)', IOSTAT=ios) nR_ws
      ALLOCATE(iRlist_ws(3, nR_ws))
      DO iR = 1, nR_ws
        READ(iun, '(3I8)', IOSTAT=ios) iRlist_ws(:, iR)
      ENDDO
      !
      CLOSE(iun, STATUS='keep')
    ENDIF
    !
    CALL mp_bcast(nR_ws, ionode_id, intra_image_comm)
    IF (.NOT. ionode) ALLOCATE(iRlist_ws(3, nR_ws))
    CALL mp_bcast(iRlist_ws, ionode_id, intra_image_comm)
    !
    ! Print the Wigner-Seitz grid points
    !
    WRITE(stdout, '(/5x, A, I8)') "Number of Wigner-Seitz grid points = ", nR_ws
    DO iR = 1, nR_ws
      WRITE(stdout, '(5x, A, I6, A, 3I8)') "R-vector ", iR, " = ", iRlist_ws(:, iR)
    ENDDO
    !
    WRITE(stdout, '(/5x, A)') "Reading collected Wannier functions, writing them to buffer"
    !
    ! Open a file to write the Wannier functions in real space
    !
    iuwan_r = 100
    lrwan_r = dffts%nnr * npol * num_wann
    CALL open_buffer(iuwan_r, 'wan_r', lrwan_r, io_level, exst_mem, exst, tmp_dir_buffer)
    !
    ! Read wan_Rr wavefunction from file and write to the buffer
    !
    ALLOCATE(wan_Rr(dffts%nnr, npol, num_wann))
    !
    DO iR = 1, nR_ws
      CALL w90_readwrite_wan_Rr(iRlist_ws(:, iR), wan_Rr, -1, tmp_dir_collected)
      CALL save_buffer(wan_Rr, lrwan_r, iuwan_r, iR)
    ENDDO ! iR
    !
    DEALLOCATE(wan_Rr)
    !
    CALL stop_clock("w90_buf_save")
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE w90_wan_Rr_save_buffer
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_wan_Rr_get_buffer(iR, wan_Rr)
    !----------------------------------------------------------------------------
    !! Read real-space Wannier functions from the distrubuted buffer.
    !! The buffer is created by w90_write_wan_Rr_from_collected_to_buffer.
    !----------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE buffers,          ONLY : get_buffer
    USE fft_base,         ONLY : dffts
    USE noncollin_module, ONLY : npol
    USE w90_interface,    ONLY : num_wann
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iR
    !! Real space lattice vector index
    COMPLEX(DP), INTENT(inout) :: wan_Rr(dffts%nnr, npol, num_wann)
    !! Wavefunction in the Wannier gauge in real space
    !
    CALL start_clock("w90_buf_get")
    !
    CALL get_buffer(wan_Rr, lrwan_r, iuwan_r, iR)
    !
    CALL stop_clock("w90_buf_get")
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE w90_wan_Rr_get_buffer
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_wan_Rr_close_buffer()
    !----------------------------------------------------------------------------
    !! Close the buffer for real-space Wannier functions.
    !! The buffer is created by w90_write_wan_Rr_from_collected_to_buffer.
    !----------------------------------------------------------------------------
    USE buffers,          ONLY : close_buffer
    !
    IMPLICIT NONE
    !
    CALL close_buffer(iuwan_r, 'DELETE')
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE w90_wan_Rr_close_buffer
  !----------------------------------------------------------------------------
  !
END MODULE w90_wan_Rr_buffer
