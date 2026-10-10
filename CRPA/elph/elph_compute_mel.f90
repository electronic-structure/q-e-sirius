!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------------
SUBROUTINE elph_compute_mel()
  !----------------------------------------------------------------------------
  !! Compute the electron-phonon matrix elements
  !
  !! $$ \text{epmat}= \langle\psi(k+q) | dV_{SCF}/du^q_{i a} | \psi(k)\rangle $$
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi, rytoev
  USE mp,               ONLY : mp_sum
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE io_global,        ONLY : ionode
  USE io_files,         ONLY : restart_dir
  USE buffers,          ONLY : get_buffer
  USE fft_base,         ONLY : dffts
  USE noncollin_module, ONLY : nspin_mag
  USE klist,            ONLY : xk, igk_k, ngk
  USE wvfct,            ONLY : nbnd, npwx
  USE wavefunctions,    ONLY : evc
  USE lsda_mod,         ONLY : lsda, current_spin, isk
  USE noncollin_module, ONLY : npol
  USE uspp,             ONLY : vkb
  USE uspp_init,        ONLY : init_us_2
  USE modes,            ONLY : nmodes
  USE lrus,             ONLY : becp1
  USE phus,             ONLY : alphap
  USE qpoint,           ONLY : ikks, ikqs, nksq
  USE units_lr,         ONLY : iuwfc, lrwfc
  USE eqv,              ONLY : dvpsi, evq
  USE apply_dpot_mod,   ONLY : apply_dpot_allocate, apply_dpot_deallocate, apply_dpot_bands
  USE control_ph,       ONLY : current_iq
  USE w90_interface,    ONLY : num_wann
  USE elphcom,          ONLY : folder_elph
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: filename
  INTEGER :: ik, ikk, ikq
  !! k point index
  ! INTEGER :: ik_global
  ! !! k point index in the global index for Wannier functions
  INTEGER :: imode
  !! Mode index
  INTEGER :: ibnd, jbnd
  !! Band indices
  INTEGER :: nbnd_elph
  !! Number of bands for electron-phonon matrix elements
  !! FIXME: Currently, we calculate all nbnd, and use indices 1 to nbnd_elph.
  INTEGER :: npw, npwq
  !! Number of plane waves
  INTEGER :: iun
  !! Unit number for the file
  INTEGER :: recl
  !! Record length
  INTEGER :: ios
  !! IO status
  INTEGER :: direct_io_factor
  !! Direct IO factor
  INTEGER*8 :: unf_recl
  !! double precision to prevent integer overflow
  REAL(DP) :: dummy
  !! Dummy variable for direct_io_factor
  COMPLEX(DP) :: upert(nmodes)
  !! Temporary variable for the phonon perturbation vector
  COMPLEX(DP), ALLOCATABLE :: epmat(:, :, :)
  !! Electron-phonon matrix elements
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  !! Auxiliary wavefunction buffer
  COMPLEX(DP), ALLOCATABLE :: dvscfs(:, :, :)
  !! dvscf at the iq-th q point, in the atomic Cartesian basis, on the soft grid
  !! Size (dffts%nnr, nspin_mag, nmodes)
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  ! FIXME: Set nbnd_elph in init
  nbnd_elph = num_wann
  !
  CALL apply_dpot_allocate()
  ALLOCATE(epmat(nbnd, nbnd, nmodes))
  ALLOCATE(aux(npwx * npol, nbnd))
  ALLOCATE(dvscfs(dffts%nnr, nspin_mag, nmodes))
  !
  ! Read the self-consistent potential for the phonon perturbation
  !
  CALL elph_read_dvscf(current_iq, dvscfs)
  !
  ! Open binary file for the electron-phonon matrix elements
  !
  recl = 2 * nbnd_elph * nbnd_elph * nmodes
  filename = TRIM(folder_elph) // "elph.iq" // TRIM(int_to_char(current_iq))
  !
  INQUIRE(IOLENGTH = direct_io_factor) dummy
  unf_recl = direct_io_factor * INT(recl, KIND = KIND(unf_recl))
  IF (unf_recl <= 0) CALL errore('epw_write', 'wrong record length', 3)
  !
  IF (ionode) THEN
    iun = find_free_unit()
    OPEN(UNIT=iun, FILE = TRIM(ADJUSTL(filename)), IOSTAT = ios, FORM = 'unformatted', &
          STATUS = 'unknown', ACCESS = 'direct', RECL = unf_recl)
    IF (ios /= 0) CALL errore('elph_compute_mel', 'error opening ' // TRIM(filename), 1)
  ENDIF
  !
  DO ik = 1, nksq
    !
    ikk = ikks(ik)
    ikq = ikqs(ik)
    IF (lsda) current_spin = isk (ikk)
    npw = ngk(ikk)
    npwq= ngk(ikq)
    !
    ! compute beta functions for k-point ikq
    !
    CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb)
    !
    ! Load wavefunction for the current k point
    !
    CALL get_buffer(evc, lrwfc, iuwfc, ikk)
    CALL get_buffer(evq, lrwfc, iuwfc, ikq)
    !
    DO imode = 1, nmodes
      !
      ! Perturbation vector in the atomic Cartesian basis
      !
      upert(:) = (0.d0, 0.d0)
      upert(imode) = (1.d0, 0.d0)
      !
      ! Compute the bare phonon potential contribution
      ! (Output is stored in variable dvpsi in module eqv)
      !
      CALL dvqpsi_us(ik, upert, .FALSE., becp1, alphap)
      !
      ! Add the self-consistent potential contribution
      !
      CALL apply_dpot_bands(ik, nbnd, dvscfs(:, :, imode), evc, aux)
      dvpsi = dvpsi + aux
      !
      ! Calculate epmat(m, n) = <evc_{k+q,m} | dvpsi_{k,n}> for this perturbation
      !
      CALL ZGEMM('C', 'N', nbnd, nbnd, npwx*npol, (1.d0, 0.d0), evq, npwx*npol, &
                 dvpsi, npwx*npol, (0.d0, 0.d0), epmat(:, :, imode), nbnd)
      !
    ENDDO ! imode
    !
    CALL mp_sum(epmat, intra_bgrp_comm)
    !
    IF (ionode) CALL davcio(epmat(1:nbnd_elph, 1:nbnd_elph, :), recl, iun, ik, +1)
    !
  ENDDO ! ik
  !
  CALL apply_dpot_deallocate()
  !
  IF (ionode) CLOSE(iun, STATUS = 'keep')
  !
END SUBROUTINE elph_compute_mel
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
SUBROUTINE elph_read_dvscf(iq, dvscfs_cart)
  !----------------------------------------------------------------------------
  !! Read folder_phonon/prefix.dvscf files and rotate them from the pattern basis
  !! to the atomic Cartesian basis.
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE io_global,        ONLY : ionode
  USE io_files,         ONLY : prefix, diropn
  USE mp,               ONLY : mp_sum, mp_barrier, mp_bcast
  USE fft_interfaces,   ONLY : fft_interpolate
  USE fft_base,         ONLY : dffts, dfftp
  USE scatter_mod,      ONLY : scatter_grid, gather_grid
  USE gvecs,            ONLY : doublegrid
  USE noncollin_module, ONLY : nspin_mag
  USE control_lr,       ONLY : lgamma
  USE ph_restart,       ONLY : ph_readfile
  USE modes,            ONLY : u, nmodes
  USE control_ph,       ONLY : tmp_dir_ph
  USE elphcom,          ONLY : fildvscf, folder_phonon
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !! index of q point to read
  COMPLEX(DP), INTENT(INOUT) :: dvscfs_cart(dffts%nnr, nspin_mag, nmodes)
  !! dvscf at the iq-th q point, in the atomic Cartesian basis, on the soft grid
  !! Size (dffts%nnr, nspin_mag, nmodes)
  !
  CHARACTER(LEN=256) :: folder
  !! directory where dvscf file for iq is located
  CHARACTER(LEN=256) :: tmp_dir_ph_bak
  !! backup of tmp_dir_ph
  LOGICAL :: exst
  !! on output of diropn, exst is .True. if opened file exists
  INTEGER :: ierr
  !! error variable
  INTEGER :: iudvscf
  !! unit for reading dvscf
  INTEGER :: is
  !! index of spin
  INTEGER :: imode, jmode
  !! index of mode
  INTEGER :: lrdrho
  !! the length of the deltarho files ( = length of dvscf files)
  !! dvscf rotated from iq to iq
  COMPLEX(DP), ALLOCATABLE :: dvscfp(:, :)
  !! Phonon potential read from file in pattern basis. Size (dfftp%nnr, nspin_mag)
  COMPLEX(DP), ALLOCATABLE :: dvscfp_cart(:, :, :)
  !! Phonon potential read from file, converted to atomic Cartesian basis.
  !! Size (dfftp%nnr, nspin_mag, nmodes)
  !
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(len=6), EXTERNAL :: int_to_char
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CALL start_clock('read_dvscf')
  !
  ALLOCATE(dvscfp(dfftp%nnr, nspin_mag))
  ALLOCATE(dvscfp_cart(dfftp%nnr, nspin_mag, nmodes))
  !
  lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
  !
  ! ph_readfile uses the global variable tmp_dir_ph to read the phonon data, so we
  ! temporarily set it here. The output will be written to variable u in module modes.
  !
  tmp_dir_ph_bak = tmp_dir_ph
  tmp_dir_ph = folder_phonon
  CALL ph_readfile('data_u', iq, 0, ierr)
  tmp_dir_ph = tmp_dir_ph_bak
  !
  ! Open dvscf file
  !
  IF (iq == 1 .AND. lgamma) THEN
    folder = TRIM(folder_phonon)
  ELSE
    folder = TRIM(folder_phonon) // TRIM(prefix) // '.q_' // TRIM(int_to_char(iq))
  ENDIF
  folder = trimcheck(folder)
  !
  iudvscf = find_free_unit()
  IF (ionode) THEN
    CALL diropn(iudvscf, fildvscf, lrdrho, exst, folder)
    IF (.NOT. exst) CALL errore('read_dvscf', 'dvscf file ' // TRIM(fildvscf) // ' does not exist', 1)
  ENDIF
  !
  ! Read dvscf file ad rotate to atomic Cartesian coordinate
  ! (Adapted from dfile_star.f90)
  !
  dvscfp_cart = (0.d0, 0.d0)
  ! dvscfp = (0.d0, 0.d0)
  !
  DO imode = 1, nmodes
    !
    CALL davcio_drho(dvscfp, lrdrho, iudvscf, imode, -1)
! print*, "SUM(dvscfp_imode) = ", SUM(dvscfp), imode
    !
    DO jmode = 1, nmodes
      dvscfp_cart(:, :, jmode) = dvscfp_cart(:, :, jmode) + CONJG(u(jmode, imode)) * dvscfp(:, :)
    ENDDO
    !
  ENDDO ! imode
! print*, "SUM(dvscfp_cart) = ", SUM(dvscfp_cart)
! print*, "SUM(ABS(u)) = ", SUM(ABS(u))
  !
  IF (ionode) CLOSE(UNIT = iudvscf, STATUS = 'KEEP')
  !
  IF (doublegrid) THEN
    DO is = 1, nspin_mag
      DO imode = 1, nmodes
        CALL fft_interpolate(dfftp, dvscfp_cart(:, is, imode), dffts, dvscfs_cart(:, is, imode))
      ENDDO
    ENDDO
  ELSE
    CALL zcopy(dffts%nnr * nspin_mag * nmodes, dvscfp_cart, 1, dvscfs_cart, 1)
  ENDIF
  !
  DEALLOCATE(dvscfp)
  DEALLOCATE(dvscfp_cart)
  !
  CALL stop_clock('read_dvscf')
  !
!------------------------------------------------------------------------------
END SUBROUTINE elph_read_dvscf
!------------------------------------------------------------------------------
