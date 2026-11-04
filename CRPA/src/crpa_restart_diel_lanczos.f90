!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Checkpoint / restart of the dielectric-function Lanczos chain
! (crpa_do_lanczos_epsilon.f90), as TDDFPT's lr_write_restart / lr_restart.
! A checkpoint written at loop index iter (iter_progress) holds:
!  - coefficients 1..iter-1, in tmp_dir_ph/diel_lanczos_meta.dat (ionode);
!  - in per-process buffers: the (q,p) Lanczos vectors of step iter and the
!    previous pair ('dielevc1'), the dense-grid response potential built from
!    step iter ('dieldvscf') and, for PAW, int3_paw ('dielint3paw'). The
!    caller rebuilds dvscfs and int3 exactly as the chain loop does.
! A chain can be resumed, or extended/truncated to a new
! crpa_diel_lanczos_steps. The restart needs the same input and the same
! number of processes and pools.
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_write_restart_diel(iter_progress, dvscfp)
  !-----------------------------------------------------------------------
  !! Writes a checkpoint at loop index iter_progress. Called every
  !! crpa_diel_lanczos_restart_step iterations (if > 0) and, with
  !! crpa_diel_lanczos_save_final, once at the end of the chain.
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : ionode, ionode_id
  USE mp,                ONLY : mp_bcast
  USE mp_images,         ONLY : intra_image_comm, nproc_image
  USE mp_pools,          ONLY : npool
  USE control_ph,        ONLY : tmp_dir_ph
  USE noncollin_module,  ONLY : nspin_mag
  USE fft_base,          ONLY : dfftp
  USE crpa_diel_lanczos, ONLY : diel_beta_store, diel_gamma_store, diel_zeta_store
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iter_progress
  COMPLEX(DP), INTENT(INOUT) :: dvscfp(dfftp%nnr, nspin_mag, 3)
  !! response potential (dense grid), ready for the next iteration
  !
  INTEGER :: iunit, ios, ncoef
  CHARACTER(LEN=320) :: filename
  !
  ncoef = iter_progress - 1
  filename = TRIM(tmp_dir_ph) // 'diel_lanczos_meta.dat'
  ios = 0
  IF (ionode) THEN
     OPEN(NEWUNIT=iunit, FILE=TRIM(filename), FORM='unformatted', STATUS='replace', IOSTAT=ios)
     IF (ios == 0) THEN
        WRITE(iunit) iter_progress, ncoef, nproc_image, npool
        IF (ncoef > 0) THEN
           WRITE(iunit) diel_beta_store(:, 1:ncoef)
           WRITE(iunit) diel_gamma_store(:, 1:ncoef)
           WRITE(iunit) diel_zeta_store(:, :, 1:ncoef)
        ENDIF
        CLOSE(iunit, STATUS='keep')
     ENDIF
  ENDIF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0) CALL errore('crpa_write_restart_diel', 'cannot write '//TRIM(filename), 1)
  !
  CALL crpa_diel_restart_buffers(dvscfp, .TRUE.)
  !
END SUBROUTINE crpa_write_restart_diel
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_read_restart_diel(iter_progress, ncoef, dvscfp)
  !-----------------------------------------------------------------------
  !! Reads back a checkpoint (crpa_diel_lanczos_restart = .true.); a missing
  !! checkpoint is an error. On return ncoef is the number of coefficients it
  !! holds, the first MIN(ncoef, crpa_diel_lanczos_steps) of which are loaded.
  !! If ncoef >= crpa_diel_lanczos_steps the caller skips the loop
  !! (truncation) and the buffers are not read.
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : ionode, ionode_id
  USE mp,                ONLY : mp_bcast
  USE mp_images,         ONLY : intra_image_comm, nproc_image
  USE mp_pools,          ONLY : npool
  USE control_ph,        ONLY : tmp_dir_ph
  USE noncollin_module,  ONLY : nspin_mag
  USE fft_base,          ONLY : dfftp
  USE crpacom,           ONLY : crpa_diel_lanczos_steps
  USE crpa_diel_lanczos, ONLY : diel_beta_store, diel_gamma_store, diel_zeta_store
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: iter_progress, ncoef
  COMPLEX(DP), INTENT(OUT) :: dvscfp(dfftp%nnr, nspin_mag, 3)
  !
  CHARACTER(LEN=*), PARAMETER :: sub = 'crpa_read_restart_diel'
  INTEGER :: iunit, ierr, nload, lay(2)
  REAL(DP), ALLOCATABLE :: btmp(:,:), gtmp(:,:)
  COMPLEX(DP), ALLOCATABLE :: ztmp(:,:,:)
  LOGICAL :: exst
  CHARACTER(LEN=320) :: filename
  !
  iter_progress = 0
  ncoef = 0
  lay = 0
  dvscfp = (0.0_DP, 0.0_DP)
  !
  filename = TRIM(tmp_dir_ph) // 'diel_lanczos_meta.dat'
  ierr = 0
  IF (ionode) THEN
     INQUIRE(FILE=TRIM(filename), EXIST=exst)
     ierr = 1
     IF (exst) OPEN(NEWUNIT=iunit, FILE=TRIM(filename), FORM='unformatted', STATUS='old', IOSTAT=ierr)
     IF (ierr == 0) READ(iunit, IOSTAT=ierr) iter_progress, ncoef, lay
     IF (ierr == 0 .AND. ncoef > 0) THEN
        ALLOCATE(btmp(3, ncoef), gtmp(3, ncoef), ztmp(3, 3, ncoef))
        READ(iunit, IOSTAT=ierr) btmp
        IF (ierr == 0) READ(iunit, IOSTAT=ierr) gtmp
        IF (ierr == 0) READ(iunit, IOSTAT=ierr) ztmp
        nload = MIN(ncoef, crpa_diel_lanczos_steps)
        diel_beta_store(:, 1:nload)    = btmp(:, 1:nload)
        diel_gamma_store(:, 1:nload)   = gtmp(:, 1:nload)
        diel_zeta_store(:, :, 1:nload) = ztmp(:, :, 1:nload)
        DEALLOCATE(btmp, gtmp, ztmp)
     ENDIF
     IF (exst) CLOSE(iunit, STATUS='keep')
  ENDIF
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  IF (ierr == 1) CALL errore(sub, 'crpa_diel_lanczos_restart = .true. but no checkpoint '// &
     TRIM(filename)//' was found', 1)
  IF (ierr /= 0) CALL errore(sub, 'error reading '//TRIM(filename), ABS(ierr))
  !
  CALL mp_bcast(iter_progress, ionode_id, intra_image_comm)
  CALL mp_bcast(ncoef, ionode_id, intra_image_comm)
  CALL mp_bcast(lay, ionode_id, intra_image_comm)
  CALL mp_bcast(diel_beta_store, ionode_id, intra_image_comm)
  CALL mp_bcast(diel_gamma_store, ionode_id, intra_image_comm)
  CALL mp_bcast(diel_zeta_store, ionode_id, intra_image_comm)
  !
  IF (ncoef >= crpa_diel_lanczos_steps) RETURN
  !
  IF (lay(1) /= nproc_image .OR. lay(2) /= npool) CALL errore(sub, 'the checkpoint was '// &
     'written with a different number of MPI processes or pools', 1)
  CALL crpa_diel_restart_buffers(dvscfp, .FALSE.)
  !
END SUBROUTINE crpa_read_restart_diel
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_diel_restart_buffers(dvscfp, lwrite)
  !-----------------------------------------------------------------------
  !! Writes (lwrite) or reads the per-process buffers of the checkpoint.
  !
  USE kinds,             ONLY : DP
  USE buffers,           ONLY : open_buffer, save_buffer, get_buffer, close_buffer
  USE wvfct,             ONLY : npwx, nbnd
  USE noncollin_module,  ONLY : npol, nspin_mag
  USE qpoint,            ONLY : nksq
  USE fft_base,          ONLY : dfftp
  USE ions_base,         ONLY : nat
  USE uspp_param,        ONLY : nhm
  USE paw_variables,     ONLY : okpaw
  USE lrus,              ONLY : int3_paw
  USE crpa_diel_lanczos, ONLY : diel_evc1, diel_evc1_old
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: dvscfp(dfftp%nnr, nspin_mag, 3)
  LOGICAL, INTENT(IN) :: lwrite
  !
  INTEGER :: nword
  LOGICAL :: exst
  !
  nword = npwx * npol * nbnd * nksq * 3 * 2
  CALL open_buffer(91, 'dielevc1', nword, 1, exst)
  IF (lwrite) THEN
     CALL save_buffer(diel_evc1(1,1,1,1,1), nword, 91, 1)
     CALL save_buffer(diel_evc1_old(1,1,1,1,1), nword, 91, 2)
  ELSE
     IF (.NOT. exst) CALL errore('crpa_diel_restart_buffers', 'dielevc1 not found', 1)
     CALL get_buffer(diel_evc1(1,1,1,1,1), nword, 91, 1)
     CALL get_buffer(diel_evc1_old(1,1,1,1,1), nword, 91, 2)
  ENDIF
  CALL close_buffer(91, 'keep')
  !
  nword = dfftp%nnr * nspin_mag * 3
  CALL open_buffer(92, 'dieldvscf', nword, 1, exst)
  IF (lwrite) THEN
     CALL save_buffer(dvscfp, nword, 92, 1)
  ELSE
     IF (.NOT. exst) CALL errore('crpa_diel_restart_buffers', 'dieldvscf not found', 1)
     CALL get_buffer(dvscfp, nword, 92, 1)
  ENDIF
  CALL close_buffer(92, 'keep')
  !
  IF (okpaw) THEN
     nword = nhm * nhm * nat * nspin_mag * 3
     CALL open_buffer(93, 'dielint3paw', nword, 1, exst)
     IF (lwrite) THEN
        CALL save_buffer(int3_paw(1,1,1,1,1), nword, 93, 1)
     ELSE
        IF (.NOT. exst) CALL errore('crpa_diel_restart_buffers', 'dielint3paw not found', 1)
        CALL get_buffer(int3_paw(1,1,1,1,1), nword, 93, 1)
     ENDIF
     CALL close_buffer(93, 'keep')
  ENDIF
  !
END SUBROUTINE crpa_diel_restart_buffers
