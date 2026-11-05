!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Restart / extension of the screened-Coulomb Lanczos chains
! (crpa_do_lanczos_coulomb.f90), with crpa_coul_restart = .true., on three
! levels, since every q point runs its own NSCF and npert_tot chains:
!
!  1. q level (crpa_coul_q_from_chain_file): a q point whose chain file
!     filU.coul_lanczos_chain#iq holds crpa_coul_lanczos_steps steps is
!     rebuilt from it (truncated), as crpa_spectrum_only does.
!
!  2. perturbation level: after every finished chain the image's ionode
!     writes tmp_dir_crpa/coul_lanczos_meta_q#iq.dat with the coefficients of
!     every chain done so far; complete chains are taken from it.
!
!  3. step level (crpa_coul_lanczos_save_final, _restart_step): the Lanczos
!     vectors coul_evc1/coul_evc1_old and the dense-grid response potential
!     dvscfp (+ int3_paw for PAW) of one chain at loop index iter, with
!     coefficients 1..iter-1 in the meta file, so that the chain resumes at
!     iter+1 or is extended. Per-process files prefix.coulevc_q#iq,
!     .couldv_q#iq, .coulpaw_q#iq in tmp_dir_crpa, one record per
!     perturbation.
!
! As for TDDFPT, the restart needs the same input; the vector files also
! need the same number of processes and pools.
!
!-----------------------------------------------------------------------
MODULE crpa_coul_lanczos_restart
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  PRIVATE
  !
  PUBLIC :: crpa_coul_ckpt_init, crpa_coul_ckpt_get, crpa_coul_ckpt_done, &
            crpa_coul_ckpt_save_vec, crpa_coul_ckpt_vec_io, crpa_coul_ckpt_end
  !
  INTEGER :: ck_iq = 0
  !! q point the checkpoint state refers to
  INTEGER, ALLOCATABLE :: ck_ncoef(:)
  !! number of coefficients in coul_*_store for each perturbation
  INTEGER, ALLOCATABLE :: ck_iterp(:)
  !! loop index of each perturbation's vector checkpoint (0: none)
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  FUNCTION meta_file(iq)
    !-----------------------------------------------------------------------
    USE crpacom, ONLY : tmp_dir_crpa
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iq
    CHARACTER(LEN=320) :: meta_file
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    meta_file = TRIM(tmp_dir_crpa) // 'coul_lanczos_meta_q' // TRIM(int_to_char(iq)) // '.dat'
  END FUNCTION meta_file
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_coul_ckpt_init(iq)
    !-----------------------------------------------------------------------
    !! Called once per q point, after coul_*_store is allocated. With
    !! crpa_coul_restart, loads the meta file of this q, if any, into
    !! coul_*_store (truncated to crpa_coul_lanczos_steps).
    !
    USE io_global,         ONLY : ionode, ionode_id, stdout
    USE mp,                ONLY : mp_bcast
    USE mp_images,         ONLY : intra_image_comm, nproc_image
    USE mp_pools,          ONLY : npool
    USE crpacom,           ONLY : crpa_coul_lanczos_steps, crpa_coul_restart
    USE crpa_pert,         ONLY : npert_tot, nmels_tot
    USE crpa_coul_lanczos, ONLY : coul_beta_store, coul_gamma_store, coul_zeta_store
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iq
    !
    CHARACTER(LEN=*), PARAMETER :: sub = 'crpa_coul_ckpt_init'
    INTEGER :: iunit, ierr, npert0, nmels0, lay0(2), ipert, n, nsteps
    LOGICAL :: exst
    CHARACTER(LEN=320) :: filename
    !
    CALL crpa_coul_ckpt_end()
    ck_iq = iq
    ALLOCATE(ck_ncoef(npert_tot), ck_iterp(npert_tot))
    ck_ncoef = 0
    ck_iterp = 0
    IF (.NOT. crpa_coul_restart) RETURN
    !
    nsteps = crpa_coul_lanczos_steps
    filename = meta_file(iq)
    exst = .FALSE.
    ierr = 0
    IF (ionode) THEN
       INQUIRE(FILE=TRIM(filename), EXIST=exst)
       IF (exst) THEN
          OPEN(NEWUNIT=iunit, FILE=TRIM(filename), FORM='unformatted', STATUS='old', IOSTAT=ierr)
          IF (ierr == 0) READ(iunit, IOSTAT=ierr) npert0, nmels0, lay0
          IF (ierr == 0 .AND. (npert0 /= npert_tot .OR. nmels0 /= nmels_tot)) ierr = 2
          IF (ierr == 0) READ(iunit, IOSTAT=ierr) ck_ncoef, ck_iterp
          ! records hold ck_ncoef entries; only the first nsteps are kept
          DO ipert = 1, npert_tot
             n = MIN(ck_ncoef(ipert), nsteps)
             IF (ierr == 0 .AND. n > 0) READ(iunit, IOSTAT=ierr) coul_beta_store(ipert, 1:n)
             IF (ierr == 0 .AND. n > 0) READ(iunit, IOSTAT=ierr) coul_gamma_store(ipert, 1:n)
             IF (ierr == 0 .AND. n > 0) READ(iunit, IOSTAT=ierr) coul_zeta_store(ipert, :, 1:n)
          ENDDO
          CLOSE(iunit)
       ENDIF
    ENDIF
    CALL mp_bcast(exst, ionode_id, intra_image_comm)
    CALL mp_bcast(ierr, ionode_id, intra_image_comm)
    IF (ierr == 2) CALL errore(sub, TRIM(filename)//' does not match this calculation '// &
       '(number of perturbations or matrix elements)', 1)
    IF (ierr /= 0) CALL errore(sub, 'error reading '//TRIM(filename), ABS(ierr))
    IF (.NOT. exst) RETURN
    !
    CALL mp_bcast(lay0, ionode_id, intra_image_comm)
    CALL mp_bcast(ck_ncoef, ionode_id, intra_image_comm)
    CALL mp_bcast(ck_iterp, ionode_id, intra_image_comm)
    CALL mp_bcast(coul_beta_store, ionode_id, intra_image_comm)
    CALL mp_bcast(coul_gamma_store, ionode_id, intra_image_comm)
    CALL mp_bcast(coul_zeta_store, ionode_id, intra_image_comm)
    !
    ! a chain truncated here can no longer use its vectors
    WHERE (ck_ncoef > nsteps) ck_iterp = 0
    ck_ncoef = MIN(ck_ncoef, nsteps)
    !
    ! vector checkpoints are per-process files
    IF (lay0(1) /= nproc_image .OR. lay0(2) /= npool) THEN
       IF (ANY(ck_iterp > 0 .AND. ck_ncoef < nsteps)) CALL errore(sub, 'the Coulomb '// &
          'Lanczos checkpoint '//TRIM(filename)//' was written with a different number '// &
          'of MPI processes or pools', 1)
       ck_iterp = 0
    ENDIF
    !
    WRITE(stdout, '(/,5x,"Coulomb Lanczos checkpoint of q point #",i4," read: ",i5, &
         &" of ",i5," chains complete")') iq, COUNT(ck_ncoef >= nsteps), npert_tot
    !
  END SUBROUTINE crpa_coul_ckpt_init
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_coul_ckpt_get(ipert, iter0)
    !-----------------------------------------------------------------------
    !! First loop index of the chain of perturbation ipert:
    !!  - crpa_coul_lanczos_steps + 2: complete in coul_*_store, skip it;
    !!  - > 1: resume from the vector checkpoint (crpa_coul_ckpt_vec_io),
    !!    coefficients 1..iter0-2 in coul_*_store;
    !!  - 1: start from scratch.
    !
    USE io_global,  ONLY : stdout
    USE crpacom,    ONLY : crpa_coul_lanczos_steps
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ipert
    INTEGER, INTENT(OUT) :: iter0
    !
    iter0 = 1
    IF (ck_ncoef(ipert) >= crpa_coul_lanczos_steps) THEN
       iter0 = crpa_coul_lanczos_steps + 2
       WRITE(stdout, '(5x,"pert",i5,": chain taken from the checkpoint")') ipert
    ELSEIF (ck_iterp(ipert) > 0) THEN
       iter0 = ck_iterp(ipert) + 1
       WRITE(stdout, '(5x,"pert",i5,": chain resumed from the checkpoint after step",i8)') &
            ipert, ck_iterp(ipert) - 1
    ENDIF
    !
  END SUBROUTINE crpa_coul_ckpt_get
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_coul_ckpt_done(ipert)
    !-----------------------------------------------------------------------
    !! The chain of perturbation ipert is complete: record it in the meta file.
    !
    USE crpacom,  ONLY : crpa_coul_lanczos_steps
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ipert
    !
    ck_ncoef(ipert) = crpa_coul_lanczos_steps
    CALL write_meta()
    !
  END SUBROUTINE crpa_coul_ckpt_done
  !
  !-----------------------------------------------------------------------
  SUBROUTINE write_meta()
    !-----------------------------------------------------------------------
    !! (Re)writes the meta file of the current q point (image's ionode).
    !
    USE io_global,         ONLY : ionode, ionode_id
    USE mp,                ONLY : mp_bcast
    USE mp_images,         ONLY : intra_image_comm, nproc_image
    USE mp_pools,          ONLY : npool
    USE crpa_pert,         ONLY : npert_tot, nmels_tot
    USE crpa_coul_lanczos, ONLY : coul_beta_store, coul_gamma_store, coul_zeta_store
    !
    IMPLICIT NONE
    !
    INTEGER :: iunit, ios, ipert, n
    CHARACTER(LEN=320) :: filename
    !
    filename = meta_file(ck_iq)
    ios = 0
    IF (ionode) THEN
       OPEN(NEWUNIT=iunit, FILE=TRIM(filename), FORM='unformatted', STATUS='replace', IOSTAT=ios)
       IF (ios == 0) THEN
          WRITE(iunit) npert_tot, nmels_tot, nproc_image, npool
          WRITE(iunit) ck_ncoef, ck_iterp
          DO ipert = 1, npert_tot
             n = ck_ncoef(ipert)
             IF (n == 0) CYCLE
             WRITE(iunit) coul_beta_store(ipert, 1:n)
             WRITE(iunit) coul_gamma_store(ipert, 1:n)
             WRITE(iunit) coul_zeta_store(ipert, :, 1:n)
          ENDDO
          CLOSE(iunit, STATUS='keep')
       ENDIF
    ENDIF
    CALL mp_bcast(ios, ionode_id, intra_image_comm)
    IF (ios /= 0) CALL errore('crpa_coul_ckpt', 'cannot write '//TRIM(filename), 1)
    !
  END SUBROUTINE write_meta
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_coul_ckpt_save_vec(ipert, iter, dvscfp)
    !-----------------------------------------------------------------------
    !! Vector checkpoint of perturbation ipert at loop index iter, plus the
    !! meta file. The meta file is first rewritten without this checkpoint,
    !! so that an interruption while the vectors are being written cannot
    !! leave a meta file pointing to half-written vectors.
    !
    USE fft_base,          ONLY : dfftp
    USE noncollin_module,  ONLY : nspin_mag
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ipert, iter
    COMPLEX(DP), INTENT(INOUT) :: dvscfp(dfftp%nnr, nspin_mag)
    !
    ck_ncoef(ipert) = iter - 1
    ck_iterp(ipert) = 0
    CALL write_meta()
    CALL crpa_coul_ckpt_vec_io(ipert, dvscfp, +1)
    ck_iterp(ipert) = iter
    CALL write_meta()
    !
  END SUBROUTINE crpa_coul_ckpt_save_vec
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_coul_ckpt_vec_io(ipert, dvscfp, idir)
    !-----------------------------------------------------------------------
    !! Writes (idir = +1) or reads (idir = -1) the vectors of perturbation
    !! ipert: coul_evc1, coul_evc1_old, the dense-grid potential and, for
    !! PAW, int3_paw, in the per-process files of the current q point. After
    !! a read the caller rebuilds the smooth-grid potential and int3.
    !
    USE io_files,          ONLY : diropn
    USE wvfct,             ONLY : npwx, nbnd
    USE noncollin_module,  ONLY : npol, nspin_mag
    USE qpoint,            ONLY : nksq
    USE fft_base,          ONLY : dfftp
    USE ions_base,         ONLY : nat
    USE uspp_param,        ONLY : nhm
    USE paw_variables,     ONLY : okpaw
    USE lrus,              ONLY : int3_paw
    USE crpacom,           ONLY : tmp_dir_crpa
    USE crpa_coul_lanczos, ONLY : coul_evc1, coul_evc1_old
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ipert, idir
    COMPLEX(DP), INTENT(INOUT) :: dvscfp(dfftp%nnr, nspin_mag)
    !! INTENT(INOUT) because davcio's dummy argument is
    !
    INTEGER, EXTERNAL :: find_free_unit
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    INTEGER :: iun, nw
    LOGICAL :: exst
    !
    ! record lengths in real words: two Krylov components, complex
    nw = 4*npwx*npol*nbnd*nksq
    iun = find_free_unit()
    CALL diropn(iun, 'coulevc_q'//TRIM(int_to_char(ck_iq)), nw, exst, tmp_dir_crpa)
    CALL davcio(coul_evc1, nw, iun, 2*ipert-1, idir)
    CALL davcio(coul_evc1_old, nw, iun, 2*ipert, idir)
    CLOSE(iun, STATUS='keep')
    !
    nw = 2*dfftp%nnr*nspin_mag
    iun = find_free_unit()
    CALL diropn(iun, 'couldv_q'//TRIM(int_to_char(ck_iq)), nw, exst, tmp_dir_crpa)
    CALL davcio(dvscfp, nw, iun, ipert, idir)
    CLOSE(iun, STATUS='keep')
    !
    IF (okpaw) THEN
       nw = 2*nhm*nhm*nat*nspin_mag
       iun = find_free_unit()
       CALL diropn(iun, 'coulpaw_q'//TRIM(int_to_char(ck_iq)), nw, exst, tmp_dir_crpa)
       CALL davcio(int3_paw, nw, iun, ipert, idir)
       CLOSE(iun, STATUS='keep')
    ENDIF
    !
  END SUBROUTINE crpa_coul_ckpt_vec_io
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_coul_ckpt_end()
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    IF (ALLOCATED(ck_ncoef)) DEALLOCATE(ck_ncoef)
    IF (ALLOCATED(ck_iterp)) DEALLOCATE(ck_iterp)
    ck_iq = 0
  END SUBROUTINE crpa_coul_ckpt_end
  !
END MODULE crpa_coul_lanczos_restart
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_coul_q_from_chain_file(iq, done)
  !-----------------------------------------------------------------------
  !! q-level restart (crpa_coul_restart): if the chain file of q
  !! point iq holds at least crpa_coul_lanczos_steps steps, rebuild this q
  !! from it - truncated to that length - exactly as crpa_run_spectrum_only
  !! does, and return done = .true.; the caller then skips the NSCF and the
  !! chains. Called after crpa_prepare_q (xq, lgamma). The file is not
  !! rewritten, so a longer chain is kept.
  !
  USE io_global,         ONLY : stdout
  USE crpacom,           ONLY : crpa_coul_lanczos_steps, v_coul_bare, v_coul_scrd
  USE control_lr,        ONLY : lgamma
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  LOGICAL, INTENT(OUT) :: done
  !
  INTEGER :: nfile
  !
  done = .FALSE.
  nfile = 0
  !
  ! the perturbation lists of the current calculation (also allocates
  ! v_coul_bare)
  !
  CALL crpa_prepare_perturbations_q(iq)
  CALL crpa_read_coul_lanczos_chain(iq, done, crpa_coul_lanczos_steps, nfile)
  !
  IF (done) THEN
     WRITE(stdout, '(/,5x,"q point #",i4,": screened Coulomb taken from the chain file",&
          &" (crpa_coul_restart), no NSCF and no Lanczos chains")') iq
     CALL crpa_calc_coulomb_lanczos(iq)
     IF (lgamma) CALL crpa_q0correction()
     CALL crpa_write_coulomb_lanczos(iq)
     CALL crpa_deallocate_coul_lanczos()
  ELSEIF (nfile > 0) THEN
     WRITE(stdout, '(/,5x,"q point #",i4,": the chain file holds",i8," <",i8, &
          &" steps, the chains are continued or recomputed")') iq, nfile, crpa_coul_lanczos_steps
  ENDIF
  !
  ! v_coul_bare/v_coul_scrd are reallocated by crpa_prepare_perturbations_q
  ! (again, if the q point is computed after all)
  !
  IF (ALLOCATED(v_coul_bare)) DEALLOCATE(v_coul_bare)
  IF (ALLOCATED(v_coul_scrd)) DEALLOCATE(v_coul_scrd)
  !
END SUBROUTINE crpa_coul_q_from_chain_file
