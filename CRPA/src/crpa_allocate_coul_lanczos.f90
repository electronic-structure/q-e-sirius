!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_allocate_coul_lanczos()
  !-----------------------------------------------------------------------
  !! Allocates the arrays of the module crpa_coul_lanczos, used to build the
  !! Lanczos chains for the screened Coulomb matrix elements at the current q.
  !
  USE kinds,             ONLY : DP
  USE klist,             ONLY : ngk
  USE uspp,              ONLY : okvan, nkb
  USE wvfct,             ONLY : npwx, nbnd
  USE noncollin_module,  ONLY : noncolin, npol
  USE qpoint,            ONLY : nksq
  USE lrus,              ONLY : bbg, bbk, bbnc
  USE control_flags,     ONLY : gamma_only
  USE crpacom,           ONLY : crpa_coul_lanczos_steps, crpa_coul_lanczos_steps_ext
  USE crpa_pert,         ONLY : npert_tot, nmels_tot
  USE crpa_coul_lanczos, ONLY : coul_evc0, coul_evq0, coul_sevq0, coul_d0psi_mels, &
                                coul_evc1, coul_evc1_new, coul_evc1_old,           &
                                coul_sevc1, coul_beta_store,                       &
                                coul_gamma_store, coul_zeta_store,                 &
                                coul_beta_store_ext, coul_gamma_store_ext
  !
  IMPLICIT NONE
  !
  ALLOCATE(coul_evc0(npwx*npol, nbnd, nksq))
  ALLOCATE(coul_evq0(npwx*npol, nbnd, nksq))
  coul_evc0(:,:,:) = (0.0_DP, 0.0_DP)
  coul_evq0(:,:,:) = (0.0_DP, 0.0_DP)
  !
  IF (okvan) THEN
     ALLOCATE(coul_sevq0(npwx*npol, nbnd, nksq))
     coul_sevq0(:,:,:) = (0.0_DP, 0.0_DP)
  ELSE
     coul_sevq0 => coul_evq0
  ENDIF
  !
  ALLOCATE(coul_d0psi_mels(npwx*npol, nbnd, nksq, nmels_tot))
  coul_d0psi_mels(:,:,:,:) = (0.0_DP, 0.0_DP)
  !
  ALLOCATE(coul_evc1_old(npwx*npol, nbnd, nksq, 2))
  ALLOCATE(coul_evc1(npwx*npol, nbnd, nksq, 2))
  ALLOCATE(coul_evc1_new(npwx*npol, nbnd, nksq, 2))
  ALLOCATE(coul_sevc1(npwx*npol, nbnd, nksq))
  coul_evc1_old(:,:,:,:) = (0.0_DP, 0.0_DP)
  coul_evc1(:,:,:,:)     = (0.0_DP, 0.0_DP)
  coul_evc1_new(:,:,:,:) = (0.0_DP, 0.0_DP)
  coul_sevc1(:,:,:)      = (0.0_DP, 0.0_DP)
  !
  ALLOCATE(coul_beta_store(npert_tot, crpa_coul_lanczos_steps))
  ALLOCATE(coul_gamma_store(npert_tot, crpa_coul_lanczos_steps))
  ALLOCATE(coul_zeta_store(npert_tot, nmels_tot, crpa_coul_lanczos_steps))
  ALLOCATE(coul_beta_store_ext(npert_tot, crpa_coul_lanczos_steps_ext))
  ALLOCATE(coul_gamma_store_ext(npert_tot, crpa_coul_lanczos_steps_ext))
  !
  ! Coefficients of S^-1, needed by lr_sm1_psi (LR_Modules/lr_sm1_psi.f90).
  ! Guarded, since crpa_epsil_only runs are a separate mode and never reach
  ! this routine, but a run could in principle also have crpa_diel_llanczos
  ! active earlier in the SAME q-point loop iteration (not currently
  ! possible, since crpa_epsil_only always exits before the q-loop body
  ! that calls this routine, but guard anyway to be safe).
  !
  IF (okvan) THEN
     IF (gamma_only) THEN
        IF (.NOT. ALLOCATED(bbg)) ALLOCATE(bbg(nkb, nkb))
     ELSEIF (noncolin) THEN
        IF (.NOT. ALLOCATED(bbnc)) ALLOCATE(bbnc(nkb*npol, nkb*npol, nksq))
     ELSE
        IF (.NOT. ALLOCATED(bbk)) ALLOCATE(bbk(nkb, nkb, nksq))
     ENDIF
     CALL lr_sm1_initialize()
  ENDIF
  !
  RETURN
  !
END SUBROUTINE crpa_allocate_coul_lanczos
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_deallocate_coul_lanczos()
  !-----------------------------------------------------------------------
  !! Deallocates the per-q working arrays of the module crpa_coul_lanczos.
  !
  USE uspp,              ONLY : okvan
  USE crpa_coul_lanczos, ONLY : coul_evc0, coul_evq0, coul_sevq0, coul_d0psi_mels, &
                                coul_evc1, coul_evc1_new, coul_evc1_old,           &
                                coul_sevc1, coul_beta_store,                       &
                                coul_gamma_store, coul_zeta_store,                 &
                                coul_beta_store_ext, coul_gamma_store_ext
  !
  IMPLICIT NONE
  !
  ! coul_sevq0 is set together with coul_evq0 in crpa_allocate_coul_lanczos:
  ! nothing to release if that routine did not run for this q point (q-level
  ! restart from the chain file)
  IF (ALLOCATED(coul_evq0)) THEN
     IF (okvan) THEN
        DEALLOCATE(coul_sevq0)
     ELSE
        NULLIFY(coul_sevq0)
     ENDIF
  ENDIF
  !
  IF (ALLOCATED(coul_evc0))       DEALLOCATE(coul_evc0)
  IF (ALLOCATED(coul_evq0))       DEALLOCATE(coul_evq0)
  IF (ALLOCATED(coul_d0psi_mels)) DEALLOCATE(coul_d0psi_mels)
  IF (ALLOCATED(coul_evc1_old))   DEALLOCATE(coul_evc1_old)
  IF (ALLOCATED(coul_evc1))       DEALLOCATE(coul_evc1)
  IF (ALLOCATED(coul_evc1_new))   DEALLOCATE(coul_evc1_new)
  IF (ALLOCATED(coul_sevc1))      DEALLOCATE(coul_sevc1)
  IF (ALLOCATED(coul_beta_store))  DEALLOCATE(coul_beta_store)
  IF (ALLOCATED(coul_gamma_store)) DEALLOCATE(coul_gamma_store)
  IF (ALLOCATED(coul_zeta_store))  DEALLOCATE(coul_zeta_store)
  IF (ALLOCATED(coul_beta_store_ext))  DEALLOCATE(coul_beta_store_ext)
  IF (ALLOCATED(coul_gamma_store_ext)) DEALLOCATE(coul_gamma_store_ext)
  !
  RETURN
  !
END SUBROUTINE crpa_deallocate_coul_lanczos
