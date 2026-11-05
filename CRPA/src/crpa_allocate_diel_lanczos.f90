!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_allocate_diel_lanczos()
  !-----------------------------------------------------------------------
  !! Allocates the arrays of the module crpa_diel_lanczos, used to build the
  !! Lanczos chain for the complex macroscopic dielectric tensor at q=Gamma.
  !! rpert=3 (the three Cartesian polarizations) always, since this chain is
  !! only used for the q=Gamma electric-field response.
  !
  USE kinds,             ONLY : DP
  USE klist,             ONLY : ngk
  USE uspp,              ONLY : okvan, nkb
  USE wvfct,             ONLY : npwx, nbnd
  USE noncollin_module,  ONLY : noncolin, npol
  USE qpoint,            ONLY : nksq
  USE lrus,              ONLY : bbg, bbk, bbnc
  USE control_flags,     ONLY : gamma_only
  USE crpacom,           ONLY : crpa_diel_lanczos_steps, crpa_nfreq, epsilon_w_diel, &
                                crpa_diel_lanczos_steps_ext
  USE crpa_diel_lanczos, ONLY : diel_evc0, diel_sevq0, diel_d0psi,      &
                                diel_d0psi_left,                        &
                                diel_evc1, diel_evc1_new, diel_evc1_old, &
                                diel_sevc1, diel_size_evc1,             &
                                diel_beta_store, diel_gamma_store,      &
                                diel_zeta_store, diel_beta_store_ext,   &
                                diel_gamma_store_ext
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: rpert = 3
  !! Number of Cartesian polarizations (independent Lanczos chains)
  !
  diel_size_evc1 = npwx * npol * nbnd * nksq
  !
  ALLOCATE(diel_evc0(npwx*npol, nbnd, nksq))
  diel_evc0(:,:,:) = (0.0_DP, 0.0_DP)
  !
  IF (okvan) THEN
     ALLOCATE(diel_sevq0(npwx*npol, nbnd, nksq))
     diel_sevq0(:,:,:) = (0.0_DP, 0.0_DP)
  ELSE
     diel_sevq0 => diel_evc0
  ENDIF
  !
  ALLOCATE(diel_d0psi(npwx*npol, nbnd, nksq, rpert))
  diel_d0psi(:,:,:,:) = (0.0_DP, 0.0_DP)
  ALLOCATE(diel_d0psi_left(npwx*npol, nbnd, nksq, rpert))
  diel_d0psi_left(:,:,:,:) = (0.0_DP, 0.0_DP)
  !
  ! 3 non-Hermitian chains advanced in lock-step (see comments in
  ! crpa_diel_lanczos.f90)
  !
  ALLOCATE(diel_evc1_old(npwx*npol, nbnd, nksq, rpert, 2))
  ALLOCATE(diel_evc1(npwx*npol, nbnd, nksq, rpert, 2))
  ALLOCATE(diel_evc1_new(npwx*npol, nbnd, nksq, rpert, 2))
  ALLOCATE(diel_sevc1(npwx*npol, nbnd, nksq, rpert))
  diel_evc1_old(:,:,:,:,:) = (0.0_DP, 0.0_DP)
  diel_evc1(:,:,:,:,:)     = (0.0_DP, 0.0_DP)
  diel_evc1_new(:,:,:,:,:) = (0.0_DP, 0.0_DP)
  diel_sevc1(:,:,:,:)      = (0.0_DP, 0.0_DP)
  !
  ALLOCATE(diel_beta_store(rpert, crpa_diel_lanczos_steps))
  ALLOCATE(diel_gamma_store(rpert, crpa_diel_lanczos_steps))
  ALLOCATE(diel_zeta_store(rpert, rpert, crpa_diel_lanczos_steps))
  ALLOCATE(diel_beta_store_ext(rpert, crpa_diel_lanczos_steps_ext))
  ALLOCATE(diel_gamma_store_ext(rpert, crpa_diel_lanczos_steps_ext))
  !
  IF (ALLOCATED(epsilon_w_diel)) DEALLOCATE(epsilon_w_diel)
  ALLOCATE(epsilon_w_diel(3, 3, crpa_nfreq))
  epsilon_w_diel(:,:,:) = (0.0_DP, 0.0_DP)
  !
  ! Coefficients of S^-1, needed by lr_sm1_psi (LR_Modules/lr_sm1_psi.f90)
  !
  IF (okvan) THEN
     IF (gamma_only) THEN
        ALLOCATE(bbg(nkb, nkb))
     ELSEIF (noncolin) THEN
        ALLOCATE(bbnc(nkb*npol, nkb*npol, nksq))
     ELSE
        ALLOCATE(bbk(nkb, nkb, nksq))
     ENDIF
     CALL lr_sm1_initialize()
  ENDIF
  !
  RETURN
  !
END SUBROUTINE crpa_allocate_diel_lanczos
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_deallocate_diel_lanczos()
  !-----------------------------------------------------------------------
  !! Deallocates the arrays of the module crpa_diel_lanczos.
  !
  USE uspp,              ONLY : okvan
  USE lrus,              ONLY : bbg, bbk, bbnc
  USE crpa_diel_lanczos, ONLY : diel_evc0, diel_sevq0, diel_d0psi,      &
                                diel_d0psi_left,                        &
                                diel_evc1, diel_evc1_new, diel_evc1_old, &
                                diel_sevc1, diel_beta_store,            &
                                diel_gamma_store, diel_zeta_store,      &
                                diel_beta_store_ext, diel_gamma_store_ext
  !
  IMPLICIT NONE
  !
  IF (okvan) THEN
     IF (ASSOCIATED(diel_sevq0)) DEALLOCATE(diel_sevq0)
  ELSE
     NULLIFY(diel_sevq0)
  ENDIF
  !
  IF (ALLOCATED(diel_evc0))     DEALLOCATE(diel_evc0)
  IF (ALLOCATED(diel_d0psi))    DEALLOCATE(diel_d0psi)
  IF (ALLOCATED(diel_d0psi_left)) DEALLOCATE(diel_d0psi_left)
  IF (ALLOCATED(diel_evc1_old)) DEALLOCATE(diel_evc1_old)
  IF (ALLOCATED(diel_evc1))     DEALLOCATE(diel_evc1)
  IF (ALLOCATED(diel_evc1_new)) DEALLOCATE(diel_evc1_new)
  IF (ALLOCATED(diel_sevc1))    DEALLOCATE(diel_sevc1)
  IF (ALLOCATED(diel_beta_store))  DEALLOCATE(diel_beta_store)
  IF (ALLOCATED(diel_gamma_store)) DEALLOCATE(diel_gamma_store)
  IF (ALLOCATED(diel_zeta_store))  DEALLOCATE(diel_zeta_store)
  IF (ALLOCATED(diel_beta_store_ext))  DEALLOCATE(diel_beta_store_ext)
  IF (ALLOCATED(diel_gamma_store_ext)) DEALLOCATE(diel_gamma_store_ext)
  !
  IF (ALLOCATED(bbg))  DEALLOCATE(bbg)
  IF (ALLOCATED(bbk))  DEALLOCATE(bbk)
  IF (ALLOCATED(bbnc)) DEALLOCATE(bbnc)
  !
  RETURN
  !
END SUBROUTINE crpa_deallocate_diel_lanczos
