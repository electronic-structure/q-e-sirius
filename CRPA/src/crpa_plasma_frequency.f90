!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_plasma_frequency()
  !-----------------------------------------------------------------------
  !! Computes the intraband (Drude) plasma-frequency tensor omega_p^2(3,3),
  !! for metals only, and stores it in crpacom::plasma_freq2. This is the
  !! q=Gamma, same-band contribution that the interband CG/Lanczos machinery
  !! (crpa_solve_e_w/crpa_do_lanczos_epsilon, via dvpsi_e -> orthogonalize)
  !! deliberately excludes: for lgauss/ltetra, orthogonalize() subtracts the
  !! full fractional-occupation-weighted projection of [H,r]|psi_nk> onto
  !! EVERY band in nbnd_occ(ikk), including the band itself, so the diagonal
  !! (intraband) velocity matrix element never reaches the Sternheimer/
  !! Lanczos response. It must be added back separately, analytically:
  !!
  !!   epsilon_head(w) += - omega_p^2 / ( w * (w + i*eta) )
  !!
  !! following the standard metallic-RPA split (e.g. RESPACK's chiqw.F,
  !! wpl/wd construction), evaluated here as the Fermi-surface integral
  !!
  !!   omega_p^2_ab = (fpi*e2/omega) * sum_k wk(k) * sum_n
  !!                  ( -df/dE )|_{n,k} * v_a(n,k) * v_b(n,k)
  !!
  !! with the smearing derivative (-df/dE) given by w0gauss(...)/degauss for
  !! lgauss, or dfpt_tetra_delta(n,k) for ltetra (same convention used by
  !! LR_Modules/efermi_shift.f90). v_a(n,k) = i * <psi_nk|[H,r.at(:,a)]|psi_nk>
  !! is the diagonal velocity matrix element in crystal axis "a", obtained
  !! from the same commutator_Hx_psi() used by dvpsi_e.f90 (no Sternheimer
  !! solve needed: this is a single non-iterative pass over bands/k-points).
  !!
  !! The (fpi*e2/omega) prefactor and crystal-axis convention match
  !! crpa_calc_chi_epsilon.f90's "fact", so that plasma_freq2 can be added to
  !! epsilon_w(_diel)'s head directly, before crys_to_cart/symmatrix, exactly
  !! like epsilon_w itself. NOTE: the overall prefactor should be validated
  !! against a free-electron-gas reference before trusting production
  !! results (see e.g. RESPACK's own "factor of 2" note in chiqw_sub.F).
  !!
  !! CAVEAT: an unshifted k-mesh whose points are ALL time-reversal-invariant
  !! momenta (TRIM: every component 0 or 0.5) gives EXACTLY zero group
  !! velocity for every band at every k-point in a centrosymmetric crystal
  !! (v(k) = -v(-k) by time reversal, v(-k) = v(k) by inversion, k = -k+G at
  !! a TRIM => v(k) = 0), hence plasma_freq2 = 0 - not a bug, but a mesh
  !! choice that cannot see any intraband/Drude weight. A shifted mesh (or
  !! any mesh with non-TRIM points) is required to get a nonzero result.
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : fpi, e2, rytoev
  USE io_global,         ONLY : stdout
  USE cell_base,         ONLY : omega, at
  USE klist,             ONLY : xk, wk, ngk, igk_k, lgauss, ltetra, degauss, ngauss
  USE ener,              ONLY : ef
  USE wvfct,             ONLY : npwx, nbnd, et
  USE noncollin_module,  ONLY : npol
  USE buffers,           ONLY : get_buffer
  USE wavefunctions,     ONLY : evc
  USE uspp,              ONLY : vkb, nkb
  USE uspp_init,         ONLY : init_us_2
  USE lsda_mod,          ONLY : lsda, current_spin, isk
  USE becmod,            ONLY : bec_type, allocate_bec_type, deallocate_bec_type
  USE lrus,              ONLY : becp1
  USE qpoint,            ONLY : nksq, ikks
  USE units_lr,          ONLY : iuwfc, lrwfc
  USE control_lr,        ONLY : nbnd_occ
  USE symme,             ONLY : symmatrix, crys_to_cart
  USE dfpt_tetra_mod,    ONLY : dfpt_tetra_delta
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE mp,                ONLY : mp_sum
  USE crpacom,           ONLY : plasma_freq2, crpa_metal, crpa_diel_flat, crpa_diel_deg_thr
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ikk, npw, ibnd, jbnd, ipol, jpol, nocc
  REAL(DP) :: fact, weight_ef
  REAL(DP) :: velocity(3, nbnd)
  ! diagonal velocity matrix element v_a(ibnd,ik), crystal axis a=1,2,3
  COMPLEX(DP), ALLOCATABLE :: dpsi_c(:,:)
  COMPLEX(DP), ALLOCATABLE :: vmat(:,:,:)
  ! crpa_diel_flat only: full velocity matrix vmat(m,n,a) = i<psi_m|[H,r.a]|psi_n>,
  ! needed for the trace over degenerate blocks |e_m - e_n| < crpa_diel_deg_thr
  TYPE(bec_type) :: becp2
  COMPLEX(DP), EXTERNAL :: zdotc
  REAL(DP), EXTERNAL :: w0gauss
  !
  crpa_metal = lgauss .OR. ltetra
  plasma_freq2(:,:) = 0.0_DP
  IF (.NOT. crpa_metal) RETURN
  !
  CALL start_clock('crpa_plasmafreq')
  !
  ALLOCATE(dpsi_c(npwx*npol, nbnd))
  IF (crpa_diel_flat) ALLOCATE(vmat(nbnd, nbnd, 3))
  CALL allocate_bec_type(nkb, nbnd, becp2)
  !
  DO ik = 1, nksq
     ikk = ikks(ik)
     npw = ngk(ikk)
     IF (lsda) current_spin = isk(ikk)
     IF (nksq > 1) CALL get_buffer(evc, lrwfc, iuwfc, ikk)
     CALL init_us_2(npw, igk_k(1,ikk), xk(1,ikk), vkb, .true.)
     !$acc update host(vkb)   ! init_us_2 fills vkb on the device; calbec/zaxpy/commutator below use the host copy
     !
     velocity(:,:) = 0.0_DP
     nocc = nbnd_occ(ikk)
     IF (crpa_diel_flat) vmat(:,:,:) = (0.0_DP, 0.0_DP)
     DO ipol = 1, 3
        dpsi_c(:,:) = (0.0_DP, 0.0_DP)
        CALL commutator_Hx_psi(ikk, nbnd_occ(ikk), at(:,ipol), becp1(ik), becp2, &
                                dpsi_c(:, 1:nbnd_occ(ikk)))
        DO ibnd = 1, nbnd_occ(ikk)
           ! v = i * <psi|[H,r.vpol]|psi>; the commutator is anti-Hermitian,
           ! so this expectation value is purely imaginary and i*(...) is real
           velocity(ipol, ibnd) = &
              DBLE(CMPLX(0.0_DP, 1.0_DP, KIND=DP) * &
                   zdotc(npwx*npol, evc(1,ibnd), 1, dpsi_c(1,ibnd), 1))
        ENDDO
        IF (crpa_diel_flat) CALL zgemm('C', 'N', nocc, nocc, npwx*npol, &
             (0.0_DP, 1.0_DP), evc, npwx*npol, dpsi_c, npwx*npol, &
             (0.0_DP, 0.0_DP), vmat(1,1,ipol), nbnd)
     ENDDO
     !
     ! velocity(ipol,ibnd) = i*<psi|[H,r]|psi> is a plane-wave inner product
     ! (zdotc above), summed only over the G-vectors local to this rank under
     ! intra_bgrp_comm's plane-wave domain decomposition (npwx is the per-rank
     ! local count, not the global one) - it must be reduced to the true,
     ! full-G-vector value here, before being squared/multiplied into
     ! plasma_freq2 below.
     !
     CALL mp_sum(velocity, intra_bgrp_comm)
     IF (crpa_diel_flat) CALL mp_sum(vmat, intra_bgrp_comm)
     !
     DO ibnd = 1, nbnd_occ(ikk)
        IF (ltetra) THEN
           weight_ef = dfpt_tetra_delta(ibnd, ikk)
        ELSE
           weight_ef = w0gauss((ef - et(ibnd,ikk)) / degauss, ngauss) / degauss
        ENDIF
        IF (ABS(weight_ef) < 1.0d-12) CYCLE
        IF (crpa_diel_flat) THEN
           !
           ! trace over the degenerate block of ibnd (jbnd = ibnd is the usual
           ! diagonal term): the pairs removed from the interband part by the
           ! occupation-free projector (same threshold) are counted here
           !
           DO jbnd = 1, nocc
              IF (ABS(et(jbnd,ikk) - et(ibnd,ikk)) >= crpa_diel_deg_thr) CYCLE
              DO ipol = 1, 3
                 DO jpol = 1, 3
                    plasma_freq2(ipol,jpol) = plasma_freq2(ipol,jpol) + wk(ikk) * weight_ef * &
                       DBLE(vmat(ibnd,jbnd,ipol) * CONJG(vmat(ibnd,jbnd,jpol)))
                 ENDDO
              ENDDO
           ENDDO
           CYCLE
        ENDIF
        DO ipol = 1, 3
           DO jpol = 1, 3
              plasma_freq2(ipol,jpol) = plasma_freq2(ipol,jpol) + wk(ikk) * weight_ef * &
                 velocity(ipol,ibnd) * velocity(jpol,ibnd)
           ENDDO
        ENDDO
     ENDDO
     !
  ENDDO ! ik
  !
  IF (nkb > 0) CALL deallocate_bec_type(becp2)
  DEALLOCATE(dpsi_c)
  IF (crpa_diel_flat) DEALLOCATE(vmat)
  !
  CALL mp_sum_plasma_freq2()
  !
  fact = fpi * e2 / omega
  plasma_freq2(:,:) = plasma_freq2(:,:) * fact
  !
  CALL crys_to_cart(plasma_freq2)
  CALL symmatrix(plasma_freq2)
  !
  ! Intraband plasma frequency, cartesian axis: omega_p = sqrt(plasma_freq2).
  !
  WRITE(stdout,'(/,5x,"Intraband (Drude) plasma frequency, cartesian axis")')
  DO ipol = 1, 3
     WRITE(stdout,'(5x,"omega_p^2 (Ry^2) =",3f14.8)') plasma_freq2(ipol,:)
  ENDDO
  WRITE(stdout,'(5x,"omega_p   (eV)   =",3f14.8)') &
     (SQRT(MAX(plasma_freq2(ipol,ipol), 0.0_DP)) * rytoev, ipol = 1, 3)
  !
  CALL stop_clock('crpa_plasmafreq')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE mp_sum_plasma_freq2()
    !! Only inter_pool_comm is reduced here: the k-point loop above is split
    !! across pools (nksq/ikks), so each pool genuinely holds a distinct
    !! partial sum that must be combined.
    USE mp_pools, ONLY : inter_pool_comm
    USE mp,       ONLY : mp_sum
    IMPLICIT NONE
    CALL mp_sum(plasma_freq2, inter_pool_comm)
  END SUBROUTINE mp_sum_plasma_freq2
  !
END SUBROUTINE crpa_plasma_frequency
