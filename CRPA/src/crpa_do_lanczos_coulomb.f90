!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_do_lanczos_coulomb(iq)
  !-----------------------------------------------------------------------
  !! Driver that replaces crpa_solve_linear_system() when crpa_llanczos =
  !! .TRUE. Builds, for each of the npert_tot bare (Wannier-pair)
  !! perturbations at the iq-th q point, an independent non-Hermitian
  !! Lanczos chain for the frequency-independent operator A = S^-1 H,
  !! following the same construction as crpa_do_lanczos_epsilon.f90, but:
  !!
  !!  - the perturbation is a bare Wannier-pair potential at general q
  !!    (crpa_compute_perturbation_q_wannier), not the q=Gamma E-field;
  !!  - the npert_tot chains are run SEQUENTIALLY, not in lock-step: unlike
  !!    the 3 Cartesian E-field polarizations (genuinely mixed into each
  !!    other by the crystal point group, see crpa_do_lanczos_epsilon.f90),
  !!    distinct CRPA perturbations are not mixed by any symmetry operation
  !!    exploited by this code (crpa_set_upert gives every perturbation a
  !!    trivial, unmixed representation - see crpa_symmetry.f90). No
  !!    per-iteration density symmetrization is performed either: the
  !!    calculation is assumed to run without exploiting crystal symmetry.
  !!  - instead of a scalar dielectric tensor, the observable is the
  !!    screened Coulomb matrix element between the bare perturbation and
  !!    every entry of the "mels" (matrix-element) list. This is obtained,
  !!    WITHOUT ever reconstructing the real-space induced potential, from
  !!    the same zeta-projection trick used for the dielectric tensor:
  !!
  !!      v_coul_scrd(j,i,w) = <drho_bare,mels(j) | dvscf_i(w)>
  !!                         = <drho_i(w) | dvbare,mels(j)>            (kernel symmetry)
  !!                         = <d0psi_mels(j) | X_i(w)>                (bilinear density identity)
  !!
  !!    where X_i(w) = (w-A)^-1 d0psi_i is exactly the quantity the Lanczos
  !!    chain gives access to via its continued fraction, and d0psi_mels(j)
  !!    is built from mels(j)'s own bare potential exactly like d0psi_i is
  !!    built from perturbation i's bare potential (both via
  !!    apply_dpot_bands + orthogonalize, as in crpa_do_lanczos_epsilon.f90).
  !!    d0psi_mels is built ONCE per q (nmels_tot vectors), shared by every
  !!    perturbation's chain, and zeta_i(j,n) = <d0psi_mels(j)|q_n> is
  !!    accumulated every iteration for all nmels_tot targets at once.
  !!    See crpa_calc_coulomb_lanczos.f90 for the final continued-fraction
  !!    reconstruction of v_coul_scrd from beta/gamma/zeta.
  !!
  !! v_coul_bare is computed exactly as in crpa_solve_linear_system.f90
  !! (non-iterative, untouched).
  !!
  !! Restrictions of this first implementation:
  !! - pert_basis = 'wannier' only
  !! - no DFPT+U, no noncollinear magnetism (domag), no metals
  !!
  !! Restart / extension (crpa_coul_restart, crpa_coul_lanczos_save_final,
  !! crpa_coul_lanczos_restart_step): see crpa_restart_coul_lanczos.f90.
  !
  USE kinds,                     ONLY : DP
  USE io_global,                 ONLY : stdout
  USE ions_base,                 ONLY : nat
  USE klist,                     ONLY : xk, wk, ngk, igk_k, ltetra
  USE dfpt_tetra_mod,            ONLY : dfpt_tetra_beta, dfpt_tetra_ttheta
  USE gvecs,                     ONLY : doublegrid
  USE fft_base,                  ONLY : dffts, dfftp
  USE fft_interfaces,            ONLY : fft_interpolate
  USE lsda_mod,                  ONLY : lsda, nspin, current_spin, isk
  USE control_flags,             ONLY : gamma_only
  USE noncollin_module,          ONLY : noncolin, npol, nspin_mag, domag
  USE wvfct,                     ONLY : nbnd, npwx, et
  USE wavefunctions,             ONLY : evc
  USE eqv,                       ONLY : dpsi, dvpsi, evq
  USE becmod,                    ONLY : becp, calbec, allocate_bec_type_acc, deallocate_bec_type_acc
  USE uspp,                      ONLY : okvan, vkb, nkb
  USE uspp_param,                ONLY : nhm
  USE paw_variables,             ONLY : okpaw
  USE paw_onecenter,             ONLY : PAW_dpotential
  USE lrus,                      ONLY : int3, int3_paw, int3_nc
  USE scf,                       ONLY : rho
  USE ldaU,                      ONLY : lda_plus_u
  USE control_lr,                ONLY : alpha_pv, nbnd_occ, lgamma, lrpa, lnoloc
  USE units_lr,                  ONLY : iuwfc, lrwfc
  USE buffers,                   ONLY : get_buffer
  USE qpoint,                    ONLY : nksq, ikks, ikqs
  USE dv_of_drho_lr,             ONLY : dv_of_drho
  USE apply_dpot_mod,            ONLY : apply_dpot_bands, apply_dpot_allocate, apply_dpot_deallocate
  USE response_kernels,          ONLY : sternheimer_postprocess
  USE dfpt_type,                 ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data, &
                                         dfpt_dvscfp_to_dvscfs
  USE uspp_init,                 ONLY : init_us_2
  USE crpacom,                   ONLY : crpa_coul_lanczos_steps, dvbare, dvbares, v_coul_bare, &
                                       crpa_xcmode, crpa_coul_lanczos_save_final, &
                                       crpa_coul_lanczos_restart_step
  USE crpa_coul_lanczos_restart, ONLY : crpa_coul_ckpt_init, crpa_coul_ckpt_get, crpa_coul_ckpt_done, &
                                         crpa_coul_ckpt_save_vec, crpa_coul_ckpt_vec_io, crpa_coul_ckpt_end
  USE crpa_pert,                 ONLY : pert_basis, npert_tot, nmels_tot
  USE crpa_coul_lanczos,         ONLY : coul_evc0, coul_evq0, coul_sevq0, coul_d0psi_mels,     &
                                       coul_evc1, coul_evc1_new, coul_evc1_old, coul_sevc1,  &
                                       coul_beta_store, coul_gamma_store, coul_zeta_store
  USE incdrhoscf_mod,            ONLY : incdrhoscf, incdrhoscf_nc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !! index of the q point
  !
  EXTERNAL :: lanczos_nonhermitian
  !
  TYPE(dfpt_data_type) :: dfpt_data
  COMPLEX(DP), ALLOCATABLE :: dbecsum_nc(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: aux2(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvscf_tmp(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvscf_tmps(:,:)
  !
  INTEGER :: iter, ik, ikk, ikq, npwq, jmels, ipert, is, ich, inch
  INTEGER :: iter0
  !! first loop index of the current chain (> 1 when resumed from a checkpoint)
  INTEGER :: nnr
  REAL(DP) :: alpha_pv0
  REAL(DP), ALLOCATABLE :: tetra_beta0(:,:,:)
  REAL(DP) :: alpha, beta, gamma_l
  COMPLEX(DP), ALLOCATABLE :: zeta(:)
  !
  IF (pert_basis /= 'wannier') CALL errore('crpa_do_lanczos_coulomb', &
     'crpa_llanczos = .true. is currently only implemented for pert_basis = wannier', 1)
  IF (lda_plus_u) CALL errore('crpa_do_lanczos_coulomb', &
     'DFPT+U is not yet supported with crpa_llanczos = .true.', 1)
  IF (noncolin .AND. domag) CALL errore('crpa_do_lanczos_coulomb', &
     'noncolin with domag is not yet supported with crpa_llanczos = .true.', 1)
  !
  ! Metals (lgauss/ltetra) are allowed here for general q (q=0 is handled
  ! separately by crpa_do_lanczos_epsilon.f90/crpa_plasma_frequency.f90, not
  ! this routine): dvpsi_e/apply_dpot_bands + orthogonalize() already build
  ! the fractional-occupation-weighted response generically for any q, and
  ! there is no q->0 Drude-type divergence to handle at finite q.
  !
  !
  !
  CALL start_clock('crpa_coul_lncz')
  !
  WRITE(stdout,*) "     =--------------------------------------------="
  WRITE(stdout, '(13x,"SOLVE THE LINEAR SYSTEM (Lanczos, screened Coulomb)")')
  WRITE(stdout,*) "     =--------------------------------------------="
  !
  CALL crpa_allocate_coul_lanczos()
  !
  CALL allocate_dfpt_data(dfpt_data, 1)
  ALLOCATE(dvbare(dfftp%nnr, nspin_mag))
  ALLOCATE(dvbares(dffts%nnr, nspin_mag))
  ALLOCATE(dvscf_tmp(dfftp%nnr, nspin_mag))
  ALLOCATE(dvscf_tmps(dffts%nnr, nspin_mag))
  !
  ! USPP/PAW augmentation integrals, npert = 1 (one chain at a time). Not
  ! allocated anywhere else on this path - see crpa_solve_linear_system.f90.
  !
  IF (okvan) THEN
     ALLOCATE (int3 (nhm, nhm, nat, nspin_mag, 1))
     IF (okpaw)    ALLOCATE (int3_paw (nhm, nhm, nat, nspin_mag, 1))
     IF (noncolin) ALLOCATE (int3_nc  (nhm, nhm, nat, nspin, 1))
  ENDIF
  IF (noncolin) ALLOCATE(dbecsum_nc(nhm, nhm, nat, nspin, 1))
  ALLOCATE(aux2(npwx*npol, nbnd))
  ALLOCATE(zeta(nmels_tot))
  !
  nnr = dffts%nnr
  !$acc enter data create(dfpt_data, dfpt_data%dvscfs(1:nnr, 1:nspin_mag, 1:1))
  !$acc enter data create(dvbare(1:dfftp%nnr, 1:nspin_mag))
  ! smooth-grid copies: apply_dpot_bands (here and in crpa_update_coulomb) only
  ! updates its dv argument on the device, so they must already be present
  !$acc enter data create(dvbares(1:dffts%nnr, 1:nspin_mag))
  !$acc enter data create(dvscf_tmps(1:dffts%nnr, 1:nspin_mag))
  !$acc enter data create(dvscf_tmp(1:dfftp%nnr, 1:nspin_mag))
  !$acc enter data create(aux2)
  !
  CALL allocate_bec_type_acc(nkb, nbnd, becp)
  !
  CALL crpa_set_upert()
  !
  alpha_pv0 = alpha_pv
  !
  ! alpha_pv=0 for the REST of this routine: unlike the CG path
  ! (crpa_solve_linear_system.f90 -> dfpt_kernel -> cgsolve_all), this file
  ! never solves a linear system - ch_psi_all is only ever applied directly,
  ! to build the Krylov chain. LR_Modules/orthogonalize.f90's metallic branch
  ! multiplies alpha_pv into its occupation-difference weight wwg (the
  ! (wgp-wg1)/deltae and w0g terms); that term is only exactly cancelled by
  ! ch_psi_all's alpha_pv*P_v shift when the SAME alpha_pv is used in both
  ! places INSIDE an iterative solve (the standard metal-DFPT alpha_pv
  ! cancellation, Rev. Mod. Phys. 73, 515 (2001)) - with no solve here (every
  ! orthogonalize() call below, including the d0psi_mels/starting-vector
  ! construction at line ~199/245, not just the one inside the chain loop),
  ! there is nothing for alpha_pv to cancel against, so it must be 0
  ! everywhere in this routine, not just around ch_psi_all. Left nonzero,
  ! this leaves an UNCANCELLED alpha_pv contribution in wwg (of order
  ! alpha_pv = emax-emin, i.e. tens of Ry - LR_Modules/setup_alpha_pv.f90)
  ! for a metal (insulators never hit this: orthogonalize's lgauss/ltetra
  ! branch, the only one alpha_pv enters, is not taken there). Same fix
  ! already applied to the cdfpt-constrained variant of this same routine
  ! (LR_Modules/constrained_dfpt.f90) for the same reason - see the
  ! crpa_cdfpt_lanczos_coulomb memory note - and to
  ! crpa_do_lanczos_epsilon.f90's chain loop (that file's d0psi construction
  ! uses dvpsi_e's real cgsolve_all, so is correctly left at alpha_pv0 there).
  !
  alpha_pv = 0.0_DP
  !
  ! With tetrahedra alpha_pv cannot simply be set to 0: dfpt_tetra_main() bakes it
  ! into the weight array once, at setup (LR_Modules/dfpt_tetra_mod.f90:178,
  ! dfpt_tetra_beta = dfpt_tetra_ttheta + alpha_pv*dfpt_tetra_beta), whereas the
  ! lgauss branch of orthogonalize() builds its alpha_pv terms at call time. So the
  ! switch above is paired with a swap to the alpha_pv-free dfpt_tetra_ttheta, the
  ! exact tetrahedron counterpart of the Gaussian wwg. Restored before returning.
  ! Same fix as crpa_do_lanczos_epsilon.f90. lgauss/insulators are unaffected.
  !
  IF (ltetra) THEN
     ALLOCATE(tetra_beta0, SOURCE=dfpt_tetra_beta)
     dfpt_tetra_beta = dfpt_tetra_ttheta
  ENDIF
  !
  ! Cache the unperturbed wavefunctions at k and k+q (and S*evq, for USPP)
  ! for the whole calculation at this q: they do not change during the
  ! Lanczos chains.
  !
  DO ik = 1, nksq
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npwq = ngk(ikq)
     IF (lsda) current_spin = isk(ikk)
     CALL get_buffer(evc, lrwfc, iuwfc, ikk)
     !$acc update device(evc)
     coul_evc0(:,:,ik) = evc(:,:)
     IF (.NOT. lgamma) THEN
        CALL get_buffer(evq, lrwfc, iuwfc, ikq)
        !$acc update device(evq)
        coul_evq0(:,:,ik) = evq(:,:)
     ELSE
        coul_evq0(:,:,ik) = evc(:,:)
     ENDIF
     CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb, .true.)
     !$acc update host(vkb)
     IF (okvan) THEN
        CALL calbec(npwq, vkb, coul_evq0(:,:,ik), becp, nbnd)
        ! host calbec; with OpenACC the s_psi below reads becp on the device
        ! (present_or_copyin), so refresh the device copy if one exists
        IF (gamma_only) THEN
           !$acc update device(becp%r) if_present
        ELSEIF (noncolin) THEN
           !$acc update device(becp%nc) if_present
        ELSE
           !$acc update device(becp%k) if_present
        ENDIF
        CALL s_psi(npwx, npwq, nbnd, coul_evq0(:,:,ik), coul_sevq0(1,1,ik))
     ENDIF
  ENDDO
  !
  ! Checkpoint state of this q point (crpa_restart_coul_lanczos.f90)
  !
  CALL crpa_coul_ckpt_init(iq)
  !
  ! Build the nmels_tot LEFT vectors d0psi_mels = dV_mels(j) |psi_k>, one per
  ! entry of the matrix-element list, shared by every perturbation's chain
  ! (see the module-level comment above for why this gives the screened
  ! Coulomb matrix element without ever reconstructing the real-space induced
  ! potential).
  !
  ! These are deliberately NOT run through orthogonalize(), unlike the chain
  ! SEED below (which must be, since it is the RHS of the Sternheimer
  ! equation). The identity these left vectors implement is
  !
  !   U_ind,ij = <rho_i|v_ext|drho_j> = <dV_i|drho_j>
  !            = sum_nk wk <dV_i psi_nk | dpsi_j,nk+q>
  !
  ! with dV_i = v_ext rho_i (dvscf_tmp here, via crpa_dv_bare: v for 'TCTC',
  ! v+f_xc for 'TETE') applied to the UNPROJECTED |psi>.
  ! Projecting the left vector too applies the response weight a SECOND time.
  ! In the KS eigenbasis, orthogonalize() (LR_Modules/orthogonalize.f90)
  ! multiplies by W_nm = theta_F,n - wwg(m,n) = Theta_nm (f_n - f_m) (times
  ! the 1+alpha_pv/de factor), and dpsi_j already carries exactly one factor
  ! of W_nm from its own RHS. So a projected left vector yields sum_nm W_nm^2
  ! /(de-w) instead of the correct sum_nm W_nm/(de-w).
  !
  ! For an INSULATOR W_nm is 0 or 1, so W^2 = W and both forms agree exactly
  ! - which is why the CG and Lanczos paths matched to ~1e-7 there and this
  ! was never noticed. For a METAL, pairs straddling E_F have 0 < W_nm < 1,
  ! W^2 < W, the induced (screening) term comes out too small, and U is
  ! systematically OVERestimated (measured on SrVO3, q#2, w=0: CG 1.11 eV vs
  ! projected-left-vector Lanczos 4.2-7.8 eV, worsening with chain length
  ! because it converges to the WRONG operator - not a truncation error).
  !
  ! The unprojected form above is the general one and is what the CG path
  ! effectively evaluates (crpa_update_coulomb.f90:crpa_compute_mel_wannier
  ! contracts the converged real-space dV_scf with the Wannier functions,
  ! so it never involves P_c^+ at all).
  !
  ! NOTE (USPP): for okvan the density also has an augmentation part that is
  ! accumulated through dbecsum, not through this plane-wave bilinear, so the
  ! okvan case needs a separate check before being trusted for metals.
  !
  WRITE(stdout, '(/5x, A)') 'Building matrix-element (mels) starting vectors'
  !
  DO jmels = 1, nmels_tot
     !
     CALL crpa_compute_mels_density_wannier(jmels, dvscf_tmp)
     ! Induced external leg of U_ij: v*rho_j for crpa_xcmode='TCTC' (default),
     ! (v+f_xc)*rho_j for 'TETE'.
     CALL crpa_dv_bare(dvscf_tmp)
     !
     ! Smooth-grid copy for apply_dpot_bands, and int3 for the adddvscf below.
     !
     ! USPP/PAW: the zeta bra must carry the AUGMENTED perturbation
     ! V^_j = dV_j + sum_nm |beta_n> int3_nm <beta_m|, exactly like the chain
     ! seed (turboEELS does the same, lr_addus_dvpsi before its dvpsi2 copy).
     !
     CALL crpa_interpolate_dense_to_smooth(dvscf_tmp, dvscf_tmps)
     !
     ! PAW: int3_paw is not yet built for this zeta bra- zero it
     ! so newdq's "IF (okpaw) int3 = int3 + int3_paw" does not pick up stale
     ! state left over from the previous jmels/ipert's converged chain.
     !
     IF (okpaw) int3_paw = (0.0_DP, 0.0_DP)
     IF (okvan) CALL newdq(dvscf_tmp, 1)
     !
     CALL apply_dpot_allocate()
     !$acc enter data create(dvpsi)
     !
     DO ik = 1, nksq
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npwq = ngk(ikq)
        IF (lsda) current_spin = isk(ikk)
        CALL apply_dpot_bands(ik, nbnd, dvscf_tmps, coul_evc0(:,:,ik), dvpsi)
        IF (okvan) THEN
           ! adddvscf needs vkb = beta(k+q); it does not build it itself.
           CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb, .true.)
           !$acc update host(vkb)
           CALL adddvscf(1, ik, .FALSE.)
        ENDIF
        ! no orthogonalize() here - see the long comment above. The sign is
        ! +dvpsi (not -dvpsi): orthogonalize() applied -P_c^+, so the previous
        ! "-dvpsi" was recovering +P_c^+ dV|psi>; without it, +dV|psi> is
        ! already the wanted left vector.
        coul_d0psi_mels(:,:,ik,jmels) = dvpsi(:,:)
     ENDDO
     !
     !$acc exit data delete(dvpsi)
     CALL apply_dpot_deallocate()
     !
  ENDDO
  !
  ! Bare Coulomb matrix elements: computed exactly as in
  ! crpa_solve_linear_system.f90 (non-iterative, untouched).
  !
  WRITE(stdout, '(/,5x,"Non-Hermitian Lanczos chains for the screened", &
                 &" Coulomb matrix elements (",i6," perturbations, sequential)")') npert_tot
  !
  DO ipert = 1, npert_tot
     !
     WRITE(stdout,'(/6x,"q point #",i4,3x,"Perturbation #",i3)') iq, ipert
     !
     CALL crpa_compute_perturbation_q_wannier(ipert, dvbare)
     ! Chain perturbing leg of U_ij: same crpa_dv_bare kernel
     ! as the mels leg (v for 'TCTC', v+f_xc for 'TETE').
     CALL crpa_dv_bare(dvbare)
     !
     ! Smooth-grid copy, and int3 built from dvbare: both the bare matrix
     ! elements (crpa_update_coulomb just below, called with dvscfs = 0) and
     ! the chain seed apply dvbare alone, so this single newdq serves both.
     ! The chain loop later overwrites int3 via newdq(dvscfp) each iteration,
     ! which is what its own screening term needs.
     !
     CALL crpa_interpolate_dense_to_smooth(dvbare, dvbares)
     !
     ! PAW: this bare/seed newdq must not inherit the previous perturbation's
     ! converged int3_paw.
     IF (okpaw) int3_paw = (0.0_DP, 0.0_DP)
     IF (okvan) CALL newdq(dvbare, 1)
     !
     ! v_coul_bare must be computed with dvscfs = 0 (pure bare term): reset
     ! it here, right before the call, since by the end of the PREVIOUS
     ! perturbation's chain (if any) dvscfs holds that chain's converged
     ! self-consistent screening potential, not zero.
     !
     dfpt_data%dvscfs(:,:,:) = (0.0_DP, 0.0_DP)
     CALL crpa_update_coulomb(iq, ipert, v_coul_bare, dfpt_data%dvscfs)
     !
     ! Build the starting vector d0psi_i for this perturbation's chain, and
     ! reset the chain's working arrays (the npert_tot chains run
     ! sequentially, so these are reused, not lock-stepped).
     !
     dfpt_data%dvscfp(:,:,:) = (0.0_DP, 0.0_DP)
     coul_evc1_old(:,:,:,:) = (0.0_DP, 0.0_DP)
     !
     CALL apply_dpot_allocate()
     !$acc enter data create(dvpsi)
     !
     DO ik = 1, nksq
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npwq = ngk(ikq)
        IF (lsda) current_spin = isk(ikk)
        CALL apply_dpot_bands(ik, nbnd, dvbares, coul_evc0(:,:,ik), dvpsi)
        IF (okvan) THEN
           ! adddvscf needs vkb = beta(k+q); it does not build it itself.
           CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb, .true.)
           !$acc update host(vkb)
           CALL adddvscf(1, ik, .FALSE.)
           ! adddvscf works on the host copy; dvpsi is device-resident here
           ! (enter data create above), and orthogonalize works on the device
           !$acc update device(dvpsi)
        ENDIF
        CALL orthogonalize(dvpsi, coul_evq0(:,:,ik), ikk, ikq, coul_sevq0(:,:,ik), npwq, .TRUE.)
        ! dvpsi is present, so orthogonalize's present_or_copy does not copy the
        ! projected vector back: bring it to the host before it is used below
        !$acc update self(dvpsi)
        !
        ! USPP/PAW: the chain seed must be S^-1 P_c^+ dV|psi>, not P_c^+ dV|psi>.
        ! The operator this chain tridiagonalizes is A = S^-1 (H - eS), and S^-1
        ! is duly applied to coul_evc1_new at every iteration below - but it was
        ! missing HERE, on the initial vector, so for okvan the chain started
        ! from a vector inconsistent with its own recursion. Both references do
        ! apply it: thermo_pw/qe/do_lanczos.f90 ("At the first iteration evc1 is
        ! not known and is set to S^-1 P_c^+ V_ext u_kv") and turboEELS
        ! TDDFPT/src/lr_dvpsi_eels.f90 (lr_sm1_psi on dvpsi1 when
        ! calculator=='lanczos'). Our own dielectric chain
        ! (crpa_do_lanczos_epsilon.f90, iter==1 branch) already did it.
        ! NOTE the left vector coul_d0psi_mels must NOT get S^-1 - upstream
        ! keeps dvpsi2/d0psi free of it, and it is built before this loop.
        ! No effect for norm-conserving PPs (okvan=.FALSE.).
        !
        IF (okvan) THEN
           dpsi(:,:) = (0.0_DP, 0.0_DP)
           CALL lr_sm1_psi(ik, npwx, npwq, nbnd_occ(ikk), dvpsi, dpsi)
           dvpsi(:,:) = dpsi(:,:)
        ENDIF
        coul_evc1(:,:,ik,1) = -dvpsi(:,:)
        coul_evc1(:,:,ik,2) = -dvpsi(:,:)
     ENDDO
     !
     !$acc exit data delete(dvpsi)
     CALL apply_dpot_deallocate()
     !
     ! Restart: a chain already complete in the checkpoint is skipped (its
     ! bare matrix elements were recomputed above); one with a vector
     ! checkpoint resumes after it, with the response potential and int3
     ! rebuilt exactly as at the end of every loop iteration below.
     !
     CALL crpa_coul_ckpt_get(ipert, iter0)
     IF (iter0 > crpa_coul_lanczos_steps + 1) CYCLE
     IF (iter0 > 1) THEN
        CALL crpa_coul_ckpt_vec_io(ipert, dfpt_data%dvscfp(:,:,1), -1)
        CALL dfpt_dvscfp_to_dvscfs(dfpt_data)
        IF (okvan) CALL newdq(dfpt_data%dvscfp, 1)
     ENDIF
     !
     DO iter = iter0, crpa_coul_lanczos_steps + 1
        !
        IF (MOD(iter,2) == 0) THEN
           ich = 2; inch = 1
        ELSE
           ich = 1; inch = 2
        ENDIF
        !
        ! The induced density built below accumulates (incdrhoscf): reset it
        ! every iteration, including iter = 1, where the previous chain's last
        ! density build (crpa_coul_lanczos_save_final) may still be there.
        !
        dfpt_data%drhos(:,:,:) = (0.0_DP, 0.0_DP)
        dfpt_data%drhop(:,:,:) = (0.0_DP, 0.0_DP)
        dfpt_data%dbecsum(:,:,:,:) = (0.0_DP, 0.0_DP)
        IF (noncolin) dbecsum_nc(:,:,:,:,:) = (0.0_DP, 0.0_DP)
        !
        IF (iter > 1) THEN
           !
           DO ik = 1, nksq
              !
              ikk = ikks(ik)
              ikq = ikqs(ik)
              npwq = ngk(ikq)
              IF (lsda) current_spin = isk(ikk)
              evc(:,:) = coul_evc0(:,:,ik)
              !$acc update device(evc)
              IF (.NOT. lgamma) THEN
                 evq(:,:) = coul_evq0(:,:,ik)
                 !$acc update device(evq)
              ENDIF
              CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb, .true.)
              !$acc update host(vkb)
              !
              ! S * p_(iter-1), needed for the bi-orthogonal metric
              ! beta = <q|S|p>. Computed from the CURRENT (not yet
              ! updated) vector, before coul_evc1 is overwritten below.
              !
              IF (okvan) THEN
                 CALL calbec(npwq, vkb, coul_evc1(:,:,ik,2), becp, nbnd_occ(ikk))
                 ! host calbec; with OpenACC the s_psi below reads becp on the device
                 ! (present_or_copyin), so refresh the device copy if one exists
                 IF (gamma_only) THEN
                    !$acc update device(becp%r) if_present
                 ELSEIF (noncolin) THEN
                    !$acc update device(becp%nc) if_present
                 ELSE
                    !$acc update device(becp%k) if_present
                 ENDIF
                 CALL s_psi(npwx, npwq, nbnd_occ(ikk), coul_evc1(:,:,ik,2), coul_sevc1(:,:,ik))
              ELSE
                 coul_sevc1(:,:,ik) = coul_evc1(:,:,ik,2)
              ENDIF
              !
              ! Apply (H - e S), the frequency-independent part of the
              ! Liouvillian, to BOTH Krylov vectors.
              !
              CALL g2_kin(ikq)
              ! ch_psi_all requires its h/ah arguments to already be
              ! !$acc present: coul_evc1/coul_evc1_new are otherwise kept
              ! host-resident throughout this routine (lanczos_nonhermitian
              ! updates them via host BLAS), so register them on device
              ! only for the duration of this call, rather than keeping
              ! them permanently device-resident (which would make lr_dot's
              ! lenient present_or_copyin silently pick up stale device
              ! copies inside lanczos_nonhermitian).
              !$acc data copyin(coul_evc1(:,:,ik,1)) copyout(coul_evc1_new(:,:,ik,1))
              CALL ch_psi_all(npwq, coul_evc1(:,:,ik,1), coul_evc1_new(:,:,ik,1), &
                              et(:,ikk), ik, nbnd_occ(ikk))
              !$acc end data
              !$acc data copyin(coul_evc1(:,:,ik,2)) copyout(coul_evc1_new(:,:,ik,2))
              CALL ch_psi_all(npwq, coul_evc1(:,:,ik,2), coul_evc1_new(:,:,ik,2), &
                              et(:,ikk), ik, nbnd_occ(ikk))
              !$acc end data
              !
              ! Add the self-consistent screening term P_c^+ dV_Hxc q to
              ! the "ich" component only, every iteration.
              !
              dvpsi(:,:) = (0.0_DP, 0.0_DP)
              CALL apply_dpot_bands(ik, nbnd_occ(ikk), dfpt_data%dvscfs(:,:,1), evc, aux2)
              dvpsi(:,:) = dvpsi(:,:) + aux2(:,:)
              CALL adddvscf(1, ik, .FALSE.)
              CALL orthogonalize(dvpsi, coul_evq0(:,:,ik), ikk, ikq, coul_sevq0(:,:,ik), npwq, .TRUE.)
              coul_evc1_new(:,:,ik,ich) = coul_evc1_new(:,:,ik,ich) - dvpsi(:,:)
              !
              IF (okvan) THEN
                 dpsi(:,:) = (0.0_DP, 0.0_DP)
                 CALL lr_sm1_psi(ik, npwx, npwq, nbnd_occ(ikk), coul_evc1_new(:,:,ik,1), dpsi)
                 coul_evc1_new(:,:,ik,1) = dpsi(:,:)
                 dpsi(:,:) = (0.0_DP, 0.0_DP)
                 CALL lr_sm1_psi(ik, npwx, npwq, nbnd_occ(ikk), coul_evc1_new(:,:,ik,2), dpsi)
                 coul_evc1_new(:,:,ik,2) = dpsi(:,:)
              ENDIF
              !
           ENDDO ! ik
           !
        ENDIF ! iter > 1
        !
        ! Lanczos step of iteration (iter-1): normalizes q_(iter-1)/p_(iter-1),
        ! computes beta/gamma and the zeta projections onto the nmels_tot
        ! matrix-element target vectors, and turns evc1_new into (q_iter, p_iter).
        !
        IF (iter > 1) THEN
           !
           CALL lanczos_nonhermitian(iter-1, npwx*npol, nbnd, nksq,          &
                    coul_evc1, coul_evc1_new, coul_sevc1, coul_evc1_old,     &
                    nmels_tot, coul_d0psi_mels, alpha, beta, gamma_l, zeta)
           !
           coul_beta_store(ipert, iter-1)  = beta
           coul_gamma_store(ipert, iter-1) = gamma_l
           coul_zeta_store(ipert, :, iter-1) = zeta(:)
           !
           WRITE(stdout,'(5x,"pert",i5,": beta (",i8.8,")=",f15.8,5x,"gamma(",i8.8,")=",f15.8)') &
                ipert, iter-1, beta, iter-1, gamma_l
           !
        ENDIF
        !
        ! Last coefficient done. With crpa_coul_lanczos_save_final the response
        ! potential for the next step is still built (below), so that the
        ! chain can be extended later; otherwise stop here.
        !
        IF (iter == crpa_coul_lanczos_steps + 1 .AND. .NOT. crpa_coul_lanczos_save_final) EXIT
        !
        ! Compute the induced charge density from the just-updated "inch"
        ! Krylov vector, ready for the next iteration's screening term.
        !
        DO ik = 1, nksq
           ikk = ikks(ik)
           ikq = ikqs(ik)
           IF (lsda) current_spin = isk(ikk)
           IF (okvan) THEN
              npwq = ngk(ikq)
              CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb, .true.)
              !$acc update host(vkb)
           ENDIF
           evc(:,:) = coul_evc0(:,:,ik)
           !$acc update device(evc)
           dpsi(:,:) = coul_evc1(:,:,ik,inch)
           IF (noncolin) THEN
              CALL incdrhoscf_nc(dfpt_data%drhos(1,1,1), wk(ikk), ik, &
                                  dbecsum_nc(1,1,1,1,1), dpsi, 1.0_DP)
           ELSE
              CALL incdrhoscf(dfpt_data%drhos(1,current_spin,1), wk(ikk), ik, &
                               dfpt_data%dbecsum(1,1,current_spin,1), dpsi)
           ENDIF
        ENDDO
        !
        IF (noncolin) THEN
           CALL sternheimer_postprocess(1, dfpt_data, dbecsum_nc)
        ELSE
           CALL sternheimer_postprocess(1, dfpt_data)
        ENDIF
        !
        IF (okpaw) dfpt_data%dbecsum(:,:,:,:) = 2.0_DP * dfpt_data%dbecsum(:,:,:,:)
        !
        ! NOTE: no density symmetrization here (unlike
        ! crpa_do_lanczos_epsilon.f90): the calculation is assumed to run
        ! without exploiting crystal symmetry for these perturbations.
        !
        IF (lnoloc) THEN
           ! No local-field effects, as in dfpt_kernels.f90 and the epsilon chain.
           dfpt_data%dvscfp(:,:,1) = (0.0_DP, 0.0_DP)
        ELSE
           CALL zcopy(dfftp%nnr*nspin_mag, dfpt_data%drhop(1,1,1), 1, dfpt_data%dvscfp(1,1,1), 1)
           CALL dv_of_drho(dfpt_data%dvscfp(:,:,1))
        ENDIF
        !
        IF (doublegrid) THEN
           DO is = 1, nspin_mag
              CALL fft_interpolate(dfftp, dfpt_data%dvscfp(:,is,1), &
                                   dffts, dfpt_data%dvscfs(:,is,1))
           ENDDO
        ELSE
           dfpt_data%dvscfs(:,:,1) = dfpt_data%dvscfp(:,:,1)
        ENDIF
        !
        IF (okvan) THEN
           IF (okpaw) CALL PAW_dpotential(dfpt_data%dbecsum, rho%bec, int3_paw, 1)
           CALL newdq(dfpt_data%dvscfp, 1)
        ENDIF
        !
        IF (iter == crpa_coul_lanczos_steps + 1) THEN
           CALL crpa_coul_ckpt_save_vec(ipert, iter, dfpt_data%dvscfp(:,:,1))
           EXIT
        ENDIF
        IF (crpa_coul_lanczos_restart_step > 0) THEN
           IF (MOD(iter, crpa_coul_lanczos_restart_step) == 0) &
              CALL crpa_coul_ckpt_save_vec(ipert, iter, dfpt_data%dvscfp(:,:,1))
        ENDIF
        !
     ENDDO ! iter
     !
     ! chain complete: record it in the checkpoint of this q (always written)
     !
     CALL crpa_coul_ckpt_done(ipert)
     !
  ENDDO ! ipert
  !
  CALL crpa_coul_ckpt_end()
  !
  CALL crpa_deallocate_upert()
  CALL deallocate_bec_type_acc(becp)
  !
  !$acc exit data delete(aux2)
  DEALLOCATE(aux2)
  DEALLOCATE(zeta)
  !
  ! restore the alpha_pv form of the tetrahedron weights for whatever runs next
  !
  IF (ltetra) THEN
     dfpt_tetra_beta = tetra_beta0
     DEALLOCATE(tetra_beta0)
  ENDIF
  IF (noncolin) DEALLOCATE(dbecsum_nc)
  !
  !$acc exit data delete(dvbare)
  !$acc exit data delete(dvbares)
  !$acc exit data delete(dvscf_tmps)
  !$acc exit data delete(dvscf_tmp)
  !$acc exit data delete(dfpt_data, dfpt_data%dvscfs)
  !
  CALL deallocate_dfpt_data(dfpt_data)
  DEALLOCATE(dvbare)
  DEALLOCATE(dvbares)
  DEALLOCATE(dvscf_tmps)
  IF (okvan) THEN
     DEALLOCATE (int3)
     IF (okpaw)    DEALLOCATE (int3_paw)
     IF (noncolin) DEALLOCATE (int3_nc)
  ENDIF
  DEALLOCATE(dvscf_tmp)
  !
  ! Save the raw chains (before extrapolation), for crpa_spectrum_only and
  ! the q-level restart.
  !
  CALL crpa_write_coul_lanczos_chain(iq)
  CALL crpa_calc_coulomb_lanczos(iq)
  !
  CALL crpa_deallocate_coul_lanczos()
  !
  WRITE( stdout,*) "     "
  WRITE( stdout,*) "     =--------------------------------------------="
  WRITE( stdout,   '(13x,"CONVERGENCE HAS BEEN REACHED")')
  WRITE( stdout,*) "     =--------------------------------------------="
  !
  CALL stop_clock('crpa_coul_lncz')
  !
  RETURN
  !
END SUBROUTINE crpa_do_lanczos_coulomb
