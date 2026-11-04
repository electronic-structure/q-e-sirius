!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_do_lanczos_epsilon()
  !-----------------------------------------------------------------------
  !! Driver that replaces crpa_solve_e_w() when crpa_diel_llanczos = .TRUE.
  !! Builds ONE non-Hermitian Lanczos chain, spanning the 3 Cartesian
  !! polarizations at q=Gamma jointly, for the frequency-independent operator
  !! A = S^-1 H (two coupled Krylov vectors "q" and "p", since a genuinely
  !! complex frequency breaks the cheaper pseudo-Hermitian shortcut). The
  !! chain does not depend on frequency: the resulting beta/gamma/zeta
  !! coefficients are stored in crpa_diel_lanczos and later used by
  !! crpa_calc_chi_epsilon() to build epsilon_w_diel over the whole frequency
  !! grid from that single chain.
  !!
  !! The 3 polarization blocks are ONE vector of the direct-sum space, carrying
  !! a single normalization (crpa_nh_lanczos_step, at the end of this file):
  !! they are NOT 3 independent chains. The per-iteration symmetrization of the
  !! induced density rotates the 3 Cartesian components into each other, so the
  !! Liouvillian is not block diagonal in the polarization index and a per-block
  !! beta would not tridiagonalize any single operator (it breaks the
  !! reciprocity of the response tensor, chi_ij /= chi_ji, and puts epsilon a
  !! few % off the Sternheimer/CG result). This is the same joint chain that
  !! thermo_pw's do_lanczos.f90 builds by stacking the rpert perturbations.
  !!
  !! The blocks are therefore advanced in LOCK-STEP, one Lanczos iteration at a
  !! time for all 3 together (as in thermo_pw's do_lanczos.f90), rather than one
  !! full chain at a time. This is also required so that the
  !! induced density can be symmetrized every iteration (crpa_psym_diel.f90,
  !! a CRPA-local routine built from LR_Modules/lr_sym_mod.f90 building
  !! blocks - see there for why), which mixes the 3 Cartesian components
  !! under the crystal's point group and therefore needs all 3 to be
  !! simultaneously up to date. Without this, when the calculation runs on a
  !! symmetry-reduced (irreducible) k-point set, the induced density built
  !! from incdrhoscf is only a partial (unsymmetrized) representation of the
  !! true crystal-symmetric density: the missing star-related components are
  !! never reconstructed, so the gap compounds every iteration until the
  !! chain's bi-orthogonal metric becomes ill-conditioned and beta/gamma
  !! blow up. (With a non-symmetry-reduced, full k-point set there is
  !! nothing to reconstruct, which is why that case was not affected.)
  !!
  !! Restrictions of this first implementation:
  !! - no DFPT+U, no noncollinear magnetism (domag)
  !! Metals (lgauss/ltetra) are allowed: dvpsi_e + orthogonalize() already
  !! build the fractional-occupation-weighted d0psi/response generically (see
  !! LR_Modules/orthogonalize.f90), and the intraband/Drude contribution to
  !! the head of epsilon_w_diel is added separately in crpa_calc_chi_epsilon()
  !! from crpa_plasma_frequency(). No Fermi-energy shift is applied: with
  !! alpha_pv=0 the chain has Delta N(G=0)=0 by construction, as CG at w/=0.
  !! PAW is supported (non-magnetic and plain-noncollinear cases), following
  !! the same dbecsum "*2" / PAW_dpotential / newdq sequence used by the CG
  !! path (LR_Modules/dfpt_kernels.f90) and by thermo_pw's do_lanczos.f90.
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout
  USE ions_base,              ONLY : nat
  USE klist,                  ONLY : xk, wk, ngk, igk_k, ltetra
  USE dfpt_tetra_mod,         ONLY : dfpt_tetra_beta, dfpt_tetra_ttheta
  USE gvecs,                  ONLY : doublegrid
  USE fft_base,                ONLY : dffts, dfftp
  USE fft_interfaces,          ONLY : fft_interpolate
  USE lsda_mod,                ONLY : lsda, nspin, current_spin, isk
  USE control_flags,           ONLY : gamma_only
  USE noncollin_module,        ONLY : noncolin, npol, nspin_mag, domag
  USE wvfct,                   ONLY : nbnd, npwx, et
  USE wavefunctions,           ONLY : evc
  USE eqv,                     ONLY : dpsi, dvpsi
  USE becmod,                  ONLY : becp, calbec
  USE uspp,                    ONLY : okvan, vkb
  USE uspp_param,              ONLY : nhm
  USE paw_variables,           ONLY : okpaw
  USE paw_onecenter,           ONLY : PAW_dpotential
  USE lrus,                    ONLY : int3_paw
  USE scf,                     ONLY : rho
  USE ldaU,                    ONLY : lda_plus_u
  USE control_lr,              ONLY : alpha_pv, nbnd_occ, lgamma_gamma, lnoloc, &
                                      lflat_proj, flat_deg_thr
  USE units_lr,                ONLY : iuwfc, lrwfc
  USE buffers,                 ONLY : get_buffer
  USE qpoint,                  ONLY : nksq, ikks
  USE dv_of_drho_lr,           ONLY : dv_of_drho
  USE apply_dpot_mod,          ONLY : apply_dpot_bands
  USE response_kernels,        ONLY : sternheimer_postprocess
  USE dfpt_type,               ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data, &
                                       dfpt_dvscfp_to_dvscfs
  USE uspp_init,               ONLY : init_us_2
  USE crpacom,                 ONLY : crpa_diel_flat, crpa_diel_deg_thr
  USE crpacom,                 ONLY : crpa_diel_lanczos_steps, crpa_diel_lanczos_restart_step, &
                                      crpa_diel_lanczos_restart, crpa_diel_lanczos_save_final
  USE crpa_diel_lanczos,       ONLY : diel_evc0, diel_sevq0, diel_d0psi,        &
                                       diel_d0psi_left,                         &
                                       diel_evc1, diel_evc1_new, diel_evc1_old, &
                                       diel_sevc1, diel_beta_store,             &
                                       diel_gamma_store, diel_zeta_store
  USE incdrhoscf_mod,          ONLY : incdrhoscf, incdrhoscf_nc
  !
  IMPLICIT NONE
  !
  EXTERNAL :: crpa_nh_lanczos_step
  !
  TYPE(dfpt_data_type) :: dfpt_data
  COMPLEX(DP), ALLOCATABLE :: dbecsum_nc(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: aux2(:,:)
  !
  INTEGER :: iter, iter0, ik, ikk, npw, ipol, is, ich, inch
  INTEGER :: nnr
  !! dffts%nnr, used for acc
  INTEGER :: iter_progress
  !! restart bookkeeping: see crpa_restart_diel_lanczos.f90
  INTEGER :: ncoef_restart
  !! number of chain coefficients held by the checkpoint
  REAL(DP) :: alpha_pv0
  REAL(DP), ALLOCATABLE :: tetra_beta0(:,:,:)
  REAL(DP) :: beta, gamma_l
  COMPLEX(DP) :: zeta_m(3,3)
  !
  IF (lda_plus_u) CALL errore('crpa_do_lanczos_epsilon', &
     'DFPT+U is not yet supported with crpa_diel_llanczos = .true.', 1)
  IF (noncolin .AND. domag) CALL errore('crpa_do_lanczos_epsilon', &
     'noncolin with domag is not yet supported with crpa_diel_llanczos = .true.', 1)
  !
  CALL start_clock('crpa_diel_lncz')
  !
  CALL crpa_allocate_diel_lanczos()
  !
  CALL allocate_dfpt_data(dfpt_data, 3)
  IF (noncolin) ALLOCATE(dbecsum_nc(nhm, nhm, nat, nspin, 3))
  ALLOCATE(aux2(npwx*npol, nbnd))
  !
  nnr = dffts%nnr
  !$acc enter data create(dfpt_data, dfpt_data%dvscfs(1:nnr, 1:nspin_mag, 1:3))
  !$acc enter data create(aux2)
  !
  ! Fresh start; a checkpoint (crpa_diel_lanczos_restart) is restored below,
  ! after the starting vectors are built.
  !
  dfpt_data%dvscfs(:,:,:) = (0.0_DP, 0.0_DP)
  dfpt_data%dvscfp(:,:,:) = (0.0_DP, 0.0_DP)
  diel_evc1_old(:,:,:,:,:) = (0.0_DP, 0.0_DP)
  iter0 = 1
  !
  alpha_pv0 = alpha_pv
  !
  ! With tetrahedra alpha_pv cannot simply be set to 0: dfpt_tetra_main() bakes it
  ! into the weight array once, at setup (LR_Modules/dfpt_tetra_mod.f90:178,
  ! dfpt_tetra_beta = dfpt_tetra_ttheta + alpha_pv*dfpt_tetra_beta), whereas the
  ! lgauss branch of orthogonalize() builds its alpha_pv terms at call time. So
  ! every alpha_pv switch below is paired with a swap to the alpha_pv-free
  ! dfpt_tetra_ttheta, which is the exact tetrahedron counterpart of the Gaussian
  ! wwg = wg1*(1-theta) + wgp*theta. Without this, empty states enter the chain
  ! with weight -alpha_pv*beta ~ -1 instead of 0.
  !
  IF (ltetra) ALLOCATE(tetra_beta0, SOURCE=dfpt_tetra_beta)
  !
  ! Cache the unperturbed wavefunctions (and S*evc, for USPP) for the whole
  ! calculation: they do not change during the Lanczos chains.
  !
  DO ik = 1, nksq
     ikk = ikks(ik)
     npw = ngk(ikk)
     IF (lsda) current_spin = isk(ikk)
     IF (nksq > 1) THEN
        CALL get_buffer(evc, lrwfc, iuwfc, ikk)
        !$acc update device(evc)
     ENDIF
     diel_evc0(:,:,ik) = evc(:,:)
     CALL init_us_2(npw, igk_k(1,ikk), xk(1,ikk), vkb, .true.)
     !$acc update host(vkb)
     IF (okvan) THEN
        CALL calbec(npw, vkb, evc, becp, nbnd)
        ! host calbec; with OpenACC the s_psi below reads becp on the device
        ! (present_or_copyin), so refresh the device copy if one exists
        IF (gamma_only) THEN
           !$acc update device(becp%r) if_present
        ELSEIF (noncolin) THEN
           !$acc update device(becp%nc) if_present
        ELSE
           !$acc update device(becp%k) if_present
        ENDIF
        CALL s_psi(npwx, npw, nbnd, evc, diel_sevq0(1,1,ik))
     ENDIF
  ENDDO
  !
  ! Build the 3 starting vectors d0psi = P_c^+ V_ext|psi_k> (one per Cartesian
  ! polarization): lanczos_nonhermitian computes, at every iteration, the
  ! projection of the current Lanczos vector onto ALL 3 d0psi components
  ! (needed for the off-diagonal epsilon_w_diel(ipol,jpol), ipol /= jpol).
  !
  DO ipol = 1, 3
     DO ik = 1, nksq
        ikk = ikks(ik)
        npw = ngk(ikk)
        IF (lsda) current_spin = isk(ikk)
        evc(:,:) = diel_evc0(:,:,ik)
        !$acc update device(evc)
        CALL init_us_2(npw, igk_k(1,ikk), xk(1,ikk), vkb, .true.)
        !$acc update host(vkb)
        !
        ! crpa_diel_flat: occupation-free projector Pbar_n inside dvpsi_e only
        ! (see crpa_solve_e_w.f90); the seed orthogonalize() below keeps W.
        !
        lflat_proj   = crpa_diel_flat
        flat_deg_thr = crpa_diel_deg_thr
        CALL dvpsi_e(ik, ipol)
        lflat_proj   = .FALSE.
        !
        ! LEFT vector for the zeta projections: the RAW dvpsi_e() output, with
        ! NO further orthogonalize(). This is exactly what the CG path uses as
        ! its left vector (crpa_dielec.f90 reads the same quantity straight
        ! back from buffer iuebar and contracts it with dpsi, applying no
        ! projector of its own), so the two paths now extract epsilon with the
        ! same operator.
        !
        ! Why it matters (metals only): dvpsi_e() ALREADY carries one factor of
        ! the response weight. Its internal orthogonalize() (PHonon/PH/dvpsi_e.f90:91)
        ! multiplies [H,x]|psi_n> by W_nm, and the cgsolve_all that follows
        ! divides by the same (eps_m - eps_n) the commutator introduced, leaving
        ! dvpsi_e = sum_m |psi_m> W_nm <psi_m|x|psi_n>. Applying orthogonalize()
        ! a SECOND time here (as the previous code did) squares that weight.
        ! With W_nm = Theta_nm (f_n - f_m), an insulator has W in {0,1} so
        ! W^2 = W and the extra call was a harmless no-op - which is why CG and
        ! Lanczos agreed to ~1e-7 there. For a metal, Fermi-surface pairs have
        ! 0 < W < 1, W^2 < W, and the two paths diverge. Same root cause as the
        ! zeta left-vector bug fixed in crpa_do_lanczos_coulomb.f90.
        !
        diel_d0psi_left(:,:,ik,ipol) = dvpsi(:,:)
        !
        ! Chain SEED: keeps the orthogonalize(), because the seed must be the
        ! RHS of the Sternheimer equation, and that is what CG feeds its solver
        ! too (LR_Modules/response_kernels.f90:245 orthogonalizes the buffer
        ! contents + dV_scf|psi> before cgsolve_all). Do NOT "simplify" these
        ! two into one array again: they are deliberately different operators
        ! for metals, and identical only for insulators.
        !
        ! alpha_pv must be 0 here even though dvpsi_e above needed alpha_pv0:
        ! for metals orthogonalize() multiplies pairs v' <= nbnd_occ by
        ! (1 + alpha_pv/de), which is only cancelled by a subsequent solve with
        ! (H - e + alpha_pv P_v). The chain applies (H - e) with alpha_pv = 0,
        ! so a seed built with alpha_pv0 over-weights Fermi-surface pairs.
        alpha_pv = 0.0_DP
        IF (ltetra) dfpt_tetra_beta = dfpt_tetra_ttheta
        CALL orthogonalize(dvpsi, evc, ikk, ikk, diel_sevq0(:,:,ik), npw, .TRUE.)
        alpha_pv = alpha_pv0
        IF (ltetra) dfpt_tetra_beta = tetra_beta0
        diel_d0psi(:,:,ik,ipol) = -dvpsi(:,:)
     ENDDO
  ENDDO
  !
  ! alpha_pv=0 for the REST of this routine (chain-building loop below), not
  ! just around ch_psi_all: for a metal, LR_Modules/orthogonalize.f90's
  ! metallic branch multiplies alpha_pv into its occupation-difference weight
  ! wwg (the (wgp-wg1)/deltae and w0g terms), which is only exactly cancelled
  ! by ch_psi_all's alpha_pv*P_v shift when the SAME alpha_pv is used in both
  ! places at once (the standard metal-DFPT alpha_pv cancellation, see e.g.
  ! Rev. Mod. Phys. 73, 515 (2001)). The loop below sets alpha_pv=0 only
  ! around ch_psi_all and restores alpha_pv0 (nonzero) before the very next
  ! call, orthogonalize() for the self-consistent screening term: for a
  ! metal this leaves an UNCANCELLED alpha_pv contribution in wwg (of order
  ! alpha_pv = emax-emin, i.e. tens of Ry - LR_Modules/setup_alpha_pv.f90),
  ! corrupting every iteration's screening increment. Insulators never hit
  ! this: orthogonalize's lgauss/ltetra branch (the only one alpha_pv enters)
  ! is not taken there. This is the same fix already applied to the
  ! cdfpt-constrained Coulomb chain (LR_Modules/constrained_dfpt.f90) for the
  ! same reason - see crpa_cdfpt_lanczos_coulomb memory note.
  !
  alpha_pv = 0.0_DP
  IF (ltetra) dfpt_tetra_beta = dfpt_tetra_ttheta
  !
  ! Intraband (Drude) plasma-frequency tensor, metals only (crpa_metal=.false.
  ! and plasma_freq2=0 otherwise): a single non-iterative pass over
  ! bands/k-points, independent of the chain itself, using the same
  ! commutator_Hx_psi() machinery as dvpsi_e above. Computed here (rather
  ! than in crpa_calc_chi_epsilon, which also runs with no live wavefunctions
  ! in crpa_spectrum_only mode) and saved alongside the chain in
  ! crpa_write_diel_lanczos_chain(), see crpa_plasma_frequency.f90.
  !
  CALL crpa_plasma_frequency()
  !
  ! Checkpoint / restart (crpa_restart_diel_lanczos.f90, as TDDFPT's
  ! lr_restart). The starting vectors above are always rebuilt (diel_d0psi_left
  ! is needed for the zeta projections). On a restart the chain resumes at
  ! iter_progress+1 with the saved Lanczos vectors and response potential, or,
  ! if the checkpoint already holds crpa_diel_lanczos_steps coefficients, is
  ! truncated and the loop skipped.
  !
  IF (crpa_diel_lanczos_restart) THEN
     CALL crpa_read_restart_diel(iter_progress, ncoef_restart, dfpt_data%dvscfp)
     IF (ncoef_restart >= crpa_diel_lanczos_steps) THEN
        iter0 = crpa_diel_lanczos_steps + 2
        WRITE(stdout, '(/,5x,"Dielectric Lanczos checkpoint holds",i8," steps: truncated to", &
             &i8,", no new steps")') ncoef_restart, crpa_diel_lanczos_steps
     ELSE
        ! dense-grid potential from the checkpoint; smooth-grid copy and int3
        ! rebuilt exactly as at the end of every loop iteration
        CALL dfpt_dvscfp_to_dvscfs(dfpt_data)
        IF (okvan) CALL newdq(dfpt_data%dvscfp, 3)
        iter0 = iter_progress + 1
        WRITE(stdout, '(/,5x,"Restarting the dielectric Lanczos chain after step",i8, &
             &", continuing to",i8)') ncoef_restart, crpa_diel_lanczos_steps
     ENDIF
  ENDIF
  !
  WRITE(stdout, '(/,5x,"Non-Hermitian Lanczos chains for the dielectric", &
                 &" function (3 polarizations, lock-step)")')
  !
  DO iter = iter0, crpa_diel_lanczos_steps + 1
     !
     ! ich/inch select, this iteration, which of the two coupled Krylov
     ! vectors (1="q", 2="p") receives the self-consistent dV_Hxc term
     ! (ich) and which one is used right after to build the charge
     ! density for the NEXT iteration's dV_Hxc (inch); they swap every
     ! iteration, following thermo_pw's do_lanczos.f90 convention.
     !
     IF (MOD(iter,2) == 0) THEN
        ich = 2; inch = 1
     ELSE
        ich = 1; inch = 2
     ENDIF
     !
     dfpt_data%drhos(:,:,:) = (0.0_DP, 0.0_DP)
     dfpt_data%drhop(:,:,:) = (0.0_DP, 0.0_DP)
     dfpt_data%dbecsum(:,:,:,:) = (0.0_DP, 0.0_DP)
     IF (noncolin) dbecsum_nc(:,:,:,:,:) = (0.0_DP, 0.0_DP)
     !
     DO ipol = 1, 3
        DO ik = 1, nksq
           !
           ikk = ikks(ik)
           npw = ngk(ikk)
           IF (lsda) current_spin = isk(ikk)
           evc(:,:) = diel_evc0(:,:,ik)
           !$acc update device(evc)
           CALL init_us_2(npw, igk_k(1,ikk), xk(1,ikk), vkb, .true.)
           !$acc update host(vkb)
           !
           IF (iter == 1) THEN
              !
              ! Starting vector: d0psi = P_c^+ V_ext |psi_k>, already built
              ! for all 3 polarizations above. Both Krylov vectors start
              ! from the same vector (thermo_pw's do_lanczos.f90 convention).
              !
              dvpsi(:,:) = -diel_d0psi(:,:,ik,ipol)
              !
              IF (okvan) THEN
                 dpsi(:,:) = (0.0_DP, 0.0_DP)
                 CALL lr_sm1_psi(ik, npwx, npw, nbnd_occ(ikk), dvpsi, dpsi)
                 dvpsi(:,:) = dpsi(:,:)
              ENDIF
              !
              diel_evc1(:,:,ik,ipol,1) = -dvpsi(:,:)
              diel_evc1(:,:,ik,ipol,2) = -dvpsi(:,:)
              !
           ELSE
              !
              ! S * p_(iter-1), needed by lanczos_nonhermitian for the
              ! bi-orthogonal metric beta = <q|S|p>. Must be computed from
              ! the CURRENT (not yet updated) vector, before diel_evc1 is
              ! overwritten by the lanczos_nonhermitian call below.
              !
              IF (okvan) THEN
                 CALL calbec(npw, vkb, diel_evc1(:,:,ik,ipol,2), becp, nbnd_occ(ikk))
                 ! host calbec; with OpenACC the s_psi below reads becp on the device
                 ! (present_or_copyin), so refresh the device copy if one exists
                 IF (gamma_only) THEN
                    !$acc update device(becp%r) if_present
                 ELSEIF (noncolin) THEN
                    !$acc update device(becp%nc) if_present
                 ELSE
                    !$acc update device(becp%k) if_present
                 ENDIF
                 CALL s_psi(npwx, npw, nbnd_occ(ikk), diel_evc1(:,:,ik,ipol,2), diel_sevc1(:,:,ik,ipol))
              ELSE
                 diel_sevc1(:,:,ik,ipol) = diel_evc1(:,:,ik,ipol,2)
              ENDIF
              !
              ! Apply (H - e S), the frequency-independent part of the
              ! Liouvillian, to BOTH Krylov vectors.
              !
              CALL g2_kin(ikk)
              ! ch_psi_all requires its h/ah arguments to already be
              ! !$acc present: diel_evc1/diel_evc1_new are otherwise kept
              ! host-resident throughout this routine (lanczos_nonhermitian
              ! updates them via host BLAS), so register them on device
              ! only for the duration of this call, rather than keeping
              ! them permanently device-resident (which would make lr_dot's
              ! lenient present_or_copyin silently pick up stale device
              ! copies inside lanczos_nonhermitian).
              !$acc data copyin(diel_evc1(:,:,ik,ipol,1)) copyout(diel_evc1_new(:,:,ik,ipol,1))
              CALL ch_psi_all(npw, diel_evc1(:,:,ik,ipol,1), diel_evc1_new(:,:,ik,ipol,1), &
                              et(:,ikk), ik, nbnd_occ(ikk))
              !$acc end data
              !$acc data copyin(diel_evc1(:,:,ik,ipol,2)) copyout(diel_evc1_new(:,:,ik,ipol,2))
              CALL ch_psi_all(npw, diel_evc1(:,:,ik,ipol,2), diel_evc1_new(:,:,ik,ipol,2), &
                              et(:,ikk), ik, nbnd_occ(ikk))
              !$acc end data
              !
              ! Add the self-consistent screening term P_c^+ dV_Hxc q to the
              ! "ich" component only (every iteration, unlike the
              ! pseudo-Hermitian case where this happens every other step).
              !
              dvpsi(:,:) = (0.0_DP, 0.0_DP)
              CALL apply_dpot_bands(ik, nbnd_occ(ikk), dfpt_data%dvscfs(:,:,ipol), evc, aux2)
              dvpsi(:,:) = dvpsi(:,:) + aux2(:,:)
              CALL adddvscf(ipol, ik, .FALSE.)
              CALL orthogonalize(dvpsi, evc, ikk, ikk, diel_sevq0(:,:,ik), npw, .TRUE.)
              diel_evc1_new(:,:,ik,ipol,ich) = diel_evc1_new(:,:,ik,ipol,ich) - dvpsi(:,:)
              !
              IF (okvan) THEN
                 dpsi(:,:) = (0.0_DP, 0.0_DP)
                 CALL lr_sm1_psi(ik, npwx, npw, nbnd_occ(ikk), diel_evc1_new(:,:,ik,ipol,1), dpsi)
                 diel_evc1_new(:,:,ik,ipol,1) = dpsi(:,:)
                 dpsi(:,:) = (0.0_DP, 0.0_DP)
                 CALL lr_sm1_psi(ik, npwx, npw, nbnd_occ(ikk), diel_evc1_new(:,:,ik,ipol,2), dpsi)
                 diel_evc1_new(:,:,ik,ipol,2) = dpsi(:,:)
              ENDIF
              !
           ENDIF
           !
        ENDDO ! ik
     ENDDO ! ipol
     !
     ! Lanczos step of iteration (iter-1), for each of the 3 polarizations:
     ! normalizes q_(iter-1)/p_(iter-1), computes beta/gamma and the zeta
     ! projections onto the 3 starting vectors, and turns evc1_new into the
     ! new (q_iter, p_iter).
     !
     IF (iter > 1) THEN
        !
        ! ONE chain spanning all 3 polarizations, not 3 independent ones.
        ! The per-iteration symmetrization below (crpa_psym_diel) rotates the
        ! induced density with upert(jpert,ipert,isym) = s(ipert,jpert,isym),
        ! which genuinely MIXES the 3 Cartesian components, so the Liouvillian
        ! is not block diagonal in the polarization index. Three chains, each
        ! renormalized by its own beta, therefore do not tridiagonalize any
        ! single operator: the cross terms lose reciprocity (chi_ij /= chi_ji)
        ! and epsilon comes out a few % off the Sternheimer/CG result. Using a
        ! single beta/gamma for all 3 blocks makes this the joint chain over
        ! the direct-sum space, exactly as thermo_pw's do_lanczos.f90 does
        ! (its evc1 stacks the rpert perturbations into one vector of size
        ! npwx*npol*nbnd*nksq*rpert, normalized by a single beta_store(iter)).
        !
        ! The size argument here must match the actual allocated shape of
        ! diel_evc1 etc. (nbnd, from crpa_allocate_diel_lanczos), not
        ! nbnd_occ: crpa_nh_lanczos_step declares its qj/Aqj/... dummy
        ! arguments as explicit-shape arrays sized by this argument, so
        ! passing anything other than the true second dimension would
        ! misalign the per-k memory blocks whenever nksq > 1 (exactly the
        ! convention TDDFPT/src/lr_lanczos.f90 follows via size(evc1,2)).
        !
        CALL crpa_nh_lanczos_step(iter-1, npwx*npol, nbnd, nksq, 3,      &
                 diel_evc1, diel_evc1_new, diel_sevc1, diel_evc1_old,    &
                 diel_d0psi_left, beta, gamma_l, zeta_m)
        !
        ! beta/gamma are chain properties, not per-polarization ones: store
        ! the same value in all 3 rows so that the per-ipol tridiagonal solve
        ! of crpa_calc_chi_epsilon() reduces to the single joint solve (the
        ! sub/super-diagonals, and hence r, become identical for all ipol).
        !
        DO ipol = 1, 3
           diel_beta_store(ipol, iter-1)  = beta
           diel_gamma_store(ipol, iter-1) = gamma_l
           diel_zeta_store(ipol, :, iter-1) = zeta_m(ipol, :)
        ENDDO
        !
        WRITE(stdout,'(5x,"beta (",i8.8,")=",f15.8,5x,"gamma(",i8.8,")=",f15.8)') &
             iter-1, beta, iter-1, gamma_l
        !
     ENDIF
     !
     ! Last coefficient done. With crpa_diel_lanczos_save_final the response
     ! potential for the next step is still built (below), so that the final
     ! checkpoint can be extended later; otherwise stop here.
     !
     IF (iter == crpa_diel_lanczos_steps + 1 .AND. .NOT. crpa_diel_lanczos_save_final) EXIT
     !
     ! Compute the induced charge density from the just-updated "inch"
     ! Krylov vector, ready for the next iteration's screening term
     ! (built every iteration, unlike the pseudo-Hermitian case).
     !
     DO ipol = 1, 3
        DO ik = 1, nksq
           ikk = ikks(ik)
           IF (lsda) current_spin = isk(ikk)
           IF (okvan) THEN
              npw = ngk(ikk)
              CALL init_us_2(npw, igk_k(1,ikk), xk(1,ikk), vkb, .true.)
              !$acc update host(vkb)
           ENDIF
           !
           ! incdrhoscf builds drho from CONJG(evc)*dpsi and takes the
           ! unperturbed wavefunctions from the wavefunctions module (evc),
           ! not from an argument, so evc must be reset to this k point here
           ! (as thermo_pw's do_lanczos.f90 does before its own incdrhoscf
           ! call). Without this, evc still holds the last k point left over
           ! from the ch_psi_all loop above and every k but that one enters
           ! the induced density with the wrong |psi>: harmless only when a
           ! pool owns a single k point, wrong (and hence dependent on how
           ! the k points are distributed over pools) as soon as nksq > 1.
           !
           evc(:,:) = diel_evc0(:,:,ik)
           !$acc update device(evc)
           dpsi(:,:) = diel_evc1(:,:,ik,ipol,inch)
           IF (noncolin) THEN
              CALL incdrhoscf_nc(dfpt_data%drhos(1,1,ipol), wk(ikk), ik, &
                                  dbecsum_nc(1,1,1,1,ipol), dpsi, 1.0_DP)
           ELSE
              CALL incdrhoscf(dfpt_data%drhos(1,current_spin,ipol), wk(ikk), ik, &
                               dfpt_data%dbecsum(1,1,current_spin,ipol), dpsi)
           ENDIF
        ENDDO
     ENDDO
     !
     IF (noncolin) THEN
        CALL sternheimer_postprocess(1, dfpt_data, dbecsum_nc)
     ELSE
        CALL sternheimer_postprocess(1, dfpt_data)
     ENDIF
     !
     ! PAW: the c.c. term in the dbecsum formula gives a factor of 2 (as
     ! in dfpt_kernels.f90); applied after sternheimer_postprocess (which
     ! already used the un-doubled dbecsum for the USPP-like augmentation
     ! charge via lr_addusddens), before PAW_dpotential below.
     !
     IF (okpaw) dfpt_data%dbecsum(:,:,:,:) = 2.0_DP * dfpt_data%dbecsum(:,:,:,:)
     !
     ! Symmetrize the induced density every iteration, exactly as the CG
     ! path (LR_Modules/dfpt_kernels.f90) and thermo_pw's do_lanczos.f90
     ! both do: with a symmetry-reduced (irreducible) k-point set, the
     ! density built above from incdrhoscf is only a partial representation
     ! of the true crystal-symmetric density, and this reconstructs the
     ! missing star-related components. The rotation-and-accumulation
     ! genuinely mixes the 3 Cartesian polarizations under the point group
     ! (e.g. a 3-fold axis permutes x/y/z), which is why all 3 must be up
     ! to date here.
     !
     ! crpa_psym_diel (CRPA-local, crpa_psym_diel.f90) is used here instead
     ! of LR_Modules/psymdvscf.f90, which computes the same result but
     ! redundantly loops over the FULL FFT grid on every MPI rank rather
     ! than distributing the work across ranks (its own per-call cost does
     ! not decrease with more processes). That is negligible for the
     ! ~10-30 calls of the CG dielectric solve, but dominates the runtime
     ! when called thousands of times here.
     !
     IF (.NOT. lgamma_gamma) THEN
        ! drhos is used by sternheimer_postprocess above and never read again, so symmetrizing it is dead work
        !CALL crpa_psym_diel(dfpt_data%drhos, dffts)
        CALL crpa_psym_diel(dfpt_data%drhop, dfftp)
        IF (okpaw) CALL PAW_dsymmetrize(dfpt_data%dbecsum)
     ENDIF
     !
     ! dv_of_drho works in place (drho in -> dV out): the induced density
     ! built above by incdrhoscf/sternheimer_postprocess lives in drhop,
     ! not dvscfp, so it must be copied in first (see the dvscftmp copy
     ! in LR_Modules/dfpt_kernels.f90 before its own dv_of_drho call).
     !
     IF (lnoloc) THEN
        ! No local field effect: set dvscf to 0 (as LR_Modules/dfpt_kernels.f90
        ! does for the CG path, and do_lanczos.f90 for thermo_pw's chain).
        dfpt_data%dvscfp(:,:,:) = (0.0_DP, 0.0_DP)
     ELSE
        DO ipol = 1, 3
           CALL zcopy(dfftp%nnr*nspin_mag, dfpt_data%drhop(1,1,ipol), 1, dfpt_data%dvscfp(1,1,ipol), 1)
           CALL dv_of_drho(dfpt_data%dvscfp(1,1,ipol))
        ENDDO
     ENDIF
     !
     IF (doublegrid) THEN
        DO ipol = 1, 3
           DO is = 1, nspin_mag
              CALL fft_interpolate(dfftp, dfpt_data%dvscfp(:,is,ipol), &
                                   dffts, dfpt_data%dvscfs(:,is,ipol))
           ENDDO
        ENDDO
     ELSE
        dfpt_data%dvscfs(:,:,:) = dfpt_data%dvscfp(:,:,:)
     ENDIF
     !
     ! NOTE: no explicit "update device(dvscfs)" is needed here:
     ! apply_dpot_bands (LR_Modules/apply_dpot_mod.f90) already does
     ! "!$acc update device(dv(...))" on its own dv argument on entry.
     !
     IF (okvan) THEN
        IF (okpaw) CALL PAW_dpotential(dfpt_data%dbecsum, rho%bec, int3_paw, 3)
        CALL newdq(dfpt_data%dvscfp, 3)
     ENDIF
     !
     IF (iter == crpa_diel_lanczos_steps + 1) THEN
        CALL crpa_write_restart_diel(iter, dfpt_data%dvscfp)
        EXIT
     ENDIF
     IF (crpa_diel_lanczos_restart_step > 0) THEN
        IF (MOD(iter, crpa_diel_lanczos_restart_step) == 0) &
           CALL crpa_write_restart_diel(iter, dfpt_data%dvscfp)
     ENDIF
     !
  ENDDO ! iter
  !
  !$acc exit data delete(aux2)
  DEALLOCATE(aux2)
  IF (noncolin) DEALLOCATE(dbecsum_nc)
  !
  ! restore the alpha_pv form of the tetrahedron weights for whatever runs next
  !
  IF (ltetra) THEN
     dfpt_tetra_beta = tetra_beta0
     DEALLOCATE(tetra_beta0)
  ENDIF
  !
  !$acc exit data delete(dfpt_data, dfpt_data%dvscfs)
  CALL deallocate_dfpt_data(dfpt_data)
  !
  CALL stop_clock('crpa_diel_lncz')
  !
  RETURN
  !
END SUBROUTINE crpa_do_lanczos_epsilon
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_nh_lanczos_step(j, npwx_npol, nbnd_l, nksq_l, npert, &
                                qj, Aqj, Sqj, qjold, u, beta, gamma, zeta)
  !-----------------------------------------------------------------------
  !! One step of the bi-orthogonal (non-Hermitian) Lanczos algorithm for the
  !! q=Gamma electric-field chain, for all npert Cartesian polarizations at
  !! once. Same recursion as LR_Modules/lanczos_nonhermitian.f90, with the
  !! npert blocks treated as ONE vector of the direct-sum space: beta/gamma
  !! are formed from the inner product summed over all blocks and the same
  !! beta/gamma renormalizes every block.
  !!
  !! That is required because the per-iteration symmetrization of the induced
  !! density (crpa_psym_diel.f90) rotates the 3 Cartesian components into each
  !! other with upert(jpert,ipert,isym) = s(ipert,jpert,isym), so the
  !! Liouvillian L is NOT block diagonal in the polarization index. Writing
  !! D_n = diag(beta_1,beta_2,beta_3) for a per-block normalization, three
  !! separately renormalized chains tridiagonalize D_n^-1 L D_(n-1), which is
  !! similar to L only when D_n is a multiple of the identity. With a single
  !! beta, D_n = beta*I and the standard two-sided recursion applies to L
  !! itself, so one tridiagonal matrix serves all npert^2 response components
  !! and the reciprocity chi_ij = chi_ji is recovered.
  !!
  !! Array layout: qj/Aqj/qjold are (..., nksq, npert, 2), i.e. the npert
  !! blocks are CONTIGUOUS within each Krylov component, so a component spans
  !! nksq*npert blocks in memory. This is the same layout thermo_pw uses
  !! (evc1(npwx*npol,nbnd,nksq*rpert,2)) and the one a block-aware
  !! lanczos_nonhermitian would expect, so the BLAS below are single calls
  !! over the whole direct-sum vector rather than per-block loops.
  !!
  !! The two quantities that still have to be assembled here, rather than
  !! taken from LR_Modules/lanczos_nonhermitian.f90 directly, are exactly the
  !! two that LR_Modules/lr_dot.f90 cannot express: it is hardwired to the
  !! qpoint module's nksq (dummy shape and loop bound), so it can neither sum
  !! an inner product over the whole direct-sum vector (needed for the joint
  !! beta) nor return the block-resolved zeta(ipert,jpert). Everything else -
  !! the scaling, the three-term recurrence and the vector shuffle - is the
  !! shared routine's body verbatim, and lr_dot is still used for every
  !! individual inner product.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: j
  !! iteration index (zeta is only accumulated on even steps)
  INTEGER, INTENT(IN) :: npwx_npol
  !! first dimension of the Krylov vectors, npwx*npol
  INTEGER, INTENT(IN) :: nbnd_l
  !! second dimension, nbnd (the ALLOCATED band dimension, not nbnd_occ)
  INTEGER, INTENT(IN) :: nksq_l
  !! third dimension, nksq
  INTEGER, INTENT(IN) :: npert
  !! number of perturbations (3 Cartesian polarizations at q=Gamma)
  COMPLEX(DP), INTENT(INOUT) :: qj(npwx_npol, nbnd_l, nksq_l, npert, 2)
  !! (q_j, p_j); becomes (q_j+1, p_j+1) on output
  COMPLEX(DP), INTENT(INOUT) :: Aqj(npwx_npol, nbnd_l, nksq_l, npert, 2)
  !! (A q_j, A^T p_j) on input; used as workspace
  COMPLEX(DP), INTENT(IN) :: Sqj(npwx_npol, nbnd_l, nksq_l, npert)
  !! S p_j (a copy of p_j when okvan = .FALSE.)
  COMPLEX(DP), INTENT(INOUT) :: qjold(npwx_npol, nbnd_l, nksq_l, npert, 2)
  !! (q_j-1, p_j-1); becomes (q_j, p_j) on output
  COMPLEX(DP), INTENT(IN) :: u(npwx_npol, nbnd_l, nksq_l, npert)
  !! the npert starting vectors d0psi, for the zeta projections
  REAL(DP), INTENT(OUT) :: beta
  !! sub-diagonal coefficient of the tridiagonal matrix (whole chain)
  REAL(DP), INTENT(OUT) :: gamma
  !! super-diagonal coefficient of the tridiagonal matrix (whole chain)
  COMPLEX(DP), INTENT(OUT) :: zeta(npert, npert)
  !! zeta(ipert,jpert) = <u_jpert | q_j restricted to block ipert>
  !
  COMPLEX(DP), EXTERNAL :: lr_dot
  INTEGER :: size_all, ipert, jpert
  !
  ! one Krylov component = the npert blocks, contiguous
  !
  size_all = npwx_npol * nbnd_l * nksq_l * npert
  !
  ! Orthogonality requirement <p|S|q> = 1 over the whole direct-sum vector.
  ! lr_dot can only contract one nksq block at a time, so the blocks are
  ! summed here.
  !
  beta = 0.0_DP
  DO ipert = 1, npert
     beta = beta + DBLE(lr_dot(qj(1,1,1,ipert,1), Sqj(1,1,1,ipert)))
  ENDDO
  !
  IF (beta < 0.0_DP) THEN
     beta = SQRT(-beta)
     gamma = -beta
  ELSEIF (beta > 0.0_DP) THEN
     beta = SQRT(beta)
     gamma = beta
  ENDIF
  !
  ! Renormalize q(j), Aq(j) and p(j), Ap(j) - the whole direct-sum vector by
  ! the SAME factor, in one BLAS call per component.
  !
  CALL zscal(size_all, CMPLX(1.0_DP/beta,  0.0_DP, KIND=DP), qj(1,1,1,1,1),  1)
  CALL zscal(size_all, CMPLX(1.0_DP/beta,  0.0_DP, KIND=DP), Aqj(1,1,1,1,1), 1)
  CALL zscal(size_all, CMPLX(1.0_DP/gamma, 0.0_DP, KIND=DP), qj(1,1,1,1,2),  1)
  CALL zscal(size_all, CMPLX(1.0_DP/gamma, 0.0_DP, KIND=DP), Aqj(1,1,1,1,2), 1)
  !
  ! zeta coefficients, Eq. (35) of Malcioglu et al., Comput. Phys. Commun.
  ! 182, 1744 (2011): the block-resolved projections of the (single) current
  ! Lanczos vector onto the npert starting vectors. Nonzero only on even
  ! steps, as in lanczos_nonhermitian.f90 / nh_lanczos_step.f90.
  !
  IF (MOD(j,2) == 0) THEN
     DO ipert = 1, npert
        DO jpert = 1, npert
           zeta(ipert,jpert) = lr_dot(u(1,1,1,jpert), qj(1,1,1,ipert,1))
        ENDDO
     ENDDO
  ELSE
     zeta(:,:) = (0.0_DP, 0.0_DP)
  ENDIF
  !
  ! q(j+1) = A q(j) - gamma(j) q(j-1), p(j+1) = A^T p(j) - beta(j) p(j-1);
  ! then shift (q(j-1), p(j-1)) <- (q(j), p(j)) <- (q(j+1), p(j+1)).
  !
  CALL zaxpy(size_all, -CMPLX(gamma,0.0_DP,KIND=DP), qjold(1,1,1,1,1), 1, Aqj(1,1,1,1,1), 1)
  CALL zaxpy(size_all, -CMPLX(beta, 0.0_DP,KIND=DP), qjold(1,1,1,1,2), 1, Aqj(1,1,1,1,2), 1)
  !
  CALL zcopy(size_all, qj(1,1,1,1,1),  1, qjold(1,1,1,1,1), 1)
  CALL zcopy(size_all, Aqj(1,1,1,1,1), 1, qj(1,1,1,1,1),    1)
  CALL zcopy(size_all, qj(1,1,1,1,2),  1, qjold(1,1,1,1,2), 1)
  CALL zcopy(size_all, Aqj(1,1,1,1,2), 1, qj(1,1,1,1,2),    1)
  !
  RETURN
  !
END SUBROUTINE crpa_nh_lanczos_step
