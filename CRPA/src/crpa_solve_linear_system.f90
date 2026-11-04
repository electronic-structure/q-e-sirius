!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_solve_linear_system (iq)
   !-----------------------------------------------------------------------
   !
   ! This is a driver routine for the solution of the linear-response Kohn-Sham
   ! equations. The solution defines the change of Kohn-Sham wavefunctions due
   ! to the perturbation.
   !
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE fft_base,             ONLY : dfftp, dffts
   USE buffers,              ONLY : save_buffer, get_buffer
   USE ions_base,            ONLY : nat
   USE wavefunctions,        ONLY : evc
   USE klist,                ONLY : ngk, lgauss, ltetra, igk_k, xk
   USE gvecs,                ONLY : doublegrid
   USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
   USE wvfct,                ONLY : nbnd, npwx
   USE becmod,               ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, becp
   USE uspp_param,           ONLY : nhm
   USE uspp,                 ONLY : nkb, okvan, vkb
   USE uspp_init,            ONLY : init_us_2
   USE uspp_param,           ONLY : nhm
   USE paw_variables,        ONLY : okpaw
   USE lrus,                 ONLY : int3, int3_paw, int3_nc
   USE fft_interfaces,       ONLY : fft_interpolate
   USE noncollin_module,     ONLY : npol, nspin_mag, noncolin, domag
   USE qpoint,               ONLY : nksq, ikks, ikqs, xq, ikmks
   USE control_lr,           ONLY : rec_code, where_rec, convt, lrpa, lgamma
   USE units_lr,             ONLY : iuwfc, lrwfc, lrdvwfc, iudvwfc
   USE dv_of_drho_lr,        ONLY : dv_of_drho
   USE apply_dpot_mod,       ONLY : apply_dpot_allocate, apply_dpot_deallocate, &
                                    apply_dpot_bands
   USE eqv,                  ONLY : dvpsi
   USE crpacom,              ONLY : code, dvbare, dvbares, v_coul_bare, v_coul_scrd, &
                                    crpa_xcmode, crpa_nfreq, crpa_freq, current_w
   USE crpa_pert,            ONLY : npert_tot
   USE dfpt_type,            ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data, &
                                    dfpt_dvscfp_to_dvscfs
   USE dfpt_kernels,         ONLY : dfpt_kernel
   USE constrained_dfpt,     ONLY : cdfpt, cdfpt_setup_q
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iq
   !! number of the q point
   !
   INTEGER :: ipert
   !! Index for perturbation
   INTEGER :: ifreq
   !! Index for frequency
   !
   REAL(DP) :: dr2
   !! self-consistency error
   !
   TYPE(dfpt_data_type) :: dfpt_data
   !! Data that describes linear response quantities
   !
   COMPLEX(DP), ALLOCATABLE :: dvtot(:,:)
   !! dvbare + dvscfp on the dense grid, USPP/PAW only (see newdq call below)
   !
   INTEGER :: ik, ikk,    & ! counter on k points
      ikq,        & ! index of k+q
      ndim,       &
      is,         & ! counter on spin polarizations
      npw,        & ! number of plane waves at k
      npwq,       & ! number of plane waves at k+q
      nsolv,      & ! number of linear systems
      isolv,      & ! counter on linear systems
      ikmk,       & ! index of mk
      nrec          ! the record number for dvpsi
   !
   CALL start_clock ('crpa_solve_linear_system')
   !
   WRITE( stdout,*) "     =--------------------------------------------="
   WRITE( stdout, '(13x,"    SOLVE THE LINEAR SYSTEM")')
   WRITE( stdout,*) "     =--------------------------------------------="
   !
   ! Allocate arrays for the SCF density/potential
   !
  !
   ! USPP/PAW now RUNS here (int3 allocation, dense->smooth interpolation and the
   ! augmentation terms were added 2026-09-02), but it is NOT VALIDATED: on the
   ! Si USPP test (crpa-test-Si/lanczos_runs/frozen_gauge_test/uspp_vs_nc) the CG
   ! and Lanczos solvers still disagree by 4-10%% on the screened U, and even the
   ! BARE elements differ by 4e-2 eV - while the norm-conserving case agrees to
   ! 5e-5 eV. Something in the USPP perturbation/density construction is still
   ! inconsistent between the two drivers. Treat USPP results as UNRELIABLE.
   !
   IF (okvan) THEN
      WRITE(stdout,'(/,5x,"*** WARNING: USPP/PAW support in the CRPA screened-Coulomb")')
      WRITE(stdout,'(5x,  "*** calculation is INCOMPLETE and UNVALIDATED. CG and Lanczos")')
      WRITE(stdout,'(5x,  "*** disagree by 4-10%% on Si (norm-conserving agrees to 5e-5 eV).")')
      WRITE(stdout,'(5x,  "*** Do not use these numbers for production. ***",/)')
   ENDIF
   !
   CALL allocate_dfpt_data(dfpt_data, 1)
   ALLOCATE (dvbare(dfftp%nnr, nspin_mag))
   ALLOCATE (dvbares(dffts%nnr, nspin_mag))
   !
   ! USPP/PAW: int3_nm = \int dV(r) Q_nm(r) dr, the augmentation term that
   ! adddvscf() adds to dV|psi>. NOT allocated anywhere else on this path
   ! (crpa_phescf_w.f90 allocates it only for the DIELECTRIC calculation), so
   ! adddvscf used to dereference an unallocated array and segfault.
   ! npert = 1 here: this path solves one perturbation at a time.
   !
   IF (okvan) THEN
      ALLOCATE (int3 (nhm, nhm, nat, nspin_mag, 1))
      IF (okpaw)    ALLOCATE (int3_paw (nhm, nhm, nat, nspin_mag, 1))
      IF (noncolin) ALLOCATE (int3_nc  (nhm, nhm, nat, nspin, 1))
   ENDIF
   !$acc enter data create(dfpt_data, dfpt_data%dvscfs(1:dffts%nnr, 1:nspin_mag, 1:1))
   !$acc enter data create(dvbare(1:dfftp%nnr, 1:nspin_mag))
   ! smooth-grid copy: apply_dpot_bands (here and in crpa_update_coulomb) only
   ! updates its dv argument on the device, so it must already be present
   !$acc enter data create(dvbares(1:dffts%nnr, 1:nspin_mag))
   !
   ! USPP-specific allocations
   !
   CALL allocate_bec_type_acc(nkb, nbnd, becp)
   !
   nsolv = 1
   IF (noncolin .AND. domag) nsolv = 2
   !
   CALL crpa_set_upert()
   !
   IF (cdfpt) CALL cdfpt_setup_q()
   !
   DO ipert = 1, npert_tot
      !
      WRITE(stdout,'(/6x,"q point #",i4,3x,"Perturbation #",i3)') iq, ipert
      !
      CALL crpa_compute_perturbation_q(ipert, dvbare, dfpt_data)
      !
      ! Perturbing leg of U_ij: crpa_xcmode selects v for 'TCTC', v+f_xc for 'TETE'.
      CALL crpa_dv_bare(dvbare)
      !
      ! dvbare lives on the DENSE grid; apply_dpot_bands wants the SMOOTH one.
      ! Identical arrays when doublegrid = .FALSE.
      !
      CALL crpa_interpolate_dense_to_smooth(dvbare, dvbares)
      !
      ! USPP/PAW: int3 must correspond to the potential actually being applied.
      ! Here (and for the bare matrix elements just below, called with
      ! dvscfs = 0) that potential is dvbare alone.
      !
      ! PAW: int3_paw is not yet built for this perturbation at this point -
      ! PAW_dpotential (which fills it) only runs later, inside dfpt_kernel,
      ! from the converged response density. Has to be zeroed or int3_paw is
      ! uninitialized or from the previous perturbation
      IF (okpaw) int3_paw = (0.0_DP, 0.0_DP)
      IF (okvan) CALL newdq(dvbare, 1)
      !
      ! Compute dV_bare * psi and write to buffer iubar
      !
      CALL apply_dpot_allocate()
      !$acc enter data create(dvpsi)
      !
      DO ik = 1, nksq
         !
         ikk  = ikks(ik)
         npw  = ngk(ikk)
         !
         IF (lsda) current_spin = isk(ikk)
         !
         DO isolv = 1, nsolv
            !
            IF (isolv == 1) THEN
               ikmk = ikks(ik)
            ELSE
               ikmk = ikmks(ik)
            ENDIF
            !
            ! Read unperturbed KS wavefuctions psi(k) from buffer
            !
            ! IF (nksq > 1 .OR. nsolv == 2)
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
            !$acc update device(evc)
            !
            ! Compute the action of the bare perturbation potential on the unperturbed KS
            ! wavefunctions: |dvpsi> = dV_pert * |evc>. Save the result to buffer iudvwfc.
            !
            nrec = ik + (isolv - 1) * nksq
            CALL apply_dpot_bands(ik, nbnd, dvbares, evc, dvpsi)
            !
            ! USPP/PAW: the response must be driven by the AUGMENTED
            ! perturbation V^ = dV(r) + sum_nm |beta_n> int3_nm <beta_m| (int3
            ! from the newdq(dvbare) above), not by the local part alone.
            ! adddvscf() does NOT build its own vkb: it zaxpy's int3-weighted
            ! vkb(:,k+q) into dvpsi, so vkb must first be set to the beta
            ! functions at k+q (npwq of them). Without this, adddvscf used
            ! whatever vkb crpa_init_q left behind (beta at the LAST k-point of
            ! this pool), which made the bare RHS - and hence the whole
            ! response - depend on the k-point/pool partition.
            !
            IF (okvan) THEN
               ikq  = ikqs(ik)
               npwq = ngk(ikq)
               CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb, .true.)
               !$acc update host(vkb)
               CALL adddvscf(1, ik, .FALSE.)
            ENDIF
            CALL save_buffer(dvpsi, lrdvwfc, iudvwfc, nrec)
            !
         ENDDO
         !
      ENDDO
      !
      !$acc exit data delete(dvpsi)
      !
      CALL apply_dpot_deallocate()
      !
      ! Compute bare Coulomb matrix elements
      dfpt_data%dvscfs = (0.0_DP, 0.0_DP)
      CALL crpa_update_coulomb(iq, ipert, v_coul_bare, dfpt_data%dvscfs)
      !
      ! Set records for restart
      !
      rec_code = 10
      where_rec = 'crpa_solve'
      !
      ! One self-consistent solve per frequency. Everything above (the bare
      ! perturbation, dV_bare*psi in iudvwfc, the bare matrix elements) does
      ! not depend on the frequency and is done once per perturbation.
      !
      DO ifreq = 1, crpa_nfreq
         !
         current_w = crpa_freq(ifreq)
         IF (crpa_nfreq > 1) WRITE(stdout, '(/6x,"Frequency #",i4," of",i4,": w = (",f10.6, &
                                          &",",f10.6,") Ry")') ifreq, crpa_nfreq, current_w
         !
         ! Start from a zero induced potential: dfpt_kernel takes dvscfp as the
         ! input of its first mixing step, otherwise the converged one of the
         ! previous solve
         !
         dfpt_data%dvscfp = (0.0_DP, 0.0_DP)
         IF (okpaw) dfpt_data%dbecsum = (0.0_DP, 0.0_DP)
         !
         ! alpha_pv stays on for every occupation scheme: with tetrahedra the static
         ! solve uses the alpha_pv form of dfpt_tetra_beta (as ph.x), and the
         ! complex-frequency projector builds its alpha_pv term at call time with the
         ! (deltae - w) denominator (orthogonalize_omega), as in the Gaussian case.
         !
         CALL dfpt_kernel(code, 1, 0, lrdvwfc, iudvwfc, dr2, dfpt_data, 0, 0, w_freq = current_w)
         !
         IF (.NOT. convt) THEN
            WRITE(stdout, '(/,5x, "DFPT iteration did not converge")')
            CALL crpa_stop_smoothly(.FALSE.)
         ENDIF
         !
         ! TCTC screened Coulomb: the induced external leg of U_ij carries the
         ! bare Coulomb potential of the converged response density, v*drho, not
         ! the full V_Hxc = (v+f_xc)*drho that dfpt_kernel leaves in dvscfp/dvscfs.
         ! The induced potential dvscfp is rebuilt from the converged drhop
         ! with the bare kernel.
         !
         IF (TRIM(crpa_xcmode) == 'TCTC' .AND. .NOT. lrpa) THEN
            dfpt_data%dvscfp(:, :, 1) = dfpt_data%drhop(:, :, 1)
            CALL crpa_dv_bare(dfpt_data%dvscfp(:, :, 1))
            IF ((lgauss .OR. ltetra) .AND. lgamma) THEN
               dfpt_data%dvscfp(:, 1, 1) = dfpt_data%dvscfp(:, 1, 1) - dfpt_data%def(1)
               IF (nspin_mag == 2) &
                  dfpt_data%dvscfp(:, 2, 1) = dfpt_data%dvscfp(:, 2, 1) - dfpt_data%def(1)
            ENDIF
            !
            CALL dfpt_dvscfp_to_dvscfs(dfpt_data)
            !$acc update device(dfpt_data%dvscfs(1:dffts%nnr, 1:nspin_mag, 1:1))
         ENDIF
         !
         ! Compute screened Coulomb matrix elements.
         !
         ! USPP/PAW: crpa_compute_mel_wannier applies (dvbare + dvscf) locally and
         ! then calls adddvscf ONCE, so int3 must be built from that same SUM.
         ! dfpt_kernel left int3 = \int dvscf Q (its own newdq), which is only the
         ! second half - recompute it on the total, on the DENSE grid.
         !
         IF (okvan) THEN
            ALLOCATE (dvtot(dfftp%nnr, nspin_mag))
            dvtot(:,:) = dvbare(:,:) + dfpt_data%dvscfp(:,:,1)
            CALL newdq(dvtot, 1)
            DEALLOCATE (dvtot)
         ENDIF
         !
         CALL crpa_update_coulomb(iq, ipert, v_coul_scrd(:,:,ifreq), dfpt_data%dvscfs)
         !
      ENDDO ! ifreq
      !
   ENDDO ! ipert
   !
   CALL crpa_deallocate_upert()
   !
   !$acc exit data delete(dvbare)
   !$acc exit data delete(dvbares)
   !$acc exit data delete(dfpt_data, dfpt_data%dvscfs)
   !
   CALL deallocate_dfpt_data(dfpt_data)
   CALL deallocate_bec_type_acc (becp)
   DEALLOCATE(dvbare)
   DEALLOCATE(dvbares)
   IF (okvan) THEN
      DEALLOCATE (int3)
      IF (okpaw)    DEALLOCATE (int3_paw)
      IF (noncolin) DEALLOCATE (int3_nc)
   ENDIF
   !
   WRITE( stdout,*) "     "
   WRITE( stdout,*) "     =--------------------------------------------="
   WRITE( stdout,   '(13x,"CONVERGENCE HAS BEEN REACHED")')
   WRITE( stdout,*) "     =--------------------------------------------="
   !
   CALL stop_clock ('crpa_solve_linear_system')
   !
END SUBROUTINE crpa_solve_linear_system
