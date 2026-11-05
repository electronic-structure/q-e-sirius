!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE crpa_update_coulomb(iq, ipert, v_coul, dvscfs)
!------------------------------------------------------------
   !
   ! Update v_coul using the Coulomb matrix elements computed for the ipert-th perturbation
   ! at the iq-th q-point.
   !
   USE kinds,                ONLY : DP
   USE mp,                   ONLY : mp_sum
   USE mp_pools,             ONLY : inter_pool_comm, npool
   USE constants,            ONLY : tpi
   USE fft_base,             ONLY : dffts
   USE cell_base,            ONLY : at
   USE klist,                ONLY : xk
   USE qpoint,               ONLY : nksq, nksqtot, ikks, ikqs
   USE lsda_mod,             ONLY : nspin
   USE noncollin_module,     ONLY : nspin_mag
   USE wvfct,                ONLY : nbnd
   USE units_lr,             ONLY : lrdvwfc, iudvwfc
   USE crpacom,              ONLY : ik_to_ik_orig
   USE crpa_pert,            ONLY : pert_basis, npert_tot, pert_nbnd, pert_ibnd_min, &
                                    nmels_tot, mels_iwlist, mels_jwlist, mels_Rlist
   USE w90_interface,        ONLY : num_wann, u_matrix
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iq
   !! index of the q point
   INTEGER, INTENT(IN) :: ipert
   !! index of the perturbation
   COMPLEX(DP), INTENT(INOUT) :: v_coul(nmels_tot, npert_tot)
   !! Coulomb matrix elements for the ipert-th perturbation at the iq-th q-point
   COMPLEX(DP), INTENT(IN) :: dvscfs(dffts%nnr, nspin_mag, 1)
  !! Induced potenttial real space
   !
   INTEGER :: ib, jb, ib_, jb_
   !! Counters for bands
   INTEGER :: ik, ikk
   !! Counter for k points
   INTEGER :: ik_global
   !! Global index of k point
   INTEGER :: iw, jw
   !! Counter for Wannier functions
   INTEGER :: jpert
   !! Counter for perturbations
   REAL(DP) :: arg
   !! k dot r
   REAL(DP) :: xkk(3)
   !! k point in Cartesian coordinates
   REAL(DP) :: R(3)
   !! lattice vector in Cartesian coordinates
   COMPLEX(DP) :: phase
   !! exp(ikR) phase
   COMPLEX(DP) :: val
   !! Coulomn matrix in the Wannier basis
   COMPLEX(DP), ALLOCATABLE :: tmp(:, :)
   !! Temporary array
   COMPLEX(DP), ALLOCATABLE :: coulomb_matrix(:, :, :)
   !! Coulomb matrix elements in the band basis
   !
   CALL start_clock('crpa_upd_coul')
   !
   ! FIXME: Pool parallelization
   !
   IF (pert_basis == 'bands') THEN
      IF (npool > 1) CALL errore("crpa_update_coulomb", "npool > 1 not implemented for pert_basis = bands", 1)
      !
      ALLOCATE(coulomb_matrix(nbnd, nbnd, nksq))
      CALL lr_compute_matrix_elements(1, lrdvwfc, iudvwfc, dvscfs(1, 1, 1), coulomb_matrix)
      !
      DO ik = 1, nksq
         ik_global = ik  ! FIXME: Pool parallelization
         DO jb = 1, pert_nbnd
            DO ib = 1, pert_nbnd
               ib_ = ib + pert_ibnd_min - 1
               jb_ = jb + pert_ibnd_min - 1
               !
               jpert = (ik_global - 1) * pert_nbnd * pert_nbnd + (jb - 1) * pert_nbnd + ib
               !
               v_coul(jpert, ipert) = v_coul(jpert, ipert) + coulomb_matrix(ib_, jb_, ik)
            ENDDO
         ENDDO
      ENDDO
      !
      DEALLOCATE(coulomb_matrix)
      !
   ELSEIF (pert_basis == 'wannier') THEN
      !
      ! v_coul += U(ikq)^\dagger * coulomb_matrix * U(k) * exp(-ikR)
      !
      ALLOCATE(coulomb_matrix(num_wann, num_wann, nksq))
      CALL crpa_compute_mel_wannier(dvscfs(1, 1, 1), coulomb_matrix)
      !
      ALLOCATE(tmp(num_wann, nbnd))
      !
      DO jpert = 1, nmels_tot
         !
         R = REAL(mels_Rlist(:, jpert), DP)
         CALL cryst_to_cart(1, R, at, +1)
         !
         iw = mels_iwlist(jpert)
         jw = mels_jwlist(jpert)
         !
         val = (0.d0, 0.d0)
         !
         DO ik = 1, nksq
            !
            ikk = ikks(ik)
            xkk = xk(:, ikk)
            !
            arg = SUM(xkk * R) * tpi
            phase = CMPLX( COS(arg), - SIN(arg), KIND = DP)
            !
            val = val + coulomb_matrix(iw, jw, ik) * phase
            !
         ENDDO
         !
         CALL mp_sum(val, inter_pool_comm)
         val = val / nksqtot
         !
         v_coul(jpert, ipert) = v_coul(jpert, ipert) + val
         !
      ENDDO ! jpert
      !
      DEALLOCATE(tmp)
      DEALLOCATE(coulomb_matrix)
      !
   ENDIF
   !
   CALL stop_clock('crpa_upd_coul')
   !
END SUBROUTINE crpa_update_coulomb



SUBROUTINE crpa_compute_mel_wannier(dvscf, mel)
   !----------------------------------------------------------------------
   !
   ! Compute the matrix element for total potential perturbation.
   ! mel = < psi | dvpsi > (bare perturbation) + < psi | dvscf | psi > (induced perturbation)
   !
   ! TODO: Merge with lr_compute_matrix_elements ?
   ! (To do so, need to save WFs to buffer and pass buffer id to lr_compute_matrix_elements)
   !
   USE kinds,                ONLY : DP
   USE mp,                   ONLY : mp_sum
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE fft_base,             ONLY : dffts
   USE buffers,              ONLY : get_buffer
   USE lsda_mod,             ONLY : lsda, current_spin, isk
   USE wvfct,                ONLY : nbnd, npwx
   USE ldaU,                 ONLY : lda_plus_u
   USE noncollin_module,     ONLY : nspin_mag, npol
   USE eqv,                  ONLY : dvpsi, dpsi
   USE qpoint,               ONLY : ikks, ikqs, nksq
   USE apply_dpot_mod,       ONLY : apply_dpot_bands
   USE w90_interface,        ONLY : num_wann
   USE crpacom,              ONLY : dvbares, iuwan, lrwan
   USE uspp,                 ONLY : okvan, nkb, vkb
   USE uspp_param,           ONLY : upf, nh
   USE ions_base,            ONLY : nat, ityp, ntyp => nsp
   USE lrus,                 ONLY : int3
   USE klist,                ONLY : ngk, igk_k, xk
   USE uspp_init,            ONLY : init_us_2
   USE becmod,               ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN) :: dvscf(dffts%nnr, nspin_mag)
   !! Wavefunction pair density in real space
   COMPLEX(DP), INTENT(INOUT) :: mel(num_wann, num_wann, nksq)
   !! Matrix elements of the perturbation potential
   !
   LOGICAL :: time_reversed
   !! If solving time reversed Sternheimer equation
   INTEGER :: ik, ikk, ikq, nrec, npert
   INTEGER :: na, nt, ih, jh, ikb, jkb, ijkb0, ibnd, npw, npwq
   COMPLEX(DP) :: sum_us
   TYPE(bec_type) :: becp_wan
   !! <beta|w>, the WANNIER projections (not becp1, which is <beta|psi_KS>)
   !! indices
   COMPLEX(DP), ALLOCATABLE :: wan_k(:, :), wan_kq(:, :)
   !! Wannier functions at k and k+q, read from the iuwan buffer. These must NOT
   !! be the global evc/evq: those hold the KS states that the Sternheimer solver
   !! called right after this routine (crpa_solve_linear_system -> dfpt_kernel)
   !! relies on, and when nksq == 1 that solver does not re-read them
   !! (see the guard in LR_Modules/response_kernels.f90). Overwriting them here
   !! made the DFPT response contract against Wannier functions instead of KS
   !! states, silently corrupting the screened Coulomb matrix elements.
   !
   CALL start_clock('crpa_mel_wan')
   !
   time_reversed = .FALSE.
   !
   mel = (0.d0, 0.d0)
   !
   ALLOCATE(wan_k(npwx*npol, num_wann))
   IF (okvan) CALL allocate_bec_type(nkb, num_wann, becp_wan)
   ALLOCATE(wan_kq(npwx*npol, num_wann))
   !
   !$acc enter data create(dvpsi, dpsi)
   !
   DO ik = 1, nksq
      !
      ikk = ikks(ik)
      ikq = ikqs(ik)
      npw  = ngk(ikk)
      npwq = ngk(ikq)
      !
      ! Read Wannier function from buffer
      !
      CALL get_buffer(wan_k, lrwan, iuwan, ikk)
      CALL get_buffer(wan_kq, lrwan, iuwan, ikq)
      !
      ! Compute dvpsi = dV_tot * psi_k
      !
      ! (1) Bare term dV_bare * psi_k is stored in iudvwfc
      !
      CALL apply_dpot_bands(ik, num_wann, dvbares, wan_k, dvpsi)   ! smooth grid (see crpa_interpolate_dense_to_smooth)
      !
      ! (2) Induced term dvscf * psi_k
      !
      ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
      ! dvscf_q from previous iteration (mix_potential)
      !
      CALL apply_dpot_bands(ik, num_wann, dvscf, wan_k, dpsi)
      dvpsi = dvpsi + dpsi
      !
      ! (3) Additional contribution for USPP
      !
      !  In the case of US pseudopotentials there is an additional
      !  selfconsist term which comes from the dependence of D on
      !  V_{eff} on the bare change of the potential
      !
      ! USPP/PAW augmentation of this matrix element:
      !
      !    mel(j,i) += sum_nm int3_nm <w_{k+q,j}|beta_n^{k+q}> <beta_m^k|w_{k,i}>
      !
      ! implemented by adding  sum_n |beta_n^{k+q}> sum_m int3_nm <beta_m^k|w_i>
      ! to dvpsi, so that the ZGEMM below (<wan_kq|dvpsi>) picks it up. Body
      ! mirrors LR_Modules/adddvscf.f90, but adddvscf itself CANNOT be used
      ! here: it hardcodes becp1(ik) from lrus, i.e. <beta|psi_KS>, whereas both
      ! projections here must be onto the WANNIER functions (num_wann columns,
      ! not nbnd). Calling it mixed two different band-index spaces.
      ! int3 is built by the caller's newdq() from the SAME potential that
      ! apply_dpot_bands applied above (dvbare, or dvbare + dvscf).
      !
      IF (okvan) THEN
         ! KET projector at k, BRA at k+q - as in adddvscf, where
         ! becp1(ik) = <beta^k|psi_k> while the zaxpy uses vkb at k+q. Building
         ! the ket against the k+q vkb (and with npwq) cost ~1.8 % on the bare
         ! U, and vanished at q=Gamma where k+q = k.
         CALL init_us_2(npw, igk_k(1,ikk), xk(1,ikk), vkb, .true.)
         !$acc update host(vkb)   ! init_us_2 fills vkb on the device; calbec/zaxpy/commutator below use the host copy
         CALL calbec(npw, vkb, wan_k, becp_wan, num_wann)
         CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb, .true.)
         !$acc update host(vkb)   ! init_us_2 fills vkb on the device; calbec/zaxpy/commutator below use the host copy
         ijkb0 = 0
         DO nt = 1, ntyp
            IF (upf(nt)%tvanp) THEN
               DO na = 1, nat
                  IF (ityp(na) == nt) THEN
                     DO ibnd = 1, num_wann
                        DO ih = 1, nh(nt)
                           ikb = ijkb0 + ih
                           sum_us = (0.0_DP, 0.0_DP)
                           DO jh = 1, nh(nt)
                              jkb = ijkb0 + jh
                              sum_us = sum_us + int3(ih, jh, na, current_spin, 1) * &
                                                becp_wan%k(jkb, ibnd)
                           ENDDO
                           CALL zaxpy(npwq, sum_us, vkb(1,ikb), 1, dvpsi(1,ibnd), 1)
                        ENDDO
                     ENDDO
                     ijkb0 = ijkb0 + nh(nt)
                  ENDIF
               ENDDO
            ELSE
               DO na = 1, nat
                  IF (ityp(na) == nt) ijkb0 = ijkb0 + nh(nt)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      !
      ! IF (lda_plus_u) CALL adddvhubscf(1, ik)   ! see the comment above
      !
      ! Compute matrix element
      ! mel(j, i) = <psi_{k+q,j} | dV_tot * psi_{k,i}>
      !
      CALL ZGEMM('C', 'N', num_wann, num_wann, npwx*npol, (1.d0, 0.d0), &
                 wan_kq, npwx*npol, dvpsi, npwx*npol, (0.d0, 0.d0), &
                 mel(:, :, ik), num_wann)
      !
   ENDDO ! ik
   !
   DEALLOCATE(wan_k, wan_kq)
   !
   !$acc exit data delete(dvpsi, dpsi)
   !
   CALL mp_sum(mel, intra_bgrp_comm)
   !
   CALL stop_clock('crpa_mel_wan')
   !
END SUBROUTINE crpa_compute_mel_wannier

!----------------------------------------------------------------------
SUBROUTINE crpa_interpolate_dense_to_smooth(dv_dense, dv_smooth)
   !----------------------------------------------------------------------
   !
   ! Copy a perturbing potential from the DENSE (dfftp) grid to the SMOOTH
   ! (dffts) one, per spin component.
   !
   ! Needed because dv_of_drho() and newdq() work on dfftp while
   ! apply_dpot_bands() declares its dummy as dv(dffts%nnr, nspin_mag). The two
   ! grids coincide when doublegrid = .FALSE. (norm-conserving with the default
   ! ecutrho = 4*ecutwfc), which is why every NC test passed while handing
   ! apply_dpot_bands a dfftp-sized array; with USPP (ecutrho >> 4*ecutwfc,
   ! e.g. dense 45^3 vs smooth 32^3) that silently reinterpreted dense data on
   ! the smooth grid.
   !
   USE kinds,           ONLY : DP
   USE fft_base,        ONLY : dfftp, dffts
   USE fft_interfaces,  ONLY : fft_interpolate
   USE gvecs,           ONLY : doublegrid
   USE noncollin_module,ONLY : nspin_mag
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN)  :: dv_dense (dfftp%nnr, nspin_mag)
   COMPLEX(DP), INTENT(OUT) :: dv_smooth(dffts%nnr, nspin_mag)
   !
   INTEGER :: is
   !
   IF (doublegrid) THEN
      DO is = 1, nspin_mag
         CALL fft_interpolate(dfftp, dv_dense(:,is), dffts, dv_smooth(:,is))
      ENDDO
   ELSE
      dv_smooth(:,:) = dv_dense(:,:)
   ENDIF
   !
   RETURN
   !
END SUBROUTINE crpa_interpolate_dense_to_smooth
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_dv_bare(drho)
   !-----------------------------------------------------------------------
   !! dv_of_drho() for an external leg of the screened Coulomb interaction W.
   !! The kernel it applies is selected by crpa_xcmode:
   !!   'TCTC' (default): force the bare Coulomb kernel (lrpa = .true.)
   !!   'TETE'          : leave lrpa as given (must be .false.), so the leg
   !!                     carries v + f_xc
   !! The screening kernel is applied by dv_of_drho and always follows lrpa.
   !
   USE kinds,         ONLY : DP
   USE fft_base,      ONLY : dfftp
   USE noncollin_module, ONLY : nspin_mag
   USE control_lr,    ONLY : lrpa
   USE dv_of_drho_lr, ONLY : dv_of_drho
   USE crpacom,       ONLY : crpa_xcmode
   !
   IMPLICIT NONE
   COMPLEX(DP), INTENT(INOUT) :: drho(dfftp%nnr, nspin_mag)
   LOGICAL :: lrpa_save
   !
   lrpa_save = lrpa
   IF (TRIM(crpa_xcmode) == 'TCTC') lrpa = .TRUE.
   ! drho is filled on the host, but may also be on the device (dvbare):
   ! dv_of_drho then works on the device copy (present_or_copy), so refresh
   ! it before and bring the result back after
   !$acc update device(drho) if_present
   CALL dv_of_drho(drho)
   !$acc update self(drho) if_present
   lrpa      = lrpa_save
   !
   RETURN
END SUBROUTINE crpa_dv_bare
