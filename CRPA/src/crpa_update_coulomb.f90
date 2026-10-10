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
   USE wavefunctions,        ONLY : evc
   USE wvfct,                ONLY : nbnd, npwx
   USE ldaU,                 ONLY : lda_plus_u
   USE noncollin_module,     ONLY : nspin_mag, npol
   USE eqv,                  ONLY : dvpsi, evq, dpsi
   USE qpoint,               ONLY : ikks, ikqs, nksq
   USE apply_dpot_mod,       ONLY : apply_dpot_bands
   USE w90_interface,        ONLY : num_wann
   USE crpacom,              ONLY : dvbare, iuwan, lrwan
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
   INTEGER :: ik, ikk, ikq, ipert, nrec, npert
   !! indices
   !
   time_reversed = .FALSE.
   !
   mel = (0.d0, 0.d0)
   !
   DO ik = 1, nksq
      !
      ikk = ikks(ik)
      ikq = ikqs(ik)
      !
      ! Read Wannier function from buffer
      !
      CALL get_buffer(evc, lrwan, iuwan, ikk)
      CALL get_buffer(evq, lrwan, iuwan, ikq)
      !
      ! Compute dvpsi = dV_tot * psi_k
      !
      ! (1) Bare term dV_bare * psi_k is stored in iudvwfc
      !
      CALL apply_dpot_bands(ik, num_wann, dvbare, evc, dvpsi)
      !
      ! (2) Induced term dvscf * psi_k
      !
      ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
      ! dvscf_q from previous iteration (mix_potential)
      !
      CALL apply_dpot_bands(ik, num_wann, dvscf, evc, dpsi)
      dvpsi = dvpsi + dpsi
      !
      ! (3) Additional contribution for USPP
      !
      !  In the case of US pseudopotentials there is an additional
      !  selfconsist term which comes from the dependence of D on
      !  V_{eff} on the bare change of the potential
      !
      IF (time_reversed) THEN
         CALL adddvscf_ph_mag(ipert, ik)
      ELSE
         CALL adddvscf(ipert, ik)
      ENDIF
      !
      ! (4) Additional contribution for DFPT + U
      !
      ! DFPT+U: add to dvpsi the scf part of the response
      ! Hubbard potential dV_hub
      !
      IF (lda_plus_u) CALL adddvhubscf(ipert, ik)
      !
      ! Compute matrix element
      ! mel(j, i) = <psi_{k+q,j} | dV_tot * psi_{k,i}>
      !
      CALL ZGEMM('C', 'N', num_wann, num_wann, npwx*npol, (1.d0, 0.d0), &
                 evq, npwx*npol, dvpsi, npwx*npol, (0.d0, 0.d0), &
                 mel(:, :, ik), num_wann)
      !
   ENDDO ! ik
   !
   CALL mp_sum(mel, intra_bgrp_comm)
   !
END SUBROUTINE crpa_compute_mel_wannier
