!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE lr_compute_matrix_elements(npert, lrdvpsi, iudvpsi, dvscf, mel)
  !----------------------------------------------------------------------
  !
  ! Compute the matrix element for total potential perturbation.
  ! mel = < psi | dvpsi > (bare perturbation) + < psi | dvscf | psi > (induced perturbation)
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
  USE units_lr,             ONLY : iuwfc, lrwfc
  USE eqv,                  ONLY : dvpsi, evq, dpsi
  USE qpoint,               ONLY : ikks, ikqs, nksq
  USE apply_dpot_mod,       ONLY : apply_dpot_bands
  USE control_lr,           ONLY : lgamma
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npert
  !! number of perturbations
  INTEGER, INTENT(IN) :: lrdvpsi
  !! record length for the buffer storing dV_bare * psi
  INTEGER, INTENT(IN) :: iudvpsi
  !! unit for the buffer storing dV_bare * psi
  COMPLEX(DP), INTENT(IN) :: dvscf(dffts%nnr, nspin_mag, npert)
  !! Wavefunction pair density in real space
  COMPLEX(DP), INTENT(INOUT) :: mel(nbnd, nbnd, nksq, npert)
  !! Matrix elements of the perturbation potential
  !
  LOGICAL :: time_reversed
  !! If solving time reversed Sternheimer equation
  INTEGER :: ik, ikk, ikq, ipert, nrec
  !! indices
  COMPLEX(DP), ALLOCATABLE :: aux(:, :)
  !! Auxiliary wavefunction
  !
  time_reversed = .FALSE.
  ! FIXME: time_reversed = .TRUE. for noncollinear magnetism. What is the proper matrix element?
  !
  CALL start_clock ('lr_matrix_elements')
  !
  ALLOCATE(aux(npwx*npol, nbnd))
  !
  mel = (0.d0, 0.d0)
  !
  DO ipert = 1, npert
    !
    DO ik = 1, nksq
      !
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      !
      IF (lsda) current_spin = isk(ikk)
      !
      CALL get_buffer(evc, lrwfc, iuwfc, ikks(ik))
      IF (lgamma) CALL get_buffer(evq, lrwfc, iuwfc, ikqs(ik))
      !
      ! Compute dvpsi = dV_tot * psi_k
      !
      ! (1) Bare term dV_bare * psi_k is stored in iudvpsi
      !
      nrec = (ipert - 1) * nksq + ik
      IF (time_reversed) nrec = nrec + npert * nksq
      !
      CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, nrec)
      !
      ! (2) Induced term dvscf * psi_k
      !
      ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
      ! dvscf_q from previous iteration (mix_potential)
      !
      CALL apply_dpot_bands(ik, nbnd, dvscf, evc, aux)
      dvpsi = dvpsi + aux
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
      CALL ZGEMM('C', 'N', nbnd, nbnd, npwx*npol, (1.d0, 0.d0), &
                evq, npwx*npol, dvpsi, npwx*npol, (0.d0, 0.d0), &
                mel(:, :, ik, ipert), nbnd)
      !
    ENDDO ! ik
    !
  ENDDO ! ipert
  !
  CALL mp_sum(mel, intra_bgrp_comm)
  !
  DEALLOCATE(aux)
  !
  CALL stop_clock ('lr_matrix_elements')
  !
END SUBROUTINE lr_compute_matrix_elements
