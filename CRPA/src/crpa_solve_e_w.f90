!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine crpa_solve_e_w
  !-----------------------------------------------------------------------
  !! This routine is a driver for the solution of the linear system which
  !! defines the change of the wavefunction due to an electric field.
  !! It performs the following tasks:
  !! a) computes the bare potential term times \(|\psi\rangle \);
  !! b) adds to it the screening term \(\Delta V_\text{SCF}|psi\rangle\).
  !!    If \(\text{lda_plus_u}=\text{TRUE}\) compute also the SCF part
  !!    of the response Hubbard potential;
  !! c) applies \(P_c^+\) (orthogonalization to valence states);
  !! d) calls \(\texttt{cgsolve_all}\) to solve the linear system;
  !! e) computes \(\Delta \rho\), \(\Delta V_\text{SCF}|psi\rangle\) and
  !!    symmetrizes them;
  !! f) if \(\text{lda_plus_u}=\text{TRUE}\) compute also the response
  !!    occupation matrices dnsscf.
  !! Step b, c, d are done in \(\text{sternheimer_kernel}\).
  !
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : ionode
  USE io_files,              ONLY : diropn
  USE klist,                 ONLY : xk, ngk, igk_k
  USE fft_base,              ONLY : dffts
  USE lsda_mod,              ONLY : lsda, current_spin, isk
  USE buffers,               ONLY : get_buffer
  USE wavefunctions,         ONLY : evc
  USE uspp,                  ONLY : vkb
  USE noncollin_module,      ONLY : nspin_mag
  USE ldaU,                  ONLY : lda_plus_u
  USE units_ph,              ONLY : lrdrho, iudrho, lrebar, iuebar
  USE units_lr,              ONLY : iuwfc, lrwfc
  USE output,                ONLY : fildrho
  USE qpoint,                ONLY : nksq, ikks
  USE control_lr,            ONLY : lgamma, convt, rec_code, where_rec, &
                                    lflat_proj, flat_deg_thr
  USE uspp_init,             ONLY : init_us_2
  USE dfpt_type,             ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data
  USE dfpt_kernels,          ONLY : dfpt_kernel

  USE crpacom,               ONLY : current_w, crpa_diel_flat, crpa_diel_deg_thr
  !
  IMPLICIT NONE
  !
  EXTERNAL :: stop_smoothly_ph
  !
  LOGICAL :: exst
  !!
  INTEGER :: ikk, npw, iter0, ipol, ik
  !! counters
  INTEGER :: nnr
  !! dffts%nnr, defined used for acc
  REAL(DP) :: dr2
  !! self-consistency error
  TYPE(dfpt_data_type) :: dfpt_data
  !! Data that describes linear response quantities
  !
  call start_clock ('crpa_solve_e_w')
  !
  !  This routine is task group aware
  !
  CALL allocate_dfpt_data(dfpt_data, 3)
  !
  nnr = dffts%nnr
  !$acc enter data create(dfpt_data, dfpt_data%dvscfs(1:nnr, 1:nspin_mag, 1:3))
  !
  convt = .false.
  iter0 = 0
  dr2 = 0.d0
  !
  ! q=0 for a metal is allowed: dfpt_kernel below detects lgauss/ltetra and
  ! applies the Fermi-energy-shift correction (LR_Modules/efermi_shift.f90)
  ! internally. Only the q/=0 case (no CRPA use for solve_e) is disallowed.
  !
  if ( .not. lgamma ) call errore ('solve_e', &
       'called in the wrong case', 1)
  !
  ! Compute P_c^+ x psi for all polarization and k points and store in buffer.
  ! With crpa_diel_flat the projector inside dvpsi_e is the occupation-free
  ! Pbar_n (LR_Modules/orthogonalize.f90, lflat_proj; Gaussian and tetrahedra): the buffer holds the
  ! unweighted auxiliary function, and the single metallic weight is applied
  ! later by dfpt_kernel. Only around dvpsi_e: every other orthogonalize call
  ! must keep the metallic weights.
  !
  lflat_proj   = crpa_diel_flat
  flat_deg_thr = crpa_diel_deg_thr
  !
  DO ik = 1, nksq
     DO ipol = 1, 3
        ikk = ikks(ik)
        npw = ngk(ikk)
        IF (lsda) current_spin = isk(ikk)
        !
        ! reads unperturbed wavefunctions psi_k in G_space, for all bands
        !
        IF (nksq > 1) THEN
           CALL get_buffer(evc, lrwfc, iuwfc, ikk)
           !$acc update device(evc)
        ENDIF
        !
        CALL init_us_2(npw, igk_k(1, ikk), xk(1, ikk), vkb, .true.)
        !$acc update host(vkb)
        !
        ! computes P_c^+ x psi_kpoint, written to buffer iuebar.
        !
        CALL dvpsi_e(ik, ipol)
        !
     ENDDO ! ipol
  ENDDO ! ik
  !
  lflat_proj = .FALSE.
  !
  ! Set records for restart
  !
  rec_code = -20  ! Electric field
  where_rec = 'solve_e...'
  !
  !   Solve DFPT fixed-point equation
  !
  ! alpha_pv stays on for every occupation scheme. With tetrahedra the static
  ! solve (w = 0) uses the alpha_pv form of dfpt_tetra_beta, as ph.x does, and the
  ! complex-frequency projector (orthogonalize_omega) builds its alpha_pv term at
  ! call time with the (deltae - w) denominator, so it cancels the alpha_pv*P_v
  ! shift of ch_psi_all_complex exactly, as in the Gaussian case.
  !
  CALL dfpt_kernel('CRPA', 3, iter0, lrebar, iuebar, dr2, dfpt_data, 1, 0, &
                   stop_callback = stop_smoothly_ph, w_freq = current_w)
  !
  IF (lda_plus_u) CALL dnsq_store(3, 0)
  !
  IF ( fildrho /= ' ') THEN
     IF ( ionode ) THEN
        INQUIRE (UNIT = iudrho, OPENED = exst)
        IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
        CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
     END IF
     DO ipol=1,3
        CALL davcio_drho(dfpt_data%drhop(1,1,ipol),lrdrho, iudrho,ipol,+1)
     END DO
  END IF
  !
  !$acc exit data delete(dfpt_data, dfpt_data%dvscfs)
  CALL deallocate_dfpt_data(dfpt_data)
  !
  call stop_clock ('crpa_solve_e_w')
  return
end subroutine crpa_solve_e_w
