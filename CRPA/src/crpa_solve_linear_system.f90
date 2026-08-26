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
   USE klist,                ONLY : ngk
   USE gvecs,                ONLY : doublegrid
   USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
   USE wvfct,                ONLY : nbnd, npwx
   USE becmod,               ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, becp
   USE uspp_param,           ONLY : nhm
   USE uspp,                 ONLY : nkb
   USE noncollin_module,     ONLY : npol, nspin_mag, noncolin, domag
   USE qpoint,               ONLY : nksq, ikks, xq, ikmks
   USE control_lr,           ONLY : rec_code, where_rec, convt
   USE units_lr,             ONLY : iuwfc, lrwfc, lrdvwfc, iudvwfc
   USE dv_of_drho_lr,        ONLY : dv_of_drho
   USE apply_dpot_mod,       ONLY : apply_dpot_allocate, apply_dpot_deallocate, &
                                    apply_dpot_bands
   USE eqv,                  ONLY : dvpsi
   USE crpacom,              ONLY : code, dvbare, v_coul_bare, v_coul_scrd, w_freq
   USE crpa_pert,            ONLY : npert_tot
   USE dfpt_type,            ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data
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
   !
   REAL(DP) :: dr2
   !! self-consistency error
   !
   TYPE(dfpt_data_type) :: dfpt_data
   !! Data that describes linear response quantities
   !
   INTEGER :: ik, ikk,    & ! counter on k points
      ndim,       &
      is,         & ! counter on spin polarizations
      npw,        & ! number of plane waves at k
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
   CALL allocate_dfpt_data(dfpt_data, 1)
   ALLOCATE (dvbare(dfftp%nnr, nspin_mag))
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
      CALL dv_of_drho(dvbare)
      !
      ! Compute dV_bare * psi and write to buffer iubar
      !
      CALL apply_dpot_allocate()
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
            !
            ! Compute the action of the bare perturbation potential on the unperturbed KS
            ! wavefunctions: |dvpsi> = dV_pert * |evc>. Save the result to buffer iudvwfc.
            !
            nrec = ik + (isolv - 1) * nksq
            CALL apply_dpot_bands(ik, nbnd, dvbare, evc, dvpsi)
            CALL save_buffer(dvpsi, lrdvwfc, iudvwfc, nrec)
            !
         ENDDO
         !
      ENDDO
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
      CALL dfpt_kernel(code, 1, 0, lrdvwfc, iudvwfc, dr2, dfpt_data, 0, 0, w_freq = w_freq)
      !
      IF (.NOT. convt) THEN
         WRITE(stdout, '(/,5x, "DFPT iteration did not converge")')
         CALL crpa_stop_smoothly(.FALSE.)
      ENDIF
      !
      ! Compute screened Coulomb matrix elements
      !
      CALL crpa_update_coulomb(iq, ipert, v_coul_scrd, dfpt_data%dvscfs)
      !
   ENDDO ! ipert
   !
   CALL crpa_deallocate_upert()
   CALL deallocate_dfpt_data(dfpt_data)
   CALL deallocate_bec_type_acc (becp)
   DEALLOCATE(dvbare)
   !
   WRITE( stdout,*) "     "
   WRITE( stdout,*) "     =--------------------------------------------="
   WRITE( stdout,   '(13x,"CONVERGENCE HAS BEEN REACHED")')
   WRITE( stdout,*) "     =--------------------------------------------="
   !
   CALL stop_clock ('crpa_solve_linear_system')
   !
END SUBROUTINE crpa_solve_linear_system
