!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_run()
  !-----------------------------------------------------------------------
  !
  ! Driver of the cRPA calculation (hp.x with method = 'cRPA').
  !
  ! Input has already been read (hp_readin -> crpa_readin) and the
  ! stop condition initialized in hp_main. MPI startup, the preamble,
  ! citations and the environment are also handled by hp_main; this
  ! routine performs the cRPA solver's own initialization, the q-point
  ! loop, and its clock printing.
  !
  USE io_global,         ONLY : stdout
  USE qpoint,            ONLY : nqs
  USE w90_wan_Rr_buffer, ONLY : w90_wan_Rr_close_buffer
  USE crpa_pert,         ONLY : pert_basis
  !
  IMPLICIT NONE
  !
  INTEGER :: iq
  LOGICAL :: do_iq, setup_pw
  !
  ! Allocate and initialize various arrays needed for the cRPA calculation
  !
  CALL crpa_init()
  !
  ! Print a summary
  !
  CALL crpa_summary()
  !
  ! Initialize the q mesh (q_points) and R mesh (R_points)
  !
  CALL crpa_generate_grids()
  !
  CALL crpa_setup_wannier()
  !
  ! Loop over the q points
  !
  DO iq = 1, nqs
     !
     ! Setup a calculation for a specific q point
     !
     CALL crpa_prepare_q(iq, do_iq, setup_pw)
     !
     !  If this q is not done in this run, cycle
     !
     IF (.NOT. do_iq) CYCLE
     !
     ! If necessary the bands are recalculated
     !
     IF (setup_pw) CALL crpa_run_nscf(.TRUE.)
     !
     ! Initialize the quantities which do not depend on
     ! the linear response of the system
     !
     CALL crpa_load_q()
     !
     CALL crpa_prepare_perturbations_q(iq)
     !
     ! SCF solution of the linearized Kohn-Sham equation
     !
     CALL crpa_solve_linear_system(iq)
     !
     ! Write the Coulomb matrix for the current q to file
     !
     CALL crpa_postprocess_coulomb(iq)
     !
     ! Cleanup of the variables for the next q point
     !
     CALL crpa_clean_q(.TRUE.)
     !
  ENDDO
  !
  IF (pert_basis == 'wannier') THEN
     ! Close the buffer
     CALL w90_wan_Rr_close_buffer()
  ENDIF
  !
  ! Print clocks
  !
  WRITE( stdout, * )
  WRITE( stdout, * )  '    PRINTING TIMING FROM PWSCF ROUTINES: '
  CALL print_clock_pw()
  CALL crpa_print_clock()
  !
  RETURN
  !
END SUBROUTINE crpa_run
