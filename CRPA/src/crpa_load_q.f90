!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_load_q()
  !-----------------------------------------------------------------------
  !
  ! This is a driver to the CRPA initialization routines.
  !
  USE io_global,        ONLY : stdout
  USE qpoint,           ONLY : qpoint_setup_k_plus_q_indices
  USE crpacom,          ONLY : code
  !
  IMPLICIT NONE
  !
  ! Setup k and k+q point indices (also for -k and -k-q for magnetic calculations)
  !
  CALL qpoint_setup_k_plus_q_indices()
  !
  ! Allocate various arrays
  !
  CALL crpa_allocate_q()
  !
  ! Setup various control variables
  !
  CALL crpa_setup_q()
  !
  ! Output summary of the main variables
  !
  CALL crpa_summary_q()
  !
  ! Open all necessary files
  !
  CALL crpa_openfil_q()
  !
  ! Initialize all quantities which do not depend on the
  ! linear response to the perturbation
  !
  CALL crpa_init_q()
  !
  WRITE( stdout, '(/5x,"Total time spent up to now is:")')
  !
  CALL print_clock (code)
  !
  RETURN
  !
END SUBROUTINE crpa_load_q
