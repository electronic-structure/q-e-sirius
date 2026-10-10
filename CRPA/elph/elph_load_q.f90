!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE elph_load_q()
  !-----------------------------------------------------------------------
  !! This is a driver to the phonon initialization routines.
  !! Same as PHonon/PH/initialize_ph, but
  !! - call elph_setup_wfc_from_w90 between openfilq and phq_init
  !
  USE qpoint, ONLY : qpoint_setup_k_plus_q_indices
  !
  IMPLICIT NONE
  !
  ! Setup k and k+q point indices (also for -k and -k-q for magnetic calculations)
  !
  CALL qpoint_setup_k_plus_q_indices()
  !
  !  Allocate the phonon variables
  !
  CALL allocate_phq()
  !
  !  Set the main control variable of the phonon code
  !
  CALL phq_setup()
  !
  !  Recover the status if available
  !
  CALL phq_recover()
  !
  !  Output summary of the main variables of the phonon code
  !
  CALL phq_summary()
  !
  !  Open the files of the ELPH code
  !
  CALL openfilq()
  !
  ! Compute the Wannier functions at k and k+q points
  !
  CALL elph_setup_wfc_from_w90()
  !
  !  Initialize all quantities which do not depend on the
  !  linear response to the perturbation
  !
  CALL phq_init()
  !
  CALL print_clock( 'ELPH' )
  !
  RETURN
  !
END SUBROUTINE elph_load_q


SUBROUTINE elph_setup_wfc_from_w90()
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : tmp_dir
  USE control_flags,        ONLY : io_level
  USE buffers,              ONLY : save_buffer, open_buffer, close_buffer
  USE klist,                ONLY : nks
  USE wvfct,                ONLY : nbnd, npwx
  USE noncollin_module,     ONLY : npol
  USE units_lr,             ONLY : lrwfc, iuwfc
  USE w90_interpolate,      ONLY : w90_interpolate_wfc
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !! k point index
  LOGICAL :: exst
  !! logical variable to check file exists
  LOGICAL :: exst_mem
  !! logical variable to check file exists in memory
  COMPLEX(DP), ALLOCATABLE :: wan(:, :)
  !! Wannier functions
  !
  WRITE( stdout, '(/5x, A, I0, A)') "Computing Wannier functions at k and k+q, ", nks, " points"
  !
  ! We have num_wann Wannier functions, which are less than nbnd.
  ! But, we want to store them in the buffer iuwfc which is of size nbnd.
  ! The reason is because (i) we want the elph code to be generic for all types of
  ! wavefunctions (band eigenstates or Wannier functions), and (ii) phq_init assumes
  ! wavefunctions are stored in the buffer iuwfc and has size nbnd.
  ! Therefore, we allocate wan to have size nbnd, but we will only use the first
  ! num_wann columns.
  !
  ALLOCATE(wan(npwx*npol, nbnd))
  wan(:, :) = (0.0_DP, 0.0_DP)
  !
  ! Close and reopen wavefunction buffer for robustness.
  ! (See PH/openfilq.f90: If io_level = 0, iuwfc = 10 is used. This may cause differences
  ! in lrwfc due to change in the number of plane waves.)
  !
  ! IF (iuwfc == 10) THEN
    CALL close_buffer(iuwfc, 'keep')
    ! iuwfc = 20
    CALL open_buffer(iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
  ! ENDIF
  !
  DO ik = 1, nks
    CALL w90_interpolate_wfc(ik, wan)
    CALL save_buffer(wan, lrwfc, iuwfc, ik)
  ENDDO ! ik
  !
END SUBROUTINE elph_setup_wfc_from_w90
