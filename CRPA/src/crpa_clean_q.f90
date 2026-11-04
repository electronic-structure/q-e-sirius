!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_clean_q (flag)
  !-----------------------------------------------------------------------
  !
  ! This routine deallocates the variables of PWscf and of the
  ! CRPA code, and resets the same variables as after reading input in
  ! crpa_readin, so that it is possible to start a calculation at a new q.
  !
  USE lr_symm_base,    ONLY : nsymq
  !
  IMPLICIT NONE
  LOGICAL :: flag
  !
  CALL clean_pw(.FALSE.)
  !
  ! Deallocate the arrays
  !
  CALL crpa_dealloc_q()
  !
  nsymq = 0
  !
  ! Close the files
  !
  CALL crpa_close_q (flag)
  !
  RETURN
  !
END SUBROUTINE crpa_clean_q
