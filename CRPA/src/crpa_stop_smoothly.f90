!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE crpa_stop_smoothly(flag)
  !----------------------------------------------------------------------
  !
  ! Deallocate dynamical arrays, close files, and stop.
  !
  USE environment,    ONLY : environment_end
  USE mp_global,      ONLY : mp_global_end
  USE crpacom,        ONLY : code
  !
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: flag
  !
  CALL crpa_clean_q (.FALSE.)
  !
  ! Deallocate some arrays
  !
  CALL crpa_dealloc_1()
  CALL crpa_dealloc_2()
  !
  ! Print clocks
  !
  IF (flag) THEN
     CALL print_clock_pw()
     CALL crpa_print_clock()
  ENDIF
  !
  CALL environment_end(code)
  !
  CALL mp_global_end()
  !
  STOP 1
  !
END SUBROUTINE crpa_stop_smoothly
