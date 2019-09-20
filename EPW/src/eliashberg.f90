  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE eliashberg_eqs()
  !-----------------------------------------------------------------------
  !!
  !! This is the main driver for solving the Eliashberg equations 
  !!
  USE io_global,         ONLY : stdout 
  USE epwcom,            ONLY : liso, fila2f, gap_edge, lreal, limag, laniso 
  USE eliashbergcom,     ONLY : gap0
  USE superconductivity, ONLY : eliashberg_init, estimate_tc_gap, deallocate_eliashberg
  USE io_eliashberg,     ONLY : read_a2f, read_frequencies, read_eigenvalues, read_ephmat, &
                                read_kqmap  
  USE superconductivity_iso,   ONLY : eliashberg_iso_iaxis, eliashberg_iso_raxis
  USE superconductivity_aniso, ONLY : eliashberg_aniso_iaxis, evaluate_a2f_lambda
  !
  IMPLICIT NONE
  !
  CALL start_clock('ELIASHBERG')
  !
  IF (liso) THEN
    WRITE(stdout, '(/5x,a)') REPEAT('=',67)
    WRITE(stdout, '(5x,"Solve isotropic Eliashberg equations")')
    WRITE(stdout, '(5x,a/)') REPEAT('=',67)
    IF (fila2f /= ' ') THEN
      CALL read_a2f
      CALL eliashberg_init
    ELSE
      CALL read_frequencies()
      CALL read_eigenvalues()
      CALL read_kqmap()
      CALL read_ephmat()
      CALL eliashberg_init()
      CALL evaluate_a2f_lambda()
    ENDIF
    ! 
    CALL estimate_tc_gap()
    IF (gap_edge > 0.d0) THEN
      gap0 = gap_edge
    ENDIF
    IF (lreal) CALL eliashberg_iso_raxis()
    IF (limag) CALL eliashberg_iso_iaxis()
  ENDIF
  !
  IF (laniso) THEN
    WRITE(stdout, '(/5x,a)') REPEAT('=',67)
    WRITE(stdout, '(5x,"Solve anisotropic Eliashberg equations")')
    WRITE(stdout, '(5x,a/)') REPEAT('=',67)
    CALL read_frequencies()
    CALL read_eigenvalues()
    CALL read_kqmap()
    CALL read_ephmat()
    CALL eliashberg_init()
    CALL evaluate_a2f_lambda()
    CALL estimate_tc_gap()
    IF (gap_edge > 0.d0) THEN 
      gap0 = gap_edge
    ENDIF
    IF (limag) CALL eliashberg_aniso_iaxis()
  ENDIF
  !
  IF (.NOT. liso .AND. .NOT. laniso) THEN 
    WRITE(stdout, '(/5x,a)') REPEAT('=',67)
    WRITE(stdout, '(5x,"Calculate Eliashberg spectral function")')
    WRITE(stdout, '(5x,a/)') REPEAT('=',67)
    CALL read_frequencies()
    CALL read_eigenvalues()
    CALL read_kqmap()
    CALL read_ephmat()
    !
    CALL eliashberg_init()
    CALL evaluate_a2f_lambda()
    CALL estimate_tc_gap() 
  ENDIF
  !
  CALL deallocate_eliashberg()
  !
  CALL stop_clock('ELIASHBERG')
  !
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE eliashberg_eqs
  !-----------------------------------------------------------------------
