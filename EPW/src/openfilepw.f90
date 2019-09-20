  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/openfilq - Quantum-ESPRESSO group                
  !-----------------------------------------------------------------------
  SUBROUTINE openfilepw()
  !-----------------------------------------------------------------------
  !!
  !! This SUBROUTINE opens all the files necessary for the EPW
  !! calculation.
  !!
  !! RM - Nov/Dec 2014
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !-----------------------------------------------------------------------
  USE io_files,         ONLY : prefix, diropn, seqopn
  USE units_lr,         ONLY : iuwfc, lrwfc
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : npol, nspin_mag
  USE units_ph,         ONLY : lrdrho
  USE fft_base,         ONLY : dfftp
  USE uspp,             ONLY : okvan
  !
  IMPLICIT NONE
  !
  ! Local variables
  LOGICAL :: exst
  !! logical variable to check file existe
  !
  IF (len_TRIM(prefix) == 0) CALL errore('openfilepw', 'wrong prefix', 1)
  !
  ! The file with the wavefunctions
  !
  iuwfc = 20 
  lrwfc = 2 * nbnd * npwx * npol 
  CALL diropn(iuwfc, 'wfc', lrwfc, exst) 
  IF (.NOT. exst) CALL errore ('openfilepw', 'file ' // TRIM(prefix) // '.wfc' // ' not found', 1)
  !
  ! file for setting unitary gauges of eigenstates
  !
  lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin_mag
  !
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE openfilepw
  !-----------------------------------------------------------------------
