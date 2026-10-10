!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE crpa_close_q ( flag )
  !----------------------------------------------------------------------------
  !
  ! This subroutine closes all files.
  ! Called at the end of the run with flag=.TRUE. (removes 'recover')
  ! or during execution with flag=.FALSE. (does not remove 'recover')
  !
  USE buffers,        ONLY : close_buffer
  USE io_files,       ONLY : iunhub
  USE units_lr,       ONLY : iuwfc, iuatswfc, iudwf, iudvwfc
  USE control_lr,     ONLY : lgamma
  USE ldaU,           ONLY : lda_plus_u
  USE crpacom,        ONLY : lrwan, iuwan
  USE crpa_pert,      ONLY : pert_basis
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  LOGICAL :: opnd
  !
  CALL close_buffer(iuwfc,'delete')
  !
  IF (flag) THEN
     CALL close_buffer(iudwf,'delete')
     CALL close_buffer(iudvwfc,'delete')
  ELSE
     CALL close_buffer(iudwf,'keep')
     CALL close_buffer(iudvwfc,'keep')
  ENDIF
  !
  IF (lda_plus_u) THEN
     CALL close_buffer(iuatswfc,'delete')
  ENDIF
  IF (lgamma) CALL close_buffer(iunhub,'delete')
  !
  IF (pert_basis == 'wannier') CALL close_buffer(iuwan,'delete')
  !
END SUBROUTINE crpa_close_q
