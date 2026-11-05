!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------------------
SUBROUTINE crpa_set_upert()
   !--------------------------------------------------------------------------------------
   !! Set lr_npert, upert, and upert_mp for a single perturbation without symmetry.
   !--------------------------------------------------------------------------------------
   !
   USE symm_base,    ONLY : s
   USE lr_symm_base, ONLY : nsymq, minus_q, lr_npert, upert, upert_mq
   !
   IMPLICIT NONE
   !
   INTEGER :: ipol, jpol
   !! Counter on perturbations
   INTEGER :: isym
   !! Counter on symmetries
   !
   ! Set symmetry representation in lr_symm_base
   !
   lr_npert = 1
   !
   ALLOCATE(upert(lr_npert, lr_npert, nsymq))
   !
   DO isym = 1, nsymq
      upert(1, 1, isym) = (1.d0, 0.d0)
   ENDDO
   !
   IF (minus_q) THEN
      !
      ! upert_mq is the rotation matrix for symmetry S such that T * S * q = q + G.
      !
      ALLOCATE(upert_mq(lr_npert, lr_npert))
      upert_mq(1, 1) = (1.d0, 0.d0)
   ENDIF ! minus_q
   !
END SUBROUTINE crpa_set_upert
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
SUBROUTINE crpa_deallocate_upert()
!----------------------------------------------------------------------------------------
   USE lr_symm_base, ONLY : nsymq, minus_q, lr_npert, upert, upert_mq
   IMPLICIT NONE
   DEALLOCATE(upert)
   IF (minus_q) DEALLOCATE(upert_mq)
!----------------------------------------------------------------------------------------
END SUBROUTINE crpa_deallocate_upert
!----------------------------------------------------------------------------------------
