!
! Copyright (C) 2009-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE compute_deff( deff, et )
  !-----------------------------------------------------------------------
  !!  This routine computes the effective value of the D-eS coefficients
  !!  which appear often in many expressions in the US or PAW case. 
  !!  This routine is for the collinear case.
  !
  USE kinds,       ONLY: DP
  USE ions_base,   ONLY: nsp, nat, ityp
  USE uspp,        ONLY: deeq, qq_at, okvan
  USE uspp_param,  ONLY: nhm
  USE lsda_mod,    ONLY: current_spin
  !
  USE uspp_gpum, ONLY : using_deeq, using_qq_at
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: et
  !! The eigenvalues of the hamiltonian
  REAL(DP), INTENT(OUT) :: deff(nhm,nhm,nat)
  !! Effective values of the D-eS coefficients
  !
  ! ... local variables
  !
  INTEGER :: nt, na, is
  !
  CALL using_deeq(0)
  CALL using_qq_at(0)
  !
  deff(:,:,:) = deeq(:,:,:,current_spin)
  !
  IF (okvan) THEN
     !
     DO nt = 1, nsp
        DO na = 1, nat
           !
           IF ( ityp(na) == nt ) THEN
              deff(:,:,na) = deff(:,:,na) - et*qq_at(:,:,na)
           ENDIF
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE compute_deff
!
!
!---------------------------------------------------------------------------
SUBROUTINE compute_deff_nc( deff, et )
  !-------------------------------------------------------------------------
  !! This routine computes the effective value of the D-eS coefficients
  !! which appears often in many expressions. This routine is for the
  !! noncollinear case.
  !
  USE kinds,            ONLY: DP
  USE ions_base,        ONLY: nsp, nat, ityp
  USE spin_orb,         ONLY: lspinorb
  USE noncollin_module, ONLY: npol
  USE uspp,             ONLY: deeq_nc, qq_at, qq_so, okvan
  USE uspp_param,       ONLY: nhm
  USE lsda_mod,         ONLY: nspin
  !
  USE uspp_gpum, ONLY : using_deeq_nc, using_qq_at, using_qq_so
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: et
  !! The eigenvalues of the hamiltonian
  COMPLEX(DP), INTENT(OUT) :: deff(nhm,nhm,nat,nspin) 
  !! Effective values of the D-eS coefficients
  !
  ! ... local variables
  !
  INTEGER :: nt, na, is, js, ijs
  !
  CALL using_deeq_nc(0)
  IF (.not. lspinorb) CALL using_qq_at(0)
  IF (lspinorb) CALL using_qq_so(0)
  !
  deff=deeq_nc
  IF (okvan) THEN
     !
     DO nt = 1, nsp
        DO na = 1, nat
           !
           IF ( ityp(na) == nt ) THEN
              IF (lspinorb) THEN
                 deff(:,:,na,:) = deff(:,:,na,:) - et * qq_so(:,:,:,nt)
              ELSE
                 ijs=0
                 !
                 DO is=1,npol
                    DO js=1,npol
                       !
                       ijs=ijs+1
                       IF (is==js) deff(:,:,na,ijs)=deff(:,:,na,ijs)-et*qq_at(:,:,na)
                       !
                    ENDDO
                 ENDDO
                 !
              ENDIF
           ENDIF
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE compute_deff_nc
