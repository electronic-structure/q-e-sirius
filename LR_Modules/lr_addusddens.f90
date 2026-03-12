!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
SUBROUTINE lr_addusddens (npert, dbecsum, drhop)
  !---------------------------------------------------------------------------
  !
  ! Calculate the additional charge in reciprocal space due to US PP's
  ! See Eq.(36) in B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007)
  ! Then sum up the normal and ultrasoft charges.
  ! It assumes that the array dbecsum has already been computed.
  ! Inspired by PH/addusddens.f90
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, g, eigts1, eigts2, eigts3, mill
  USE noncollin_module,     ONLY : nspin_mag
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE qpoint,               ONLY : xq, eigqts
  USE mod_sirius
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npert
  ! input : number of perturbations
  COMPLEX(DP), INTENT(in)    :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, npert)
  ! input : the ultrasoft term
  COMPLEX(DP), INTENT(inout) :: drhop(dfftp%nnr, nspin_mag, npert)
  ! input/output : change of the charge density
  !
  ! the local variables
  !
  INTEGER :: ig, na, nt, ih, jh, is, ijh, ir, na1, ld
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on r vectors
  ! counter on spin
  ! counter on combined beta functions
  INTEGER :: ipert
  !! counter on perturbations
  !
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:, :, :), aux_r(:), dm1(:,:), tmp1(:,:)
  ! the structure factor
  ! q_lm(G)
  ! auxiliary variable for drho(G)
  ! auxiliary variable for drho(r)
  !
  COMPLEX(DP) :: z1
  !
  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('lr_addusddens')
  !
  ALLOCATE (aux(ngm, nspin_mag, npert))
  ALLOCATE (aux_r(dfftp%nnr))
  !
  aux(:, :, :) = (0.d0, 0.d0)

  DO nt = 1, ntyp
     IF (upf(nt)%tvanp) THEN
        ! get number of atoms of the given type
        na1 = 0
        DO na = 1, nat
           IF (ityp (na) .EQ. nt) THEN
              na1 = na1 + 1
           ENDIF
        ENDDO
        ld = nh(nt) * (nh(nt) + 1) / 2
        ALLOCATE(dm1(ld, na1))
        ALLOCATE(tmp1(na1, ngm))
        DO ipert = 1, npert
           DO is = 1, nspin_mag
              na1 = 0
              DO na = 1, nat
                 IF (ityp (na) .EQ. nt) THEN
                    na1 = na1 + 1
                    dm1(1:ld, na1) = dbecsum(1:ld, na, is, ipert)
                 ENDIF
              ENDDO !ia
              CALL zgemm('T', 'N', na1, ngm, ld, cmplx(1.d0, 0.d0, kind=kind(0.0d0)), dm1, ld,&
                         atom_type(nt)%qpw_t, ld, cmplx(0.d0, 0.d0, kind=kind(0.0d0)), tmp1, na1)
              na1 = 0
              DO na = 1, nat
                 IF (ityp (na) .eq.nt) THEN
                    na1 = na1 + 1
!$omp parallel do default(shared) private(z1)
                    DO ig = 1, ngm
                       z1 = eigts1(mill(1,ig),na) * &
                            eigts2(mill(2,ig),na) * &
                            eigts3(mill(3,ig),na) * &
                            eigqts(na)
                        aux(ig, is, ipert) = aux(ig, is, ipert) + 2.0d0 * tmp1(na1, ig) * z1
                    ENDDO !ig
!$omp end parallel do
                 ENDIF
              ENDDO !na
           ENDDO !is
        ENDDO !ipert
        DEALLOCATE(dm1)
        DEALLOCATE(tmp1)
     ENDIF
  ENDDO
  !
  !
  ! Convert aux to real space, and add to the charge density.
  !
  DO ipert = 1, npert
     DO is = 1, nspin_mag
         !
         aux_r(:) = (0.d0, 0.d0)
         !
         DO ig = 1, ngm
            aux_r(dfftp%nl(ig)) = aux(ig, is, ipert)
         ENDDO
         !
         CALL invfft('Rho', aux_r, dfftp)
         !
         DO ir = 1, dfftp%nnr
            drhop(ir, is, ipert) = drhop(ir, is, ipert) + aux_r(ir)
         ENDDO
         !
     ENDDO
  ENDDO ! ipert
  !
  DEALLOCATE (aux)
  DEALLOCATE (aux_r)
  !
  CALL stop_clock ('lr_addusddens')
  !
  RETURN
  !
END SUBROUTINE lr_addusddens
