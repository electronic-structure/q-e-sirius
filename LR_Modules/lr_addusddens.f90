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
  INTEGER :: ig, na, nt, ih, jh, is, ijh, ir
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
  REAL(DP), ALLOCATABLE :: qmod(:), qpg(:,:), ylmk0(:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics
  !
  COMPLEX(DP), ALLOCATABLE :: sk(:), qgm(:), aux(:, :, :), aux_r(:)
  ! the structure factor
  ! q_lm(G)
  ! auxiliary variable for drho(G)
  ! auxiliary variable for drho(r)
  !
  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('lr_addusddens')
  !
  ALLOCATE (aux(ngm, nspin_mag, npert))
  ALLOCATE (aux_r(dfftp%nnr))
  ALLOCATE (sk(ngm))
  ALLOCATE (ylmk0(ngm,lmaxq * lmaxq))
  ALLOCATE (qgm(ngm))
  ALLOCATE (qmod(ngm))
  ALLOCATE (qpg(3,ngm))
  !
  aux(:, :, :) = (0.d0, 0.d0)
  !
  ! Calculate the q+G vector, its modulus, and the spherical harmonics.
  !
  CALL setqmod (ngm, xq, g, qmod, qpg)
  !
  CALL ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
  !
  DO ig = 1, ngm
     qmod(ig) = sqrt(qmod(ig)) * tpiba
  ENDDO
  !
  ! TODO: check performance. There was a bottleneck here which we fixed with the following call
  !#if defined(__SIRIUS)
  !DO nt = 1, ntyp
  !   IF (upf(nt)%tvanp) THEN
  !      CALL sirius_generate_rhoaug_q(gs_handler, nt, nat, ngm, nspin_mag, atom_type(nt)%qpw_t, &
  !          & nh(nt) * (nh(nt) + 1) / 2, eigqts, mill, dbecsum, nhm * (nhm + 1) / 2, aux)
  !   ENDIF
  !ENDDO
  !#else
  !
  ! in the newer version becsum and aux dimensions have changed
  ! logic of the code is clear: compute Q_aug once and then add to all pertubations, but the benefit is not
  ! clear to me as Q_aug depends on the perturbation index and all of them have to be generated

  DO nt = 1, ntyp
     IF (upf(nt)%tvanp) THEN
        ijh = 0
        DO ih = 1, nh (nt)
           DO jh = ih, nh (nt)
              !
              ! Calculate the Fourier transform of the Q functions,
              ! and put the result in qgm.
              !
              CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              !
              ijh = ijh + 1
              DO na = 1, nat
                 IF (ityp (na) .eq.nt) THEN
                    !
                    ! Calculate the second term in Eq.(36) of the ultrasoft paper.
                    !
                    ! TODO: add omp statements here
!                       !$omp parallel default(shared) private(is, z1)
!                       DO is = 1, nspin_mag
!   !$omp do
!                          DO ig = 1, ngm
                    DO ipert = 1, npert
                       DO is = 1, nspin_mag
                          DO ig = 1, ngm
                             !
                             ! Calculate the structure factor
                             !
                             sk(ig) = eigts1(mill(1,ig),na) * &
                                      eigts2(mill(2,ig),na) * &
                                      eigts3(mill(3,ig),na) * &
                                      eigqts(na)
                             !
                             aux(ig, is, ipert) = aux(ig, is, ipert) &
                                + 2.0d0 * qgm(ig) * sk(ig) * dbecsum(ijh,na,is, ipert)
                             !
                          ENDDO
                       ENDDO
                    ENDDO ! ipert
                    !
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
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
  DEALLOCATE (qpg)
  DEALLOCATE (qmod)
  DEALLOCATE (qgm)
  DEALLOCATE (ylmk0)
  DEALLOCATE (sk)
  DEALLOCATE (aux)
  DEALLOCATE (aux_r)
  !
  CALL stop_clock ('lr_addusddens')
  !
  RETURN
  !
END SUBROUTINE lr_addusddens
