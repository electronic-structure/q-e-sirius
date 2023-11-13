MODULE mod_lr_addons

CONTAINS

  SUBROUTINE generate_qpw()
    USE mod_sirius
    USE kinds,      ONLY : DP
    USE ions_base,  ONLY : nsp
    USE gvect,      ONLY : ngm, g, gg
    USE uspp_param, ONLY : nh, lmaxq
    USE control_lr, ONLY : lgamma
    USE qpoint,     ONLY : xq
    USE cell_base,  ONLY : tpiba
    !
    IMPLICIT NONE
    INTEGER :: iat, ih, jh, ijh, ig
    REAL(DP), ALLOCATABLE :: qmod (:), qg (:,:), ylmk0 (:,:)

    ALLOCATE (ylmk0(ngm , lmaxq * lmaxq))
    ALLOCATE (qmod (ngm))

    IF (.NOT.lgamma) THEN
      ALLOCATE (qg (3,  ngm))
      CALL setqmod (ngm, xq, g, qmod, qg)
      CALL ylmr2 (lmaxq * lmaxq, ngm, qg, qmod, ylmk0)
!$omp parallel do default(shared)
      DO ig = 1, ngm
        qmod (ig) = SQRT (qmod (ig) ) * tpiba
      ENDDO
!$omp end parallel do
      DEALLOCATE(qg)
    ELSE
      CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
!$omp parallel do default(shared)
      DO ig = 1, ngm
        qmod (ig) = SQRT (gg (ig) ) * tpiba
      ENDDO
!$omp end parallel do
    ENDIF

    DO iat = 1, nsp
      IF (ALLOCATED(atom_type(iat)%qpw)) DEALLOCATE(atom_type(iat)%qpw)
      ALLOCATE(atom_type(iat)%qpw( ngm, nh(iat) * (1 + nh(iat)) / 2 ))
      ijh = 0
      DO ih = 1, nh (iat)
        DO jh = ih, nh (iat)
          ijh = ijh + 1
          CALL qvan2 (ngm, ih, jh, iat, qmod, atom_type(iat)%qpw(:, ijh), ylmk0)
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(ylmk0, qmod)

  END SUBROUTINE


END MODULE
