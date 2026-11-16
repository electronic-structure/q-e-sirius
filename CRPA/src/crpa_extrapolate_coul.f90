!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_extrapolate_coul()
  !-----------------------------------------------------------------------
  !! Extends coul_beta_store/coul_gamma_store (length crpa_coul_lanczos_steps)
  !! into coul_beta_store_ext/coul_gamma_store_ext (length
  !! crpa_coul_lanczos_steps_ext), independently for each of the npert_tot
  !! bare-perturbation chains. beta/gamma are assumed to settle into a
  !! stable plateau value that may alternate in a two-step (odd/even)
  !! pattern; a running average/amplitude of that plateau is computed over
  !! most of the computed chain, rejecting  any coefficient (together with
  !! its very next partner, to keep the odd/even pairing balanced) that
  !! departs from the running average by more than outlier_cutoff. That
  !! plateau (or just its average, if crpa_coul_extrapolation = 'constant')
  !! is then extrapolated for the remaining, cheaper, "virtual" steps
  !! instead of paying for more actual (expensive) Lanczos iterations.
  !! When crpa_coul_extrapolation = 'no', coul_beta/gamma_store_ext are
  !! simply a copy of coul_beta/gamma_store (crpa_coul_lanczos_steps_ext =
  !! crpa_coul_lanczos_steps is enforced in crpa_readin.f90 in that case).
  !
  USE kinds,             ONLY : DP
  USE crpacom,           ONLY : crpa_coul_lanczos_steps, crpa_coul_lanczos_steps_ext, &
                                 crpa_coul_extrapolation
  USE crpa_pert,         ONLY : npert_tot
  USE crpa_coul_lanczos, ONLY : coul_beta_store, coul_gamma_store, &
                                 coul_beta_store_ext, coul_gamma_store_ext
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: outlier_cutoff = 2.0_DP
  !! Same unit as beta: Ry
  !
  INTEGER  :: ipert, i, istart, counter
  REAL(DP) :: average, av_amplitude
  LOGICAL  :: skip
  !
  DO ipert = 1, npert_tot
     !
     ! Shifted by one: see crpa_extrapolate_diel.f90's comment for why
     ! (matches thermo_pw's extrapolate() convention).
     !
     coul_beta_store_ext(ipert, 1:crpa_coul_lanczos_steps-1)  = coul_beta_store(ipert, 2:crpa_coul_lanczos_steps)
     coul_gamma_store_ext(ipert, 1:crpa_coul_lanczos_steps-1) = coul_gamma_store(ipert, 2:crpa_coul_lanczos_steps)
     coul_beta_store_ext(ipert, crpa_coul_lanczos_steps)  = coul_beta_store(ipert, crpa_coul_lanczos_steps)
     coul_gamma_store_ext(ipert, crpa_coul_lanczos_steps) = coul_gamma_store(ipert, crpa_coul_lanczos_steps)
     !
     IF (TRIM(crpa_coul_extrapolation) == 'no') CYCLE
     !
     ! Running average/amplitude of the plateau (odd steps contribute +beta,
     ! even steps -beta, capturing the odd/even oscillation), over most of
     ! the computed chain, with outlier rejection.
     !
     istart = MAX(1, MIN(151, crpa_coul_lanczos_steps/4))
     average = 0.0_DP; av_amplitude = 0.0_DP; counter = 0
     skip = .FALSE.
     !
     DO i = istart, crpa_coul_lanczos_steps
        !
        IF (skip) THEN
           skip = .FALSE.
           CYCLE
        ENDIF
        !
        IF (i /= istart .AND. &
            ABS(coul_beta_store(ipert,i) - average/counter) > outlier_cutoff) THEN
           skip = .TRUE.   ! reject this point and its odd/even partner
           CYCLE
        ENDIF
        !
        average = average + coul_beta_store(ipert,i)
        IF (MOD(i,2) == 1) THEN
           av_amplitude = av_amplitude + coul_beta_store(ipert,i)
        ELSE
           av_amplitude = av_amplitude - coul_beta_store(ipert,i)
        ENDIF
        counter = counter + 1
        !
     ENDDO
     !
     IF (counter > 0) THEN
        average      = average      / counter
        av_amplitude = av_amplitude / counter
     ENDIF
     IF (TRIM(crpa_coul_extrapolation) == 'constant') av_amplitude = 0.0_DP
     !
     DO i = crpa_coul_lanczos_steps+1, crpa_coul_lanczos_steps_ext
        IF (MOD(i,2) == 1) THEN
           coul_beta_store_ext(ipert,i)  = average + av_amplitude
           coul_gamma_store_ext(ipert,i) = average + av_amplitude
        ELSE
           coul_beta_store_ext(ipert,i)  = average - av_amplitude
           coul_gamma_store_ext(ipert,i) = average - av_amplitude
        ENDIF
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE crpa_extrapolate_coul
