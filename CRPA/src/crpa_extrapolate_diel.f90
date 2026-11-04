!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_extrapolate_diel()
  !-----------------------------------------------------------------------
  !! Extends diel_beta_store/diel_gamma_store (length crpa_diel_lanczos_steps)
  !! into diel_beta_store_ext/diel_gamma_store_ext (length
  !! crpa_diel_lanczos_steps_ext), independently for each of the 3
  !! polarizations, following TDDFPT/src/turbo_spectrum.f90): beta/gamma are assumed to
  !! settle into a stable plateau value that may alternate in a two-step
  !! (odd/even) pattern; a running average/amplitude of that plateau is
  !! computed over most of the computed chain, rejecting any coefficient
  !! (together with its very next partner, to keep the odd/even pairing
  !! balanced) that departs from the running average by more than
  !! outlier_cutoff. That plateau (or just its average, if
  !! crpa_diel_extrapolation = 'constant') is then extrapolated for the
  !! remaining, cheaper, "virtual" steps instead of paying for more actual
  !! (expensive) Lanczos iterations.
  !! When crpa_diel_extrapolation = 'no', diel_beta/gamma_store_ext are
  !! simply a copy of diel_beta/gamma_store (crpa_diel_lanczos_steps_ext =
  !! crpa_diel_lanczos_steps is enforced in crpa_readin.f90 in that case).
  !
  USE kinds,             ONLY : DP
  USE crpacom,           ONLY : crpa_diel_lanczos_steps, crpa_diel_lanczos_steps_ext, &
                                 crpa_diel_extrapolation
  USE crpa_diel_lanczos, ONLY : diel_beta_store, diel_gamma_store, &
                                 diel_beta_store_ext, diel_gamma_store_ext
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: outlier_cutoff = 2.0_DP
  !! Same unit as beta: Ry
  !
  INTEGER  :: ipol, i, istart, counter
  REAL(DP) :: average, av_amplitude
  LOGICAL  :: skip
  !
  DO ipol = 1, 3
     !
     ! Shifted by one: the tridiagonal solve in crpa_calc_chi_epsilon.f90
     ! places diel_beta/gamma_store_ext(i) at the (i+1,i)/(i,i+1) matrix
     ! entries, i.e. it needs beta_{i+1}/gamma_{i+1} (not beta_i/gamma_i)
     ! at position i, since A q_j = beta_{j+1} q_{j+1} + gamma_j q_{j-1}
     ! (see LR_Modules/lanczos_nonhermitian.f90's recursion). This matches
     ! thermo_pw's own extrapolate() convention:
     !   beta_store_ext(1:lanczos_steps-1) = beta_store(2:lanczos_steps)
     !   beta_store_ext(lanczos_steps)     = beta_store(lanczos_steps)
     !
     diel_beta_store_ext(ipol, 1:crpa_diel_lanczos_steps-1)  = diel_beta_store(ipol, 2:crpa_diel_lanczos_steps)
     diel_gamma_store_ext(ipol, 1:crpa_diel_lanczos_steps-1) = diel_gamma_store(ipol, 2:crpa_diel_lanczos_steps)
     diel_beta_store_ext(ipol, crpa_diel_lanczos_steps)  = diel_beta_store(ipol, crpa_diel_lanczos_steps)
     diel_gamma_store_ext(ipol, crpa_diel_lanczos_steps) = diel_gamma_store(ipol, crpa_diel_lanczos_steps)
     !
     IF (TRIM(crpa_diel_extrapolation) == 'no') CYCLE
     !
     ! Running average/amplitude of the plateau (odd steps contribute +beta,
     ! even steps -beta, capturing the odd/even oscillation), over most of
     ! the computed chain, with outlier rejection.
     !
     istart = MAX(1, MIN(151, crpa_diel_lanczos_steps/4))
     average = 0.0_DP; av_amplitude = 0.0_DP; counter = 0
     skip = .FALSE.
     !
     DO i = istart, crpa_diel_lanczos_steps
        !
        IF (skip) THEN
           skip = .FALSE.
           CYCLE
        ENDIF
        !
        IF (i /= istart .AND. &
            ABS(diel_beta_store(ipol,i) - average/counter) > outlier_cutoff) THEN
           skip = .TRUE.   ! reject this point and its odd/even partner
           CYCLE
        ENDIF
        !
        average = average + diel_beta_store(ipol,i)
        IF (MOD(i,2) == 1) THEN
           av_amplitude = av_amplitude + diel_beta_store(ipol,i)
        ELSE
           av_amplitude = av_amplitude - diel_beta_store(ipol,i)
        ENDIF
        counter = counter + 1
        !
     ENDDO
     !
     IF (counter > 0) THEN
        average      = average      / counter
        av_amplitude = av_amplitude / counter
     ENDIF
     IF (TRIM(crpa_diel_extrapolation) == 'constant') av_amplitude = 0.0_DP
     !
     DO i = crpa_diel_lanczos_steps+1, crpa_diel_lanczos_steps_ext
        IF (MOD(i,2) == 1) THEN
           diel_beta_store_ext(ipol,i)  = average + av_amplitude
           diel_gamma_store_ext(ipol,i) = average + av_amplitude
        ELSE
           diel_beta_store_ext(ipol,i)  = average - av_amplitude
           diel_gamma_store_ext(ipol,i) = average - av_amplitude
        ENDIF
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE crpa_extrapolate_diel
