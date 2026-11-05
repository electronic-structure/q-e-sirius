!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_write_diel_lanczos_chain()
  !-----------------------------------------------------------------------
  !! Writes the RAW (pre-extrapolation) non-Hermitian Lanczos chain
  !! coefficients built by crpa_do_lanczos_epsilon() to a formatted file
  !! ('diel_lanczos_chain.dat', in the run's working directory), so that
  !! epsilon(w) can later be reconstructed at a different frequency grid or
  !! extrapolation scheme (crpa_spectrum_only = .true., see
  !! crpa_run_spectrum_only.f90) without redoing the chain itself: only the
  !! final continued-fraction solve in crpa_calc_chi_epsilon.f90 depends on
  !! frequency. Mirrors the role of TDDFPT's lr_write_restart.f90 (which
  !! this format is loosely modeled on), but self-contained in a single
  !! file with no separate wavefunction/projection files.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode, stdout
  USE lsda_mod,   ONLY : nspin
  USE crpacom,    ONLY : crpa_diel_lanczos_steps, plasma_freq2, crpa_metal
  USE crpa_diel_lanczos, ONLY : diel_beta_store, diel_gamma_store, diel_zeta_store
  !
  IMPLICIT NONE
  !
  INTEGER :: iun, ios, ipol, jpol
  REAL(DP) :: degspin
  !
  IF (.NOT. ionode) RETURN
  !
  degspin = 2.0_DP
  IF (nspin /= 1) degspin = 1.0_DP
  !
  OPEN(NEWUNIT=iun, FILE='diel_lanczos_chain.dat', STATUS='replace', IOSTAT=ios)
  IF (ios /= 0) CALL errore('crpa_write_diel_lanczos_chain', &
     'cannot open diel_lanczos_chain.dat', 1)
  !
  WRITE(iun, '(a)') "# CRPA dielectric Lanczos chain (3 lock-step polarizations, q=Gamma)"
  WRITE(iun, *) crpa_diel_lanczos_steps
  WRITE(iun, *) degspin
  !
  DO ipol = 1, 3
     WRITE(iun, *) diel_beta_store(ipol, :)
  ENDDO
  DO ipol = 1, 3
     WRITE(iun, *) diel_gamma_store(ipol, :)
  ENDDO
  DO ipol = 1, 3
     DO jpol = 1, 3
        WRITE(iun, *) diel_zeta_store(ipol, jpol, :)
     ENDDO
  ENDDO
  !
  ! Intraband (Drude) plasma-frequency tensor (crpa_plasma_frequency.f90):
  ! computed once, from wavefunctions, while building the chain
  ! (crpa_do_lanczos_epsilon), and saved here so crpa_spectrum_only mode
  ! (crpa_run_spectrum_only.f90 -> crpa_calc_chi_epsilon, no live
  ! wavefunctions) can still add the metallic head correction.
  !
  WRITE(iun, *) crpa_metal
  DO ipol = 1, 3
     WRITE(iun, *) plasma_freq2(ipol, :)
  ENDDO
  !
  CLOSE(iun)
  !
  WRITE(stdout, '(5X,A)') "INFO: Dielectric Lanczos chain written to diel_lanczos_chain.dat"
  !
  RETURN
  !
END SUBROUTINE crpa_write_diel_lanczos_chain
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_read_diel_lanczos_chain()
  !-----------------------------------------------------------------------
  !! Reads back the chain written by crpa_write_diel_lanczos_chain(),
  !! allocating diel_beta_store/diel_gamma_store/diel_zeta_store (and the
  !! _ext arrays used by crpa_extrapolate_diel/crpa_calc_chi_epsilon) from
  !! the file's own chain length, overriding crpa_diel_lanczos_steps (and,
  !! if needed, crpa_diel_lanczos_steps_ext) from the current namelist:
  !! only the frequencies and extrapolation scheme are taken fresh from
  !! the current input, exactly as intended for crpa_spectrum_only mode.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE crpacom,    ONLY : crpa_diel_lanczos_steps, crpa_diel_lanczos_steps_ext, &
                          crpa_diel_extrapolation, crpa_nfreq, epsilon_w_diel, &
                          plasma_freq2, crpa_metal
  USE crpa_diel_lanczos, ONLY : diel_beta_store, diel_gamma_store, diel_zeta_store, &
                                 diel_beta_store_ext, diel_gamma_store_ext
  !
  IMPLICIT NONE
  !
  INTEGER :: iun, ios, ipol, jpol, nsteps
  REAL(DP) :: degspin
  CHARACTER(LEN=256) :: line
  LOGICAL :: exst
  !
  INQUIRE(FILE='diel_lanczos_chain.dat', EXIST=exst)
  IF (.NOT. exst) CALL errore('crpa_read_diel_lanczos_chain', &
     'diel_lanczos_chain.dat not found', 1)
  !
  OPEN(NEWUNIT=iun, FILE='diel_lanczos_chain.dat', STATUS='old', IOSTAT=ios)
  IF (ios /= 0) CALL errore('crpa_read_diel_lanczos_chain', &
     'cannot open diel_lanczos_chain.dat', 1)
  !
  READ(iun, '(A)') line
  READ(iun, *) nsteps
  READ(iun, *) degspin
  !
  crpa_diel_lanczos_steps = nsteps
  IF (TRIM(crpa_diel_extrapolation) == 'no') THEN
     crpa_diel_lanczos_steps_ext = crpa_diel_lanczos_steps
  ELSE IF (crpa_diel_lanczos_steps_ext < crpa_diel_lanczos_steps) THEN
     CALL errore('crpa_read_diel_lanczos_chain', &
        'crpa_diel_lanczos_steps_ext must be >= the chain length stored on file', 1)
  ENDIF
  !
  IF (ALLOCATED(diel_beta_store))  DEALLOCATE(diel_beta_store)
  IF (ALLOCATED(diel_gamma_store)) DEALLOCATE(diel_gamma_store)
  IF (ALLOCATED(diel_zeta_store))  DEALLOCATE(diel_zeta_store)
  IF (ALLOCATED(diel_beta_store_ext))  DEALLOCATE(diel_beta_store_ext)
  IF (ALLOCATED(diel_gamma_store_ext)) DEALLOCATE(diel_gamma_store_ext)
  IF (ALLOCATED(epsilon_w_diel))       DEALLOCATE(epsilon_w_diel)
  !
  ! Normally allocated by crpa_allocate_diel_lanczos(), not called in
  ! crpa_spectrum_only mode: the frequencies are those of the CURRENT input
  ! (may differ from the run that produced the saved chain).
  !
  ALLOCATE(epsilon_w_diel(3, 3, crpa_nfreq))
  epsilon_w_diel(:,:,:) = (0.0_DP, 0.0_DP)
  !
  ALLOCATE(diel_beta_store(3, nsteps))
  ALLOCATE(diel_gamma_store(3, nsteps))
  ALLOCATE(diel_zeta_store(3, 3, nsteps))
  ALLOCATE(diel_beta_store_ext(3, crpa_diel_lanczos_steps_ext))
  ALLOCATE(diel_gamma_store_ext(3, crpa_diel_lanczos_steps_ext))
  !
  DO ipol = 1, 3
     READ(iun, *) diel_beta_store(ipol, :)
  ENDDO
  DO ipol = 1, 3
     READ(iun, *) diel_gamma_store(ipol, :)
  ENDDO
  DO ipol = 1, 3
     DO jpol = 1, 3
        READ(iun, *) diel_zeta_store(ipol, jpol, :)
     ENDDO
  ENDDO
  !
  ! Intraband (Drude) plasma-frequency tensor: see crpa_write_diel_lanczos_chain.
  ! Older chain files without this record will fail here; regenerate the
  ! chain once with the updated code to add it.
  !
  READ(iun, *) crpa_metal
  DO ipol = 1, 3
     READ(iun, *) plasma_freq2(ipol, :)
  ENDDO
  !
  CLOSE(iun)
  !
  WRITE(stdout, '(5X,A,I8,A)') "INFO: Dielectric Lanczos chain read from diel_lanczos_chain.dat (", &
     nsteps, " steps)"
  !
  RETURN
  !
END SUBROUTINE crpa_read_diel_lanczos_chain
