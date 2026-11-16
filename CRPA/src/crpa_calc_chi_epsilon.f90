!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_calc_chi_epsilon()
  !-----------------------------------------------------------------------
  !! Reconstructs the complex macroscopic dielectric tensor
  !! epsilon_w_diel(3,3,crpa_nfreq) from the non-Hermitian Lanczos
  !! chain coefficients built by crpa_do_lanczos_epsilon(), by solving the
  !! tridiagonal system representing the continued fraction of the chain,
  !! evaluated at every frequency of crpacom%crpa_freq (w_freq or the
  !! crpa_freq_* grid). Only this final solve depends on frequency: the
  !! chain itself (alpha/beta/gamma/zeta) is built once and reused for the
  !! whole grid, exactly as in TDDFPT's calculate_spectrum/calc_chi.
  !!
  !! If crpa_diel_extrapolation /= 'no', crpa_extrapolate_diel() first
  !! extends the beta/gamma chain out to crpa_diel_lanczos_steps_ext (see
  !! that routine); otherwise the tridiagonal system below is solved with
  !! the chain exactly as computed (steps_ext = steps).
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE constants,  ONLY : fpi, e2
  USE cell_base,  ONLY : omega
  USE lsda_mod,   ONLY : nspin
  USE symme,      ONLY : symmatrix, crys_to_cart
  USE crpacom,    ONLY : crpa_done_epsil, crpa_diel_lanczos_steps, &
                          crpa_diel_lanczos_steps_ext, crpa_nfreq, crpa_freq, &
                          epsilon_w_diel, crpa_spectrum_only, &
                          plasma_freq2, crpa_metal, crpa_diel_intraband, &
                          w_drude_min
  USE crpa_diel_lanczos, ONLY : diel_beta_store, diel_beta_store_ext, &
                                diel_gamma_store_ext, diel_zeta_store
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: a(:), b(:), c(:), r(:), b0(:), c0(:)
  COMPLEX(DP), EXTERNAL :: zdotc
  INTEGER :: ipol, jpol, ifreq, i, info
  REAL(DP) :: degspin, repsilon_w(3,3), iepsilon_w(3,3), fact, wdr
  !
  CALL start_clock('crpa_calc_chi')
  !
  ! Save the raw chain before extrapolation, so it can be reused later at a
  ! different frequency grid/extrapolation scheme (crpa_spectrum_only,
  ! crpa_run_spectrum_only.f90). Not repeated when the chain itself was just
  ! read back from that same file (crpa_spectrum_only = .true.).
  !
  IF (.NOT. crpa_spectrum_only) CALL crpa_write_diel_lanczos_chain()
  !
  CALL crpa_extrapolate_diel()
  !
  degspin = 2.0_DP
  IF (nspin /= 1) degspin = 1.0_DP
  !
  ALLOCATE(a(crpa_diel_lanczos_steps_ext))
  ALLOCATE(b(crpa_diel_lanczos_steps_ext))
  ALLOCATE(c(crpa_diel_lanczos_steps_ext))
  ALLOCATE(r(crpa_diel_lanczos_steps_ext))
  ALLOCATE(b0(crpa_diel_lanczos_steps_ext))
  ALLOCATE(c0(crpa_diel_lanczos_steps_ext))
  !
  epsilon_w_diel(:,:,:) = (0.0_DP, 0.0_DP)
  !
  DO ipol = 1, 3
     !
     ! -beta, -gamma (the sub-/super-diagonal of the tridiagonal Lanczos
     ! matrix) do not depend on frequency: built once per chain (ipol) and
     ! copied fresh into b/c before each ZGTSV call below, since ZGTSV
     ! overwrites its DL/D/DU/B arguments in place.
     !
     b0(:) = (0.0_DP, 0.0_DP)
     c0(:) = (0.0_DP, 0.0_DP)
     DO i = 1, crpa_diel_lanczos_steps_ext - 1
        b0(i) = CMPLX(-diel_beta_store_ext(ipol,i), 0.0_DP, DP)
        c0(i) = CMPLX(-diel_gamma_store_ext(ipol,i), 0.0_DP, DP)
     ENDDO
     !
     DO ifreq = 1, crpa_nfreq
        !
        a(:) = crpa_freq(ifreq)
        b(:) = b0(:)
        c(:) = c0(:)
        r(:) = (0.0_DP, 0.0_DP)
        r(1) = (1.0_DP, 0.0_DP)
        !
        ! |w_t> = (w - L)^-1 |1,0,0,...,0>
        !
        CALL ZGTSV(crpa_diel_lanczos_steps_ext, 1, b, a, c, r, &
                   crpa_diel_lanczos_steps_ext, info)
        !
        IF (info /= 0) CALL errore('crpa_calc_chi_epsilon', &
           'unable to solve the tridiagonal Lanczos system', 1)
        !
        ! chi = -<zeta|w_t>, renormalized by the norm of the first chain
        ! vector. See Malcioglu et al., Comput. Phys. Commun. 182, 1744
        ! (2011), Eq. (35), and calc_chi.f90 (thermo_pw) for the same formula
        ! (valid for both the pseudo-Hermitian and non-Hermitian chains).
        ! Only the first crpa_diel_lanczos_steps entries of r enter here:
        ! zeta is only defined (from the actual chain) up to that point,
        ! even though r was obtained from the longer, possibly extrapolated,
        ! tridiagonal system.
        !
        DO jpol = 1, 3
           epsilon_w_diel(ipol,jpol,ifreq) = zdotc(crpa_diel_lanczos_steps, &
                diel_zeta_store(ipol,jpol,:), 1, r, 1)
           epsilon_w_diel(ipol,jpol,ifreq) = epsilon_w_diel(ipol,jpol,ifreq) &
                * CMPLX(diel_beta_store(ipol,1), 0.0_DP, DP)
           epsilon_w_diel(ipol,jpol,ifreq) = epsilon_w_diel(ipol,jpol,ifreq) &
                * CMPLX(-2.0_DP*degspin, 0.0_DP, DP)
        ENDDO
        !
     ENDDO ! ifreq
     !
  ENDDO ! ipol
  !
  DEALLOCATE(a, b, c, r, b0, c0)
  !
  ! Symmetrize (as in crpa_dielec.f90) and add the diagonal part, frequency
  ! by frequency.
  !
  fact = fpi * e2 / omega
  DO ifreq = 1, crpa_nfreq
     repsilon_w = DBLE(epsilon_w_diel(:,:,ifreq))
     CALL crys_to_cart(repsilon_w)
     CALL symmatrix(repsilon_w)
     iepsilon_w = AIMAG(epsilon_w_diel(:,:,ifreq))
     CALL crys_to_cart(iepsilon_w)
     CALL symmatrix(iepsilon_w)
     epsilon_w_diel(:,:,ifreq) = CMPLX(repsilon_w, iepsilon_w, KIND=DP) * fact
     !
     ! Intraband (Drude) contribution, metals only: see crpa_dielec.f90 for
     ! the same term on the CG path. plasma_freq2 (crpa_plasma_frequency.f90)
     ! is already in cartesian axis and symmetrized.
     !
     IF (crpa_metal .AND. crpa_diel_intraband) THEN
        ! -omega_p^2 / ( w * (w+i*eta) ), w = Re(crpa_freq(ifreq))
        ! Re(w) = 0 (static, with or without broadening) is floored
        ! exactly as on the CG path; see w_drude_min.
        wdr = DBLE(crpa_freq(ifreq))
        IF (ABS(wdr) < w_drude_min) THEN
           wdr = w_drude_min
           WRITE(stdout,'(5x,A,ES9.2,A)') 'Drude term: Re(w) = 0; denominator floored to Re(w) = ', w_drude_min,' Ry.'
        ENDIF
        DO ipol = 1, 3
           DO jpol = 1, 3
              epsilon_w_diel(ipol,jpol,ifreq) = epsilon_w_diel(ipol,jpol,ifreq) &
                 - CMPLX(plasma_freq2(ipol,jpol), 0.d0, KIND=DP) &
                   / ( CMPLX(wdr, 0.d0, KIND=DP) * CMPLX(wdr, AIMAG(crpa_freq(ifreq)), KIND=DP) )
           ENDDO
        ENDDO
     ENDIF
     !
     DO ipol = 1, 3
        epsilon_w_diel(ipol,ipol,ifreq) = epsilon_w_diel(ipol,ipol,ifreq) + 1.0_DP
     ENDDO
  ENDDO
  !
  crpa_done_epsil = .TRUE.
  CALL crpa_summarize_epsilon_diel()
  CALL write_eps_file('epsilon_w.dat', crpa_nfreq, crpa_freq, epsilon_w_diel)
  !
  ! The Lanczos chain arrays (allocated in crpa_do_lanczos_epsilon) are no
  ! longer needed once epsilon_w_diel has been built from them.
  !
  CALL crpa_deallocate_diel_lanczos()
  !
  CALL stop_clock('crpa_calc_chi')
  !
  RETURN
  !
END SUBROUTINE crpa_calc_chi_epsilon
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_summarize_epsilon_diel()
  !-----------------------------------------------------------------------
  !! Print the complex macroscopic dielectric tensor epsilon_w_diel(3,3,:)
  !! to stdout, at every frequency of the grid.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE crpacom,    ONLY : crpa_nfreq, crpa_freq, epsilon_w_diel
  !
  IMPLICIT NONE
  INTEGER :: ipol, jpol, ifreq
  !
  DO ifreq = 1, crpa_nfreq
     WRITE(stdout, '(/,10x,"Dielectric tensor in cartesian axis at w = (", &
          &f12.6,",",f12.6,") Ry",/)') DBLE(crpa_freq(ifreq)), AIMAG(crpa_freq(ifreq))
     DO ipol = 1, 3
        WRITE(stdout, '(10x,3("(",f15.6,",",f15.6,")",2x))') &
             (DBLE(epsilon_w_diel(ipol,jpol,ifreq)), &
              AIMAG(epsilon_w_diel(ipol,jpol,ifreq)), jpol=1,3)
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE crpa_summarize_epsilon_diel
