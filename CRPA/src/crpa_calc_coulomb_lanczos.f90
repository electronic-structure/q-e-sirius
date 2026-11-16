!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_calc_coulomb_lanczos(iq)
  !-----------------------------------------------------------------------
  !! Reconstructs the frequency-resolved screened Coulomb matrix elements
  !! v_coul_scrd(nmels_tot, npert_tot, crpa_nfreq) from the
  !! non-Hermitian Lanczos chains built by crpa_do_lanczos_coulomb(), by
  !! solving the tridiagonal system representing the continued fraction of
  !! each chain, evaluated at every frequency of crpacom%crpa_freq (the same
  !! frequencies as crpa_calc_chi_epsilon). Only this final solve depends on frequency:
  !! each chain (beta/gamma/zeta) is built once and reused for the whole
  !! grid, exactly as in crpa_calc_chi_epsilon.f90.
  !!
  !! See the module-level comment in crpa_do_lanczos_coulomb.f90 for why
  !! the INDUCED part of v_coul_scrd(j,i,w) is
  !! <coul_zeta_store(i,j,:) | r(w)> * coul_beta_store(i,1) (no real-space
  !! reconstruction needed). v_coul_bare(j,i) (frequency-independent) is
  !! added on top of that below, since v_coul_scrd/crpa_compute_mel_wannier
  !! (the CG path) defines its "screened" result as bare+induced combined
  !! (dvpsi = dvbare*psi + dvscf*psi), not the induced correction alone.
  !!
  !! If crpa_coul_extrapolation /= 'no', crpa_extrapolate_coul() first
  !! extends each chain's beta/gamma out to crpa_coul_lanczos_steps_ext (see
  !! that routine); otherwise the tridiagonal system below is solved with
  !! the chain exactly as computed (steps_ext = steps). Only the first
  !! crpa_coul_lanczos_steps entries of r enter the zeta dot product below:
  !! zeta is only defined (from the actual chain) up to that point, even
  !! though r was obtained from the longer, possibly extrapolated,
  !! tridiagonal system (same convention as crpa_calc_chi_epsilon.f90).
  !!
  !! NOTE: this uses v_coul_bare as computed by crpa_do_lanczos_coulomb,
  !! i.e. BEFORE crpa_q0correction() (crpa_main.f90) has (potentially)
  !! shifted its q=0 diagonal entries: crpa_llanczos does not yet apply
  !! that q=0 correction to its own screened result (only the CG path's
  !! v_coul_bare/v_coul_scrd get corrected). Deliberately kept this way for
  !! now (extending the q0 correction to the Lanczos path is future work):
  !! summing the UNCORRECTED bare term with the chain's own SCF/induced
  !! contribution is the physically well-defined quantity this path
  !! currently produces.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout, ionode
  USE lsda_mod,   ONLY : nspin
  USE crpacom,    ONLY : crpa_coul_lanczos_steps, crpa_coul_lanczos_steps_ext, &
                          crpa_nfreq, crpa_freq, filU, v_coul_bare, v_coul_scrd
  USE klist,      ONLY : lgauss, ltetra
  USE qpoint,     ONLY : xq
  USE crpa_pert,  ONLY : npert_tot, nmels_tot
  USE crpa_coul_lanczos, ONLY : coul_beta_store, coul_beta_store_ext, &
                                 coul_gamma_store_ext, coul_zeta_store
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !! index of the q point
  !
  COMPLEX(DP), ALLOCATABLE :: a(:), b(:), c(:), r(:), b0(:), c0(:)
  COMPLEX(DP), EXTERNAL :: zdotc
  INTEGER :: ipert, ifreq, i, info
  REAL(DP) :: degspin
  !
  CALL start_clock('crpa_calc_coul_lncz')
  !
  ! The raw chain (before extrapolation) is saved by the caller,
  ! crpa_do_lanczos_coulomb (crpa_write_coul_lanczos_chain), not here: this
  ! routine also runs on chains read back from that file (crpa_spectrum_only,
  ! and the q-level restart, which may have truncated them).
  !
  CALL crpa_extrapolate_coul()
  !
  degspin = 2.0_DP
  IF (nspin /= 1) degspin = 1.0_DP
  !
  ALLOCATE(a(crpa_coul_lanczos_steps_ext))
  ALLOCATE(b(crpa_coul_lanczos_steps_ext))
  ALLOCATE(c(crpa_coul_lanczos_steps_ext))
  ALLOCATE(r(crpa_coul_lanczos_steps_ext))
  ALLOCATE(b0(crpa_coul_lanczos_steps_ext))
  ALLOCATE(c0(crpa_coul_lanczos_steps_ext))
  !
  v_coul_scrd(:,:,:) = (0.0_DP, 0.0_DP)
  !
  ! The Fermi-energy shift for metals exists only at q = Gamma  and is
  ! absent by construction in the Lanczos chain. For a metal the w = (0,0) value
  ! returned below is therefore the w -> 0 limit of the dynamical response, NOT
  ! the adiabatic static value the CG solve gives.
  !
  IF ((lgauss .OR. ltetra) .AND. ALL(xq(:) == 0.0_DP) .AND. &
      ANY(ABS(crpa_freq(1:crpa_nfreq)) < 1.0d-10)) &
     WRITE(stdout,'(5x,A)') &
        'Warning: metal at w = (0,0). The Lanczos chain has no static branch,&
        & no Fermi-energy shift is applied'
  !
  DO ipert = 1, npert_tot
     !
     ! -beta, -gamma (sub-/super-diagonal of the tridiagonal Lanczos
     ! matrix) do not depend on frequency: built once per chain (using the
     ! possibly-extrapolated _ext arrays) and copied fresh into b/c before
     ! each ZGTSV call, since ZGTSV overwrites its DL/D/DU/B arguments in
     ! place.
     !
     b0(:) = (0.0_DP, 0.0_DP)
     c0(:) = (0.0_DP, 0.0_DP)
     DO i = 1, crpa_coul_lanczos_steps_ext - 1
        b0(i) = CMPLX(-coul_beta_store_ext(ipert,i), 0.0_DP, DP)
        c0(i) = CMPLX(-coul_gamma_store_ext(ipert,i), 0.0_DP, DP)
     ENDDO
     !
     DO ifreq = 1, crpa_nfreq
        !
        a(:) = crpa_freq(ifreq)
        b(:) = b0(:)
        c(:) = c0(:)
        !
        ! |w_t> = (w - L)^-1 |1,0,0,...,0>
        !
        r(:) = (0.0_DP, 0.0_DP)
        r(1) = (1.0_DP, 0.0_DP)
        !
        CALL ZGTSV(crpa_coul_lanczos_steps_ext, 1, b, a, c, r, &
                   crpa_coul_lanczos_steps_ext, info)
        !
        IF (info /= 0) CALL errore('crpa_calc_coulomb_lanczos', &
           'unable to solve the tridiagonal Lanczos system', 1)
        !
        ! coul_zeta_store(ipert,:,:) has shape (nmels_tot, steps): for each
        ! jmels this is the zdotc(steps, zeta(jmels,:), 1, r, 1) reduction
        ! used (per polarization pair) in crpa_calc_chi_epsilon, done here
        ! for all nmels_tot targets at once. Only the first
        ! crpa_coul_lanczos_steps entries of r enter here (see module-level
        ! comment above).
        !
        ! NOTE the sign here is +2*degspin, not the dielectric case's
        ! -2*degspin (crpa_calc_chi_epsilon.f90): that formula is a direct
        ! resolvent sandwich between two source vectors, chi0_ij =
        ! <d0psi_i|(w-A)^-1|d0psi_j>. Ours instead comes from converting a
        ! real-space density integral (integral drho_bare,mels_j *
        ! dvscf_i(w)) into a bra-ket projection via the bilinear identity
        ! integral B(r) drho_i(r) dr = sum 2Re<psi|B|X_i>, combined with the
        ! d0psi := -P_c^+ (potential)|psi> convention appearing on BOTH the
        ! seed (pert i) and target (mels j) sides - that combination flips
        ! an extra sign relative to the direct-sandwich case (confirmed
        ! empirically: without this flip, every general-q induced term came
        ! out as the exact negative, ratio -1.0000, of the CG reference).
        !
        v_coul_scrd(:, ipert, ifreq) = v_coul_bare(:, ipert) + &
           MATMUL(CONJG(coul_zeta_store(ipert,:,:)), r(1:crpa_coul_lanczos_steps)) &
           * CMPLX(coul_beta_store(ipert,1), 0.0_DP, DP) &
           * CMPLX(2.0_DP*degspin, 0.0_DP, DP)
        !
     ENDDO ! ifreq
     !
  ENDDO ! ipert
  !
  DEALLOCATE(a, b, c, r, b0, c0)
  !
  ! NOTE: crpa_write_coulomb_lanczos(iq) is NOT called here: it must run
  ! after crpa_q0correction() (crpa_main.f90), exactly like the CG path
  ! writes v_coul_bare only after that correction has been applied to it.
  !
  CALL stop_clock('crpa_calc_coul_lncz')
  !
  RETURN
  !
END SUBROUTINE crpa_calc_coulomb_lanczos
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_write_coulomb_lanczos(iq)
  !-----------------------------------------------------------------------
  !! Write filU.bare#iq (v_coul_bare) and filU.scrd#iq (v_coul_scrd), one
  !! line per (mels, pert, frequency): iw1 jw1 R1 iw2 jw2 R2 w0 eta (Ry)
  !! Re(V) Im(V) (eV), for every Wannier-basis run (CG, at one or more
  !! frequencies, and Lanczos).
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode, stdout
  USE constants,  ONLY : RYTOEV
  USE crpacom,    ONLY : filU, v_coul_bare, v_coul_scrd, crpa_nfreq, crpa_freq
  USE crpa_pert,  ONLY : npert_tot, nmels_tot, pert_iwlist, pert_jwlist, pert_Rlist, &
                          mels_iwlist, mels_jwlist, mels_Rlist
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !! index of the q point
  !
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=6) :: int_to_char
  INTEGER :: iun, ios, ipert, jmels, ifreq
  INTEGER :: iw1, jw1, R1(3), iw2, jw2, R2(3)
  !
  IF (.NOT. ionode) RETURN
  !
  WRITE(stdout, '(/5x, a)') "Writing bare Coulomb matrix elements to file"
  !
  filename = TRIM(filU) // ".bare" // TRIM(int_to_char(iq))
  OPEN (NEWUNIT = iun, FILE = TRIM(filename), ERR = 100, IOSTAT = ios)
  !
  WRITE(iun, '(a)') "# iw1 jw1 R1a R1b R1c iw2 jw2 R2a R2b R2c      Re V (eV)      Im V (eV)"
  !
  DO ipert = 1, npert_tot
     !
     iw2 = pert_iwlist(ipert)
     jw2 = pert_jwlist(ipert)
     R2  = pert_Rlist(:, ipert)
     !
     DO jmels = 1, nmels_tot
        !
        iw1 = mels_iwlist(jmels)
        jw1 = mels_jwlist(jmels)
        R1  = mels_Rlist(:, jmels)
        !
        WRITE(iun, '(10I4, 2F20.12)') iw1, jw1, R1, iw2, jw2, R2, &
           v_coul_bare(jmels, ipert) * RYTOEV
        !
     ENDDO
     !
  ENDDO
  !
  CLOSE(iun, STATUS = "keep")
  !
  WRITE(stdout, '(/5x, a)') "Writing frequency-resolved screened Coulomb matrix elements to file"
  !
  filename = TRIM(filU) // ".scrd" // TRIM(int_to_char(iq))
  OPEN (NEWUNIT = iun, FILE = TRIM(filename), ERR = 100, IOSTAT = ios)
  !
  WRITE(iun, '(a)') "# iw1 jw1 R1a R1b R1c iw2 jw2 R2a R2b R2c   w0 (Ry)   eta (Ry)   &
                     &Re V (eV)      Im V (eV)"
  !
  DO ipert = 1, npert_tot
     !
     iw2 = pert_iwlist(ipert)
     jw2 = pert_jwlist(ipert)
     R2  = pert_Rlist(:, ipert)
     !
     DO jmels = 1, nmels_tot
        !
        iw1 = mels_iwlist(jmels)
        jw1 = mels_jwlist(jmels)
        R1  = mels_Rlist(:, jmels)
        !
        DO ifreq = 1, crpa_nfreq
           WRITE(iun, '(10I4,2F12.6,2F20.12)') iw1, jw1, R1, iw2, jw2, R2, &
              DBLE(crpa_freq(ifreq)), AIMAG(crpa_freq(ifreq)), &
              v_coul_scrd(jmels, ipert, ifreq) * RYTOEV
        ENDDO
        !
     ENDDO
     !
  ENDDO
  !
  CLOSE(iun, STATUS = "keep")
  !
  WRITE(stdout, '(5x, a/)') "Done writing frequency-resolved screened Coulomb matrix elements to file"
  !
  RETURN
100 CALL errore('crpa_write_coulomb_lanczos', 'opening file' // TRIM(filename), ABS(ios))
  !
END SUBROUTINE crpa_write_coulomb_lanczos
