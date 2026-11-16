!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_run_spectrum_only()
  !-----------------------------------------------------------------------
  !! Driver for crpa_spectrum_only = .true.: reconstructs epsilon(w) and/or
  !! the screened Coulomb matrix elements from previously saved, RAW
  !! (pre-extrapolation) Lanczos chain coefficients, at whatever frequency
  !! grid and extrapolation scheme are given in the CURRENT input - without
  !! rerunning SCF/NSCF/Wannier90 or any DFPT at all. Called from
  !! crpa_main.f90 in place of crpa_generate_grids/crpa_setup_wannier/the
  !! q-point loop.
  !!
  !! Mirrors turbo_spectrum.x's role (TDDFPT/src/turbo_spectrum.f90):
  !! a saved chain plus a cheap continued-fraction solve is all that is
  !! needed to change the frequency grid/broadening/extrapolation scheme
  !! after the fact.
  !!
  !! The q=0 divergence correction (crpa_q0correction) IS reapplied here for
  !! the q=Gamma chain, exactly as in the normal flow (crpa_main.f90:
  !! calc_coulomb_lanczos -> q0correction -> write_coulomb_lanczos): this is
  !! safe without any live DFPT because crpa_q0correction reads epsilon_w.dat
  !! from disk if present, and the dielectric-chain reconstruction just above
  !! (if crpa_epsil) has already written a fresh copy of that file, at the
  !! CURRENT (possibly new) frequency grid. If crpa_epsil = .false., it falls
  !! back to its normal bare-only correction, also with no DFPT involved.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE crpacom,    ONLY : crpa_epsil, crpa_llanczos, v_coul_bare, v_coul_scrd
  USE crpa_pert,  ONLY : pert_iwlist, pert_jwlist, pert_Rlist, &
                          mels_iwlist, mels_jwlist, mels_Rlist
  USE crpa_coul_lanczos, ONLY : coul_beta_store, coul_gamma_store, coul_zeta_store, &
                                 coul_beta_store_ext, coul_gamma_store_ext
  USE qpoint,     ONLY : xq
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  INTEGER :: iq, nfile
  !
  WRITE(stdout, '(/,5X,"=-------------------------------------------------------------=")')
  WRITE(stdout, '(/,5X,"crpa_spectrum_only: reconstructing from saved Lanczos chains")')
  WRITE(stdout, '(/,5X,"=-------------------------------------------------------------=")')
  !
  IF (crpa_epsil) THEN
     !
     CALL crpa_read_diel_lanczos_chain()
     CALL crpa_calc_chi_epsilon()
     !
  ENDIF
  !
  IF (crpa_llanczos) THEN
     !
     iq = 1
     DO
        CALL crpa_read_coul_lanczos_chain(iq, exst, 0, nfile)
        IF (.NOT. exst) EXIT
        !
        CALL crpa_calc_coulomb_lanczos(iq)
        IF (ALL(xq(:) == 0.0_DP)) CALL crpa_q0correction()
        CALL crpa_write_coulomb_lanczos(iq)
        !
        ! Targeted cleanup of only what crpa_read_coul_lanczos_chain
        ! allocated (not crpa_deallocate_coul_lanczos: that routine also
        ! handles coul_evc0/coul_sevq0/coul_d0psi_mels/coul_evc1* etc.,
        ! none of which are touched in this mode, since no wavefunctions
        ! are ever read back).
        !
        DEALLOCATE(pert_iwlist, pert_jwlist, pert_Rlist)
        DEALLOCATE(mels_iwlist, mels_jwlist, mels_Rlist)
        DEALLOCATE(v_coul_bare, v_coul_scrd)
        DEALLOCATE(coul_beta_store, coul_gamma_store, coul_zeta_store)
        DEALLOCATE(coul_beta_store_ext, coul_gamma_store_ext)
        !
        iq = iq + 1
        !
     ENDDO
     !
     IF (iq == 1) CALL errore('crpa_run_spectrum_only', &
        'crpa_llanczos = .true. but no filU.coul_lanczos_chain#iq file was found', 1)
     !
     CALL crpa_average_w()
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE crpa_run_spectrum_only
