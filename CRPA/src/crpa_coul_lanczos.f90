!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE crpa_coul_lanczos
  !----------------------------------------------------------------------------
  !! Variables for the non-Hermitian Lanczos-chain evaluation of the screened
  !! Coulomb matrix elements (pert_basis = 'wannier' only), at the complex
  !! frequencies crpacom%crpa_freq, as the crpa_diel_llanczos dielectric chain.
  !!
  !! This is an alternative to the self-consistent Sternheimer solve
  !! (crpa_solve_linear_system.f90 + dfpt_kernel) used when
  !! crpa_llanczos = .TRUE. Unlike that CG path, which solves one
  !! self-consistent linear system per perturbation at a single frequency
  !! (w_freq), this evaluates the screened Coulomb matrix over the whole
  !! frequency grid from a single Lanczos chain per bare perturbation
  !! (npert_tot chains, run sequentially - NOT lock-stepped, since distinct
  !! CRPA perturbations are not mixed into each other by any symmetry
  !! operation exploited by this code, unlike the 3 Cartesian polarizations
  !! of the dielectric chain; see crpa_set_upert.f90's trivial
  !! representation). No per-iteration density symmetrization either: the
  !! calculation is assumed to run without exploiting crystal symmetry.
  !!
  !! The result goes into crpacom%v_coul_scrd(:,:,ifreq), as for the CG path.
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE, TARGET :: coul_evc0(:,:,:)
  !! Ground-state wavefunctions at k, cached for the duration of the
  !! calculation at this q. Dimensions: (npwx*npol, nbnd, nksq)
  COMPLEX(DP), ALLOCATABLE, TARGET :: coul_evq0(:,:,:)
  !! Ground-state wavefunctions at k+q, cached likewise.
  !! Dimensions: (npwx*npol, nbnd, nksq)
  COMPLEX(DP), POINTER :: coul_sevq0(:,:,:)
  !! S * coul_evq0 (points to coul_evq0 itself when okvan = .FALSE.)
  !
  COMPLEX(DP), ALLOCATABLE :: coul_d0psi_mels(:,:,:,:)
  !! Starting vectors -P_c^+ dV_bare,mels_j |psi_k>, one per entry of the
  !! matrix-element (mels_*) list, built once per q (shared by every bare
  !! perturbation's chain). Dimensions: (npwx*npol, nbnd, nksq, nmels_tot)
  !
  ! coul_evc1, coul_evc1_new, coul_evc1_old hold ONE non-Hermitian Lanczos
  ! chain at a time: unlike the dielectric module (3 polarizations advanced
  ! in lock-step), the npert_tot bare-perturbation chains here are
  ! independent and run sequentially, so these arrays are reset and reused
  ! for each perturbation in turn. The size-2 dimension holds the two
  ! coupled Krylov vectors of the bi-orthogonal algorithm ("q" and "p" in
  ! the notation of LR_Modules/lanczos_nonhermitian.f90).
  !
  COMPLEX(DP), ALLOCATABLE :: coul_evc1(:,:,:,:)
  !! Current Lanczos vectors (q_j, p_j). Dimensions: (npwx*npol, nbnd, nksq, 2)
  COMPLEX(DP), ALLOCATABLE :: coul_evc1_new(:,:,:,:)
  !! On input to the Lanczos step: (A q_j, A^T p_j). On output: (q_j+1, p_j+1)
  COMPLEX(DP), ALLOCATABLE :: coul_evc1_old(:,:,:,:)
  !! Previous Lanczos vectors (q_j-1, p_j-1)
  COMPLEX(DP), ALLOCATABLE :: coul_sevc1(:,:,:)
  !! S * coul_evc1(:,:,:,2) (p_j), used only when okvan = .TRUE., otherwise a copy
  !
  REAL(DP), ALLOCATABLE :: coul_beta_store(:,:)
  !! Sub-diagonal coefficients of the tridiagonal Lanczos matrix.
  !! Dimensions: (npert_tot, crpa_coul_lanczos_steps)
  REAL(DP), ALLOCATABLE :: coul_gamma_store(:,:)
  !! Super-diagonal coefficients of the tridiagonal Lanczos matrix.
  !! Dimensions: (npert_tot, crpa_coul_lanczos_steps)
  COMPLEX(DP), ALLOCATABLE :: coul_zeta_store(:,:,:)
  !! zeta_store(ipert,jmels,iter) = <d0psi_mels(jmels)|q_iter> projection of
  !! the current Lanczos vector of the ipert-th chain onto the jmels-th
  !! matrix-element target vector. Dimensions: (npert_tot, nmels_tot,
  !! crpa_coul_lanczos_steps)
  !
  REAL(DP), ALLOCATABLE :: coul_beta_store_ext(:,:)
  !! coul_beta_store, extended out to crpa_coul_lanczos_steps_ext by
  !! crpa_extrapolate_coul() when crpa_coul_extrapolation /= 'no' (otherwise
  !! just a copy of coul_beta_store). Used in place of coul_beta_store in
  !! the continued-fraction solve (crpa_calc_coulomb_lanczos.f90).
  !! Dimensions: (npert_tot, crpa_coul_lanczos_steps_ext)
  REAL(DP), ALLOCATABLE :: coul_gamma_store_ext(:,:)
  !! coul_gamma_store, extended (see coul_beta_store_ext above).
  !
END MODULE crpa_coul_lanczos
