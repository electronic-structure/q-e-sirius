!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE crpa_diel_lanczos
  !----------------------------------------------------------------------------
  !! Variables for the non-Hermitian Lanczos-chain evaluation of the complex
  !! macroscopic dielectric tensor epsilon_w_diel(3,3,nfreq) at q=Gamma
  !! (electric-field response), at the complex frequencies crpacom%crpa_freq.
  !!
  !! This is an alternative to the self-consistent Sternheimer solve
  !! (crpa_solve_e_w + dfpt_kernel) used when crpa_diel_llanczos = .TRUE.
  !! The chain itself (alpha/beta/gamma/zeta below) does not depend on
  !! frequency: only the final continued-fraction evaluation in
  !! crpa_calc_chi_epsilon.f90 does, so a single chain per polarization
  !! gives the dielectric function at arbitrarily many frequencies cheaply.
  !!
  !! NOTE: This module is specific to the dielectric-function calculation.
  !! It is unrelated to the (future) Lanczos chain that will replace the
  !! self-consistent solve in crpa_solve_linear_system.f90 for the general
  !! cRPA screened-U calculation (different perturbations, different q).
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  ! NOTE: the control flags/parameters for this calculation (whether to use
  ! the Lanczos path at all, number of steps, restart interval, frequency
  ! grid) are read from the input namelist and live in crpacom, not here:
  ! this module only holds the chain's working arrays.
  !
  INTEGER :: diel_size_evc1
  !! Size of the diel_evc1 (and _new, _old) arrays, used in BLAS calls
  !
  COMPLEX(DP), ALLOCATABLE, TARGET :: diel_evc0(:,:,:)
  !! Ground-state wavefunctions at k, cached for the duration of the chain.
  !! Dimensions: (npwx*npol, nbnd, nksq)
  COMPLEX(DP), POINTER :: diel_sevq0(:,:,:)
  !! S * diel_evc0 (points to diel_evc0 itself when okvan = .FALSE.)
  COMPLEX(DP), ALLOCATABLE :: diel_d0psi(:,:,:,:)
  !! Starting vectors (chain SEED) P_c^+ V_ext|psi>, one per Cartesian
  !! polarization. Dimensions: (npwx*npol, nbnd, nksq, 3)
  COMPLEX(DP), ALLOCATABLE :: diel_d0psi_left(:,:,:,:)
  !! LEFT vectors for the zeta projections, one per Cartesian polarization:
  !! the raw dvpsi_e() output, WITHOUT the extra orthogonalize() that
  !! diel_d0psi carries. Same dimensions as diel_d0psi.
  !! For an insulator these two are identical (P_c^+ is idempotent there),
  !! so this array only matters for metals - see the long comment in
  !! crpa_do_lanczos_epsilon.f90 where it is filled.
  !
  ! diel_evc1, diel_evc1_new, diel_evc1_old hold ONE non-Hermitian Lanczos
  ! chain spanning the 3 Cartesian polarizations, advanced in LOCK-STEP, one
  ! Lanczos iteration at a time for all 3 together (as in thermo_pw's
  ! do_lanczos.f90). Lock-step is required because the induced density is
  ! symmetrized every iteration via crpa_psym_diel.f90, which mixes the 3
  ! Cartesian components under the crystal's point group and therefore needs
  ! all 3 to be simultaneously up to date. That mixing also means the 3 blocks
  ! are NOT independent chains: they form one vector of the direct-sum space
  ! and carry a single normalization (see crpa_nh_lanczos_step, at the end of
  ! crpa_do_lanczos_epsilon.f90).
  !
  ! Index order is (npwx*npol, nbnd, nksq, 3, 2): the polarization dimension
  ! sits INSIDE the Krylov-component one, so the 3 blocks of a component are
  ! contiguous and span nksq*3 k-blocks in memory. That is the layout
  ! thermo_pw uses (evc1(npwx*npol,nbnd,nksq*rpert,2)) and the one a
  ! block-aware LR_Modules/lanczos_nonhermitian.f90 would expect, so the chain
  ! could be handed to such a routine with no data movement. The size-2
  ! dimension holds the two coupled Krylov vectors of the bi-orthogonal
  ! (non-Hermitian) algorithm - "q" and "p" in the notation of
  ! LR_Modules/lanczos_nonhermitian.f90 - needed because a genuinely complex
  ! frequency breaks the pseudo-Hermitian shortcut used when this driver only
  ! had to support w_freq = 0.
  !
  COMPLEX(DP), ALLOCATABLE :: diel_evc1(:,:,:,:,:)
  !! Current Lanczos vectors (q_j, p_j). Dimensions: (npwx*npol, nbnd, nksq, 3, 2)
  COMPLEX(DP), ALLOCATABLE :: diel_evc1_new(:,:,:,:,:)
  !! On input to the Lanczos step: (A q_j, A^T p_j). On output: (q_j+1, p_j+1)
  COMPLEX(DP), ALLOCATABLE :: diel_evc1_old(:,:,:,:,:)
  !! Previous Lanczos vectors (q_j-1, p_j-1)
  COMPLEX(DP), ALLOCATABLE :: diel_sevc1(:,:,:,:)
  !! S * diel_evc1(:,:,:,2,:) (p_j), used only when okvan = .TRUE., otherwise a copy
  !
  REAL(DP), ALLOCATABLE :: diel_beta_store(:,:)
  !! Sub-diagonal coefficients of the tridiagonal Lanczos matrix.
  !! Dimensions: (3, crpa_diel_lanczos_steps). There is a SINGLE chain over
  !! the 3 polarizations, so beta is one scalar per iteration: the same value
  !! is stored in all 3 rows, purely so that the per-ipol tridiagonal solve in
  !! crpa_calc_chi_epsilon.f90 reduces to the single joint solve.
  REAL(DP), ALLOCATABLE :: diel_gamma_store(:,:)
  !! Super-diagonal coefficients of the tridiagonal Lanczos matrix (in
  !! general /= diel_beta_store for the non-Hermitian algorithm).
  !! Dimensions: (3, crpa_diel_lanczos_steps)
  COMPLEX(DP), ALLOCATABLE :: diel_zeta_store(:,:,:)
  !! zeta_store(ipol,jpol,iter) = <d0psi(jpol)|q_j> projection of the
  !! current Lanczos vector of the ipol-th chain onto the jpol-th starting
  !! vector. Dimensions: (3, 3, crpa_diel_lanczos_steps)
  !
  REAL(DP), ALLOCATABLE :: diel_beta_store_ext(:,:)
  !! beta_store, extended out to crpa_diel_lanczos_steps_ext by
  !! crpa_extrapolate_diel() when crpa_diel_extrapolation /= 'no' (otherwise
  !! just a copy of diel_beta_store). Used in place of diel_beta_store in
  !! the continued-fraction solve (crpa_calc_chi_epsilon.f90).
  !! Dimensions: (3, crpa_diel_lanczos_steps_ext)
  REAL(DP), ALLOCATABLE :: diel_gamma_store_ext(:,:)
  !! gamma_store, extended (see diel_beta_store_ext above).
  !
END MODULE crpa_diel_lanczos
