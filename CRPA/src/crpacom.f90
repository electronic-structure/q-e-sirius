!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE crpacom
  !
  ! Common variables for the HP program
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  !
  SAVE
  !
  LOGICAL :: write_wan_Rr
  !! If .true. the Wannier functions in real space will be written to folder_wan_Rr.
  !! If .false. the Wannier functions in real space are read from folder_wan_Rr.
  LOGICAL :: recalc_sym
  !! If .true. we recalculate the number of symmetries of the unperturbed lattice due to
  !! the change of the atomic type of one of the atoms
  !            compute_hp,              &     ! If .true. collects all pieces of chi0 and chi
  !            sum_pertq,               &     ! If .true. collects dns0 and dnsscf for all q points
  !                                           ! (for the specific perturbed atom) and computes their
  !                                           ! sum with the phase factor
  !            determine_q_mesh_only,   &     ! If .true. determine the q mesh for a given perturbed atom and exit
  !            determine_num_pert_only, &     ! If .true. determine only which atoms must be perterbed
  LOGICAL :: skip_equivalence_q
  !! If .true. the full frid of q points will be used
  LOGICAL :: crpa_reduce_q
  !! If .true. compute only the irreducible q points (from the crystal
  !! symmetry, recovered even when the SCF/NSCF were run with nosym) and
  !! write their star weights to filU.q. The DFPT itself still runs without
  !! symmetry - the crystal group is used ONLY to classify the q grid.
  LOGICAL :: crpa_w_rspace
  !! If .true. crpa_average_w also writes the Fourier transform over the q
  !! grid, W_ij,kl(R1,R2)(R3), to filU.dir-WR (needs the full q grid)
  !            disable_type_analysis,   &     ! If .true. disable the algorithm which detects whether
  !                                           ! there are atoms of the same type but with different occupations
  !            skip_atom(500)                 ! If .true. no LR calculation will be performed
  !                                           ! for a selected atomic site.
  !                                           ! skip_atom(i), where i runs over atoms. If skip_atom(i)=.true.
  !                                           ! then no linear-response calculation will be performed for the
  !                                           ! i-th atom. This keyword cannot be used when find_atpert=1.
  !                                           ! Warning: Make sure you know what you are doing! This option might
  !                                           ! be useful in several cases:
  !                                           ! - Debugging purposes;
  !                                           ! - You know that the atom which you do not want to perturb
  !                                           !   is equivalent to some other atom (but the code does not recognizes
  !                                           !   this from the symmetry analysis). In this case check that there is
  !                                           !   at least one atom of the same type which was perturbed (this can
  !                                           !   happen only when find_atpert=3), otherwise the post-processing
  !                                           !   calculation of U will fail.
  ! !
  ! LOGICAL, ALLOCATABLE :: todo_atom(:),              & ! Which atoms must be perturbed
  !                         perturbed_atom(:),         & ! Controls which atom is perturbed in the HP
  !                                                      ! calculation
  ! !
  ! INTEGER :: nath,            &             ! Number of (real) atoms in the primitive cell
  !                                           ! with Hubbard_U \= 0
  !            nath_sc,         &             ! Total number of real+virtual Hubbard atoms
  !                                           ! in the virtual supercell
  INTEGER :: nqsh
  !! Number of q points in the grid (without symmetry reduction)
  !! = number of primitive cells in the virtual supercell
  INTEGER :: nah_pert
  !! Site number of the perturbed Hubbard atom
  !            nath_pert,       &             ! Number of actual perturbed Hubbard atoms in the primitive cell
  !            find_atpert,     &             ! Method of searching for atoms to be perturbed
  !            ntyp_new,        &             ! Maximum number of different types detected in the calculation
  !                                           ! (used only when find_atpert=3)
  !            num_neigh,       &             ! Used in the postprocessing: number of nearest neighbors of every atom
  !                                           ! which will be written to the file parameters.out
  !                                           ! (can be used only with lda_plus_u_kind = 2)
  !            lmin,            &             ! Used in the postprocessing: minimum value of the orbital quantum number
  !                                           ! of the Hubbard atoms (Hubbard_l) starting from which (and up to the maximum
  !                                           ! Hubbard_l in the system) Hubbard V will be written to the file parameters.out
  ! !
  ! INTEGER :: equiv_type(ntypx)              ! equiv_type(i)=j, will merge type i to type j
  !                                           ! (useful when nspin=2)
  ! !
  ! CHARACTER(LEN=16)  :: background          ! Background correction
  CHARACTER(LEN=256) :: tmp_dir_save
  !! Temporary directory
  CHARACTER(LEN=256) :: tmp_dir_crpa
  !! Temporary directory
  CHARACTER(LEN=256) :: folder_wan_Rr
  !! Folder with Wannier functions in real space
  CHARACTER(LEN=256) :: wannier_seedname
  !! Seedname for Wannier90 calculation
  CHARACTER(LEN=256) :: filU
  !! Output filename prefix for writing the U tensor.
  !! Bare and screened U tensors for each q are written to filU.bare#iq and filU.scrd#iq.
  CHARACTER(LEN=4)   :: code = 'CRPA'
  !! Name of the code
  !
  REAL(DP) :: conv_thr_chi
  !! Threshold for the calculation of chi
  !             docc_thr,          &          ! Threshold for the comparison of the unperturbed
  !                                           ! occupations (used only with find_atpert=1 for
  !                                           ! determination of atoms which must be perturbed)
  !             rmax,              &          ! Maximum distance (in Bohr) between two atoms
  !                                           ! to search for neighbors (used only at the
  !                                           ! postprocessing step when lda_plus_u_kind = 2).
  ! !
  !
  COMPLEX(DP) :: w_freq
  !! Single frequency (Ry), used when no grid is set. (Default: 0 (static DFPT))
  REAL(DP) :: crpa_freq_start, crpa_freq_end, crpa_freq_broadening
  !! Frequency grid (Ry): w_i = crpa_freq_start + (i-1)*step + i*crpa_freq_broadening
  INTEGER :: crpa_nfreq
  !! Input: number of grid points, 0 (default) for no grid, i.e. w_freq.
  !! After crpa_readin: number of frequencies, 1 for w_freq
  COMPLEX(DP), ALLOCATABLE :: crpa_freq(:)
  !! The frequencies of every path (CG, Lanczos, dielectric), crpa_freq(crpa_nfreq)
  COMPLEX(DP) :: current_w
  !! CG path: frequency of the Sternheimer solve being done (one of crpa_freq)
  !
  !
  ! Definition of perturbation basis (for CRPA)
  ! NOTE: For constrained DFPT active space, see LR_Modules/constrained_dfpt.f90
  !
  INTEGER :: active_bands_max
  !! For pert_basis == 'bands': Maximum band index for perturbations
  INTEGER :: active_bands_min
  !! For pert_basis == 'bands': Minimum band index for perturbations
  !
  REAL(DP), ALLOCATABLE :: ns(:)
  !! Trace of unperturbed occupations (spin up + spin down)
  REAL(DP), ALLOCATABLE :: magn(:)
  !! Unperturbed magnetization
  REAL(DP), ALLOCATABLE :: Rvect(:, :)
  !! Radius-vector of the primitive cell
  REAL(DP), ALLOCATABLE :: chi0(:,:)
  !! Bare response function (from 1st iteration)
  REAL(DP), ALLOCATABLE :: chi(:,:)
  !! SCF response function
  ! !
  COMPLEX(DP), ALLOCATABLE :: dns0(:,:,:,:,:)
  !! Bare response occupation matrix (from 1st iteration)
  COMPLEX(DP), ALLOCATABLE :: dns0_tot(:,:,:,:,:)
  !! Total bare response occupation matrix (summed over q)
  COMPLEX(DP), ALLOCATABLE :: dnsscf_tot(:,:,:,:,:)
  !! Total SCF  response occupation matrix (summed over q)
  COMPLEX(DP), ALLOCATABLE :: trace_dns_tot_old(:)
  !! Trace of the response occupation matrix (for a convergence test)
  ! !
  INTEGER, ALLOCATABLE :: ityp_new(:)
  !! Types of atoms
  !
  ! Variables related to the k point mesh
  !
  INTEGER :: nks_orig
  ! Size of the global k point mesh for the underlying SCF calculation
  REAL(DP), ALLOCATABLE :: xk_orig(:, :)
  ! Original global k point mesh for the underlying SCF calculation. (3, nks_orig)
  INTEGER, ALLOCATABLE :: ik_to_ik_orig(:)
  ! Index of current k and k+q points modulo G in the xk_orig.
  ! xk(:, ik) = xk_orig(:, ik_to_ik_orig(ik)) + G for ik = 1, ..., nks. (nks = 2 * nksq)
  ! Here, xk(:, ik) includes both k and k+q points.
  !
  ! Coulumb matrix elements
  COMPLEX(DP), ALLOCATABLE :: v_coul_bare(:, :)
  !! Bare coulomb matrix elements. (ipert, jpert).
  COMPLEX(DP), ALLOCATABLE :: v_coul_scrd(:, :, :)
  !! Screened coulomb matrix elements. (ipert, jpert, ifreq), ifreq over crpa_freq.
  !!
  COMPLEX(DP), ALLOCATABLE :: dvbare(:, :)
  !! Bare perturbing potential on the DENSE grid (dfftp%nnr, nspin_mag).
  !! Use this for dv_of_drho() and newdq(), which both expect dfftp.
  COMPLEX(DP), ALLOCATABLE :: dvbares(:, :)
  !! Same potential interpolated onto the SMOOTH grid (dffts%nnr, nspin_mag).
  !! Use this for apply_dpot_bands(), whose dummy is dv(dffts%nnr, nspin_mag).
  !! Identical to dvbare when doublegrid = .FALSE. (norm-conserving, ecutrho =
  !! 4*ecutwfc), which is why passing dvbare there worked for every NC test and
  !! silently read the wrong array for USPP (dense 45^3 vs smooth 32^3 on the
  !! Si USPP test).
  !! Bare perturbing potential
  !
  INTEGER :: iuwan
  !! Unit for reading and writing the Wannier functions
  INTEGER :: lrwan
  !! Record length of the Wannier function buffer
  !
  REAL(DP) :: dist_thr
  !! Cutoff distance for the Wannier functions in bohr units. For perturbation.
  REAL(DP) :: dist_thr_large
  !! Larger cutoff distance for the Wannier functions in bohr units. For matrix element calculation.
  !
  CHARACTER(32) :: q0div_treatment
  ! Correction scheme for the q->0 divergence of the Coulomb kernel
  CHARACTER(32) :: crpa_mode
  ! Calculation mode of the Uijkl(R) elements
  CHARACTER(LEN=4) :: crpa_xcmode = 'TCTC'
  ! Kernel carried by the external legs of the screened Coulomb interaction W.
  ! The screening kernel itself is selected by lrpa.
  !   'TCTC' (default): bare Coulomb v on both external legs,
  !                     W_ij = <rho_j|v|rho_i> + <rho_j| v chi_r v |rho_i>.
  !   'TETE'          : v+f_xc on both external legs and in the screening
  !                     Requires lrpa=.false.
  LOGICAL :: crpa_epsil
  ! Logical for the evaluation of the dielectric tensor
  LOGICAL :: crpa_done_epsil
  ! Logical that checks if the dielectric tensor calculation is finished
  LOGICAL :: crpa_epsil_only
  ! If true, only the dielectric tensor is computed; the q-loop is skipped
  COMPLEX(DP) :: epsilon_w(3,3)
  ! Complex valued macroscopic dielectric tensor at finite frequency

  REAL(DP), PARAMETER :: w_drude_min = 2.0d-5
  ! Floor applied to Re(w) in the Drude denominator (crpa_dielec.f90 and
  ! crpa_calc_chi_epsilon.f90). Matches RESPACK's own floor (chiqw_sub.F:
  ! IF (w < 1.0d-5 [Ha]) w = 1.0d-5), converted to this code's Ry convention
  ! (1 Ha = 2 Ry) so the two floors coincide in eV rather than differing by
  ! a factor of 2 as before.
  REAL(DP) :: plasma_freq2(3,3)
  ! Intraband (Drude) plasma-frequency tensor omega_p^2 (Ry^2), computed by
  ! crpa_plasma_frequency() for metals (lgauss .or. ltetra) at q=Gamma. Built
  ! in the same crystal-axis / fpi*e2/omega convention as epsilon_w and
  ! epsilon_w_diel, so that both CG and Lanczos heads can add
  ! -plasma_freq2(ipol,ipol) / (w * (w + i*eta)) directly. Zero (and unused)
  ! for insulators.
  LOGICAL :: crpa_metal
  ! .true. if lgauss .or. ltetra (set once in crpa_plasma_frequency): whether
  ! the system is a metal (nonzero plasma_freq2 is possible).

  LOGICAL :: crpa_diel_intraband
  ! Namelist-controlled switch (default .true.): whether the intraband/Drude
  ! term (plasma_freq2) is actually added to the head of epsilon_w(_diel).
  ! Sets crpa_metal .and. crpa_diel_intraband is the real "add the Drude
  ! term" condition in crpa_dielec.f90/crpa_calc_chi_epsilon.f90; kept
  ! separate from crpa_metal itself so a metal can be run with the term
  ! deliberately switched off, e.g. to isolate the interband-only result.

  LOGICAL :: crpa_diel_flat
  ! Namelist switch (default .false.), metals with lgauss only: build the
  ! electric-field auxiliary function (dvpsi_e buffer) with the occupation-free
  ! projector Pbar_n (control_lr%lflat_proj in LR_Modules/orthogonalize.f90),
  ! so that the dielectric response carries a single metallic weight W instead
  ! of W^3, and use the block-trace plasma frequency (crpa_plasma_frequency.f90)
  ! over the same degenerate blocks. See Notes/metal_efield_sternheimer.tex.

  REAL(DP) :: crpa_diel_deg_thr
  ! Degeneracy threshold (Ry, default 1.0d-5) of crpa_diel_flat: pairs with
  ! |e_m - e_n| < crpa_diel_deg_thr are removed from the interband part and
  ! counted in the plasma frequency instead.

  LOGICAL :: crpa_diel_llanczos
  ! If .true. epsilon_w(:,:,:) is computed by building a non-Hermitian
  ! Lanczos chain (see module crpa_diel_lanczos) instead of solving the
  ! Sternheimer equation self-consistently (crpa_solve_e_w + dfpt_kernel)
  ! at every frequency. It evaluates the complex dielectric tensor at all the
  ! frequencies crpa_freq from a single Lanczos chain per polarization.
  INTEGER :: crpa_diel_lanczos_steps
  ! Number of steps of the Lanczos chain used to compute epsilon_w_diel
  INTEGER :: crpa_diel_lanczos_restart_step
  ! Interval (in Lanczos steps) for writing intermediate restart points. 0 means
  ! none (a final checkpoint is still written, see crpa_diel_lanczos_save_final).
  LOGICAL :: crpa_diel_lanczos_restart
  ! .true.: resume/extend the dielectric chain from the checkpoint in tmp_dir_ph
  ! (crpa_restart_diel_lanczos.f90). .false.: start fresh, ignore any checkpoint.
  LOGICAL :: crpa_diel_lanczos_save_final
  ! .true.: write a checkpoint at the end of the chain, so that it can be
  ! extended later with crpa_diel_lanczos_restart = .true.
  COMPLEX(DP), ALLOCATABLE :: epsilon_w_diel(:,:,:)
  ! Complex macroscopic dielectric tensor epsilon_w_diel(3,3,crpa_nfreq)
  ! at the frequencies crpa_freq.

  CHARACTER(LEN=16) :: crpa_diel_extrapolation
  ! 'no' (default): use the beta/gamma chain as computed, length
  ! crpa_diel_lanczos_steps, in the continued-fraction evaluation.
  ! 'constant'/'oscillate': extrapolate beta/gamma (assumed to have settled
  ! into a stable, possibly two-periodic pattern) out to
  ! crpa_diel_lanczos_steps_ext steps before the continued-fraction solve,
  ! as in thermo_pw's calc_chi.f90/extrapolate(). 'constant' extrapolates to
  ! the average value; 'oscillate' keeps the average oscillation amplitude.
  INTEGER :: crpa_diel_lanczos_steps_ext
  ! Total chain length (>= crpa_diel_lanczos_steps) used in the
  ! continued-fraction solve when crpa_diel_extrapolation /= 'no'.

  LOGICAL :: crpa_llanczos
  ! If .true. the screened Coulomb matrix elements (crpa_coul_lanczos module)
  ! are computed by building one non-Hermitian Lanczos chain per bare
  ! perturbation (see crpa_do_lanczos_coulomb.f90), instead of solving the
  ! Sternheimer equation self-consistently for every perturbation
  ! (crpa_solve_linear_system.f90). Unlike that CG path, this evaluates the
  ! screened Coulomb matrix at all the frequencies crpa_freq, from a single
  ! chain per bare perturbation.
  ! This path is currently only implemented for pert_basis = 'wannier'.
  INTEGER :: crpa_coul_lanczos_steps
  ! Number of steps of the Lanczos chain used to compute the screened
  ! Coulomb matrix elements when crpa_llanczos = .true.

  CHARACTER(LEN=16) :: crpa_coul_extrapolation
  ! Same extrapolation scheme as crpa_diel_extrapolation ('no', 'constant',
  ! 'oscillate'), applied independently to each of the npert_tot screened
  ! Coulomb Lanczos chains: see crpa_extrapolate_coul.f90.
  INTEGER :: crpa_coul_lanczos_steps_ext
  ! Total chain length (>= crpa_coul_lanczos_steps) used in the
  ! continued-fraction solve when crpa_coul_extrapolation /= 'no'.
  LOGICAL :: crpa_coul_restart
  ! .true.: reuse the screened-Coulomb results of a previous run
  ! (crpa_restart_coul_lanczos.f90, crpa_restart_coul_cg.f90)
  LOGICAL :: crpa_coul_lanczos_save_final
  ! .true.: write the Lanczos vectors at the end of every chain, so that the
  ! chains can be extended later with crpa_coul_restart = .true.
  INTEGER :: crpa_coul_lanczos_restart_step
  ! Interval (in Lanczos steps) for writing intermediate vector checkpoints of
  ! the screened-Coulomb chains. 0 means none.

  LOGICAL :: crpa_spectrum_only
  ! If .true., crpa_main.f90 skips SCF/NSCF/Wannier setup and the q-point
  ! loop entirely, and instead reads back the RAW (pre-extrapolation)
  ! Lanczos chain coefficients previously saved by crpa_write_diel_lanczos_chain
  ! (diel_lanczos_chain.dat, if crpa_epsil) and/or crpa_write_coul_lanczos_chain
  ! (filU.coul_lanczos_chain#iq, if crpa_llanczos) - see
  ! crpa_run_spectrum_only.f90 - and just redoes the frequency-dependent
  ! continued-fraction solve, at whatever (possibly different)
  ! frequencies (w_freq or crpa_freq_*) and extrapolation scheme are
  ! given in the current input. Requires crpa_epsil and/or crpa_llanczos.

END MODULE crpacom


MODULE crpa_pert
  ! Variables for defining the perturbations
  USE kinds, ONLY : DP
  !
  CHARACTER(LEN=256) :: pert_basis
  !! = 'bands' : use product of bands with indices from active_bands_min to active_bands_max
  !! = 'wannier' : use product of Wannier functions
  !
  INTEGER :: npert_tot
  !! Total number of perturbations to be computed for the given q point
  INTEGER :: nmels_tot
  !! Total number of basis states to calculate matrix elements for the given q point
  !
  ! Perturbations defined by the band index pairs
  INTEGER :: pert_nbnd
  !! Total number of bands to be considered for perturbation.
  !! pert_nbnd = pert_ibnd_max - pert_ibnd_min + 1
  !! npert_tot = pert_nbnd * pert_nbnd * nksqtot
  INTEGER :: pert_ibnd_min
  !! Minimum band index to be considered for perturbation
  INTEGER :: pert_ibnd_max
  !! Maximum band index to be considered for perturbation
  !
  INTEGER :: pert_nk
  ! If pert_basis == 'bands': Number of k points in this core (nksq)
  ! If pert_basis == 'bands': Number of total R vectors to perturb
  !
  INTEGER, ALLOCATABLE :: pert_iwlist(:)
  ! If pert_basis == 'wannier': List of WF index i to perturb (R = 0)
  INTEGER, ALLOCATABLE :: pert_jwlist(:)
  ! If pert_basis == 'wannier': List of WF index j to perturb (R = pert_Rlist)
  INTEGER, ALLOCATABLE :: pert_Rlist(:, :)
  ! If pert_basis == 'wannier': List of R vectors to perturb in crystal coordinates
  !
  INTEGER, ALLOCATABLE :: mels_iwlist(:)
  ! If pert_basis == 'wannier': List of WF index i to calculate matrix elements (R = 0)
  INTEGER, ALLOCATABLE :: mels_jwlist(:)
  ! If pert_basis == 'wannier': List of WF index j to calculate matrix elements (R = pert_Rlist)
  INTEGER, ALLOCATABLE :: mels_Rlist(:, :)
  ! If pert_basis == 'wannier': List of R vectors to calculate matrix elements in crystal coordinates
  !
END MODULE crpa_pert

MODULE crpa_qpoints
  LOGICAL :: qplot
  !! if TRUE the q are read from input
END MODULE crpa_qpoints
