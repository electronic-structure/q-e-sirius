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
  ! Common variables for the CRPA calculation (hp.x with method = 'CRPA').
  !
  ! Some variables here (ns, magn, chi0, chi, dns0, ...) have counterparts in
  ! ldaU_hp, but they cannot be shared: in CRPA the perturbations are products
  ! of bands or Wannier functions, not Hubbard-atom occupations, so the arrays
  ! have different shapes and lifecycles than the LRT (Hubbard) ones.
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
  LOGICAL :: skip_equivalence_q
  !! If .true. the full grid of q points will be used
  !
  INTEGER :: nqsh
  !! Number of q points in the grid (without symmetry reduction)
  !! = number of primitive cells in the virtual supercell
  INTEGER :: nah_pert
  !! Site number of the perturbed Hubbard atom
  !
  CHARACTER(LEN=256) :: folder_wan_Rr
  !! Folder with Wannier functions in real space
  CHARACTER(LEN=256) :: wannier_seedname
  !! Seedname for Wannier90 calculation
  !
  REAL(DP) :: conv_thr_chi
  !! Threshold for the calculation of chi
  !
  COMPLEX(DP) :: w_freq
  !! Frequency value for finite-frequency DFPT. (Default: 0 (static DFPT))
  !
  ! Definition of perturbation basis (for CRPA)
  ! NOTE: For constrained DFPT active space, see LR_Modules/constrained_dfpt.f90
  !
  INTEGER :: pert_bands_max
  !! For pert_basis == 'bands': Maximum band index for perturbations
  INTEGER :: pert_bands_min
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
  !
  COMPLEX(DP), ALLOCATABLE :: dns0(:,:,:,:,:)
  !! Bare response occupation matrix (from 1st iteration)
  COMPLEX(DP), ALLOCATABLE :: dns0_tot(:,:,:,:,:)
  !! Total bare response occupation matrix (summed over q)
  COMPLEX(DP), ALLOCATABLE :: dnsscf_tot(:,:,:,:,:)
  !! Total SCF  response occupation matrix (summed over q)
  COMPLEX(DP), ALLOCATABLE :: trace_dns_tot_old(:)
  !! Trace of the response occupation matrix (for a convergence test)
  !
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
  COMPLEX(DP), ALLOCATABLE :: v_coul_scrd(:, :)
  !! Screened coulomb matrix elements. (ipert, jpert).
  !!
  COMPLEX(DP), ALLOCATABLE :: dvbare(:, :)
  !! Bare perturbing potential
  !
  INTEGER :: iuwan
  !! Unit for reading and writing the Wannier functions
  INTEGER :: lrwan
  !! Record length of the Wannier function buffer
  !
  REAL(DP) :: wan_dist_thr
  !! Cutoff distance for the Wannier functions in bohr units. For perturbation.
  REAL(DP) :: wan_dist_thr_mels
  !! Larger cutoff distance for the Wannier functions in bohr units. For matrix element calculation.
  !
END MODULE crpacom


MODULE crpa_pert
  ! Variables for defining the perturbations
  USE kinds, ONLY : DP
  !
  CHARACTER(LEN=256) :: pert_basis
  !! = 'bands' : use product of bands with indices from pert_bands_min to pert_bands_max
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
