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
  !! Frequency value for finite-frequency DFPT. (Default: 0 (static DFPT))
  !
  !
  ! Definition of active space
  !
  CHARACTER(LEN=256) :: active_space
  !! Method of defining the active space
  !! = 'bands' : use band indices from active_bands_min to active_bands_max
  INTEGER :: active_bands_max
  !! For active_space == 'bands': Maximum band index in the active space
  INTEGER :: active_bands_min
  !! For active_space == 'bands': Minimum band index in the active space
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
  REAL(DP) :: dist_thr
  !! Cutoff distance for the Wannier functions in bohr units. For perturbation.
  REAL(DP) :: dist_thr_large
  !! Larger cutoff distance for the Wannier functions in bohr units. For matrix element calculation.
  !
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
