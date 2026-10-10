!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE elphcom
  !
  ! Common variables for the HP program
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  !
  SAVE
  !
  ! ! LOGICAL :: skip_type(ntypx),        &     ! If .true. skip the calculation for a specific type
  ! !                                           ! (e.g. Ni_up can be used for Ni_down with some spin
  ! !                                           ! considerations)
  ! !            perturb_only_atom(500),  &     ! If perturb_only_atom(i)=.true. perterb only i-th atom
  ! LOGICAL :: recalc_sym
  ! !! If .true. we recalculate the number of symmetries of the unperturbed lattice due to
  ! !! the change of the atomic type of one of the atoms
  ! !            compute_hp,              &     ! If .true. collects all pieces of chi0 and chi
  ! !            sum_pertq,               &     ! If .true. collects dns0 and dnsscf for all q points
  ! !                                           ! (for the specific perturbed atom) and computes their
  ! !                                           ! sum with the phase factor
  ! !            determine_q_mesh_only,   &     ! If .true. determine the q mesh for a given perturbed atom and exit
  ! !            determine_num_pert_only, &     ! If .true. determine only which atoms must be perterbed
  ! LOGICAL :: skip_equivalence_q
  ! !! If .true. the full frid of q points will be used
  ! !            disable_type_analysis,   &     ! If .true. disable the algorithm which detects whether
  ! !                                           ! there are atoms of the same type but with different occupations
  ! !            skip_atom(500)                 ! If .true. no LR calculation will be performed
  ! !                                           ! for a selected atomic site.
  ! !                                           ! skip_atom(i), where i runs over atoms. If skip_atom(i)=.true.
  ! !                                           ! then no linear-response calculation will be performed for the
  ! !                                           ! i-th atom. This keyword cannot be used when find_atpert=1.
  ! !                                           ! Warning: Make sure you know what you are doing! This option might
  ! !                                           ! be useful in several cases:
  ! !                                           ! - Debugging purposes;
  ! !                                           ! - You know that the atom which you do not want to perturb
  ! !                                           !   is equivalent to some other atom (but the code does not recognizes
  ! !                                           !   this from the symmetry analysis). In this case check that there is
  ! !                                           !   at least one atom of the same type which was perturbed (this can
  ! !                                           !   happen only when find_atpert=3), otherwise the post-processing
  ! !                                           !   calculation of U will fail.
  ! ! !
  ! ! LOGICAL, ALLOCATABLE :: todo_atom(:),              & ! Which atoms must be perturbed
  ! !                         perturbed_atom(:),         & ! Controls which atom is perturbed in the HP
  ! !                                                      ! calculation
  ! ! !
  ! ! INTEGER :: nath,            &             ! Number of (real) atoms in the primitive cell
  ! !                                           ! with Hubbard_U \= 0
  ! !            nath_sc,         &             ! Total number of real+virtual Hubbard atoms
  ! !                                           ! in the virtual supercell
  ! INTEGER :: nqsh
  ! !! Number of q points in the grid (without symmetry reduction)
  ! !! = number of primitive cells in the virtual supercell
  ! INTEGER :: nah_pert
  ! !! Site number of the perturbed Hubbard atom
  ! !            nath_pert,       &             ! Number of actual perturbed Hubbard atoms in the primitive cell
  ! !            find_atpert,     &             ! Method of searching for atoms to be perturbed
  ! !            ntyp_new,        &             ! Maximum number of different types detected in the calculation
  ! !                                           ! (used only when find_atpert=3)
  ! !            num_neigh,       &             ! Used in the postprocessing: number of nearest neighbors of every atom
  ! !                                           ! which will be written to the file parameters.out
  ! !                                           ! (can be used only with lda_plus_u_kind = 2)
  ! !            lmin,            &             ! Used in the postprocessing: minimum value of the orbital quantum number
  ! !                                           ! of the Hubbard atoms (Hubbard_l) starting from which (and up to the maximum
  ! !                                           ! Hubbard_l in the system) Hubbard V will be written to the file parameters.out
  ! ! !
  ! ! INTEGER :: equiv_type(ntypx)              ! equiv_type(i)=j, will merge type i to type j
  ! !                                           ! (useful when nspin=2)
  ! ! !
  ! ! CHARACTER(LEN=16)  :: background          ! Background correction
  !
  LOGICAL :: write_wan_Rr
  !! If .true. write the collected real-space Wannier functions to file.
  !! If .false. use the existing files in folder_wan_Rr.
  LOGICAL :: qplot
  !! if TRUE the q are read from input
  !
  CHARACTER(LEN=256) :: tmp_dir_save
  !! Temporary directory
  CHARACTER(LEN=256) :: folder_wan_Rr
  !! Folder with Wannier functions in real space
  CHARACTER(LEN=256) :: folder_elph
  !! Output folder where electron-phonon matrix elements will be written
  CHARACTER(LEN=256) :: folder_phonon
  !! Folder where the phonon potential are written.
  CHARACTER(LEN=256) :: wannier_seedname
  !! Seedname for Wannier90 calculation
  CHARACTER(LEN=256) :: fildvscf
  !! Name of the dvscf (phonon potential) file
  !
  ! CHARACTER(LEN=256) :: filU
  ! !! Output filename prefix for writing the U tensor.
  ! !! Bare and screened U tensors for each q are written to filU.bare#iq and filU.scrd#iq.
  ! !
  ! REAL(DP) :: conv_thr_chi
  ! !! Threshold for the calculation of chi
  ! !             docc_thr,          &          ! Threshold for the comparison of the unperturbed
  ! !                                           ! occupations (used only with find_atpert=1 for
  ! !                                           ! determination of atoms which must be perturbed)
  ! !             rmax,              &          ! Maximum distance (in Bohr) between two atoms
  ! !                                           ! to search for neighbors (used only at the
  ! !                                           ! postprocessing step when lda_plus_u_kind = 2).
  ! ! !
  ! !
  ! !
  ! ! Definition of active space
  ! !
  ! CHARACTER(LEN=256) :: active_space
  ! !! Method of defining the active space
  ! !! = 'bands' : use band indices from active_bands_min to active_bands_max
  ! INTEGER :: active_bands_max
  ! !! For active_space == 'bands': Maximum band index in the active space
  ! INTEGER :: active_bands_min
  ! !! For active_space == 'bands': Minimum band index in the active space
  ! !
  ! REAL(DP), ALLOCATABLE :: ns(:)
  ! !! Trace of unperturbed occupations (spin up + spin down)
  ! REAL(DP), ALLOCATABLE :: magn(:)
  ! !! Unperturbed magnetization
  ! REAL(DP), ALLOCATABLE :: Rvect(:, :)
  ! !! Radius-vector of the primitive cell
  ! REAL(DP), ALLOCATABLE :: chi0(:,:)
  ! !! Bare response function (from 1st iteration)
  ! REAL(DP), ALLOCATABLE :: chi(:,:)
  ! !! SCF response function
  ! ! !
  ! COMPLEX(DP), ALLOCATABLE :: dns0(:,:,:,:,:)
  ! !! Bare response occupation matrix (from 1st iteration)
  ! COMPLEX(DP), ALLOCATABLE :: dns0_tot(:,:,:,:,:)
  ! !! Total bare response occupation matrix (summed over q)
  ! COMPLEX(DP), ALLOCATABLE :: dnsscf_tot(:,:,:,:,:)
  ! !! Total SCF  response occupation matrix (summed over q)
  ! COMPLEX(DP), ALLOCATABLE :: trace_dns_tot_old(:)
  ! !! Trace of the response occupation matrix (for a convergence test)
  ! ! !
  ! INTEGER, ALLOCATABLE :: ityp_new(:)
  ! !! Types of atoms
  ! !
  ! ! Variables related to the k point mesh
  ! !
  ! INTEGER :: nks_orig
  ! ! Size of the global k point mesh for the underlying SCF calculation
  ! REAL(DP), ALLOCATABLE :: xk_orig(:, :)
  ! ! Original global k point mesh for the underlying SCF calculation. (3, nks_orig)
  ! INTEGER, ALLOCATABLE :: ik_to_ik_orig(:)
  ! ! Index of current k and k+q points modulo G in the xk_orig.
  ! ! xk(:, ik) = xk_orig(:, ik_to_ik_orig(ik)) + G for ik = 1, ..., nks. (nks = 2 * nksq)
  ! ! Here, xk(:, ik) includes both k and k+q points.
  ! !
  ! ! Coulumb matrix elements
  ! COMPLEX(DP), ALLOCATABLE :: v_coul_bare(:, :)
  ! !! Bare coulomb matrix elements. (ipert, jpert).
  ! COMPLEX(DP), ALLOCATABLE :: v_coul_scrd(:, :)
  ! !! Screened coulomb matrix elements. (ipert, jpert).
  ! !!
  ! COMPLEX(DP), ALLOCATABLE :: dvbare(:, :)
  ! !! Bare perturbing potential
  ! !
  ! INTEGER :: iuwan
  ! !! Unit for reading and writing the Wannier functions
  ! INTEGER :: lrwan
  ! !! Record length of the Wannier function buffer
  ! !
  ! REAL(DP) :: dist_thr
  ! !! Cutoff distance for the Wannier functions in bohr units. For perturbation.
  ! REAL(DP) :: dist_thr_large
  ! !! Larger cutoff distance for the Wannier functions in bohr units. For matrix element calculation.
  !
END MODULE elphcom
