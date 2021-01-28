MODULE mod_sirius_base
USE ISO_C_BINDING
IMPLICIT NONE

! use SIRIUS to solve KS equations
LOGICAL :: use_sirius_ks_solver               = .TRUE.
! use SIRIUS to generate density
LOGICAL :: use_sirius_density                 = .TRUE.
! use SIRIUS to generate effective potential; WARNING: currently must be always set to .false.
LOGICAL :: use_sirius_potential               = .FALSE.
! use SIRIUS to generate density matrix ('bec' thing in QE) WARNING: currently must be set to the value of use_sirius_density
LOGICAL :: use_sirius_density_matrix          = .TRUE.
! use SIRIUS to compute local part of pseudopotential
LOGICAL :: use_sirius_vloc                    = .TRUE.
! use SIRIUS to compute core charge density
LOGICAL :: use_sirius_rho_core                = .TRUE.
! use SIRIUS to compute plane-wave coefficients of atomic charge density
LOGICAL :: use_sirius_rho_atomic              = .TRUE.
! use SIRIUS to compute forces
LOGICAL :: use_sirius_forces                  = .TRUE.
! use SIRIUS to compute stress tensor
LOGICAL :: use_sirius_stress                  = .TRUE.
! Setup simulation context even if SIRIUS is not used (default is False)
LOGICAL :: always_setup_sirius                = .TRUE.

! inverse of the reciprocal lattice vectors matrix
REAL(8) bg_inv(3,3)
! total number of k-points
INTEGER num_kpoints
REAL(8), ALLOCATABLE :: kpoints(:,:)
REAL(8), ALLOCATABLE :: wkpoints(:)
REAL(8), ALLOCATABLE :: beta_ri_tab(:,:,:)
REAL(8), ALLOCATABLE :: aug_ri_tab(:,:,:,:)

TYPE atom_type_t
  ! atom label
  CHARACTER(len=100, kind=C_CHAR) :: label
  ! nh(iat) in the QE notation
  INTEGER                 :: num_beta_projectors
  ! lmax for beta-projectors
  INTEGER                 :: lmax
  ! plane-wave coefficients of Q-operator
  COMPLEX(8), ALLOCATABLE :: qpw(:, :)
END TYPE atom_type_t

TYPE(atom_type_t), ALLOCATABLE :: atom_type(:)

TYPE(C_PTR) :: sctx = C_NULL_PTR
TYPE(C_PTR) :: gs_handler = C_NULL_PTR
TYPE(C_PTR) :: ks_handler = C_NULL_PTR

END MODULE mod_sirius_base
