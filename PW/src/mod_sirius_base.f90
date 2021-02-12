MODULE mod_sirius_base
USE ISO_C_BINDING
IMPLICIT NONE

! Setup simulation context even if SIRIUS is not used (default is False)
LOGICAL :: always_setup_sirius                = .FALSE.

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
