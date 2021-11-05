MODULE mod_sirius_base
USE ISO_C_BINDING
USE sirius
IMPLICIT NONE

! If true, the effective potential (Vha + Vxc + Vloc) is computed by QE
LOGICAL :: use_veff_callback = .FALSE.

! Setup simulation context even if SIRIUS is not used (default is False)
LOGICAL :: always_setup_sirius = .FALSE.

LOGICAL :: sirius_pwpp = .TRUE.

! inverse of the reciprocal lattice vectors matrix
REAL(8) bg_inv(3,3)
! total number of k-points
INTEGER num_kpoints
REAL(8), ALLOCATABLE :: kpoints(:,:)
REAL(8), ALLOCATABLE :: wkpoints(:)
! store rank and local k-point index for a gloobal k-point index
INTEGER, ALLOCATABLE :: kpoint_index_map(:,:)

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

CONTAINS

SUBROUTINE put_potential_to_sirius
  USE scf,                  ONLY : v, vltot
  USE gvect,                ONLY : mill, ngm
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : nspin_mag
  USE wavefunctions,        ONLY : psic

  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE uspp,                 ONLY : deeq
  USE uspp_param,           ONLY : nhm
  USE paw_variables,        ONLY : okpaw
  USE ions_base,            ONLY : nat
  USE sirius
  !
  IMPLICIT NONE
  !
  COMPLEX(8), ALLOCATABLE :: vxcg(:)
  COMPLEX(8) :: z1, z2
  INTEGER ig, is, ir, i, ia, j
  CHARACTER(10) label
  REAL(8), ALLOCATABLE :: deeq_tmp(:,:)
  REAL(8) :: d1,d2
  !
  IF (nspin.EQ.1.OR.nspin.EQ.4) THEN
    ! add local part of the potential and transform to PW domain
    psic(:) = v%of_r(:, 1) + vltot(:)
    CALL fwfft('Rho', psic, dfftp)
    ! convert to Hartree
    DO ig = 1, ngm
      v%of_g(ig, 1) = psic(dfftp%nl(ig)) * 0.5d0
    ENDDO
    ! set effective potential
    CALL sirius_set_pw_coeffs(gs_handler, "veff", v%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm)
  ENDIF

  IF (nspin.EQ.2) THEN
    DO is = 1, 2
      ! add local part of the potential and transform to PW domain
      psic(:) = v%of_r(:, is) + vltot(:)
      CALL fwfft('Rho', psic, dfftp)
      ! convert to Hartree
      DO ig = 1, ngm
         v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
      ENDDO
    ENDDO

    DO ig = 1, ngm
      z1 = v%of_g(ig, 1)
      z2 = v%of_g(ig, 2)
      v%of_g(ig, 1) = 0.5 * (z1 + z2)
      v%of_g(ig, 2) = 0.5 * (z1 - z2)
    ENDDO
    ! set effective potential and magnetization
    CALL sirius_set_pw_coeffs(gs_handler, "veff", v%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm)
    CALL sirius_set_pw_coeffs(gs_handler, "bz",   v%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm)
  ENDIF

  IF (nspin.EQ.4) THEN
    DO is = 2, nspin_mag
      psic(:) = v%of_r(:, is)
      CALL fwfft('Rho', psic, dfftp)
      ! convert to Hartree
      DO ig = 1, ngm
        v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
      ENDDO
      IF (is.EQ.2) label="bx"
      IF (is.EQ.3) label="by"
      IF (is.EQ.4) label="bz"

      CALL sirius_set_pw_coeffs(gs_handler, label, v%of_g(:, is), .TRUE., ngm, mill, intra_bgrp_comm)
    ENDDO
  ENDIF

END SUBROUTINE put_potential_to_sirius


SUBROUTINE get_density_from_sirius
  !
  USE scf,        ONLY : rho
  USE gvect,      ONLY : mill, ngm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE lsda_mod,   ONLY : nspin
  USE ions_base,  ONLY : nat, nsp, ityp
  USE uspp_param, ONLY : nhm, nh
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE wavefunctions, ONLY : psic
  USE control_flags,        ONLY : gamma_only
  USE sirius
  !
  IMPLICIT NONE
  !
  INTEGER iat, ig, ih, jh, ijh, na, ispn
  COMPLEX(8) z1, z2

  ! get rho(G)
  CALL sirius_get_pw_coeffs(gs_handler, "rho", rho%of_g(:, 1), ngm, mill, intra_bgrp_comm)
  IF (nspin.EQ.2) THEN
    CALL sirius_get_pw_coeffs(gs_handler, "magz", rho%of_g(:, 2), ngm, mill, intra_bgrp_comm)
  ENDIF
  IF (nspin.EQ.4) THEN
    CALL sirius_get_pw_coeffs(gs_handler, "magx", rho%of_g(:, 2), ngm, mill, intra_bgrp_comm)
    CALL sirius_get_pw_coeffs(gs_handler, "magy", rho%of_g(:, 3), ngm, mill, intra_bgrp_comm)
    CALL sirius_get_pw_coeffs(gs_handler, "magz", rho%of_g(:, 4), ngm, mill, intra_bgrp_comm)
  ENDIF
  ! get density matrix
  !CALL get_density_matrix_from_sirius
  DO ispn = 1, nspin
    psic(:) = 0.d0
    psic(dfftp%nl(:)) = rho%of_g(:, ispn)
    IF (gamma_only) psic(dfftp%nlm(:)) = CONJG(rho%of_g(:, ispn))
    CALL invfft('Rho', psic, dfftp)
    rho%of_r(:,ispn) = psic(:)
  ENDDO
END SUBROUTINE get_density_from_sirius


END MODULE mod_sirius_base
