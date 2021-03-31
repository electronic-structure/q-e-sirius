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

!  ! convert Vxc to plane-wave domain
!  if (nspin.eq.1.or.nspin.eq.4) then
!    do ir = 1, dfftp%nnr
!      psic(ir) = vxc(ir, 1)
!    enddo
!  else
!    do ir = 1, dfftp%nnr
!      psic(ir) = 0.5d0 * (vxc(ir, 1) + vxc(ir, 2))
!    enddo
!  endif
!  call fwfft('Rho', psic, dfftp)
!  allocate(vxcg(ngm))
!  ! convert to Hartree
!  do ig = 1, ngm
!     vxcg(ig) = psic(dfftp%nl(ig)) * 0.5d0
!  end do
!  ! set XC potential
!  call sirius_set_pw_coeffs("vxc", vxcg(1), ngm, mill(1, 1), intra_bgrp_comm)
!  deallocate(vxcg)
!
!  ! update D-operator matrix
!  !call sirius_generate_d_operator_matrix()
!  !if (okpaw) then
!    allocate(deeq_tmp(nhm, nhm))
!  !  !! get D-operator matrix
!  !  !do ia = 1, nat
!  !  !  do is = 1, nspin
!  !  !    call sirius_get_d_operator_matrix(ia, is, deeq(1, 1, ia, is), nhm)
!  !  !  enddo
!  !  !  if (nspin.eq.2) then
!  !  !    do i = 1, nhm
!  !  !      do j = 1, nhm
!  !  !        d1 = deeq(i, j, ia, 1)
!  !  !        d2 = deeq(i, j, ia, 2)
!  !  !        deeq(i, j, ia, 1) = d1 + d2
!  !  !        deeq(i, j, ia, 2) = d1 - d2
!  !  !      enddo
!  !  !    enddo
!  !  !  endif
!  !  !  ! convert to Ry
!  !  !  deeq(:, :, ia, :) = deeq(:, :, ia, :) * 2
!  !  !enddo
!  !  !call add_paw_to_deeq(deeq)
!    do ia = 1, nat
!      do is = 1, nspin
!        if (nspin.eq.2.and.is.eq.1) then
!          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) + deeq(:, :, ia, 2)) / 2 ! convert to Ha
!        endif
!        if (nspin.eq.2.and.is.eq.2) then
!          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) - deeq(:, :, ia, 2)) / 2 ! convert to Ha
!        endif
!        if (nspin.eq.1.or.nspin.eq.4) then
!          deeq_tmp(:, :) = deeq(:, :, ia, is) / 2 ! convert to Ha
!        endif
!        call sirius_set_d_operator_matrix(ia, is, deeq_tmp(1, 1), nhm)
!      enddo
!    enddo
!    deallocate(deeq_tmp)
!  !endif

END SUBROUTINE put_potential_to_sirius

SUBROUTINE get_density_matrix_from_sirius
  !
  USE scf,        ONLY : rho
  USE ions_base,  ONLY : nat, nsp, ityp
  USE uspp_param, ONLY : nhm, nh
  USE lsda_mod,   ONLY : nspin
  !
  IMPLICIT NONE
  !
  COMPLEX(8), ALLOCATABLE :: dens_mtrx(:,:,:)
  INTEGER iat, na, ijh, ih, jh, ispn
  ! complex density matrix in SIRIUS has at maximum three components
  ALLOCATE(dens_mtrx(nhm, nhm, 3))
  DO iat = 1, nsp
    DO na = 1, nat
      IF (ityp(na).EQ.iat.AND.ALLOCATED(rho%bec)) THEN
        rho%bec(:, na, :) = 0.d0
        CALL sirius_get_density_matrix(gs_handler, na, dens_mtrx(1, 1, 1), nhm)
        ijh = 0
        DO ih = 1, nh(iat)
          DO jh = ih, nh(iat)
            ijh = ijh + 1
            IF (nspin.LE.2) THEN
              DO ispn = 1, nspin
                rho%bec(ijh, na, ispn) = dreal(dens_mtrx(ih, jh, ispn))
              ENDDO
            ENDIF
            IF (nspin.EQ.4) THEN
              rho%bec(ijh, na, 1) = dreal(dens_mtrx(ih, jh, 1) + dens_mtrx(ih, jh, 2))
              rho%bec(ijh, na, 4) = dreal(dens_mtrx(ih, jh, 1) - dens_mtrx(ih, jh, 2))
              rho%bec(ijh, na, 2) = 2.d0 * dreal(dens_mtrx(ih, jh, 3))
              rho%bec(ijh, na, 3) = -2.d0 * dimag(dens_mtrx(ih, jh, 3))
            ENDIF
            ! off-diagonal elements have a weight of 2
            IF (ih.NE.jh) THEN
              DO ispn = 1, nspin
                rho%bec(ijh, na, ispn) = rho%bec(ijh, na, ispn) * 2.d0
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(dens_mtrx)
END SUBROUTINE get_density_matrix_from_sirius



SUBROUTINE put_density_matrix_to_sirius
  !
  USE scf,        ONLY : rho
  USE ions_base,  ONLY : nat, nsp, ityp
  USE lsda_mod,   ONLY : nspin
  USE uspp_param, ONLY : nhm, nh
  USE uspp,       ONLY : becsum
  IMPLICIT NONE
  !
  INTEGER iat, na, ih, jh, ijh, ispn
  COMPLEX(8), ALLOCATABLE :: dens_mtrx(:,:,:)
  REAL(8), ALLOCATABLE :: dens_mtrx_tmp(:, :, :)
  REAL(8) fact
  ! set density matrix
  ! complex density matrix in SIRIUS has at maximum three components
  ALLOCATE(dens_mtrx_tmp(nhm * (nhm + 1) / 2, nat, nspin))
  !if (allocated(rho%bec)) then
  !  dens_mtrx_tmp = rho%bec
  !else
    dens_mtrx_tmp = becsum
  !endif

  ALLOCATE(dens_mtrx(nhm, nhm, 3))
  DO iat = 1, nsp
    DO na = 1, nat
      IF (ityp(na).EQ.iat) THEN
        dens_mtrx = (0.d0, 0.d0)
        ijh = 0
        DO ih = 1, nh(iat)
          DO jh = ih, nh(iat)
            ijh = ijh + 1
            ! off-diagonal elements have a weight of 2
            IF (ih.NE.jh) THEN
              fact = 0.5d0
            ELSE
              fact = 1.d0
            ENDIF
            IF (nspin.LE.2) THEN
              DO ispn = 1, nspin
                dens_mtrx(ih, jh, ispn) = fact * dens_mtrx_tmp(ijh, na, ispn)
                dens_mtrx(jh, ih, ispn) = fact * dens_mtrx_tmp(ijh, na, ispn)
              ENDDO
            ENDIF
            IF (nspin.EQ.4) THEN
              ! 0.5 * (rho + mz)
              dens_mtrx(ih, jh, 1) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) + dens_mtrx_tmp(ijh, na, 4))
              dens_mtrx(jh, ih, 1) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) + dens_mtrx_tmp(ijh, na, 4))
              ! 0.5 * (rho - mz)
              dens_mtrx(ih, jh, 2) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) - dens_mtrx_tmp(ijh, na, 4))
              dens_mtrx(jh, ih, 2) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) - dens_mtrx_tmp(ijh, na, 4))
              ! 0.5 * (mx - I * my)
              dens_mtrx(ih, jh, 3) = fact * 0.5 * dcmplx(dens_mtrx_tmp(ijh, na, 2), -dens_mtrx_tmp(ijh, na, 3))
              dens_mtrx(jh, ih, 3) = fact * 0.5 * dcmplx(dens_mtrx_tmp(ijh, na, 2), -dens_mtrx_tmp(ijh, na, 3))
            ENDIF
          ENDDO
        ENDDO
        CALL sirius_set_density_matrix(gs_handler, na, dens_mtrx(1, 1, 1), nhm)
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(dens_mtrx)
  DEALLOCATE(dens_mtrx_tmp)
END SUBROUTINE put_density_matrix_to_sirius


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


SUBROUTINE put_density_to_sirius
  !
  USE scf,        ONLY : rho
  USE gvect,      ONLY : mill, ngm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE lsda_mod,   ONLY : nspin
  USE ions_base,  ONLY : nat, nsp, ityp
  USE uspp_param, ONLY : nhm, nh
  USE sirius
  IMPLICIT NONE
  !
  COMPLEX(8), ALLOCATABLE :: rho_tot(:), mag(:)
  INTEGER iat, ig, ih, jh, ijh, na, ispn
  REAL(8) :: fact
  !
  !if (nspin.eq.1.or.nspin.eq.4) then
    CALL sirius_set_pw_coeffs(gs_handler, "rho", rho%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm)
  !endif

  IF (nspin.EQ.2) THEN
    !allocate(rho_tot(ngm))
    !allocate(mag(ngm))
    !do ig = 1, ngm
    !  rho_tot(ig) = rho%of_g(ig, 1) + rho%of_g(ig, 2)
    !  mag(ig) = rho%of_g(ig, 1) - rho%of_g(ig, 2)
    !enddo
    !call sirius_set_pw_coeffs(gs_handler, "rho", rho_tot(1), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    !call sirius_set_pw_coeffs(gs_handler, "magz", mag(1), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    !deallocate(rho_tot)
    !deallocate(mag)
    CALL sirius_set_pw_coeffs(gs_handler, "magz", rho%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm)
  ENDIF

  IF (nspin.EQ.4) THEN
    CALL sirius_set_pw_coeffs(gs_handler, "magx", rho%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm)
    CALL sirius_set_pw_coeffs(gs_handler, "magy", rho%of_g(:, 3), .TRUE., ngm, mill, intra_bgrp_comm)
    CALL sirius_set_pw_coeffs(gs_handler, "magz", rho%of_g(:, 4), .TRUE., ngm, mill, intra_bgrp_comm)
  ENDIF
END SUBROUTINE put_density_to_sirius


END MODULE mod_sirius_base
