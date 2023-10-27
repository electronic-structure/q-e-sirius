MODULE mod_sirius
  !! author: Anton Kozhevnikov, Simon Pintarelli, Mathieu Taillefumier
  !!
  !!
#if defined(__SIRIUS)
  !
  USE input_parameters, ONLY : sirius_cfg, use_sirius_scf, use_sirius_nlcg
  !
  USE sirius
  !
  IMPLICIT NONE
  !
  LOGICAL :: use_veff_callback = .FALSE.
  !! If true, the effective potential (Vha + Vxc + Vloc) is computed by QE
  !
  LOGICAL :: always_setup_sirius = .FALSE.
  !! Setup simulation context even if SIRIUS is not used (default is false)
  !
  LOGICAL :: sirius_pwpp = .TRUE.
  !! Use SIRIUS in plane-wave pseudopotential mode (true) or FP-LAPW (false)
  !
  REAL(8) bg_inv(3,3)
  !! inverse of the reciprocal lattice vectors matrix
  !
  INTEGER num_kpoints
  !! total number of k-points (without x2 lsda case)
  !
  REAL(8), ALLOCATABLE :: kpoints(:,:)
  !! List of k-points in fractional coordinates
  !
  REAL(8), ALLOCATABLE :: wkpoints(:)
  !! Weights of k-points without spin factor (sum of weights = 1)
  !
  INTEGER, ALLOCATABLE :: kpoint_index_map(:,:)
  !! store rank and local k-point index for a gloobal k-point index
  !
  TYPE atom_type_t
    !! atom type parameters
    CHARACTER(len=100, kind=C_CHAR) :: label
    !! atom label
    !
    INTEGER                 :: num_beta_projectors
    !! nh(iat) in the QE notation
    !
    INTEGER                 :: lmax
    !! lmax for beta-projectors
    !
    COMPLEX(8), ALLOCATABLE :: qpw(:, :)
    !! plane-wave coefficients of Q-operator
    !
    INTEGER                 :: num_chi
    !! number of atomic wave-functions
    !
    INTEGER, ALLOCATABLE    :: l_chi(:)
    !! orbital quantum number of atomic wave-functions
    !
    INTEGER, ALLOCATABLE    :: n_chi(:)
    !! principal quantum number of atomic wave-functions
    !
    REAL(8), ALLOCATABLE    :: chi(:, :)
    !! radial part of atomic wave-functions
    !
    INTEGER, ALLOCATABLE    :: idx_chi(:,:)
    !! points to one or two radial functions in the origianl list of atomic wfs
    !
    REAL(8), ALLOCATABLE    :: occ(:)
  END TYPE atom_type_t
  !
  TYPE(atom_type_t), ALLOCATABLE :: atom_type(:)
  !! list of atom type properties
  !
  TYPE(sirius_context_handler) :: sctx
  !! SIRIUS simulation context handler
  !
  TYPE(sirius_ground_state_handler) :: gs_handler
  !! SIRIUS ground state handler
  !
  TYPE(sirius_kpoint_set_handler) :: ks_handler
  !! SIRIUS k-point set handler
  !
 CONTAINS
  !
  !--------------------------------------------------------------------
  SUBROUTINE put_potential_to_sirius()
    !------------------------------------------------------------------
    !! Put the plane-wave part of QE potential into SIRIUS.
    !
    USE scf,                  ONLY : v, vltot
    USE gvect,                ONLY : mill, ngm
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE lsda_mod,             ONLY : nspin
    USE wavefunctions,        ONLY : psic
    USE fft_base,             ONLY : dfftp
    USE fft_interfaces,       ONLY : fwfft
    !
    IMPLICIT NONE
    !
    COMPLEX(8) :: z1, z2
    INTEGER ig, is
    CHARACTER(10) label
    !
    IF (nspin.EQ.1.OR.nspin.EQ.4) THEN
      ! add local part of the potential and transform to PW domain
      psic(:) = v%of_r(:, 1) + vltot(:)
      CALL fwfft( 'Rho', psic, dfftp )
      ! convert to Hartree
      DO ig = 1, ngm
        v%of_g(ig, 1) = psic(dfftp%nl(ig)) * 0.5d0
      ENDDO
      ! set effective potential
      CALL sirius_set_pw_coeffs( gs_handler, "veff", v%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm )
    ENDIF
    !
    IF (nspin.EQ.2) THEN
      DO is = 1, 2
        ! add local part of the potential and transform to PW domain
        psic(:) = v%of_r(:, is) + vltot(:)
        CALL fwfft( 'Rho', psic, dfftp )
        ! convert to Hartree
        DO ig = 1, ngm
           v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
        ENDDO
      ENDDO
      !
      DO ig = 1, ngm
        z1 = v%of_g(ig, 1)
        z2 = v%of_g(ig, 2)
        v%of_g(ig, 1) = 0.5 * (z1 + z2)
        v%of_g(ig, 2) = 0.5 * (z1 - z2)
      ENDDO
      ! set effective potential and magnetization
      CALL sirius_set_pw_coeffs( gs_handler, "veff", v%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm )
      CALL sirius_set_pw_coeffs( gs_handler, "bz",   v%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm )
    ENDIF
    !
    IF (nspin.EQ.4) THEN
      DO is = 2, 4
        psic(:) = v%of_r(:, is)
        CALL fwfft( 'Rho', psic, dfftp )
        ! convert to Hartree
        DO ig = 1, ngm
          v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
        ENDDO
        IF (is.EQ.2) label="bx"
        IF (is.EQ.3) label="by"
        IF (is.EQ.4) label="bz"
        !
        CALL sirius_set_pw_coeffs( gs_handler, label, v%of_g(:, is), .TRUE., ngm, mill, intra_bgrp_comm )
      ENDDO
    ENDIF
    !
  END SUBROUTINE put_potential_to_sirius
  !
  !--------------------------------------------------------------------
  SUBROUTINE put_density_to_sirius()
    !------------------------------------------------------------------
    !! Put plane-wave coefficients of density to SIRIUS
    !
    USE scf,                  ONLY : rho
    USE gvect,                ONLY : mill, ngm
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE lsda_mod,             ONLY : nspin
    !
    IMPLICIT NONE
    !
    INTEGER iat, ig, ih, jh, ijh, na, ispn
    COMPLEX(8) z1, z2
    !
    ! get rho(G)
    CALL sirius_set_pw_coeffs( gs_handler, "rho", rho%of_g(:, 1), .TRUE., ngm, mill, intra_bgrp_comm )
    IF (nspin.EQ.2) THEN
      CALL sirius_set_pw_coeffs( gs_handler, "magz", rho%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm )
    ENDIF
    IF (nspin.EQ.4) THEN
      CALL sirius_set_pw_coeffs( gs_handler, "magx", rho%of_g(:, 2), .TRUE., ngm, mill, intra_bgrp_comm )
      CALL sirius_set_pw_coeffs( gs_handler, "magy", rho%of_g(:, 3), .TRUE., ngm, mill, intra_bgrp_comm )
      CALL sirius_set_pw_coeffs( gs_handler, "magz", rho%of_g(:, 4), .TRUE., ngm, mill, intra_bgrp_comm )
    ENDIF
  END SUBROUTINE put_density_to_sirius
  !
  !--------------------------------------------------------------------
  SUBROUTINE get_density_from_sirius()
    !------------------------------------------------------------------
    !! Get plane-wave coefficients of density from SIRIUS
    !
    USE scf,                  ONLY : rho
    USE gvect,                ONLY : mill, ngm
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE lsda_mod,             ONLY : nspin
    USE fft_base,             ONLY : dfftp
    USE fft_rho,              ONLY : rho_g2r
    !
    IMPLICIT NONE
    !
    INTEGER iat, ig, ih, jh, ijh, na, ispn
    COMPLEX(8) z1, z2
    !
    ! get rho(G)
    CALL sirius_get_pw_coeffs( gs_handler, "rho", rho%of_g(:, 1), ngm, mill, intra_bgrp_comm )
    IF (nspin.EQ.2) THEN
      CALL sirius_get_pw_coeffs( gs_handler, "magz", rho%of_g(:, 2), ngm, mill, intra_bgrp_comm )
    ENDIF
    IF (nspin.EQ.4) THEN
      CALL sirius_get_pw_coeffs( gs_handler, "magx", rho%of_g(:, 2), ngm, mill, intra_bgrp_comm )
      CALL sirius_get_pw_coeffs( gs_handler, "magy", rho%of_g(:, 3), ngm, mill, intra_bgrp_comm )
      CALL sirius_get_pw_coeffs( gs_handler, "magz", rho%of_g(:, 4), ngm, mill, intra_bgrp_comm )
    ENDIF
    CALL rho_g2r (dfftp, rho%of_g, rho%of_r)
    !! get density matrix
    !!CALL get_density_matrix_from_sirius
  END SUBROUTINE get_density_from_sirius
  !
  !--------------------------------------------------------------------
  SUBROUTINE put_density_matrix_to_sirius
    !------------------------------------------------------------------
    !! Put QE density matrix to SIRIUS
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
      dens_mtrx_tmp = rho%bec
    !else
    !  dens_mtrx_tmp = becsum
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
          CALL sirius_set_density_matrix(gs_handler, na, dens_mtrx, nhm)
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(dens_mtrx)
    DEALLOCATE(dens_mtrx_tmp)
  END SUBROUTINE put_density_matrix_to_sirius
  !
  !--------------------------------------------------------------------
  SUBROUTINE calc_veff() BIND(C)
    !------------------------------------------------------------------
    !! Callback function to compute effective potential by QE
    !
    USE scf,      ONLY : rho, rho_core, rhog_core, v
    USE ener,     ONLY : ehart, vtxc, etxc
    USE ldaU,     ONLY : eth
    USE extfield, ONLY : tefield, etotefield
    !
    IMPLICIT NONE
    !
    REAL(8) :: charge
    !
    CALL get_density_from_sirius()
    CALL v_of_rho( rho, rho_core, rhog_core, ehart, etxc, vtxc, eth, etotefield, charge, v )
    CALL put_potential_to_sirius()
    !
  END SUBROUTINE calc_veff
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calc_ps_rho_radial_integrals( iat, nq, q, ps_rho_ri ) BIND(C)
    !---------------------------------------------------------------------
    !! Callback function to compute radial integals of the free atomic density
    !
    USE kinds,      ONLY : DP
    USE atom,       ONLY : rgrid, msh
    USE uspp_param, ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
    REAL(KIND=c_double), INTENT(IN)  :: q(nq)
    REAL(KIND=c_double), INTENT(OUT) :: ps_rho_ri(nq)
    !
    REAL(DP) :: aux(rgrid(iat)%mesh)
    INTEGER :: iq, ir, nr
    !
    nr = msh(iat)
    !
    DO iq = 1, nq
      CALL sph_bes(nr, rgrid(iat)%r(1), q(iq), 0, aux)
      DO ir = 1, nr
        aux(ir) = aux(ir) * upf(iat)%rho_at(ir)
      ENDDO
      CALL simpson( nr, aux, rgrid(iat)%rab(1), ps_rho_ri(iq) )
      ps_rho_ri(iq) = ps_rho_ri(iq) / 12.566370614359172954d0
    ENDDO
    !
  END SUBROUTINE calc_ps_rho_radial_integrals
  !
  !--------------------------------------------------------------------
  SUBROUTINE calc_vloc_radial_integrals( iat, nq, q, vloc_ri ) BIND(C)
    !------------------------------------------------------------------
    !! Callback function to compute radial integrals of Vloc(r)
    !
    USE kinds,       ONLY : DP
    USE atom,        ONLY : rgrid, msh
    USE constants,   ONLY : e2, eps8
    USE esm,         ONLY : do_comp_esm, esm_bc
    USE Coul_cut_2D, ONLY : do_cutoff_2D, lz
    USE uspp_param,  ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
    REAL(KIND=c_double), INTENT(IN)  :: q(nq)
    REAL(KIND=c_double), INTENT(OUT) :: vloc_ri(nq)
    !
    REAL(DP) :: aux(rgrid(iat)%mesh), aux1(rgrid(iat)%mesh), r, q2
    INTEGER :: iq, ir, nr
    !
    ! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
    nr = msh(iat)
    ! q-independent radial function
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux1(ir) = r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * ERF( r )
    END DO
    ! loop over q-points
    DO iq = 1, nq
      IF (q(iq) < eps8) THEN ! q=0 case
        !
        ! first the G=0 term
        !
        IF ((do_comp_esm .AND. (esm_bc .NE. 'pbc')) .OR. do_cutoff_2D) THEN
          !
          ! ... temporarily redefine term for ESM calculation
          !
          DO ir = 1, nr
            r = rgrid(iat)%r(ir)
            aux(ir) = r * (r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * ERF( r ))
          END DO
          IF (do_cutoff_2D .AND. rgrid(iat)%r(nr) > lz) THEN
            CALL errore( 'vloc_of_g','2D cutoff is smaller than pseudo cutoff radius: &
              & increase interlayer distance (or see Modules/read_pseudo.f90)', 1 )
          END IF
        ELSE
          DO ir = 1, nr
            r = rgrid(iat)%r(ir)
            aux(ir) = r * (r * upf(iat)%vloc(ir) + upf(iat)%zp * e2)
          END DO
        END IF
        CALL simpson( nr, aux, rgrid(iat)%rab(1), vloc_ri(iq) )
      ELSE ! q > 0 case
        DO ir = 1, nr
          r = rgrid(iat)%r(ir)
          aux(ir) = aux1(ir) * SIN( q(iq) * r ) / q(iq)
        END DO
        CALL simpson( nr, aux, rgrid(iat)%rab(1), vloc_ri(iq) )
        IF ((.NOT.do_comp_esm) .OR. (esm_bc .EQ. 'pbc')) THEN
          !
          !   here we re-add the analytic fourier transform of the erf function
          !
          IF (.NOT. do_cutoff_2D) THEN
            q2 = q(iq) * q(iq)
            vloc_ri(iq) = vloc_ri(iq) - upf(iat)%zp * e2 * EXP( -q2 * 0.25d0 ) / q2
          END IF
        END IF
      END IF
      ! convert to Ha
      vloc_ri(iq) = vloc_ri(iq) / 2.0
    END DO
    !
  END SUBROUTINE calc_vloc_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_vloc_dj_radial_integrals( iat, nq, q, vloc_dj_ri ) BIND(C)
    !-----------------------------------------------------------------------
    !! A callback function to compute radial integrals of Vloc(r) with
    !! derivatives of Bessel functions
    !
    USE kinds,       ONLY : DP
    USE atom,        ONLY : rgrid, msh
    USE constants,   ONLY : e2, eps8
    USE esm,         ONLY : do_comp_esm, esm_bc
    USE Coul_cut_2D, ONLY : do_cutoff_2D
    USE uspp_param,  ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
    REAL(KIND=c_double), INTENT(IN)  :: q(nq)
    REAL(KIND=c_double), INTENT(OUT) :: vloc_dj_ri(nq)
    !
    REAL(DP) :: aux(rgrid(iat)%mesh), aux1(rgrid(iat)%mesh), r, q2
    INTEGER :: iq, ir, nr
    !
    ! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
    nr = msh(iat)
    ! q-independent radial function
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux1(ir) = r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * ERF( r )
    END DO
    ! loop over q-points
    DO iq = 1, nq
      IF (q(iq) < eps8) THEN ! q=0 case
        vloc_dj_ri(iq) = 0.d0
      ELSE ! q > 0 case
        DO ir = 1, nr
          r = rgrid(iat)%r(ir)
          aux(ir) = aux1(ir) * (SIN( q(iq) * r ) / q(iq)**2 - r * COS( q(iq) * r ) / q(iq))
        END DO
        CALL simpson( nr, aux, rgrid(iat)%rab(1), vloc_dj_ri(iq) )
        vloc_dj_ri(iq) = vloc_dj_ri(iq) / q(iq)
        IF ((.NOT.do_comp_esm)  .OR.  (esm_bc .EQ. 'pbc')) THEN
          IF (.NOT. do_cutoff_2D) THEN
            q2 = q(iq) * q(iq)
            vloc_dj_ri(iq) = vloc_dj_ri(iq) - upf(iat)%zp * e2 * &
              &EXP( -q2 * 0.25d0 ) * (q2 + 4) / 2 / q2 / q2
          END IF
        END IF
      END IF
      ! convert to Ha
      vloc_dj_ri(iq) = vloc_dj_ri(iq) / 2.0
    END DO
    !
  END SUBROUTINE calc_vloc_dj_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_rhoc_radial_integrals( iat, nq, q, rhoc_ri ) BIND(C)
    !-----------------------------------------------------------------------
    !! Ccallback function to compute radial integrals of rho_core(r)
    !
    USE kinds,       ONLY : DP
    USE atom,        ONLY : rgrid, msh
    USE constants,   ONLY : eps8
    USE uspp_param,  ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
    REAL(KIND=c_double), INTENT(IN)  :: q(nq)
    REAL(KIND=c_double), INTENT(OUT) :: rhoc_ri(nq)
    !
    REAL(DP) :: aux(rgrid(iat)%mesh), r
    INTEGER :: iq, ir, nr
    !
    ! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
    nr = msh(iat)
    ! loop over q-points
    DO iq = 1, nq
      IF (q(iq) < eps8) THEN ! q=0 case
        DO ir = 1, nr
          r = rgrid(iat)%r(ir)
          aux(ir) = r**2 * upf(iat)%rho_atc(ir)
        ENDDO
      ELSE
        CALL sph_bes( nr, rgrid(iat)%r(1), q(iq), 0, aux )
        DO ir = 1, nr
          r = rgrid(iat)%r(ir)
          aux(ir) = r**2 * upf(iat)%rho_atc(ir) * aux(ir)
        ENDDO
      END IF
      CALL simpson( nr, aux, rgrid(iat)%rab(1), rhoc_ri(iq) )
    END DO
    !
  END SUBROUTINE calc_rhoc_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_rhoc_dj_radial_integrals( iat, nq, q, rhoc_dj_ri ) BIND(C)
  !-------------------------------------------------------------------------
    !! Callback function to compute radial integrals or rho_core(r) with
    !! the derivatives of Bessel functions.
    !
    USE kinds,       ONLY : DP
    USE atom,        ONLY : rgrid, msh
    USE constants,   ONLY : eps8
    USE uspp_param,  ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
    INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
    REAL(KIND=c_double), INTENT(IN)  :: q(nq)
    REAL(KIND=c_double), INTENT(OUT) :: rhoc_dj_ri(nq)
    !
    REAL(DP) :: aux(rgrid(iat)%mesh), r
    INTEGER :: iq, ir, nr
    ! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
    nr = msh(iat)
    ! loop over q-points
    DO iq = 1, nq
      IF (q(iq) < eps8) THEN ! q=0 case
        rhoc_dj_ri(iq) = 0.d0
      ELSE
        DO ir = 1, nr
          r = rgrid(iat)%r(ir)
          aux(ir) = r * upf(iat)%rho_atc(ir) * &
            &(r * COS( q(iq) * r ) / q(iq) - SIN( q(iq) * r ) / q(iq)**2 )
        ENDDO
        CALL simpson( nr, aux, rgrid(iat)%rab(1), rhoc_dj_ri(iq) )
      END IF
    END DO
    !
  END SUBROUTINE calc_rhoc_dj_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_beta_radial_integrals( iat, q, beta_ri, ld ) BIND(C)
    !-----------------------------------------------------------------------
    !! Callback to compute radial integrals of beta-projectors
    !
    USE uspp_data,    ONLY : dq, tab, beta_ri_tab
    USE uspp_param,   ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(kind=c_int), INTENT(in), VALUE :: iat
    REAL(kind=c_double), INTENT(in), VALUE :: q
    INTEGER(kind=c_int), INTENT(in), VALUE :: ld
    REAL(kind=c_double), INTENT(out) :: beta_ri(ld)
    !
    REAL(8) :: px, ux, vx, wx
    INTEGER :: i0, i1, i2, i3, ib
    !
    IF (ld.LT.upf(iat)%nbeta) THEN
      WRITE(*,*)'not enough space to store all beta projectors, ld=',ld,' nbeta=',upf(iat)%nbeta
      STOP
    ENDIF
    IF (.NOT.ALLOCATED(tab)) THEN
      WRITE(*,*)'tab array is not allocated'
      STOP
    ENDIF
    !
    px = q / dq - INT(q / dq)
    ux = 1.d0 - px
    vx = 2.d0 - px
    wx = 3.d0 - px
    i0 = INT(q / dq) + 1
    i1 = i0 + 1
    i2 = i0 + 2
    i3 = i0 + 3
    DO ib = 1, upf(iat)%nbeta
      beta_ri(ib) = beta_ri_tab(i0, ib, iat) * ux * vx * wx / 6.d0 + &
                    beta_ri_tab(i1, ib, iat) * px * vx * wx / 2.d0 - &
                    beta_ri_tab(i2, ib, iat) * px * ux * wx / 2.d0 + &
                    beta_ri_tab(i3, ib, iat) * px * ux * vx / 6.d0
    ENDDO
    !
  END SUBROUTINE calc_beta_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_beta_dj_radial_integrals( iat, q, beta_ri, ld ) BIND(C)
    !-------------------------------------------------------------------------
    !! Callback function to compute radial integrals of beta-projectors with
    !! the derrivatives of spherical Bessel functions.
    !
    USE uspp_data,    ONLY : dq, tab, beta_ri_tab
    USE uspp_param,   ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(kind=c_int), INTENT(in), VALUE :: iat
    REAL(kind=c_double), INTENT(in), VALUE :: q
    INTEGER(kind=c_int), INTENT(in), VALUE :: ld
    REAL(kind=c_double), INTENT(out) :: beta_ri(ld)
    !
    REAL(8) :: px, ux, vx, wx
    INTEGER :: i0, i1, i2, i3, ib
    !
    IF (ld.LT.upf(iat)%nbeta) THEN
      WRITE(*,*)'not enough space to store all beta projectors, ld=',ld,' nbeta=',upf(iat)%nbeta
      STOP
    ENDIF
    IF (.NOT.ALLOCATED(tab)) THEN
      WRITE(*,*)'tab array is not allocated'
      STOP
    ENDIF
    !
    px = q / dq - INT(q / dq)
    ux = 1.d0 - px
    vx = 2.d0 - px
    wx = 3.d0 - px
    i0 = INT(q / dq) + 1
    i1 = i0 + 1
    i2 = i0 + 2
    i3 = i0 + 3
    DO ib = 1, upf(iat)%nbeta
      beta_ri(ib) = beta_ri_tab(i0, ib, iat) * (-vx*wx-ux*wx-ux*vx)/6.d0 + &
                    beta_ri_tab(i1, ib, iat) * (+vx*wx-px*wx-px*vx)/2.d0 - &
                    beta_ri_tab(i2, ib, iat) * (+ux*wx-px*wx-px*ux)/2.d0 + &
                    beta_ri_tab(i3, ib, iat) * (+ux*vx-px*vx-px*ux)/6.d0
      beta_ri(ib) = beta_ri(ib) / dq
    ENDDO
    !
  END SUBROUTINE calc_beta_dj_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_aug_radial_integrals(iat, q, aug_ri, ld1, ld2) BIND(C)
    !-----------------------------------------------------------------------
    !! Callback function to compute radial integrals of augmentation charge.
    !
    USE uspp_data,    ONLY : dq, qrad, aug_ri_tab
    USE uspp_param,   ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(kind=c_int), INTENT(in), VALUE :: iat
    REAL(kind=c_double), INTENT(in), VALUE :: q
    INTEGER(kind=c_int), INTENT(in), VALUE :: ld1
    INTEGER(kind=c_int), INTENT(in), VALUE :: ld2
    REAL(kind=c_double), INTENT(out) :: aug_ri(ld1, ld2)
    !
    REAL(8) :: px, ux, vx, wx
    INTEGER :: i0, i1, i2, i3, l, nb, mb, ijv
    !
    IF (upf(iat)%tvanp) THEN
      px = q / dq - INT(q / dq)
      ux = 1.d0 - px
      vx = 2.d0 - px
      wx = 3.d0 - px
      i0 = INT(q / dq) + 1
      i1 = i0 + 1
      i2 = i0 + 2
      i3 = i0 + 3
      DO l = 1, 2 * atom_type(iat)%lmax + 1
        DO nb = 1, upf(iat)%nbeta
          DO mb = nb, upf(iat)%nbeta
            ijv = mb * (mb-1) / 2 + nb
            aug_ri(ijv, l) = aug_ri_tab(i0, ijv, l, iat) * ux * vx * wx / 6.d0 + &
                             aug_ri_tab(i1, ijv, l, iat) * px * vx * wx / 2.d0 - &
                             aug_ri_tab(i2, ijv, l, iat) * px * ux * wx / 2.d0 + &
                             aug_ri_tab(i3, ijv, l, iat) * px * ux * vx / 6.d0
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
  END SUBROUTINE calc_aug_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_aug_dj_radial_integrals(iat, q, aug_ri, ld1, ld2) BIND(C)
    !-----------------------------------------------------------------------
    !! Callback function to compute radial integrals of augmentation charge with
    !! derivattives of spherical Bessel functions.
    !
    USE uspp_data,    ONLY : dq, qrad, aug_ri_tab
    USE uspp_param,   ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(kind=c_int), INTENT(in), VALUE :: iat
    REAL(kind=c_double), INTENT(in), VALUE :: q
    INTEGER(kind=c_int), INTENT(in), VALUE :: ld1
    INTEGER(kind=c_int), INTENT(in), VALUE :: ld2
    REAL(kind=c_double), INTENT(out) :: aug_ri(ld1, ld2)
    !
    REAL(8) :: px, ux, vx, wx
    INTEGER :: i0, i1, i2, i3, l, nb, mb, ijv
    !
    IF (upf(iat)%tvanp) THEN
      px = q / dq - INT(q / dq)
      ux = 1.d0 - px
      vx = 2.d0 - px
      wx = 3.d0 - px
      i0 = INT(q / dq) + 1
      i1 = i0 + 1
      i2 = i0 + 2
      i3 = i0 + 3
      DO l = 1, 2 * atom_type(iat)%lmax + 1
        DO nb = 1, upf(iat)%nbeta
          DO mb = nb, upf(iat)%nbeta
            ijv = mb * (mb-1) / 2 + nb
            aug_ri(ijv, l) = - aug_ri_tab(i0, ijv, l, iat) * (ux*vx + vx*wx + ux*wx) / 6.d0 &
                             + aug_ri_tab(i1, ijv, l, iat) * (wx*vx - px*wx - px*vx) * 0.5d0 &
                             - aug_ri_tab(i2, ijv, l, iat) * (wx*ux - px*wx - px*ux) * 0.5d0 &
                             + aug_ri_tab(i3, ijv, l, iat) * (ux*vx - px*ux - px*vx) / 6.d0
            aug_ri(ijv, l) = aug_ri(ijv, l) / dq
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
  END SUBROUTINE calc_aug_dj_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_atomic_wfc_radial_integrals(iat, q, wfc_ri, ld) BIND(C)
    !-----------------------------------------------------------------------
    !! Callback function to compute radial integrals of the atomic wave functions
    !
    USE iso_c_binding
    USE uspp_data,    ONLY : dq, tab, wfc_ri_tab
    USE uspp_param,   ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(kind=c_int), INTENT(in), VALUE :: iat
    REAL(kind=c_double), INTENT(in), VALUE :: q
    INTEGER(kind=c_int), INTENT(in), VALUE :: ld
    REAL(kind=c_double), INTENT(out) :: wfc_ri(ld)
    !
    REAL(8) :: px, ux, vx, wx
    REAL(8) :: wfc_ri_tmp(upf(iat)%nwfc)
    INTEGER :: i0, i1, i2, i3, ib, j
    !
    IF ( ld .LT. atom_type(iat)%num_chi ) THEN
      WRITE(*,*)'not enough space to store all atomic wave functions, ld=',ld,' nwfc=',atom_type(iat)%num_chi
      STOP
    ENDIF
    IF (.NOT.ALLOCATED(tab)) THEN
      WRITE(*,*)'tab array is not allocated'
      STOP
    ENDIF
    !
    px = q / dq - INT(q / dq)
    ux = 1.d0 - px
    vx = 2.d0 - px
    wx = 3.d0 - px
    i0 = INT(q / dq) + 1
    i1 = i0 + 1
    i2 = i0 + 2
    i3 = i0 + 3
    DO ib = 1, upf(iat)%nwfc
      wfc_ri_tmp(ib) = wfc_ri_tab(i0, ib, iat) * ux * vx * wx / 6.d0 + &
                       wfc_ri_tab(i1, ib, iat) * px * vx * wx / 2.d0 - &
                       wfc_ri_tab(i2, ib, iat) * px * ux * wx / 2.d0 + &
                       wfc_ri_tab(i3, ib, iat) * px * ux * vx / 6.d0
    ENDDO
    DO ib = 1, atom_type(iat)%num_chi
      j = atom_type(iat)%idx_chi(1, ib)
      wfc_ri(ib) = wfc_ri_tmp(j)
      IF ( atom_type(iat)%idx_chi(2, ib) .NE. -1 ) THEN
        j = atom_type(iat)%idx_chi(2, ib)
        wfc_ri(ib) = 0.5d0 * wfc_ri(ib) + 0.5d0 * wfc_ri_tmp(j)
      END IF
    ENDDO
    !
  END SUBROUTINE calc_atomic_wfc_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE calc_atomic_wfc_djl_radial_integrals(iat, q, wfc_ri, ld) BIND(C)
    !-----------------------------------------------------------------------
    !! Callback function to compute radial integrals of the atomic wave functions
    !
    USE iso_c_binding
    USE uspp_data,    ONLY : dq, tab, wfc_ri_tab
    USE uspp_param,   ONLY : upf
    !
    IMPLICIT NONE
    !
    INTEGER(kind=c_int), INTENT(in), VALUE :: iat
    REAL(kind=c_double), INTENT(in), VALUE :: q
    INTEGER(kind=c_int), INTENT(in), VALUE :: ld
    REAL(kind=c_double), INTENT(out) :: wfc_ri(ld)
    !
    REAL(8) :: px, ux, vx, wx
    REAL(8) :: wfc_ri_tmp(upf(iat)%nwfc)
    INTEGER :: i0, i1, i2, i3, ib, j
    !
    IF ( ld .LT. atom_type(iat)%num_chi ) THEN
      WRITE(*,*)'not enough space to store all atomic wave functions, ld=',ld,' nwfc=',atom_type(iat)%num_chi
      STOP
    ENDIF
    IF (.NOT.ALLOCATED(tab)) THEN
      WRITE(*,*)'tab array is not allocated'
      STOP
    ENDIF
    !
    px = q / dq - INT(q / dq)
    ux = 1.d0 - px
    vx = 2.d0 - px
    wx = 3.d0 - px
    i0 = INT(q / dq) + 1
    i1 = i0 + 1
    i2 = i0 + 2
    i3 = i0 + 3
    DO ib = 1, upf(iat)%nwfc
      wfc_ri_tmp(ib) = - wfc_ri_tab(i0, ib, iat) * (ux*vx + vx*wx + ux*wx) / 6.d0 &
                       + wfc_ri_tab(i1, ib, iat) * (wx*vx - px*wx - px*vx) * 0.5d0 &
                       - wfc_ri_tab(i2, ib, iat) * (wx*ux - px*wx - px*ux) * 0.5d0 &
                       + wfc_ri_tab(i3, ib, iat) * (ux*vx - px*ux - px*vx) / 6.d0
      wfc_ri_tmp(ib) = wfc_ri_tmp(ib) / dq
    ENDDO
    DO ib = 1, atom_type(iat)%num_chi
      j = atom_type(iat)%idx_chi(1, ib)
      wfc_ri(ib) = wfc_ri_tmp(j)
      IF ( atom_type(iat)%idx_chi(2, ib) .NE. -1 ) THEN
        j = atom_type(iat)%idx_chi(2, ib)
        wfc_ri(ib) = 0.5d0 * wfc_ri(ib) + 0.5d0 * wfc_ri_tmp(j)
      END IF
    ENDDO
    !
  END SUBROUTINE calc_atomic_wfc_djl_radial_integrals
  !
  !-------------------------------------------------------------------------
  SUBROUTINE setup_sirius()
    !-----------------------------------------------------------------------
    !! Setup SIRIUS simulation context, create k-point set and DFT ground state instance.
    !
    USE cell_base,            ONLY : alat, at, bg
    USE ions_base,            ONLY : tau, nsp, atm, zv, amass, ityp, nat
    USE uspp_param,           ONLY : upf, nhm, nh
    USE atom,                 ONLY : rgrid, msh
    USE fft_base,             ONLY : dfftp
    USE klist,                ONLY : nks, nkstot, degauss, ngauss
    USE gvect,                ONLY : ngm_g, ecutrho, ngm
    USE gvecw,                ONLY : ecutwfc
    USE control_flags,        ONLY : gamma_only, diago_full_acc, mixing_beta, nmix, iverbosity
    USE mp_world,             ONLY : mpime
    USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm, npool
    USE mp_images,            ONLY : nproc_image, intra_image_comm
    USE wvfct,                ONLY : nbnd
    USE parallel_include,     ONLY : MPI_IN_PLACE, MPI_INT, MPI_SUM
    USE input_parameters,     ONLY : diago_david_ndim
    USE noncollin_module,     ONLY : noncolin, npol, angle1, angle2, lspinorb
    USE lsda_mod,             ONLY : lsda, nspin, starting_magnetization
    USE symm_base,            ONLY : nosym, nsym
    USE ldaU,                 ONLY : lda_plus_U, Hubbard_J, Hubbard_U, Hubbard_alpha, &
                                   & Hubbard_beta, is_Hubbard, lda_plus_u_kind, &
                                   & Hubbard_J0, Hubbard_projectors, Hubbard_l, Hubbard_n, Hubbard_occ, &
                                   & ldim_u, neighood, at_sc, Hubbard_V
    USE esm,                  ONLY : do_comp_esm
    USE Coul_cut_2D,          ONLY : do_cutoff_2D
    USE constants,            ONLY : RYTOEV
    USE kinds,                ONLY : DP
    USE scf,                  ONLY : rho
    USE paw_variables,        ONLY : okpaw
    !
    IMPLICIT NONE
    !
    INTEGER :: dims(3), i, ia, iat, rank, ierr, ijv, j, l, ir, num_gvec, num_ranks_k, &
             & iwf, nmagd, viz, ia2, iat2, atom_pair(2), n_pair(2), l_pair(2)
    REAL(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3)
    REAL(8), ALLOCATABLE :: dion(:, :), vloc(:)
    INTEGER :: lmax_beta, nsymop
    CHARACTER(LEN=1024) :: conf_str
    INTEGER, EXTERNAL :: global_kpoint_index
    REAL(8), PARAMETER :: spglib_tol=1e-4
    REAL(DP), ALLOCATABLE :: r_loc(:)
    REAL(DP), ALLOCATABLE :: m_loc(:,:), initial_magn(:,:)
    INTEGER, ALLOCATABLE :: nat_of_type(:)

    CALL sirius_start_timer("setup_sirius")

    IF (do_cutoff_2D.OR.do_comp_esm) THEN
      use_veff_callback = .TRUE.
    ENDIF

    ALLOCATE(atom_type(nsp))

    DO iat = 1, nsp
      atom_type(iat)%label = TRIM(ADJUSTL(atm(iat)))

      lmax_beta = 0
      DO i = 1, upf(iat)%nbeta
        lmax_beta = MAX(lmax_beta, upf(iat)%lll(i))
      ENDDO
      atom_type(iat)%lmax = lmax_beta

      ! set the atomic radial functions
      IF (upf(iat)%has_so) THEN
        atom_type(iat)%num_chi = 0
        DO iwf = 1, upf(iat)%nwfc
          IF (upf(iat)%jchi(iwf) > upf(iat)%lchi(iwf)) THEN
            atom_type(iat)%num_chi = atom_type(iat)%num_chi + 1
          END IF
        END DO
      ELSE
        atom_type(iat)%num_chi = upf(iat)%nwfc
      END IF
      ALLOCATE(atom_type(iat)%l_chi( atom_type(iat)%num_chi) )
      ALLOCATE(atom_type(iat)%n_chi( atom_type(iat)%num_chi) )
      ALLOCATE(atom_type(iat)%chi( msh(iat), atom_type(iat)%num_chi) )
      ALLOCATE(atom_type(iat)%idx_chi( 2, atom_type(iat)%num_chi) )
      ALLOCATE(atom_type(iat)%occ( atom_type(iat)%num_chi) )
      ALLOCATE(atom_type(iat)%qpw( ngm, nh(iat) * (1 + nh(iat)) / 2 ))

      iwf = 1
      j = 1
      DO WHILE(iwf <= upf(iat)%nwfc)
        ! orbital quantum number
        atom_type(iat)%l_chi(j) = upf(iat)%lchi(iwf)
        ! principal quantum number
        IF (ALLOCATED(upf(iat)%nchi)) THEN
          atom_type(iat)%n_chi(j) = upf(iat)%nchi(iwf)
        ELSE
          atom_type(iat)%n_chi(j) = -1
        ENDIF
        ! indices of radial functions
        IF (upf(iat)%has_so) THEN
          IF (atom_type(iat)%l_chi(j) == 0) THEN
            atom_type(iat)%idx_chi(1, j) = iwf
            atom_type(iat)%idx_chi(2, j) = -1
            iwf = iwf + 1
          ELSE
            atom_type(iat)%idx_chi(1, j) = iwf
            atom_type(iat)%idx_chi(2, j) = iwf + 1
            iwf = iwf + 2
          ENDIF
        ELSE
          atom_type(iat)%idx_chi(1, j) = iwf
          atom_type(iat)%idx_chi(2, j) = -1
          iwf = iwf + 1
        END IF
        atom_type(iat)%chi(:, j) = upf(iat)%chi(1:msh(iat), atom_type(iat)%idx_chi(1, j))
        atom_type(iat)%occ(j) = upf(iat)%oc(atom_type(iat)%idx_chi(1, j))
        IF (atom_type(iat)%idx_chi(2, j) .NE. -1) THEN
          atom_type(iat)%chi(:, j) = 0.5d0 * atom_type(iat)%chi(:, j) + &
                &0.5d0 * upf(iat)%chi(1:msh(iat), atom_type(iat)%idx_chi(2, j))
          atom_type(iat)%occ(j) = atom_type(iat)%occ(j) + upf(iat)%oc(atom_type(iat)%idx_chi(2, j))
        END IF
        j = j + 1
      END DO
    ENDDO
    !
    ! create context of simulation
    CALL sirius_create_context(intra_image_comm, sctx, fcomm_k=inter_pool_comm, fcomm_band=intra_pool_comm)
    ! create initial configuration dictionary in JSON
    WRITE(conf_str, 10)diago_david_ndim, mixing_beta, nmix
    10 FORMAT('{"parameters"       : {"electronic_structure_method" : "pseudopotential", "use_scf_correction" : true}, &
               &"iterative_solver" : {"subspace_size" : ',I4,'}, &
               &"mixer"            : {"beta"        : ', F12.6, ',&
               &                      "max_history" : ', I4, ', &
               &                      "use_hartree" : true}}')
    ! set initial parameters
    CALL sirius_import_parameters(sctx, conf_str)
    ! set default verbosity
    !CALL sirius_set_parameters(sctx, verbosity=MIN(1, iverbosity))
    ! import config file
    CALL sirius_import_parameters(sctx, TRIM(ADJUSTL(sirius_cfg)))
    !
    CALL sirius_get_parameters(sctx, electronic_structure_method=conf_str)
    !
    IF (TRIM(ADJUSTL(conf_str)) .EQ. "pseudopotential") THEN
      sirius_pwpp = .TRUE.
    ELSE
      sirius_pwpp = .FALSE.
    END IF
    !
    ! derive the number of magnetic dimensions
    nmagd = 0
    IF (nspin.EQ.2) THEN
      nmagd = 1
    ENDIF
    IF (lspinorb.OR.noncolin) THEN
      nmagd = 3
    ENDIF
    ! copy dimensions of dense FFT grid
    dims(1) = dfftp%nr1
    dims(2) = dfftp%nr2
    dims(3) = dfftp%nr3
    ! set basic parameters
    ! set |G| cutoff of the dense FFT grid: convert from G^2/2 Rydbergs to |G| in [a.u.^-1]
    ! set |G+k| cutoff for the wave-functions: onvert from |G+k|^2/2 Rydbergs to |G+k| in [a.u.^-1]
    ! use symmetrization either on SIRIUS or QE side
    CALL sirius_set_parameters(sctx, num_bands=nbnd, num_mag_dims=nmagd, gamma_point=gamma_only,&
      &use_symmetry=(use_sirius_scf.OR.use_sirius_nlcg).AND..NOT.nosym, so_correction=lspinorb,&
      &pw_cutoff=SQRT(ecutrho), gk_cutoff=SQRT(ecutwfc),&
      &hubbard_correction=lda_plus_U, hubbard_correction_kind=lda_plus_u_kind,&
      &hubbard_orbitals=TRIM(ADJUSTL(Hubbard_projectors)), fft_grid_size=dims, spglib_tol=spglib_tol)

    ! degauss is converted to Ha units
    CALL sirius_set_parameters(sctx, smearing_width=degauss/2.d0)
    SELECT CASE(ngauss)
      CASE (0)
        CALL sirius_set_parameters(sctx, smearing="gaussian")
      CASE (1)
        CALL sirius_set_parameters(sctx, smearing="methfessel_paxton")
      CASE(-1)
        CALL sirius_set_parameters(sctx, smearing="cold")
      CASE(-99)
        CALL sirius_set_parameters(sctx, smearing="fermi_dirac")
    END SELECT
    !
    num_ranks_k = nproc_image / npool
    i = SQRT(DBLE(num_ranks_k) + 1d-10)
    IF (i * i .NE. num_ranks_k) THEN
      dims(1) = 1
      dims(2) = num_ranks_k
    ELSE
      dims(1) = i
      dims(2) = i
    ENDIF
    CALL sirius_set_mpi_grid_dims(sctx, 2, dims(1))
    !
    !IF (diago_full_acc) THEN
    !  CALL sirius_set_parameters(sctx, iter_solver_tol_empty=0.d0)
    !ELSE
    !  CALL sirius_set_parameters(sctx, iter_solver_tol_empty=1d-5)
    !ENDIF
    !
    CALL sirius_set_callback_function(sctx, "beta_ri", C_FUNLOC(calc_beta_radial_integrals))
    CALL sirius_set_callback_function(sctx, "beta_ri_djl", C_FUNLOC(calc_beta_dj_radial_integrals))
    CALL sirius_set_callback_function(sctx, "aug_ri", C_FUNLOC(calc_aug_radial_integrals))
    CALL sirius_set_callback_function(sctx, "aug_ri_djl", C_FUNLOC(calc_aug_dj_radial_integrals))
    CALL sirius_set_callback_function(sctx, "vloc_ri", C_FUNLOC(calc_vloc_radial_integrals))
    CALL sirius_set_callback_function(sctx, "vloc_ri_djl", C_FUNLOC(calc_vloc_dj_radial_integrals))
    CALL sirius_set_callback_function(sctx, "rhoc_ri", C_FUNLOC(calc_rhoc_radial_integrals))
    CALL sirius_set_callback_function(sctx, "rhoc_ri_djl", C_FUNLOC(calc_rhoc_dj_radial_integrals))
    CALL sirius_set_callback_function(sctx, "ps_rho_ri", C_FUNLOC(calc_ps_rho_radial_integrals))
    CALL sirius_set_callback_function(sctx, "ps_atomic_wf_ri", C_FUNLOC(calc_atomic_wfc_radial_integrals))
    CALL sirius_set_callback_function(sctx, "ps_atomic_wf_ri_djl", C_FUNLOC(calc_atomic_wfc_djl_radial_integrals))
    IF (use_veff_callback) THEN
      CALL sirius_set_callback_function(sctx, "veff", C_FUNLOC(calc_veff))
    ENDIF
    !
    ! set lattice vectors of the unit cell (length is in [a.u.])
    a1(:) = at(:, 1) * alat
    a2(:) = at(:, 2) * alat
    a3(:) = at(:, 3) * alat
    CALL sirius_set_lattice_vectors(sctx, a1(1), a2(1), a3(1))
    !
    vlat(:, 1) = a1(:)
    vlat(:, 2) = a2(:)
    vlat(:, 3) = a3(:)
    ! get the inverse of Bravais lattice vectors
    CALL invert_mtrx(vlat, vlat_inv)
    ! get the inverse of reciprocal lattice vectors
    CALL invert_mtrx(bg, bg_inv)
    ! get MPI rank associated with the distribution of k-points
    CALL mpi_comm_rank(inter_pool_comm, rank, ierr)
    !
    IF (ALLOCATED(kpoint_index_map)) DEALLOCATE(kpoint_index_map)
    ALLOCATE(kpoint_index_map(2, nkstot))
    kpoint_index_map = 0
    DO i = 1, nks
      kpoint_index_map(1, global_kpoint_index(nkstot, i)) = rank
      kpoint_index_map(2, global_kpoint_index(nkstot, i)) = i
    END DO
    !
    CALL mpi_allreduce(MPI_IN_PLACE, kpoint_index_map, 2 * nkstot, MPI_INT, MPI_SUM, inter_pool_comm, ierr)
    !
    IF (sirius_pwpp) THEN

      ! initialize atom types
      DO iat = 1, nsp

        ! add new atom type
         CALL sirius_add_atom_type(sctx, TRIM(atom_type(iat)%label), &
              & zn=NINT(zv(iat)+0.001d0), &
              & mass=amass(iat), &
              & spin_orbit=upf(iat)%has_so)

        ! set radial grid
        CALL sirius_set_atom_type_radial_grid(sctx, TRIM(atom_type(iat)%label), upf(iat)%mesh, upf(iat)%r)

        ! set beta-projectors
        DO i = 1, upf(iat)%nbeta
          l = upf(iat)%lll(i);
          IF (upf(iat)%has_so) THEN
            IF (upf(iat)%jjj(i) .LE. upf(iat)%lll(i)) THEN
              l = - upf(iat)%lll(i)
            ENDIF
          ENDIF
          CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "beta", &
               & upf(iat)%beta(1:upf(iat)%kbeta(i), i), upf(iat)%kbeta(i), l=l)
        ENDDO

        ! set the atomic radial functions
        DO j = 1, atom_type(iat)%num_chi
          CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_atomic_wf", &
               &atom_type(iat)%chi(:, j), msh(iat), l=atom_type(iat)%l_chi(j), occ=atom_type(iat)%occ(j), &
               &n=atom_type(iat)%n_chi(j))
        ENDDO

        ! QE input allow two different notations for entering the hubbard onsite interaction because 
        ! there is a bug in QE that is not fixed. 
        ! I do not set the hubbard properties right away because the Hubbard_U(iat) is not set 
        ! when the V notation is also used for onsite interaction

        IF (is_hubbard(iat)) THEN
           ! they use the second notation for onsite. I take care of this case later on
           IF (Hubbard_U(iat) .NE. 0.0) THEN
              CALL sirius_set_atom_type_hubbard(sctx, TRIM(atom_type(iat)%label), &
                   & l=Hubbard_l(iat), n=Hubbard_n(iat), occ=Hubbard_occ(iat, 1), &
                   & U=Hubbard_U(iat) / 2.0, J=Hubbard_J(1,iat) / 2.0, &
                   & alpha=Hubbard_alpha(iat) / 2.0, beta=Hubbard_beta(iat) / 2.0, &
                   & J0=Hubbard_J0(iat) / 2.0)
           ENDIF
        ENDIF
        ALLOCATE(dion(upf(iat)%nbeta, upf(iat)%nbeta))
        ! convert to hartree
        DO i = 1, upf(iat)%nbeta
          DO j = 1, upf(iat)%nbeta
            dion(i, j) = upf(iat)%dion(i, j) / 2.d0
          END DO
        END DO
        ! sed d^{ion}_{i,j}
        CALL sirius_set_atom_type_dion(sctx, TRIM(atom_type(iat)%label), upf(iat)%nbeta, dion(1, 1))
        DEALLOCATE(dion)

        ! get lmax_beta for this atom type
        lmax_beta = -1
        DO i = 1, upf(iat)%nbeta
          lmax_beta = MAX(lmax_beta, upf(iat)%lll(i))
        ENDDO

        ! set radial function of augmentation charge
        IF (upf(iat)%tvanp) THEN
          !do l = 0, upf(iat)%nqlc - 1
          DO l = 0, 2 * lmax_beta
            DO i = 1, upf(iat)%nbeta
              DO j = i, upf(iat)%nbeta
                ijv = j * (j - 1) / 2 + i
                CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "q_aug",&
                                                         &upf(iat)%qfuncl(1:upf(iat)%kkbeta, ijv, l), upf(iat)%kkbeta,&
                                                         &l=l, idxrf1=i, idxrf2=j)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        IF (upf(iat)%tpawp) THEN
          DO i = 1, upf(iat)%nbeta
            CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ae_paw_wf",&
                                                     &upf(iat)%aewfc(1:upf(iat)%paw%iraug,i), upf(iat)%paw%iraug)
            CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_paw_wf",&
                                                     &upf(iat)%pswfc(1:upf(iat)%paw%iraug,i), upf(iat)%paw%iraug)
          ENDDO
          CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ae_paw_core",&
                                                   &upf(iat)%paw%ae_rho_atc, upf(iat)%mesh)

          CALL sirius_set_atom_type_paw(sctx, TRIM(atom_type(iat)%label), upf(iat)%paw%core_energy / 2,&
                                       &upf(iat)%paw%oc, upf(iat)%nbeta)
        ENDIF

        ! set non-linear core correction
        IF (.TRUE.) THEN
          ALLOCATE(vloc(upf(iat)%mesh))
          vloc = 0.d0
          IF (ALLOCATED(upf(iat)%rho_atc)) THEN
            DO i = 1, msh(iat)
              vloc(i) = upf(iat)%rho_atc(i)
            ENDDO
          ENDIF
          CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_rho_core",&
                                                   &vloc, upf(iat)%mesh)
          DEALLOCATE(vloc)
        ENDIF

        ! set total charge density of a free atom (to compute initial rho(r))
        CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "ps_rho_total",&
                                                 &upf(iat)%rho_at, upf(iat)%mesh)

        ! the hack is done in Modules/readpp.f90
        IF (.TRUE.) THEN
          ALLOCATE(vloc(upf(iat)%mesh))
          DO i = 1, msh(iat)
            vloc(i) = upf(iat)%vloc(i)
          ENDDO
          ! convert to Hartree
          vloc = vloc / 2.d0
          ! add a correct tail
          DO i = msh(iat) + 1, upf(iat)%mesh
            vloc(i) = -zv(iat) / upf(iat)%r(i)
          ENDDO
          ! set local part of pseudo-potential
          CALL sirius_add_atom_type_radial_function(sctx, TRIM(atom_type(iat)%label), "vloc", vloc, upf(iat)%mesh)
          DEALLOCATE(vloc)
        ENDIF
      ENDDO ! iat

      IF (lda_plus_U) THEN
        DO ia = 1, nat
          !
          iat = ityp(ia)
          !
          IF (lda_plus_u_kind .EQ. 2) THEN
            !
            IF (ldim_u(iat).GT.0) THEN
              !
              DO viz = 1, neighood(ia)%num_neigh
                atom_pair(1) = ia
                ia2 = neighood(ia)%neigh(viz)
                atom_pair(2) = at_sc(ia2)%at
                iat2 = ityp(atom_pair(2))
                n_pair(1) = Hubbard_n(iat)
                n_pair(2) = Hubbard_n(iat2)
                l_pair(1) = Hubbard_l(iat)
                l_pair(2) = Hubbard_l(iat2)
                IF ((ia .EQ. ia2) .AND. (n_pair(1) .EQ. n_pair(2)) &
                & .AND. (l_pair(1) .EQ. l_pair(2)) .AND. (Hubbard_U(ia) .EQ. 0.0)) THEN
                  IF (hubbard_occ(iat,1)<0.0d0) CALL determine_hubbard_occ(iat, 1)
                  ! it is a clumsy notation as hubbard onsite correction has two different input notations.
                  CALL sirius_set_atom_type_hubbard(sctx, &
                          & TRIM(atom_type(iat)%label), &
                          & l=l_pair(1), n=n_pair(1), occ=hubbard_occ(iat, 1), &
                          & U=Hubbard_V(ia, ia2, 1) / 2.0, J=0.0D0, &
                          & alpha=0.0D0, beta=0.0D0, &
                          & J0=0.0D0)
                ELSE
                  if (ia /= ia2) then
                        ! standard-standard term in QE language
                        CALL sirius_add_hubbard_atom_pair(sctx, atom_pair(1:2), at_sc(ia2)%n(1:3), &
                                                          n_pair(1:2), l_pair(1:2), Hubbard_V(ia,ia2,1) / 2.0)
                  END IF
                ENDIF
                !
              END DO
              !
            END IF
            !
          END IF
          !=  IF (ldim_u(nt1).GT.0) THEN
          !=    DO viz = 1, neighood(iat)%num_neigh
          !=      atom_pair(1) = iat - 1
          !=      ia2 = neighood(iat)%neigh(viz)
          !=      atom_pair(2) = at_sc(ia2)%at - 1
          !=      nt2 = ityp(atom_pair(2) + 1)
          !=      !     sirius does not make any distinction
          !=      ! between orbitals contributing to the U
          !=      ! correction and orbitals with U = 0.
          !=      ! Add them to the list of interacting
          !=      ! orbitals if (V \neq 0) but (U = 0)
          !=      ! terms contributing to the standard-standard term in QE language
          !=      n_pair(1) = set_hubbard_n(upf(nt1)%psd)
          !=      n_pair(2) = set_hubbard_n(upf(nt2)%psd)
          !=      l_pair(1) = Hubbard_l(nt1)
          !=      l_pair(2) = Hubbard_l(nt2)
          !=      V_ = Hubbard_V(iat,ia2,1) * RYTOEV
          !=      if (iat /= ia2) then
          !=        call sirius_add_hubbard_atom_pair(sctx, &
          !=                atom_pair, &
          !=                at_sc(ia2)%n, &
          !=                n_pair, &
          !=                l_pair, &
          !=                V_)
          !=        ! terms contributing to the standard-background term in QE language
          !=        if (Abs(Hubbard_V(iat,ia2,2)) > 1e-8) then
          !=          n2 = set_hubbard_n(upf(nt2)%psd)
          !=          DO iwf = 1, upf(nt2)%nwfc
          !=            if (n2 /= upf(nt2)%nchi(iwf)) then
          !=              n_pair(2) = upf(nt2)%nchi(iwf)
          !=              l_pair(2) = upf(nt2)%lchi(iwf)
          !=              V_ = Hubbard_V(iat,ia2,2)  * RYTOEV
          !=              call sirius_add_hubbard_atom_pair(sctx, &
          !=                      atom_pair, &
          !=                      at_sc(ia2)%n, &
          !=                      n_pair, &
          !=                      l_pair, &
          !=                      V_)
          !=            end if
          !=          end DO
          !=         end if
          !=         ! terms contributing to the background-background  term in QE language
          !=         if (Hubbard_V(iat,ia2,3) > 1e-8) then
          !=           n1 = set_hubbard_n(upf(nt1)%psd)
          !=           n2 = set_hubbard_n(upf(nt2)%psd)
          !=           DO iwf1 = 1, upf(nt1)%nwfc
          !=             DO iwf2 = 1, upf(nt2)%nwfc
          !=               if ((n1 /= upf(nt1)%nchi(iwf)) .AND. (n2 /= upf(nt2)%nchi(iwf))) then
          !=                 n_pair(1) = upf(nt1)%nchi(iwf)
          !=                 n_pair(2) = upf(nt2)%nchi(iwf)
          !=                 l_pair(1) = upf(nt1)%lchi(iwf)
          !=                 l_pair(2) = upf(nt2)%lchi(iwf)
          !=                 V_ = Hubbard_V(iat, ia2, 3)  * RYTOEV
          !=                 call sirius_add_hubbard_atom_pair(sctx, &
          !=                         atom_pair, &
          !=                         at_sc(ia2)%n, &
          !=                         n_pair, &
          !=                         l_pair, &
          !=                         V_)
          !=               end if
          !=             end DO
          !=           end do
          !=         END if
          !=       else
          !=         CALL sirius_set_atom_type_hubbard(sctx, &
          !=                 & TRIM(atom_type(nt1)%label), &
          !=                 & Hubbard_l(nt1), &
          !=                 & set_hubbard_n(upf(nt1)%psd), &
          !=                 & hubbard_occ ( upf(nt1)%psd ), &
          !=                 & Hubbard_V(iat,ia2,1) * RYTOEV, &
          !=                 & 0.0d0, &
          !=                 & 0.0d0, &
          !=                 & 0.0d0, &
          !=                 & 0.0d0)
          !=       end if
          !=     END DO
          !=   end if
          != else
          !=   !IF (is_hubbard(iat)) THEN
          !=   !  nt1 = ityp(iat)
          !=   !  CALL sirius_set_atom_type_hubbard(sctx, &
          !=   !          & TRIM(atom_type(nt1)%label), &
          !=   !          & Hubbard_l(nt1), &
          !=   !          & set_hubbard_n(upf(nt1)%psd), &
          !=   !             & hubbard_occ ( upf(nt1)%psd ), &
          !=   !          & Hubbard_U(nt1) * RYTOEV, &
          !=   !          & Hubbard_J(1,nt1) * RYTOEV, &
          !=   !          & Hubbard_alpha(nt1) * RYTOEV, &
          !=   !          & Hubbard_beta(nt1) * RYTOEV, &
          !=   !          & Hubbard_J0(nt1) * RYTOEV)
          !=   !ENDIF
          != END IF
        END DO ! ia
      END IF
    ELSE
      ! initialize atom types with FP-LAPW species
      DO iat = 1, nsp
        ! add new atom type
         CALL sirius_add_atom_type(sctx, TRIM(atom_type(iat)%label), fname=TRIM(atom_type(iat)%label)//'.json')
      END DO
    END IF
    !
    ! compute initial magnetization for each atom type
    ALLOCATE(initial_magn(3, nsp))
    initial_magn = 0.d0
    !ALLOCATE(nat_of_type(nsp))
    !nat_of_type = 0
    !
    ! This pice of code is intended to compute integral atomic moments. The problem is that
    ! get_locals() function requires some allocated and initialized arrays. It works as expected
    ! with pw.x and doesn't work with hp.x; Attempt to allocate and initialize those arrays
    ! resulted in a crash in another place.
    !
    !!IF ( nspin .NE. 1 ) THEN
    !!  ALLOCATE( r_loc(nat), m_loc(nspin-1, nat) )
    !!  CALL get_locals( r_loc, m_loc, rho%of_r )
    !!  DO ia = 1, nat
    !!    IF (noncolin) THEN
    !!      initial_magn(:, ityp(ia)) = initial_magn(:, ityp(ia)) + m_loc(:, ia)
    !!    ELSE
    !!      initial_magn(3, ityp(ia)) = initial_magn(3, ityp(ia)) + m_loc(1, ia)
    !!    ENDIF
    !!    nat_of_type(ityp(ia)) = nat_of_type(ityp(ia)) + 1
    !!  ENDDO
    !!  DO iat = 1, nsp
    !!    initial_magn(:, iat) = initial_magn(:, iat) / nat_of_type(iat)
    !!    IF (SUM(ABS(initial_magn(:, iat))) .LT. 1e-6) THEN
    !!      initial_magn(:, iat) = 0.d0
    !!    END IF
    !!  ENDDO
    !!  DEALLOCATE(r_loc, m_loc)
    !!END IF

    ! Fallback solution: compute magentic moments on atoms in an easy way. They will be only used
    ! to determine the magentic symmetry and not as a starting guess for magnetization.
    ! Starting magnetization will be set later by the call to put_density_to_sirius() using
    ! QE values for density and magnetisation.
    DO iat = 1, nsp
      IF (noncolin) THEN
        initial_magn(1, iat) = zv(iat) * starting_magnetization(iat) * SIN(angle1(iat)) * COS(angle2(iat))
        initial_magn(2, iat) = zv(iat) * starting_magnetization(iat) * SIN(angle1(iat)) * SIN(angle2(iat))
        initial_magn(3, iat) = zv(iat) * starting_magnetization(iat) * COS(angle1(iat))
      ELSE
        initial_magn(3, iat) = zv(iat) * starting_magnetization(iat)
      ENDIF
    ENDDO
    !
    ! add atoms to the unit cell
    ! WARNING: sirius accepts only fractional coordinates
    DO ia = 1, nat
      iat = ityp(ia)
      ! Cartesian coordinates
      v1(:) = tau(:, ia) * alat
      ! fractional coordinates
      v1(:) = MATMUL(vlat_inv, v1)
      CALL sirius_add_atom(sctx, TRIM(atom_type(iat)%label), v1, initial_magn(:, iat))
    ENDDO
    !
    !DEALLOCATE(initial_magn, nat_of_type)
    DEALLOCATE(initial_magn)
    !
    CALL put_xc_functional_to_sirius()
    !
    WRITE(*,*)''
    WRITE(*,*)'=========================================='
    WRITE(*,*)'* initializing SIRIUS simulation context *'
    WRITE(*,*)'=========================================='
    !
    ! initialize global variables/indices/arrays/etc. of the simulation
    CALL sirius_initialize_context(sctx)
    !
    IF (sirius_pwpp) THEN
      DO iat = 1, nsp
        CALL sirius_get_num_beta_projectors(sctx, TRIM(atom_type(iat)%label), atom_type(iat)%num_beta_projectors)
        IF (atom_type(iat)%num_beta_projectors .NE. nh(iat)) THEN
          STOP "different number of beta-projectors"
        END IF
      ENDDO
    END IF
    !
    CALL sirius_get_parameters(sctx, num_sym_op=nsymop)
    IF (nsymop .NE. nsym) THEN
      WRITE(*,*)
      WRITE(*,'("WARNING! Different number of symmetry operations: ", I4, " QE ", I4," spglib")')nsym, nsymop
      WRITE(*,*)
    ENDIF
    !
    !IF (mpime.eq.0) THEN
    !  CALL sirius_dump_runtime_setup(sctx, "setup.json")
    !ENDIF
    !
    ! get number of g-vectors of the dense fft grid
    CALL sirius_get_num_gvec(sctx, num_gvec)
    !
    IF (.NOT.((num_gvec .EQ. ngm_g) .OR. (num_gvec * 2 - 1 .EQ. ngm_g))) THEN
      WRITE(*,*)
      WRITE(*,'("Error: wrong number of G-vectors; QE: ", I6, ", SIRIUS: ", I6)')ngm_g, num_gvec
      WRITE(*,*)
      STOP
    END IF
    !
    ! create k-point set
    ! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice and
    !          without x2 multiplication for the lsda case
    CALL sirius_create_kset(sctx, num_kpoints, kpoints, wkpoints, .FALSE., ks_handler)
    !
    ! create ground-state class
    CALL sirius_create_ground_state(ks_handler, gs_handler)
    CALL put_density_to_sirius()
    IF (okpaw) THEN
      CALL put_density_matrix_to_sirius()
      CALL sirius_generate_density(gs_handler, paw_only=.TRUE.)
    ENDIF
    CALL sirius_generate_effective_potential(gs_handler)
    !
    CALL sirius_stop_timer("setup_sirius")
    !
  END SUBROUTINE setup_sirius
  !
  !-------------------------------------------------------------------------
  SUBROUTINE update_sirius()
    !-----------------------------------------------------------------------
    !! Update parameters after the change in lattice or atomic positions
    !
    USE cell_base, ONLY : alat, at, bg
    USE ions_base, ONLY : tau, nat
    !
    IMPLICIT NONE
    !
    REAL(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), tmp
    INTEGER ia
    ! set lattice vectors of the unit cell (length is in [a.u.])
    a1(:) = at(:, 1) * alat
    a2(:) = at(:, 2) * alat
    a3(:) = at(:, 3) * alat
    CALL sirius_set_lattice_vectors(sctx, a1(1), a2(1), a3(1))
    !
    vlat(:, 1) = a1(:)
    vlat(:, 2) = a2(:)
    vlat(:, 3) = a3(:)
    ! get the inverse of Bravais lattice vectors
    CALL invert_mtrx(vlat, vlat_inv)
    ! get the inverse of reciprocal lattice vectors
    CALL invert_mtrx(bg, bg_inv)
    !
    DO ia = 1, nat
      v1(:) = tau(:, ia) * alat
      ! fractional coordinates
      v1(:) = MATMUL(vlat_inv, v1)
      CALL sirius_set_atom_position(sctx, ia, v1(1))
    ENDDO
    !
    CALL sirius_update_ground_state(gs_handler)
    !
  END SUBROUTINE update_sirius
  !
  !-------------------------------------------------------------------------
  SUBROUTINE clear_sirius()
    !-----------------------------------------------------------------------
    !! Clear all allocated objects
    !
    USE ions_base, ONLY : nsp
    !
    IMPLICIT NONE
    !
    INTEGER iat
    CALL sirius_free_handler(gs_handler)
    CALL sirius_free_handler(ks_handler)
    CALL sirius_free_handler(sctx)
    !
    IF (ALLOCATED(atom_type)) THEN
      DO iat = 1, nsp
        IF (ALLOCATED(atom_type(iat)%qpw)) DEALLOCATE(atom_type(iat)%qpw)
        IF (ALLOCATED(atom_type(iat)%l_chi)) DEALLOCATE(atom_type(iat)%l_chi)
        IF (ALLOCATED(atom_type(iat)%n_chi)) DEALLOCATE(atom_type(iat)%n_chi)
        IF (ALLOCATED(atom_type(iat)%chi)) DEALLOCATE(atom_type(iat)%chi)
        IF (ALLOCATED(atom_type(iat)%idx_chi)) DEALLOCATE(atom_type(iat)%idx_chi)
        IF (ALLOCATED(atom_type(iat)%occ)) DEALLOCATE(atom_type(iat)%occ)
      ENDDO
      DEALLOCATE(atom_type)
    ENDIF
    !
  END SUBROUTINE
  !
  !-------------------------------------------------------------------------
  SUBROUTINE invert_mtrx(vlat, vlat_inv)
    !-----------------------------------------------------------------------
    !! Auxiliary function to invert 3x3 matrix
    !
    IMPLICIT NONE
    !
    REAL(8), INTENT(in) :: vlat(3,3)
    REAL(8), INTENT(out) :: vlat_inv(3, 3)
    REAL(8) d1
    !
    d1 = vlat(1,2)*vlat(2,3)*vlat(3,1)-vlat(1,3)*vlat(2,2)*vlat(3,1)+vlat(1,3)*vlat(2,1)*vlat(3,2) &
    &   -vlat(1,1)*vlat(2,3)*vlat(3,2)+vlat(1,1)*vlat(2,2)*vlat(3,3)-vlat(1,2)*vlat(2,1)*vlat(3,3)
    d1 = 1.d0 / d1
    !
    vlat_inv(1,1)=(vlat(2,2)*vlat(3,3)-vlat(2,3)*vlat(3,2))*d1
    vlat_inv(1,2)=(vlat(1,3)*vlat(3,2)-vlat(1,2)*vlat(3,3))*d1
    vlat_inv(1,3)=(vlat(1,2)*vlat(2,3)-vlat(1,3)*vlat(2,2))*d1
    vlat_inv(2,1)=(vlat(2,3)*vlat(3,1)-vlat(2,1)*vlat(3,3))*d1
    vlat_inv(2,2)=(vlat(1,1)*vlat(3,3)-vlat(1,3)*vlat(3,1))*d1
    vlat_inv(2,3)=(vlat(1,3)*vlat(2,1)-vlat(1,1)*vlat(2,3))*d1
    vlat_inv(3,1)=(vlat(2,1)*vlat(3,2)-vlat(2,2)*vlat(3,1))*d1
    vlat_inv(3,2)=(vlat(1,2)*vlat(3,1)-vlat(1,1)*vlat(3,2))*d1
    vlat_inv(3,3)=(vlat(1,1)*vlat(2,2)-vlat(1,2)*vlat(2,1))*d1
    !
  END SUBROUTINE
  !
  !-------------------------------------------------------------------------
  SUBROUTINE get_band_energies_from_sirius()
    !-----------------------------------------------------------------------
    !! Get KS energies from SIRIUS.
    !
    USE wvfct,    ONLY : nbnd, et
    USE klist,    ONLY : nkstot, nks
    USE lsda_mod, ONLY : nspin
    !
    IMPLICIT NONE
    !
    REAL(8), ALLOCATABLE :: band_e(:,:)
    INTEGER :: ik, nk, nb, nfv
    !
    ALLOCATE(band_e(nbnd, nkstot))
    ! get band energies
    IF (nspin.NE.2) THEN
      ! non-magnetic or non-collinear case
      DO ik = 1, nkstot
        CALL sirius_get_band_energies(ks_handler, ik, 1, band_e(:, ik))
      END DO
    ELSE
      ! collinear magnetic case
      nk = nkstot / 2
      ! get band energies
      DO ik = 1, nk
        CALL sirius_get_band_energies(ks_handler, ik, 1, band_e(:, ik))
        CALL sirius_get_band_energies(ks_handler, ik, 2, band_e(:, nk + ik))
      END DO
    ENDIF
    ! convert to Ry
    DO ik = 1, nkstot
      et(:, ik) = 2.d0 * band_e(:, ik)
    ENDDO
    !
    DEALLOCATE(band_e)
    !
  END SUBROUTINE get_band_energies_from_sirius
  !
  !-------------------------------------------------------------------------
  SUBROUTINE put_band_occupancies_to_sirius
    !-------------------------------------------------------------------------
    !! Put KS occupancies to SIRIUS
    !
    USE wvfct,    ONLY : nbnd, wg
    USE klist,    ONLY : nkstot, nks, wk
    USE lsda_mod, ONLY : nspin
    USE mp_pools, ONLY : inter_pool_comm
    USE parallel_include
    !
    IMPLICIT NONE
    !
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    REAL(8), ALLOCATABLE :: bnd_occ(:, :)
    REAL(8) :: maxocc
    INTEGER :: ik, ierr, nk, nb
    ! compute occupancies
    ALLOCATE(bnd_occ(nbnd, nkstot))
    bnd_occ = 0.d0
    ! define a maximum band occupancy (2 in case of spin-unpolarized, 1 in case of spin-polarized)
    maxocc = 2.d0
    IF (nspin.GT.1) THEN
      maxocc = 1.d0
    ENDIF
    DO ik = 1, nks
      bnd_occ(:, global_kpoint_index(nkstot, ik)) = maxocc * wg(:, ik) / wk(ik)
    ENDDO
    CALL mpi_allreduce(MPI_IN_PLACE, bnd_occ(1, 1), nbnd * nkstot, MPI_DOUBLE, MPI_SUM, inter_pool_comm, ierr)
    !
    IF (nspin.NE.2) THEN
      ! set band occupancies
      DO ik = 1, nkstot
        CALL sirius_set_band_occupancies(ks_handler, ik, 1, bnd_occ(:, ik))
      ENDDO
    ELSE
      nk = nkstot / 2
      DO ik = 1, nk
        CALL sirius_set_band_occupancies(ks_handler, ik, 1, bnd_occ(:, ik))
        CALL sirius_set_band_occupancies(ks_handler, ik, 2, bnd_occ(:, ik + nk))
      ENDDO
    ENDIF
    !
    DEALLOCATE(bnd_occ)
    !
  END SUBROUTINE put_band_occupancies_to_sirius
  !
  !-------------------------------------------------------------------------
  SUBROUTINE get_band_occupancies_from_sirius
    !-------------------------------------------------------------------------
    !! Get KS occupancies from SIRIUS
    USE wvfct,    ONLY : nbnd, wg
    USE klist,    ONLY : nkstot, nks, wk
    USE lsda_mod, ONLY : nspin
    USE mp_pools, ONLY : inter_pool_comm
    USE parallel_include
    !
    IMPLICIT NONE
    !
    REAL(8), ALLOCATABLE :: bnd_occ(:, :)
    REAL(8) :: maxocc
    INTEGER :: ik, ierr, nk, n
    !
    ! compute occupancies
    ALLOCATE(bnd_occ(nbnd, nkstot))
    IF (nspin.NE.2) THEN
      ! set band occupancies
      DO ik = 1, nkstot
        CALL sirius_get_band_occupancies(ks_handler, ik, 1, bnd_occ(:, ik))
      ENDDO
    ELSE
      nk = nkstot / 2
      DO ik = 1, nk
        CALL sirius_get_band_occupancies(ks_handler, ik, 1, bnd_occ(:, ik))
        CALL sirius_get_band_occupancies(ks_handler, ik, 2, bnd_occ(:, ik + nk))
      ENDDO
    ENDIF
    ! define a maximum band occupancy (2 in case of spin-unpolarized, 1 in case of spin-polarized)
    maxocc = 2.d0
    IF (nspin.GT.1) THEN
      maxocc = 1.d0
    ENDIF
    DO ik = 1, nkstot
      wg(:, ik) = bnd_occ(:, ik) / maxocc * wk(ik)
    ENDDO
    !
    DEALLOCATE(bnd_occ)
    !
  END SUBROUTINE get_band_occupancies_from_sirius
  !
  !-------------------------------------------------------------------------
  SUBROUTINE get_wave_functions_from_sirius
    !-------------------------------------------------------------------------
    !! Get KS wave-functions.
    !
    USE klist,            ONLY : nkstot, nks, ngk, igk_k
    USE gvect,            ONLY : mill
    USE buffers,          ONLY : save_buffer, open_buffer
    USE io_files,         ONLY : iunwfc, nwordwfc
    USE control_flags,    ONLY : io_level
    USE bp,               ONLY : lelfield
    USE noncollin_module, ONLY : npol
    USE wvfct,            ONLY : npwx, nbnd
    USE wavefunctions,    ONLY : evc
    USE lsda_mod,         ONLY : isk, lsda
    USE mp_pools,         ONLY : inter_pool_comm
    USE parallel_include
    !
    IMPLICIT NONE
    !
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER, ALLOCATABLE :: vgl(:,:)
    INTEGER ig, ik, ik_, ik1, i, j, ispn, rank, ierr, nksmax, ikloc
    COMPLEX(8) z1
    LOGICAL exst_file,exst_mem
    !
    ! rank of communicator that distributes k-points
    CALL mpi_comm_rank(inter_pool_comm, rank, ierr)
    !
    ALLOCATE(vgl(3, npwx))
    !
    DO ik = 1, nkstot
      ik1 = MOD(ik - 1, num_kpoints) + 1
      ispn = 1
      IF (ik .GT. num_kpoints) THEN
        ispn = 2
      ENDIF
      IF (kpoint_index_map(1, ik) .EQ. rank) THEN
        ikloc = kpoint_index_map(2, ik)
        DO ig = 1, ngk(ikloc)
          vgl(:,ig) = mill(:, igk_k(ig, ikloc))
        ENDDO
        CALL sirius_get_wave_functions( ks_handler, vkl=kpoints(:, ik1), spin=ispn, num_gvec_loc=ngk(ikloc), &
                                      & gvec_loc=vgl, evec=evc, ld=npwx, num_spin_comp=npol )
        IF (nks > 1 .OR. lelfield) THEN
          CALL save_buffer ( evc, nwordwfc, iunwfc, ikloc )
        ENDIF
      ELSE
        CALL sirius_get_wave_functions( ks_handler )
      ENDIF
      !
      CALL mpi_barrier(inter_pool_comm, ierr)
      !
    ENDDO
    !
    !CALL mpi_allreduce(nks, nksmax, 1, MPI_INTEGER, MPI_MAX, inter_pool_comm, ierr)
    !
    !!CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
    !
    !DO ik = 1, nksmax
    !  IF (ik.LE.nks) THEN
    !    DO ig = 1, ngk(ik)
    !      gvl(:,ig) = mill(:, igk_k(ig, ik))
    !    ENDDO
    !    !
    !    ik_ = global_kpoint_index(nkstot, ik)
    !    ispn = isk(ik)
    !    IF (lsda.AND.ispn.EQ.2) THEN
    !      ik_ = ik_ - nkstot / 2
    !    ENDIF
    !    CALL sirius_get_wave_functions(ks_handler, ik_, ispn, ngk(ik), gvl(1, 1), evc(1, 1), npwx, npol)
    !    IF (nks > 1 .OR. lelfield) THEN
    !      CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
    !    ENDIF
    !  ELSE
    !    CALL sirius_get_wave_functions(ks_handler, -1, -1, -1, -1, z1, -1, -1)
    !  ENDIF
    !ENDDO
    !
    DEALLOCATE(vgl)
    !
  END SUBROUTINE get_wave_functions_from_sirius
  !
  !-------------------------------------------------------------------------
  SUBROUTINE put_xc_functional_to_sirius
    !-----------------------------------------------------------------------
    !! Put the information about XC potentials into SIRIUS.
    !
    USE xc_lib
    USE funct,                ONLY : get_dft_name
    !
    IMPLICIT NONE
    !
    INTEGER :: iexch, icorr, igcx, igcc, imeta, imetac
    !
    iexch  = xclib_get_id( 'LDA', 'EXCH' )
    icorr  = xclib_get_id( 'LDA', 'CORR' )
    igcx   = xclib_get_id( 'GGA', 'EXCH' )
    igcc   = xclib_get_id( 'GGA', 'CORR' )
    imeta  = xclib_get_id( 'MGGA','EXCH' )
    imetac = xclib_get_id( 'MGGA','CORR' )
    !WRITE(*,*)iexch,icorr,igcx,igcc,imeta,imetac
    !WRITE(*,*)trim(get_dft_name())
    !
    IF (imeta.NE.0.OR.imetac.NE.0) THEN
      STOP ("interface for meta-XC functional is not implemented")
    ENDIF
    !
    IF (iexch.NE.0.AND.igcx.EQ.0) THEN
      SELECT CASE(iexch)
      CASE(0)
      CASE(1)
        CALL sirius_add_xc_functional(sctx, "XC_LDA_X")
      CASE default
        STOP ("interface for this exchange functional is not implemented")
      END SELECT
    ENDIF
    !
    IF (iexch.NE.0.AND.igcx.NE.0) THEN
      SELECT CASE(igcx)
      CASE(0)
      CASE(2)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PW91")
      CASE(3)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE")
      CASE(101)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE")
      CASE(10)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE_SOL")
      CASE default
        WRITE(*,*)igcx
        STOP ("interface for this gradient exchange functional is not implemented")
      END SELECT
    ENDIF
    !
    IF (iexch.EQ.0.AND.igcx.NE.0) THEN
      SELECT CASE(igcx)
      CASE(101)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE")
      CASE(109)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PW91")
      CASE(116)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_X_PBE_SOL")
      CASE default
        WRITE(*,*)igcx
        STOP ("interface for this gradient exchange functional is not implemented")
      END SELECT
    ENDIF
    !
    IF (icorr.NE.0.AND.igcc.EQ.0) THEN
      SELECT CASE(icorr)
      CASE(0)
      CASE(1)
        CALL sirius_add_xc_functional(sctx, "XC_LDA_C_PZ")
      CASE(4)
        CALL sirius_add_xc_functional(sctx, "XC_LDA_C_PW")
      CASE default
        STOP ("interface for this correlation functional is not implemented")
      END SELECT
    ENDIF
    !
    IF (icorr.NE.0.AND.igcc.NE.0) THEN
      SELECT CASE(igcc)
      CASE(0)
      CASE(2)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PW91")
      CASE(4)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PBE")
      CASE(130)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PBE")
      CASE(8)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PBE_SOL")
      CASE default
        STOP ("interface for this gradient correlation functional is not implemented")
      END SELECT
    ENDIF
    !
    IF (icorr.EQ.0.AND.igcc.NE.0) THEN
      SELECT CASE(igcc)
      CASE(130)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PBE")
      CASE(133)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PBE_SOL")
      CASE(134)
        CALL sirius_add_xc_functional(sctx, "XC_GGA_C_PW91")
      CASE default
        STOP ("interface for this gradient correlation functional is not implemented")
      END SELECT
    ENDIF
    !
  END SUBROUTINE put_xc_functional_to_sirius
  !
  !-------------------------------------------------------------------------
  SUBROUTINE write_json()
    !-----------------------------------------------------------------------
    !! Auxiliary function to write the JSON dictionary at the end of the run.
    !
    USE ener
    USE force_mod
    USE ions_base
    !
    IMPLICIT NONE
    !
    INTEGER i
    !
    OPEN(200, file="output.json", action="write", form="formatted")
    WRITE(200,'("{")')
    WRITE(200,'("  ""energy"": {")')
    WRITE(200,'("    ""total"": ", G18.10)')etot
    WRITE(200, '("  },")')
    WRITE(200,'("  ""stress"": {")')
    WRITE(200,'("    ""total"": [")')
    WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') sigma(1, 1), sigma(1, 2), sigma(1, 3)
    WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') sigma(2, 1), sigma(2, 2), sigma(2, 3)
    WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "]")')  sigma(3, 1), sigma(3, 2), sigma(3, 3)
    WRITE(200,'("    ]")')
    WRITE(200,'("  },")')
    WRITE(200,'("  ""force"": {")')
    WRITE(200,'("    ""total"": [")')
    DO i = 1, nat
      IF (i.EQ.nat) THEN
        WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "]")') force(1, i), force(2, i), force(3, i)
      ELSE
        WRITE(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') force(1, i), force(2, i), force(3, i)
      ENDIF
    ENDDO
    WRITE(200,'("    ]")')
    WRITE(200,'("  }")')
    WRITE(200,'("}")')
    CLOSE(200)
    !
  END SUBROUTINE

#endif

END MODULE mod_sirius
