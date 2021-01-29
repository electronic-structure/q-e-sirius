MODULE mod_sirius_callbacks
USE mod_sirius_base
IMPLICIT NONE

CONTAINS

!! A callback function to compute band occupancies is SIRIUS using QE
!! TODO: merge to a single function with argments (energies in, occupancies out)
!SUBROUTINE calc_band_occupancies() BIND(C)
!USE ISO_C_BINDING
!IMPLICIT NONE
!!
!CALL get_band_energies_from_sirius()
!CALL weights()
!CALL put_band_occupancies_to_sirius()
!!
!END SUBROUTINE calc_band_occupancies

SUBROUTINE calc_ps_rho_radial_integrals(iat, nq, q, ps_rho_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE esm
USE Coul_cut_2D
USE uspp_param, ONLY : upf
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
    aux(ir) = aux(ir) * upf(iat)%rho_at(ir) / 12.566370614359172954d0
  ENDDO
  CALL simpson(nr, aux, rgrid(iat)%rab(1), ps_rho_ri(iq))
ENDDO
END SUBROUTINE calc_ps_rho_radial_integrals


! A callback function to compute radial integrals of Vloc(r)
SUBROUTINE calc_vloc_radial_integrals(iat, nq, q, vloc_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE esm
USE Coul_cut_2D
USE uspp_param, ONLY : upf
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
REAL(DP), EXTERNAL :: qe_erf
! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
nr = msh(iat)
! q-independent radial function
DO ir = 1, nr
  r = rgrid(iat)%r(ir)
  aux1(ir) = r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * qe_erf(r)
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
        aux(ir) = r * (r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * qe_erf(r))
      END DO
      IF (do_cutoff_2D .AND. rgrid(iat)%r(nr) > lz) THEN
        CALL errore('vloc_of_g','2D cutoff is smaller than pseudo cutoff radius: &
          & increase interlayer distance (or see Modules/read_pseudo.f90)',1)
      END IF
    ELSE
      DO ir = 1, nr
        r = rgrid(iat)%r(ir)
        aux(ir) = r * (r * upf(iat)%vloc(ir) + upf(iat)%zp * e2)
      END DO
    END IF
    CALL simpson(nr, aux, rgrid(iat)%rab(1), vloc_ri(iq))
  ELSE ! q > 0 case
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux(ir) = aux1(ir) * SIN(q(iq) * r) / q(iq)
    END DO
    CALL simpson(nr, aux, rgrid(iat)%rab(1), vloc_ri(iq))
    IF ((.NOT.do_comp_esm) .OR. (esm_bc .EQ. 'pbc')) THEN
      !
      !   here we re-add the analytic fourier transform of the erf function
      !
      IF (.NOT. do_cutoff_2D) THEN
        q2 = q(iq) * q(iq)
        vloc_ri(iq) = vloc_ri(iq) - upf(iat)%zp * e2 * EXP(-q2 * 0.25d0) / q2
      END IF
    END IF
  END IF
  ! convert to Ha
  vloc_ri(iq) = vloc_ri(iq) / 2.0
END DO

END SUBROUTINE calc_vloc_radial_integrals


! A callback function to compute radial integrals of Vloc(r) with
! derivatives of Bessel functions
SUBROUTINE calc_vloc_dj_radial_integrals(iat, nq, q, vloc_dj_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE esm
USE Coul_cut_2D
USE uspp_param, ONLY : upf
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
REAL(DP), EXTERNAL :: qe_erf
! number of points to the effective infinity (~10 a.u. hardcoded somewhere in the code)
nr = msh(iat)
! q-independent radial function
DO ir = 1, nr
  r = rgrid(iat)%r(ir)
  aux1(ir) = r * upf(iat)%vloc(ir) + upf(iat)%zp * e2 * qe_erf(r)
END DO
! loop over q-points
DO iq = 1, nq
  IF (q(iq) < eps8) THEN ! q=0 case
    vloc_dj_ri(iq) = 0.d0
  ELSE ! q > 0 case
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux(ir) = aux1(ir) * (SIN(q(iq) * r) / q(iq)**2 - r * COS(q(iq) * r) / q(iq))
    END DO
    CALL simpson(nr, aux, rgrid(iat)%rab(1), vloc_dj_ri(iq))
    vloc_dj_ri(iq) = vloc_dj_ri(iq) / q(iq)
    IF ((.NOT.do_comp_esm) .OR. (esm_bc .EQ. 'pbc')) THEN
      IF (.NOT. do_cutoff_2D) THEN
        q2 = q(iq) * q(iq)
        vloc_dj_ri(iq) = vloc_dj_ri(iq) - upf(iat)%zp * e2 * &
          &EXP(-q2 * 0.25d0) * (q2 + 4) / 2 / q2 / q2
      END IF
    END IF
  END IF
  ! convert to Ha
  vloc_dj_ri(iq) = vloc_dj_ri(iq) / 2.0
END DO

END SUBROUTINE calc_vloc_dj_radial_integrals


! A callback function to compute radial integrals of rho_core(r)
SUBROUTINE calc_rhoc_radial_integrals(iat, nq, q, rhoc_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE uspp_param, ONLY : upf
IMPLICIT NONE
!
INTEGER(KIND=c_int), INTENT(IN), VALUE :: iat
INTEGER(KIND=c_int), INTENT(IN), VALUE :: nq
REAL(KIND=c_double), INTENT(IN)  :: q(nq)
REAL(KIND=c_double), INTENT(OUT) :: rhoc_ri(nq)
!
REAL(DP) :: aux(rgrid(iat)%mesh), r
INTEGER :: iq, ir, nr
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
    CALL sph_bes(nr, rgrid(iat)%r(1), q(iq), 0, aux)
    DO ir = 1, nr
      r = rgrid(iat)%r(ir)
      aux(ir) = r**2 * upf(iat)%rho_atc(ir) * aux(ir)
    ENDDO
  END IF
  CALL simpson(nr, aux, rgrid(iat)%rab(1), rhoc_ri(iq))
END DO

END SUBROUTINE calc_rhoc_radial_integrals


! A callbacl function to compute radial integrals or rho_core(r) with
! the derivatives of Bessel functions
SUBROUTINE calc_rhoc_dj_radial_integrals(iat, nq, q, rhoc_dj_ri) BIND(C)
USE ISO_C_BINDING
USE kinds
USE atom
USE constants
USE uspp_param, ONLY : upf
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
        &(r * COS(q(iq) * r) / q(iq) - SIN(q(iq) * r) / q(iq)**2)
    ENDDO
    CALL simpson(nr, aux, rgrid(iat)%rab(1), rhoc_dj_ri(iq))
  END IF
END DO

END SUBROUTINE calc_rhoc_dj_radial_integrals


SUBROUTINE calc_beta_radial_integrals(iat, q, beta_ri, ld) BIND(C)
USE iso_c_binding
USE us,           ONLY : dq, tab
USE uspp_param,   ONLY : upf
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

END SUBROUTINE calc_beta_radial_integrals


SUBROUTINE calc_beta_dj_radial_integrals(iat, q, beta_ri, ld) BIND(C)
USE iso_c_binding
USE us,           ONLY : dq, tab
USE uspp_param,   ONLY : upf
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

END SUBROUTINE calc_beta_dj_radial_integrals


SUBROUTINE calc_aug_radial_integrals(iat, q, aug_ri, ld1, ld2) BIND(C)
USE iso_c_binding
USE us,           ONLY : dq, qrad
USE uspp_param,   ONLY : upf
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

END SUBROUTINE calc_aug_radial_integrals


SUBROUTINE calc_aug_dj_radial_integrals(iat, q, aug_ri, ld1, ld2) BIND(C)
USE iso_c_binding
USE us,           ONLY : dq, qrad
USE uspp_param,   ONLY : upf
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

END SUBROUTINE calc_aug_dj_radial_integrals

END MODULE mod_sirius_callbacks
