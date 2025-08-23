!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE w90_wigner
  !----------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ! PUBLIC :: w90_wigner_fourier_factor
  PUBLIC :: w90_find_wigner_seitz
  PUBLIC :: w90_wigner_allocate
  PUBLIC :: w90_wigner_deallocate
  PUBLIC :: ndegen
  PUBLIC :: irvec
  !
  INTEGER, PARAMETER :: NSEARCH = 3
  !! Number of Born-von Karman supercell lattice vectors to search for the Wigner-Seitz vector.
  !! [-NSEARCH, NSEARCH]^3 = (2 * NSEARCH + 1)^3 points are searched in total.
  REAL(DP), PARAMETER :: WS_DIST_TOL = 1.0E-6_DP
  !! Tolerance for finding degeneracies in the Wigner-Seitz selection.
  !! Distance differences below WS_DIST_TOL is considered equal.
  !
  INTEGER :: numT
  !! Number of supercell lattice vectors to search for the Wigner-Seitz vector. (2 * NSEARCH + 1)^3
  INTEGER, ALLOCATABLE :: irvec(:, :)
  !! Components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors.
  !! The results of w90_find_wigner_seitz are stored in (:, 1:nrr). The remaining elements are not used.
  INTEGER, ALLOCATABLE :: ndegen(:)
  !! Number of degeneracies. Same size as irvec.
  INTEGER, ALLOCATABLE :: Tvecs(:, :)
  !! Supercell lattice vectors (in fractional coordinates)
  REAL(DP), ALLOCATABLE :: Tvecs_cart(:, :)
  !! Supercell lattice vectors (in Cartesian coordinates)
  REAL(DP), ALLOCATABLE :: dist(:)
  !! Contains the distance squared |T - r|^2
  REAL(DP), ALLOCATABLE :: dist_rT(:)
  !! Contains the distance dot(T, r)
  REAL(DP), ALLOCATABLE :: dist_TT(:)
  !! Contains the distance squared |T|^2
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE w90_wigner_allocate(nc1, nc2, nc3)
  !-----------------------------------------------------------------------------
    !
    USE cell_base,     ONLY : at
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nc1
    !! size of the uniform k mesh
    INTEGER, INTENT(in) :: nc2
    !! size of the uniform k mesh
    INTEGER, INTENT(in) :: nc3
    !! size of the uniform k mesh
    !
    INTEGER :: i1, i2, i3
    !! Index for the supercell
    INTEGER :: iT
    !! Supercell lattice vector index
    !
    numT = (2 * NSEARCH + 1)**3
    !
    ALLOCATE(Tvecs(3, numT))
    ALLOCATE(Tvecs_cart(3, numT))
    ALLOCATE(dist(numT))
    ALLOCATE(dist_rT(numT))
    ALLOCATE(dist_TT(numT))
    !
    ALLOCATE(ndegen(20))
    ALLOCATE(irvec(3, 20))
    !
    iT = 0
    DO i1 = -NSEARCH, NSEARCH
      DO i2 = -NSEARCH, NSEARCH
        DO i3 = -NSEARCH, NSEARCH
          iT = iT + 1
          !
          ! T = (i1 * nc1, i2 * nc2, i3 * nc3) * at
          !
          Tvecs(1, iT) = i1 * nc1
          Tvecs(2, iT) = i2 * nc2
          Tvecs(3, iT) = i3 * nc3
          !
        ENDDO ! i3
      ENDDO ! i2
    ENDDO ! i1
    !
    Tvecs_cart = Tvecs
    CALL cryst_to_cart((2*NSEARCH+1)**3, Tvecs_cart, at, +1)
    !
    DO iT = 1, numT
      dist_TT(iT) = SUM((Tvecs_cart(:, iT))**2)
    ENDDO
    !
  END SUBROUTINE w90_wigner_allocate
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE w90_wigner_deallocate()
  !-----------------------------------------------------------------------------
    DEALLOCATE(irvec)
    DEALLOCATE(ndegen)
    DEALLOCATE(Tvecs)
    DEALLOCATE(Tvecs_cart)
    DEALLOCATE(dist)
    DEALLOCATE(dist_rT)
    DEALLOCATE(dist_TT)
  END SUBROUTINE w90_wigner_deallocate
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE w90_find_wigner_seitz(r0, nrr, mindist)
  !-----------------------------------------------------------------------------
  !!
  !! Calculates a grid of points that fall inside of (and on the surface of) the
  !! Wigner-Seitz supercell centered at r0 for the supercell lattice vectors
  !! nc1*a_1 + nc2*a_2 + nc3*a_3.
  !!
  !! The output is set in the argument nrr and module variables irvec and ndegen.
  !!
  !! Adapted from EPW/src/wigner.f90
  !!
  !-----------------------------------------------------------------------------
    USE kinds,         ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(in) :: r0(3)
    !! Center of the WS cell in real space (in Cartesian coordinates)
    INTEGER, INTENT(out) :: nrr
    !! number of Wigner-Seitz grid points
    !
    INTEGER :: iT
    !! Supercell lattice vector index
    INTEGER :: cnt
    !! Count for the Wigner-Seitz degeneracy
    REAL(DP) :: mindist
    !! Minimum distance
    REAL(DP) :: dist_rr
    !! |r0|^2
    INTEGER, ALLOCATABLE :: irvec_tmp(:, :)
    !! Temporary storage for irvec
    INTEGER, ALLOCATABLE :: ndegen_tmp(:)
    !! Temporary storage for ndegen
    !
    ! For given r, find the supercell lattice vectors inside the Wigner-Seitz cell.
    !
    ! We iterate over the supercell lattice vectors T = (i1 * nc1, i2 * nc2, i3 * nc3)
    ! and find T's such that T - r0 is inside the Wigner-Seitz cell.
    ! To do so, we compute dist(T) = |T - r0| for a set of T vectors
    ! in [-NSEARCH, NSEARCH]^3. T vectors with the smallest dist(T) are selected. If multiple T's
    ! give the same dist(T) (up to tolerance eps6) they are all selected.
    !
    ! nrr             : The number of unique T vectors.
    ! irvec(3, 1:nrr) : The set of all selected T vectors.
    ! ndegen(1:nrr)   : The number of T vectors with the same minimal dist(T).
    !
    ! Compute dist(T) for all T in [-NSEARCH, NSEARCH]^3
    ! dist(T) = |T - r0|^2 = |T|^2 - 2 * dot(T, r0) + r0^2
    ! T = (i1 * nc1, i2 * nc2, i3 * nc3) * at
    !
    CALL DGEMV('T', 3, numT, -2.d0, Tvecs_cart, 3, r0, 1, 0.d0, dist_rT, 1)
    dist_rr = SUM(r0**2)
    !
    dist(:) = SQRT(dist_TT + dist_rT + dist_rr)
    !
    ! Count the number of vectors T with the (same) smallest dist(T).
    ! This is the degeneracy.
    !
    mindist = MINVAL(dist)
    !
    cnt = 0
    DO iT = 1, numT
      IF (dist(iT) < mindist + WS_DIST_TOL) cnt = cnt + 1
    ENDDO
    !
    ! Increase the size of irvec and ndegen if needed.
    !
    IF (cnt > SIZE(ndegen)) THEN
      DEALLOCATE(ndegen)
      DEALLOCATE(irvec)
      ALLOCATE(ndegen(cnt))
      ALLOCATE(irvec(3, cnt))
    ENDIF
    !
    ! Find T's with the (same) smallest dist(T) and store them to irvec.
    !
    nrr = 0
    !
    DO iT = 1, numT
      IF (dist(iT) < mindist + WS_DIST_TOL) THEN
        nrr = nrr + 1
        ndegen(nrr) = cnt
        irvec(:, nrr) = Tvecs(:, iT)
      ENDIF
    ENDDO ! iT
    !
    ! This cannot happen, but check to be super safe.
    !
    IF (nrr == 0) CALL errore('w90_find_wigner_seitz', 'No Wigner-Seitz vectors found', 1)
    IF (nrr /= cnt) CALL errore('w90_find_wigner_seitz', 'Inconsistent number of Wigner-Seitz vectors', 1)
    !
  END SUBROUTINE w90_find_wigner_seitz
  !----------------------------------------------------------------------------
  !
  ! !-----------------------------------------------------------------------------
  ! SUBROUTINE w90_wigner_fourier_factor(xkk, R, fac_r)
  ! !-----------------------------------------------------------------------------
  ! !!
  ! !! Compute the Fourier transform factor using Wigner-Seitz minimal distance
  ! !! replica selection. Different Wigner-Seitz supercell lattice vectors are
  ! !! selected for each FFT grid point.
  ! !!
  ! !! f(r) = \sum_{T \in WS_r} exp(-i * k * (R + T))
  ! !! WS_r is the Wigner-Seitz cell centered at R - r.
  ! !!
  ! !-----------------------------------------------------------------------------
  !   USE kinds,            ONLY : DP
  !   USE constants,        ONLY : tpi
  !   USE fft_base,         ONLY : dffts
  !   USE fft_types,        ONLY : fft_index_to_3d
  !   USE cell_base,        ONLY : at
  !   !
  !   IMPLICIT NONE
  !   !
  !   REAL(DP), INTENT(in) :: R(3)
  !   !! Lattice vector in Cartesian coordinates
  !   REAL(DP), INTENT(in) :: xkk(3)
  !   !! k point in Cartesian coordinates
  !   COMPLEX(DP), INTENT(out) :: fac_r(dffts%nnr)
  !   !! Fourier factor
  !   !
  !   LOGICAL :: offrange
  !   !! Flag for out-of-range grid points
  !   INTEGER :: ir
  !   !! Real space FFT grid index
  !   INTEGER :: iT
  !   !! Supercell lattice vector index
  !   INTEGER :: i, j, k
  !   !! Direction index
  !   INTEGER :: nrr
  !   !! number of Wigner-Seitz grid points
  !   REAL(DP) :: kdotR
  !   !! Dot product of k and R
  !   REAL(DP) :: r0(3)
  !   !! Real space grid point inside unit cell
  !   REAL(DP) :: T(3)
  !   !! Spercell lattice vector
  !   REAL(DP) :: mindist
  !   !! Minimum distance (needed to call w90_find_wigner_seitz, not used here.)
  !   !
  !   CALL start_clock("w90_wigner")
  !   !
  !   fac_r(:) = (0.d0, 0.d0)
  !   !
  !   DO ir = 1, dffts%nnr
  !     !
  !     CALL fft_index_to_3d(ir, dffts, i, j, k, offrange)
  !     IF (offrange) CYCLE
  !     !
  !     r0(1) = REAL(i, DP) / REAL(dffts%nr1, DP)
  !     r0(2) = REAL(j, DP) / REAL(dffts%nr2, DP)
  !     r0(3) = REAL(k, DP) / REAL(dffts%nr3, DP)
  !     CALL cryst_to_cart(1, r0, at, +1)
  !     !
  !     CALL w90_find_wigner_seitz(r0 - R, nrr, mindist)
  !     !
  !     DO iT = 1, nrr
  !       !
  !       T = REAL(irvec(:, iT), DP)
  !       CALL cryst_to_cart(1, T, at, +1)
  !       !
  !       kdotR = tpi * SUM(xkk * (R + T))
  !       fac_r(ir) = fac_r(ir) + CMPLX(COS(kdotR), SIN(kdotR)) / REAL(ndegen(iT), DP)
  !       !
  !     ENDDO
  !     !
  !   ENDDO ! ir
  !   !
  !   CALL stop_clock("w90_wigner")
  !   !
  ! END SUBROUTINE w90_wigner_fourier_factor
  ! !----------------------------------------------------------------------------
  !
!----------------------------------------------------------------------------
END MODULE w90_wigner
!----------------------------------------------------------------------------
