!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE w90_hamiltonian
  !----------------------------------------------------------------------------
  !! Module for Wannier-interpolated Hamiltonian eigenvalues and eigenvectors.
  !!
  !! This module provides functionality to:
  !! 1. Read eigenvalues from Wannier90 .eig file
  !! 2. Convert to Wannier gauge using U matrix
  !! 3. Fourier transform to real space with optimal Wigner-Seitz degeneracy
  !! 4. Interpolate to arbitrary k-points
  !!
  !! Usage:
  !!   CALL w90_read(seedname)           ! First, read Wannier90 checkpoint
  !!   CALL w90_setup_hamiltonian(seedname)  ! Setup Hamiltonian interpolation
  !!   CALL w90_interpolate_ham(xk, eig, cuf)  ! Interpolate to any k-point
  !!   CALL w90_cleanup_hamiltonian()    ! Cleanup when done
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: w90_setup_hamiltonian
  PUBLIC :: w90_interpolate_ham
  PUBLIC :: w90_cleanup_hamiltonian
  PUBLIC :: nR_ham, iRlist_ham, inv_ndegen_ham, ham_R
  !
  INTEGER :: nR_ham
  !! Number of R vectors for Hamiltonian
  INTEGER, ALLOCATABLE :: iRlist_ham(:, :)
  !! List of R vectors in fractional coordinates (3, nR_ham)
  REAL(DP), ALLOCATABLE :: inv_ndegen_ham(:, :, :)
  !! Degeneracy of each R vector for each Wannier pair (nR_ham, num_wann, num_wann)
  !! Uses optimal WS: WS cell centered at R + r_m - r_n
  COMPLEX(DP), ALLOCATABLE :: ham_R(:, :, :)
  !! Real-space Hamiltonian H^W(R), shape (num_wann, num_wann, nR_ham)
  !
  INTEGER, PARAMETER :: NSEARCH = 3
  !! Number of Born-von Karman supercell lattice vectors to search
  REAL(DP), PARAMETER :: WS_DIST_TOL = 1.0E-6_DP
  !! Tolerance for finding degeneracies in the Wigner-Seitz selection
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_setup_hamiltonian(seedname)
    !--------------------------------------------------------------------------
    !! Setup Hamiltonian interpolation. Must be called after w90_read().
    !!
    !! This subroutine:
    !! 1. Reads eigenvalues from .eig file
    !! 2. Converts to Wannier gauge: H^W(k) = U†(k) diag(ε(k)) U(k)
    !! 3. Fourier transforms to real space: H^W(R)
    !! 4. Computes optimal Wigner-Seitz degeneracies
    !--------------------------------------------------------------------------
    !
    USE io_global,     ONLY : stdout
    USE w90_interface, ONLY : num_bands, num_kpts, num_wann
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: seedname
    !! Seedname for Wannier90 files
    !
    REAL(DP), ALLOCATABLE :: eks(:, :)
    !! Eigenvalues (num_bands, num_kpts)
    COMPLEX(DP), ALLOCATABLE :: ham_k(:, :, :)
    !! Hamiltonian in Wannier gauge at each k-point
    !
    CALL start_clock('w90_ham_setup')
    !
    WRITE(stdout, '(/5x, A)') "Setting up Hamiltonian interpolation"
    !
    ALLOCATE(eks(num_bands, num_kpts))
    ALLOCATE(ham_k(num_wann, num_wann, num_kpts))
    !
    ! 1. Read eigenvalues from .eig file
    CALL w90_read_eig(seedname, eks)
    !
    ! 2. Convert to Wannier gauge
    CALL w90_compute_ham_k(eks, ham_k)
    !
    ! 3. Fourier transform to real space (also computes WS degeneracies)
    CALL w90_compute_ham_R(ham_k)
    !
    DEALLOCATE(eks, ham_k)
    !
    WRITE(stdout, '(5x, A, I8)') "Number of R vectors: ", nR_ham
    WRITE(stdout, '(5x, A)') "Hamiltonian interpolation ready"
    !
    CALL stop_clock('w90_ham_setup')
    !
  END SUBROUTINE w90_setup_hamiltonian
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_read_eig(seedname, eks)
    !--------------------------------------------------------------------------
    !! Read eigenvalues from Wannier90 .eig file.
    !!
    !! File format (from pw2wannier90.f90):
    !!   WRITE(iun_band,'(2i10,f18.12)') ibnd, ik, eigenvalue_eV
    !--------------------------------------------------------------------------
    !
    USE constants,     ONLY : rytoev
    USE io_global,     ONLY : stdout, ionode
    USE mp,            ONLY : mp_bcast
    USE mp_images,     ONLY : intra_image_comm
    USE w90_interface, ONLY : num_bands, num_kpts
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: seedname
    !! Seedname for Wannier90 files
    REAL(DP), INTENT(OUT) :: eks(num_bands, num_kpts)
    !! Eigenvalues in Rydberg (num_bands, num_kpts)
    !
    CHARACTER(LEN=256) :: filename
    INTEGER :: iun, ibnd, ik, ibnd_read, ik_read, ios
    REAL(DP) :: eig_ev
    !
    eks = 0.d0
    !
    IF (ionode) THEN
      filename = TRIM(seedname) // '.eig'
      WRITE(stdout, '(5x, 3A)') "Reading eigenvalues from ", TRIM(filename)
      !
      OPEN(NEWUNIT=iun, FILE=filename, STATUS='old', FORM='formatted', IOSTAT=ios)
      IF (ios /= 0) CALL errore('w90_read_eig', 'Error opening '//TRIM(filename), 1)
      !
      DO ik = 1, num_kpts
        DO ibnd = 1, num_bands
          READ(iun, '(2i10, f18.12)', IOSTAT=ios) ibnd_read, ik_read, eig_ev
          IF (ios /= 0) CALL errore('w90_read_eig', 'Error reading eigenvalue', ibnd + (ik-1)*num_bands)
          !
          ! Sanity check
          IF (ibnd_read /= ibnd .OR. ik_read /= ik) THEN
            WRITE(stdout, '(5x, A, 4I8)') "Warning: index mismatch: ", ibnd, ik, ibnd_read, ik_read
          ENDIF
          !
          ! Convert eV to Rydberg
          eks(ibnd, ik) = eig_ev / rytoev
        ENDDO
      ENDDO
      !
      CLOSE(iun)
      WRITE(stdout, '(5x, A)') "Eigenvalues read successfully"
    ENDIF
    !
    CALL mp_bcast(eks, 0, intra_image_comm)
    !
  END SUBROUTINE w90_read_eig
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_compute_ham_k(eks, ham_k)
    !--------------------------------------------------------------------------
    !! Convert eigenvalues to Wannier gauge Hamiltonian.
    !!
    !! H^W_mn(k) = Σ_α U†_mα(k) ε_α(k) U_αn(k)
    !!          = [U†(k) · diag(ε(k)) · U(k)]_mn
    !--------------------------------------------------------------------------
    !
    USE w90_interface, ONLY : u_matrix, num_wann, num_kpts, num_bands
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: eks(num_bands, num_kpts)
    !! Eigenvalues in Rydberg
    COMPLEX(DP), INTENT(OUT) :: ham_k(num_wann, num_wann, num_kpts)
    !! Hamiltonian in Wannier gauge
    !
    INTEGER :: ik, ib, n
    !! Indices
    COMPLEX(DP), ALLOCATABLE :: u_ik(:, :)
    !! U matrix at k point ik, in Wannier90 band ordering (num_bands rows).
    COMPLEX(DP), ALLOCATABLE :: tmp(:, :)
    !! Temporary array for multiplication
    !
    ALLOCATE(u_ik(num_bands, num_wann))
    ALLOCATE(tmp(num_bands, num_wann))
    !
    DO ik = 1, num_kpts
      !
      ! tmp = diag(eks) * U
      ! tmp(ib, n) = eks(ib, ik) * U(ib, n, ik)
      !
      DO ib = 1, num_bands
        u_ik(ib, :) = u_matrix(ib, :, ik)
      ENDDO
      !
      DO n = 1, num_wann
        DO ib = 1, num_bands
          tmp(ib, n) = eks(ib, ik) * u_ik(ib, n)
        ENDDO
      ENDDO
      !
      ! ham_k = U† * tmp
      ! ham_k(m, n) = Σ_ib U*(ib, m) * tmp(ib, n)
      !            = Σ_ib U*(ib, m) * eks(ib) * U(ib, n)
      !
      CALL ZGEMM('C', 'N', num_wann, num_wann, num_bands, &
                 (1.d0, 0.d0), u_ik, num_bands, tmp, num_bands, &
                 (0.d0, 0.d0), ham_k(:, :, ik), num_wann)
      !
    ENDDO
    !
    DEALLOCATE(tmp)
    !
  END SUBROUTINE w90_compute_ham_k
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_compute_ham_R(ham_k)
    !--------------------------------------------------------------------------
    !! Fourier transform Hamiltonian to real space.
    !!
    !! H^W(R) = (1/Nk) Σ_k exp(-ik·R) H^W(k)
    !!
    !! Also computes optimal Wigner-Seitz degeneracies.
    !--------------------------------------------------------------------------
    !
    USE constants,     ONLY : tpi
    USE cell_base,     ONLY : at
    USE w90_interface, ONLY : num_wann, num_kpts, mp_grid, kpt_cart, wann_centers
    USE w90_wigner,    ONLY : w90_wigner_allocate, w90_wigner_deallocate, &
                              w90_find_wigner_seitz, irvec, ndegen
    USE w90_utils,     ONLY : w90_resize_in_place
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: ham_k(num_wann, num_wann, num_kpts)
    !! Hamiltonian in Wannier gauge
    !
    LOGICAL :: found
    INTEGER :: iR, iR1, iR2, iR3, ik, iT, iT2, iT_found, iw, jw
    INTEGER :: num_ws_for_R, num_ws_for_ijR
    INTEGER :: R_int(3)
    INTEGER, ALLOCATABLE :: irvec_for_R(:, :)
    REAL(DP), ALLOCATABLE :: inv_ndegen_for_R(:, :, :)
    REAL(DP) :: R(3), kdotR, mindist, r0(3)
    COMPLEX(DP) :: cfac
    !
    ! Count R vectors (product of mp_grid dimensions)
    nR_ham = mp_grid(1) * mp_grid(2) * mp_grid(3)
    !
    ! Setup Wigner-Seitz degeneracy arrays
    !
    CALL w90_wigner_allocate(mp_grid(1), mp_grid(2), mp_grid(3))
    !
    nr_ham = 0
    ALLOCATE(iRlist_ham(3, 0))
    ALLOCATE(inv_ndegen_ham(num_wann, num_wann, 0))
    !
    ! Compute Wigner-Seitz vectors and degeneracies by looping over R vectors
    !
    ! Loop over R vectors
    !
    DO iR3 = 0, mp_grid(3) - 1
      DO iR2 = 0, mp_grid(2) - 1
        DO iR1 = 0, mp_grid(1) - 1
          !
          R_int(:) = (/ iR1, iR2, iR3 /)
          !
          R = REAL(R_int, KIND=DP)
          CALL cryst_to_cart(1, R, at, +1)
          !
          ALLOCATE(irvec_for_R(3, 0))
          num_ws_for_R = 0
          !
          ! Find all the Wigner-Seitz vectors for the given R (looping over the WF pairs)
          !
          DO iw = 1, num_wann
            DO jw = 1, num_wann
              !
              r0 = wann_centers(:, iw) - wann_centers(:, jw)
              CALL w90_find_wigner_seitz(r0 - R, num_ws_for_ijR, mindist)
              !
              DO iT = 1, num_ws_for_ijR
                !
                ! Check if irvec(:, iT) is already in irvec_for_R
                !
                found = .FALSE.
                !
                DO iT2 = 1, num_ws_for_R
                  IF (ALL(irvec(:, iT) == irvec_for_R(:, iT2))) THEN
                    found = .TRUE.
                    EXIT
                  ENDIF
                ENDDO
                !
                ! If irvec(:, iT) is not found, append it to irvec_for_R
                !
                IF (.NOT. found) THEN
                  CALL w90_resize_in_place(irvec_for_R, 1)
                  irvec_for_R(:, num_ws_for_R + 1) = irvec(:, iT)
                  num_ws_for_R = num_ws_for_R + 1
                ENDIF
                !
              ENDDO
              !
            ENDDO ! jw
          ENDDO ! iw
          !
          num_ws_for_R = SIZE(irvec_for_R, 2)
          !
          ! Compute the Wigner-Seitz degeneracy for all Wigner-Seitz vectors and Wannier function pairs
          !
          ALLOCATE(inv_ndegen_for_R(num_wann, num_wann, num_ws_for_R))
          inv_ndegen_for_R = 0.d0
          !
          DO iw = 1, num_wann
            DO jw = 1, num_wann
              !
              r0 = wann_centers(:, iw) - wann_centers(:, jw)
              CALL w90_find_wigner_seitz(r0 - R, num_ws_for_ijR, mindist)
              !
              DO iT = 1, num_ws_for_ijR
                !
                ! Find index of irvec(:, iT) in irvec_for_R
                !
                iT_found = -1
                !
                DO iT2 = 1, num_ws_for_R
                  IF (ALL(irvec(:, iT) == irvec_for_R(:, iT2))) THEN
                    iT_found = iT2
                    EXIT
                  ENDIF
                ENDDO
                !
                IF (iT_found < 0) &
                  CALL errore("w90_write_real_space", "irvec(:, iT) not found in irvec_for_R", 1)
                !
                inv_ndegen_for_R(iw, jw, iT_found) = 1.d0 / REAL(ndegen(iT), DP)
                !
              ENDDO
              !
              IF (ABS(SUM(inv_ndegen_for_R(iw, jw, :)) - 1.d0) > 1.d-6) &
                CALL errore("w90_write_real_space", "Sum of Wigner-Seitz degeneracy not equal to 1", 1)
              !
            ENDDO ! jw
          ENDDO ! iw
          !
          ! Append the Wigner-Seitz vectors to iRlist_ham
          !
          CALL w90_resize_in_place(iRlist_ham, num_ws_for_R)
          CALL w90_resize_in_place(inv_ndegen_ham, num_ws_for_R)
          !
          DO iT = 1, num_ws_for_R
            iRlist_ham(:, nr_ham + iT) = R_int(:) + irvec_for_R(:, iT)
            inv_ndegen_ham(:, :, nr_ham + iT) = inv_ndegen_for_R(:, :, iT)
          ENDDO
          !
          nr_ham = nr_ham + num_ws_for_R
          !
          DEALLOCATE(irvec_for_R)
          DEALLOCATE(inv_ndegen_for_R)
          !
        ENDDO ! iR1
      ENDDO ! iR2
    ENDDO ! iR3
    !
    CALL w90_wigner_deallocate()
    !
    ! Fourier transform: H^W(R) = (1/Nk) Σ_k exp(-ik·R) H^W(k)
    !
    ALLOCATE(ham_R(num_wann, num_wann, nR_ham))
    ham_R = (0.d0, 0.d0)
    !
    DO iR = 1, nR_ham
      R(:) = REAL(iRlist_ham(:, iR), DP)
      CALL cryst_to_cart(1, R, at, +1)
      !
      DO ik = 1, num_kpts
        kdotR = tpi * SUM(kpt_cart(:, ik) * R)
        cfac = CMPLX(COS(kdotR), -SIN(kdotR), DP) / REAL(num_kpts, DP)
        ham_R(:, :, iR) = ham_R(:, :, iR) + cfac * ham_k(:, :, ik)
      ENDDO
    ENDDO
    !
    ham_R = ham_R * inv_ndegen_ham
    !
  END SUBROUTINE w90_compute_ham_R
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_interpolate_ham(xk, eig, cuf)
    !--------------------------------------------------------------------------
    !! Interpolate Hamiltonian to an arbitrary k-point and diagonalize.
    !!
    !! H^W_mn(k') = Σ_R exp(ik'·R) H^W_mn(R) / ndegen(R, m, n)
    !!
    !! Then diagonalize H^W(k') to get eigenvalues and eigenvectors.
    !--------------------------------------------------------------------------
    !
    USE constants,     ONLY : tpi
    USE cell_base,     ONLY : at
    USE w90_interface, ONLY : num_wann
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xk(3)
    !! k-point in Cartesian coordinates (units of 2π/a)
    REAL(DP), INTENT(OUT) :: eig(num_wann)
    !! Interpolated eigenvalues
    COMPLEX(DP), INTENT(OUT) :: cuf(num_wann, num_wann)
    !! Eigenvectors (rotation matrix U such that H = U diag(eig) U†)
    !
    INTEGER :: iR, m, n, info, lwork
    REAL(DP) :: R(3), kdotR
    COMPLEX(DP) :: cfac
    COMPLEX(DP), ALLOCATABLE :: ham_k(:, :), work(:)
    REAL(DP), ALLOCATABLE :: rwork(:)
    !
    ALLOCATE(ham_k(num_wann, num_wann))
    !
    ham_k = (0.d0, 0.d0)
    !
    ! Fourier transform with element-wise WS degeneracy
    ! H^W_mn(k) = Σ_R exp(ik·R) H^W_mn(R) / ndegen(R, m, n)
    !
    DO iR = 1, nR_ham
      R(:) = REAL(iRlist_ham(:, iR), DP)
      CALL cryst_to_cart(1, R, at, +1)
      kdotR = tpi * SUM(xk * R)
      cfac = CMPLX(COS(kdotR), SIN(kdotR), DP)
      !
      DO m = 1, num_wann
        DO n = 1, num_wann
          ham_k(m, n) = ham_k(m, n) + cfac * ham_R(m, n, iR)! / REAL(ndegen_ham(iR, m, n), DP)
        ENDDO
      ENDDO
    ENDDO
    !
    ! Hermitize (for numerical stability)
    DO m = 1, num_wann
      DO n = 1, m
        ham_k(m, n) = (ham_k(m, n) + CONJG(ham_k(n, m))) * 0.5d0
        ham_k(n, m) = CONJG(ham_k(m, n))
      ENDDO
    ENDDO
    !
    ! Diagonalize using ZHEEV
    ! On output, cuf contains the eigenvectors (column-wise)
    !
    cuf = ham_k
    lwork = 2 * num_wann
    ALLOCATE(work(lwork), rwork(3*num_wann - 2))
    !
    CALL ZHEEV('V', 'U', num_wann, cuf, num_wann, eig, work, lwork, rwork, info)
    !
    DEALLOCATE(work, rwork, ham_k)
    !
    IF (info /= 0) CALL errore('w90_interpolate_ham', 'ZHEEV diagonalization failed', ABS(info))
    !
  END SUBROUTINE w90_interpolate_ham
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_cleanup_hamiltonian()
    !--------------------------------------------------------------------------
    !! Deallocate Hamiltonian interpolation arrays.
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    IF (ALLOCATED(ham_R)) DEALLOCATE(ham_R)
    IF (ALLOCATED(iRlist_ham)) DEALLOCATE(iRlist_ham)
    IF (ALLOCATED(inv_ndegen_ham)) DEALLOCATE(inv_ndegen_ham)
    !
    nR_ham = 0
    !
  END SUBROUTINE w90_cleanup_hamiltonian
  !----------------------------------------------------------------------------
  !
!----------------------------------------------------------------------------
END MODULE w90_hamiltonian
!----------------------------------------------------------------------------
