!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE w90_real_space
  !----------------------------------------------------------------------------
  !! Module for writing real-space Wannier functions.
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER :: nR_ws
  !! Number of Wigner-Seitz grid points
  INTEGER, ALLOCATABLE :: iRlist_ws(:, :)
  !! List of Wigner-Seitz grid points, (3, nR_ws)
  !
  CONTAINS
  !
  SUBROUTINE w90_write_real_space(tmp_dir_wfc, tmp_dir_write)
    !----------------------------------------------------------------------------
    !! Read wavefunctions from file
    !! Rotate wavefunctions to the Wannier gauge
    !! Fourier transform to real space
    !----------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE mp,               ONLY : mp_sum, mp_bcast
    USE mp_global,        ONLY : ionode_id
    USE mp_images,        ONLY : intra_image_comm
    USE fft_interfaces,   ONLY : invfft
    USE fft_base,         ONLY : dffts
    USE fft_types,        ONLY : fft_index_to_3d
    USE klist,            ONLY : xk, nks, nkstot, igk_k, ngk
    USE io_global,        ONLY : stdout, ionode
    USE io_files,         ONLY : restart_dir, tmp_dir, prefix, create_directory
    USE wvfct,            ONLY : npwx, nbnd
    USE cell_base,        ONLY : at
    USE noncollin_module, ONLY : npol
    USE wavefunctions,    ONLY : evc
    USE pw_restart_new,   ONLY : read_collected_wfc
    USE w90_interface,    ONLY : u_matrix, num_wann, num_kpts, kpt_cart, wann_centers, mp_grid
    USE w90_wigner,       ONLY : w90_wigner_allocate, w90_wigner_deallocate, &
                                 w90_find_wigner_seitz, irvec, ndegen
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(in) :: tmp_dir_wfc
    !! Directory for where the formatted wavefunctions are stored.
    CHARACTER(LEN=256), INTENT(in) :: tmp_dir_write
    !! Directory for writing the real-space Wannier functions.
    !
    CHARACTER(LEN=256) :: tmp_dir_bak
    !! Backup of tmp_dir
    CHARACTER(LEN=256) :: filename
    !! File name for the Wigner-Seitz grid points
    LOGICAL :: offrange
    !! Flag for out-of-range grid points
    LOGICAL :: found
    !! Flag for found Wigner-Seitz vector
    INTEGER :: ik
    !! k point index
    INTEGER :: ik_global
    !! k point index in the global index for Wannier functions
    INTEGER :: iw
    !! Wannier function index
    INTEGER :: ipol
    !! Spin index
    INTEGER :: ir
    !! Real space FFT grid index
    INTEGER :: iT, iT2, iT_found
    !! Supercell lattice vector index
    INTEGER :: i, j, k
    !! Direction index
    INTEGER :: num_ws_for_r0_R
    !! number of Wigner-Seitz vectors for the given r0 and R
    INTEGER :: num_ws_for_R
    !! Number of Wigner-Seitz vectors for the given R
    INTEGER :: npw
    !! Number of plane waves
    INTEGER :: iR1, iR2, iR3
    !! Real-space lattice vector
    INTEGER :: iun
    !! Unit number for writing the Wigner-Seitz grid points
    INTEGER :: R_int(3)
    !! Real-space lattice vector in integer coordinates
    INTEGER :: RT_int(3)
    !! Real-space lattice vector in integer coordinates
    REAL(DP) :: R(3)
    !! Real-space lattice vector in Cartesian coordinates
    REAL(DP) :: xxk(3)
    !! k point in the coarse grid (in Cartesian coordinates).
    REAL(DP) :: r0(3)
    !! Real space grid point inside unit cell
    REAL(DP) :: mindist
    !! Minimum distance (needed to call w90_find_wigner_seitz, not used here.)
    COMPLEX(DP), ALLOCATABLE :: wan_G(:, :)
    !! Wavefunction in the Wannier gauge in reciprocal space
    COMPLEX(DP), ALLOCATABLE :: wan_kr(:, :, :, :)
    !! Wavefunction in the Wannier gauge in real space
    COMPLEX(DP), ALLOCATABLE :: aux(:)
    !! Auxiliary array for the wavefunctions
    !
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    INTEGER, ALLOCATABLE :: irvec_all(:, :)
    REAL(DP), ALLOCATABLE :: inv_ndegen_all(:, :)
    !
    CALL start_clock("w90_write_rs")
    !
    ALLOCATE(wan_G(npwx*npol, num_wann))
    ALLOCATE(wan_kr(dffts%nnr, npol, num_wann, nks))
    ALLOCATE(aux(dffts%nnr))
    !
    ! Compute wan_kr, the Wannier functions in the momentum space in real space grid.
    !
    IF (nkstot /= num_kpts) THEN
      WRITE(stdout, '(5x, A, I8, A, I8)') "nkstot = ", nkstot, " num_kpts = ", num_kpts
      CALL errore("w90_write_real_space", "nkstot different from num_kpts of chk file", 1)
    ENDIF
    !
    WRITE(stdout, '(/5x, A)') "Computing Wannier functions in momentum space"
    !
    ! TODO: Use buffer to save memory instead of having nks index in wan_kr
    !
    DO ik = 1, nks
      !
      xxk = xk(:, ik)
      npw = ngk(ik)
      ik_global = global_kpoint_index(nkstot, ik)
      !
      IF (ANY(ABS(xxk(1:3) - kpt_cart(1:3, ik_global)) > 1.d-6)) CALL errore("w90_write_real_space", &
        "k point mismatch between QE and Wannier90", 1)
      !
      ! restart_dir uses tmp_dir internally, so we need to change it temporarily to tmp_dir_wfc.
      !
      tmp_dir_bak = tmp_dir
      tmp_dir = tmp_dir_wfc
      CALL read_collected_wfc(restart_dir(), ik, evc)
      tmp_dir = tmp_dir_bak
      !
      ! TODO: Optimize by not computing excluded bands
      !
      CALL ZGEMM('N', 'N', npwx*npol, num_wann, nbnd, (1.d0, 0.d0), evc, npwx*npol, &
                u_matrix(1, 1, ik_global), nbnd, (0.d0, 0.d0), wan_G, npwx*npol)
      !
      DO iw = 1, num_wann
        !
        ! In real space, wan_kr(:, :, iw) is centered at wann_center(:, iw). We move the center
        ! to origin by applying shift = - wann_center(:, iw).
        !
        ! After interpolation, we shift the interpolated wavefunction back by applying
        ! shift = + wann_center(:, iw). This is done in w90_interpolate.f90.
        !
        CALL w90_shift_wfc_G(ik, -wann_centers(:, iw), wan_G(:, iw))
        !
        DO ipol = 1, npol
          !
          aux = (0.0_dp, 0.0_dp)
          aux(dffts%nl(igk_k(1:npw, ik))) = wan_G((ipol-1)*npwx+1 : (ipol-1)*npwx+npw, iw)
          CALL invfft('Wave', aux, dffts)
          !
          CALL w90_multiply_iqr(dffts, xxk, aux)
          !
          wan_kr(:, ipol, iw, ik) = aux
          !
        ENDDO
      ENDDO
      !
    ENDDO ! ik
    !
    ! Create output directory
    !
    CALL create_directory(tmp_dir_write)
    !
    WRITE(stdout, '(/5x, A)') "Computing and writing real-space Wannier functions to file"
    !
    CALL w90_wigner_allocate(mp_grid(1), mp_grid(2), mp_grid(3))
    !
    nR_ws = 0
    ALLOCATE(iRlist_ws(3, 0))
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
          ALLOCATE(irvec_all(3, 0))
          num_ws_for_R = 0
          !
          ! Find all the Wigner-Seitz vectors for the given R (looping over the grid points r)
          !
          DO ir = 1, dffts%nnr
            !
            CALL fft_index_to_3d(ir, dffts, i, j, k, offrange)
            IF (offrange) CYCLE
            !
            r0(1) = REAL(i, DP) / REAL(dffts%nr1, DP)
            r0(2) = REAL(j, DP) / REAL(dffts%nr2, DP)
            r0(3) = REAL(k, DP) / REAL(dffts%nr3, DP)
            CALL cryst_to_cart(1, r0, at, +1)
            !
            CALL w90_find_wigner_seitz(r0 - R, num_ws_for_r0_R, mindist)
            !
            DO iT = 1, num_ws_for_r0_R
              !
              ! Check if irvec(:, iT) is already in irvec_all
              !
              found = .FALSE.
              !
              DO iT2 = 1, num_ws_for_R
                IF (ALL(irvec(:, iT) == irvec_all(:, iT2))) THEN
                  found = .TRUE.
                  EXIT
                ENDIF
              ENDDO
              !
              ! If irvec(:, iT) is not found, append it to irvec_all
              !
              IF (.NOT. found) THEN
                CALL w90_resize_in_place(irvec_all, 1)
                irvec_all(:, num_ws_for_R + 1) = irvec(:, iT)
                num_ws_for_R = num_ws_for_R + 1
              ENDIF
              !
            ENDDO
            !
          ENDDO ! ir
          !
          ! Gather irvec_all over intra_pool_comm and make it unique
          !
          CALL w90_irvec_gather_and_unique(irvec_all)
          CALL mp_bcast(num_ws_for_R, ionode_id, intra_image_comm)
          num_ws_for_R = SIZE(irvec_all, 2)
          !
          ! Compute the Wigner-Seitz degeneracy for all Wigner-Seitz vectors and real-space grid points
          !
          ALLOCATE(inv_ndegen_all(dffts%nnr, num_ws_for_R))
          inv_ndegen_all = 0.d0
          !
          DO ir = 1, dffts%nnr
            !
            CALL fft_index_to_3d(ir, dffts, i, j, k, offrange)
            IF (offrange) CYCLE
            !
            r0(1) = REAL(i, DP) / REAL(dffts%nr1, DP)
            r0(2) = REAL(j, DP) / REAL(dffts%nr2, DP)
            r0(3) = REAL(k, DP) / REAL(dffts%nr3, DP)
            CALL cryst_to_cart(1, r0, at, +1)
            !
            CALL w90_find_wigner_seitz(r0 - R, num_ws_for_r0_R, mindist)
            !
            DO iT = 1, num_ws_for_r0_R
              !
              ! Find index of irvec(:, iT) in irvec_all
              !
              iT_found = -1
              !
              DO iT2 = 1, num_ws_for_R
                IF (ALL(irvec(:, iT) == irvec_all(:, iT2))) THEN
                  iT_found = iT2
                  EXIT
                ENDIF
              ENDDO
              !
              IF (iT_found < 0) &
                CALL errore("w90_write_real_space", "irvec(:, iT) not found in irvec_all", 1)
              !
              inv_ndegen_all(ir, iT_found) = 1.d0 / REAL(ndegen(iT), DP)
              !
            ENDDO
            !
            IF (ABS(SUM(inv_ndegen_all(ir, :)) - 1.d0) > 1.d-6) &
              CALL errore("w90_write_real_space", "Sum of Wigner-Seitz degeneracy not equal to 1", 1)
            !
          ENDDO ! ir
          !
          CALL w90_resize_in_place(iRlist_ws, num_ws_for_R)
          !
          DO iT = 1, num_ws_for_R
            !
            RT_int(:) = R_int(:) + irvec_all(:, iT)
            !
            ! Append the Wigner-Seitz vectors to iRlist_ws
            !
            iRlist_ws(:, nR_ws + iT) = RT_int(:)
            !
            ! Write the Wannier functions in real space to file
            !
            CALL w90_compute_and_write_wan_Rr(RT_int(:), wan_kr, inv_ndegen_all(:, iT), tmp_dir_write)
            !
          ENDDO
          !
          nR_ws = nR_ws + num_ws_for_R
          !
          DEALLOCATE(irvec_all)
          DEALLOCATE(inv_ndegen_all)
          !
        ENDDO ! iR1
      ENDDO ! iR2
    ENDDO ! iR3
    !
    CALL w90_wigner_deallocate()
    !
    ! Print the Wigner-Seitz grid points to file in tmp_dir_write
    !
    IF (ionode) THEN
      filename = TRIM(tmp_dir_write) // TRIM(prefix) // '.wan_R.dat'
      OPEN(NEWUNIT=iun, FILE=filename, STATUS='replace', FORM='formatted')
      WRITE(iun, "(I8)") nR_ws
      DO iR = 1, nR_ws
        WRITE(iun, "(3I8)") iRlist_ws(:, iR)
      ENDDO
      CLOSE(UNIT=iun, STATUS='KEEP')
    ENDIF
    !
    DEALLOCATE(iRlist_ws)
    !
    CALL stop_clock("w90_write_rs")
    !
  !------------------------------------------------------------------------------------
  END SUBROUTINE w90_write_real_space
  !------------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------------
  SUBROUTINE w90_compute_and_write_wan_Rr(iRs, wan_kr, factor_r, tmp_dir_write)
    !----------------------------------------------------------------------------
    !! Compute the Wannier function in real space and write it to file.
    !! Input wan_kr is the Wannier function in the momentum space in real space grid.
    !----------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE constants,        ONLY : tpi
    USE mp,               ONLY : mp_sum
    USE mp_global,        ONLY : ionode_id
    USE mp_pools,         ONLY : inter_pool_comm, my_pool_id
    USE fft_base,         ONLY : dffts
    USE klist,            ONLY : xk, nks
    USE cell_base,        ONLY : at
    USE noncollin_module, ONLY : npol
    USE w90_interface,    ONLY : num_wann, num_kpts
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iRs(3)
    !! Real-space lattice vector
    COMPLEX(DP), INTENT(IN) :: wan_kr(dffts%nnr, npol, num_wann, nks)
    !! Wavefunction in the Wannier gauge in real space
    REAL(DP), INTENT(IN) :: factor_r(dffts%nnr)
    !! Factor to be multiplied to the real-space Wannier functions (for Wigner-Seitz degeneracy)
    CHARACTER(LEN=256), INTENT(in) :: tmp_dir_write
    !! Directory to write the real-space Wannier functions.
    !
    INTEGER :: ik
    !! k point index
    INTEGER :: iw
    !! Wannier function index
    INTEGER :: ipol
    !! Spin index
    REAL(DP) :: R(3)
    !! Real-space lattice vector in Cartesian coordinates
    REAL(DP) :: xxk(3)
    !! k point in the coarse grid (in Cartesian coordinates).
    REAL(DP) :: kdotR
    !! Argument for the Fourier transform
    COMPLEX(DP) :: fac
    !! Factor for Fourier transformation
    COMPLEX(DP), ALLOCATABLE :: wan_Rr(:, :, :)
    !! Wannier function in real space
    !
    ALLOCATE(wan_Rr(dffts%nnr, npol, num_wann))
    !
    R(:) = REAL(iRs, KIND=DP)
    CALL cryst_to_cart(1, R, at, +1)
    !
    wan_Rr = (0.d0, 0.d0)
    !
    ! Fourier transform to real space
    ! wan_Rr = \sum_k exp(-i k R) wan_kr / num_kpts
    !
    DO ik = 1, nks
      xxk = xk(:, ik)
      kdotR = tpi * SUM(xxk * R)
      fac = CMPLX(COS(kdotR), -SIN(kdotR), DP) / REAL(num_kpts, DP)
      !
      wan_Rr(:, :, :) = wan_Rr(:, :, :) + wan_kr(:, :, :, ik) * fac
    ENDDO ! ik
    !
    CALL mp_sum(wan_Rr, inter_pool_comm)
    !
    ! Multiply by the factor_r for Wigner-Seitz degeneracy
    !
    DO iw = 1, num_wann
      DO ipol = 1, npol
        wan_Rr(:, ipol, iw) = wan_Rr(:, ipol, iw) * factor_r(:)
      ENDDO
    ENDDO
    !
    ! Write wan_Rr to file
    !
    IF (my_pool_id == ionode_id) THEN
      CALL w90_readwrite_wan_Rr(iRs, wan_Rr, +1, tmp_dir_write)
    ENDIF
    !
  !------------------------------------------------------------------------------------
  END SUBROUTINE w90_compute_and_write_wan_Rr
  !------------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------------
  SUBROUTINE w90_resize_in_place(arr, ncol_add)
    !----------------------------------------------------------------------------
    !! Resize the array in place by adding ncol_add columns.
    !! The data of the original array are kept.
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: arr(:, :)
    !! Array to be resized. Input: size (nrow, ncol), Output: size (nrow, ncol, n + ncol_add)
    INTEGER, INTENT(IN) :: ncol_add
    !! Number of columns to add
    !
    INTEGER :: nrow
    !! Number of rows in the original array
    INTEGER :: ncol
    !! Number of columns in the original array
    INTEGER, ALLOCATABLE :: arr_copy(:, :)
    !! Temporary array to hold the original data
    !
    nrow = SIZE(arr, 1)
    ncol = SIZE(arr, 2)
    ALLOCATE(arr_copy(nrow, ncol + ncol_add))
    arr_copy(:, 1:ncol) = arr(:, 1:ncol)
    !
    CALL move_alloc(arr_copy, arr)  ! Move allocation from arr_copy to arr
    !
    ! More verbose form with an explicit copy that could replace the above line :
    ! DEALLOCATE(arr)
    ! ALLOCATE(arr(nrow, ncol + ncol_add))
    ! arr(:, :) = arr_copy(:, :)
    ! DEALLOCATE(arr_copy)
    !
  !------------------------------------------------------------------------------------
  END SUBROUTINE w90_resize_in_place
  !------------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------------
  SUBROUTINE w90_irvec_gather_and_unique(arr)
    !----------------------------------------------------------------------------
    !! Gather the array over the intra_pool_comm and remove duplicates.
    !--------------------------------------------------------------------------
    !
    USE mp,               ONLY : mp_bcast, mp_gather, mp_sum
    USE mp_pools,         ONLY : nproc_pool, me_pool, intra_pool_comm, root_pool
    !
    IMPLICIT NONE
    !
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: arr(:, :)
    !! Input: Array to be gathered and made unique. Output: Resulting array with unique columns.
    !
    INTEGER :: i1, i2
    !! Column indices
    INTEGER :: ncol_unique
    !! Number of unique columns in the gathered array
    INTEGER :: iproc
    !! Process index in the intra_pool_comm
    LOGICAL, ALLOCATABLE :: keep(:)
    !! Mark unique columns
    INTEGER, ALLOCATABLE :: recvcount(:)
    !! Number of columns to receive from each process
    INTEGER, ALLOCATABLE :: displs(:)
    !! Displacements for each process in the gathered array (cumulantive sum of recvcount)
    INTEGER, ALLOCATABLE :: arr_gathered(:, :)
    !! Gathered array from all processes in the intra_pool_comm
    INTEGER, ALLOCATABLE :: arr_unique(:, :)
    !! Unique columns from the gathered array
    !
    ! Gather irvec_all over intra_pool_comm
    !
    ALLOCATE(recvcount(nproc_pool))
    recvcount(:) = 0
    recvcount(me_pool + 1) = SIZE(arr, 2)
    CALL mp_sum(recvcount, intra_pool_comm)
    !
    ALLOCATE(displs(nproc_pool))
    displs(1) = 0
    DO iproc = 2, nproc_pool
      displs(iproc) = displs(iproc - 1) + recvcount(iproc - 1)
    ENDDO
    !
    ALLOCATE(arr_gathered(SIZE(arr, 1), SUM(recvcount)))
    !
    CALL mp_gather(arr, arr_gathered, recvcount, displs, root_pool, intra_pool_comm)
    CALL mp_bcast(arr_gathered, root_pool, intra_pool_comm)
    !
    ! Remove duplicates
    !
    ncol_unique = 0
    ALLOCATE(arr_unique(SIZE(arr_gathered, 1), 0))
    ALLOCATE(keep(SIZE(arr_gathered, 2)))
    keep(:) = .TRUE.
    !
    DO i1 = 1, SIZE(arr_gathered, 2)
      !
      ! If arr_gathered(:, i1) is a duplicate of a previous column, skip it.
      !
      IF (.NOT. keep(i1)) CYCLE
      !
      ! Store the unique column in arr_unique
      !
      ncol_unique = ncol_unique + 1
      CALL w90_resize_in_place(arr_unique, 1)
      arr_unique(:, ncol_unique) = arr_gathered(:, i1)
      !
      ! Mark all duplicates with arr_gathered(:, i1)
      !
      DO i2 = i1 + 1, SIZE(arr_gathered, 2)
        IF (ALL(arr_gathered(:, i1) == arr_gathered(:, i2))) THEN
          keep(i2) = .FALSE.
        ENDIF
      ENDDO
      !
    ENDDO
    !
    CALL mp_bcast(arr_unique, root_pool, intra_pool_comm)
    !
    CALL move_alloc(arr_unique, arr)  ! Move allocation from arr_unique to arr
    !
  !------------------------------------------------------------------------------------
  END SUBROUTINE w90_irvec_gather_and_unique
  !------------------------------------------------------------------------------------
  !
END MODULE w90_real_space
