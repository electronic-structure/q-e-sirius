!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE w90_interface
  !----------------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: have_disentangled
  !
  CHARACTER(len=256) :: seedname
  !! Seedname for Wannier90 file names
  !
  INTEGER :: mp_grid(3)
  !! Coarse grid for Wannier functions
  INTEGER :: num_bands
  INTEGER :: num_exclude_bands
  INTEGER :: num_wann
  INTEGER :: num_kpts
  INTEGER, ALLOCATABLE :: exclude_bands(:)
  REAL(DP), ALLOCATABLE :: kpt_crys(:,:)
  !! Coarse k points in crystal coordinates
  REAL(DP), ALLOCATABLE :: kpt_cart(:,:)
  !! Coarse k points in crystal coordinates
  !
  REAL(DP), ALLOCATABLE :: wann_centers(:,:)
  !! Wannier function centers in Cartesian alat units. (3, num_wann)
  REAL(DP), ALLOCATABLE :: wann_spreads(:)
  !! Wannier function spreads in alat^2 units. (num_wann)
  !
  COMPLEX(DP), ALLOCATABLE :: u_matrix(:, :, :)
  !! Wannier gauge matrix (nbnd, num_wann, num_kpts).
  !! Includes both the disentanglement (if present) and maximal localization.
  COMPLEX(DP), ALLOCATABLE :: u_matrix_opt(:, :, :)
  COMPLEX(DP), ALLOCATABLE :: m_matrix(:, :, :, :)
  !
  CONTAINS
  !
  SUBROUTINE w90_read(seedname)
    !-------------------------------------------------------------------------
    !! Read Wannier90 checkpoint file and setup the final gauge matrix u_matrix
    !-------------------------------------------------------------------------
    !
    USE kinds,      ONLY : DP
    USE constants,  ONLY : BOHR_RADIUS_ANGS
    USE io_global,  ONLY : stdout
    USE wvfct,      ONLY : nbnd
    USE cell_base,  ONLY : bg, alat
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(IN) :: seedname
    !
    INTEGER :: ik
    INTEGER :: iw
    INTEGER :: ib
    INTEGER :: ib_included
    COMPLEX(DP), ALLOCATABLE :: u_matrix_tmp(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: tmp(:, :)
    !
    CALL w90_read_chkpt(seedname)
    !
    ! Convert units
    !
    ! k points: crystal to Cartesian
    !
    ALLOCATE(kpt_cart(3, num_kpts))
    kpt_cart = kpt_crys
    CALL cryst_to_cart(num_kpts, kpt_cart, bg, +1)
    !
    ! Wannier centers and spreads: Ang to alat
    wann_centers = wann_centers / (alat * BOHR_RADIUS_ANGS)
    wann_spreads = wann_spreads / (alat * BOHR_RADIUS_ANGS)**2
    !
    WRITE(stdout, '(5x, a)') "Wannier centers (alat) and spreads (alat^2)"
    DO iw = 1, num_wann
      WRITE(stdout, '(5x, 4F16.10)') wann_centers(:, iw), wann_spreads(iw)
    ENDDO
    !
    IF (nbnd /= num_bands + num_exclude_bands) CALL errore("w90_read", &
        'Number of bands in Wannier90 checkpoint file does not match QE input', 1)
    !
    ! If using disentanglement, multiply u_matrix_opt and u_matrix, store them in u_matrix.
    !
    IF (have_disentangled) THEN
        !
        ALLOCATE(u_matrix_tmp(num_bands, num_wann, num_kpts))
        DO ik = 1, num_kpts
          CALL ZGEMM('N', 'N', num_bands, num_wann, num_wann, &
              (1.d0, 0.d0), u_matrix_opt(:, :, ik), num_bands, u_matrix(:, :, ik), num_wann, &
              (0.d0, 0.d0), u_matrix_tmp(:, :, ik), num_bands)
        ENDDO
        !
        ! Now u_matrix contains the final gauge matrix. u_matrix_opt is not used anymore.
        !
        DEALLOCATE(u_matrix_opt)
        !
        ! Copy u_matrix_tmp to u_matrix
        !
        DEALLOCATE(u_matrix)
        ALLOCATE(u_matrix(num_bands, num_wann, num_kpts))
        u_matrix = u_matrix_tmp
        DEALLOCATE(u_matrix_tmp)
        !
    ENDIF
    !
    ! If some bands are excluded, fill zeros to make u_matrix
    ! have size (nbnd, num_wann, num_kpts).
    !
    ! TODO: This is not efficient. We should avoid computing stuff for the excluded bands.
    !
    IF (num_exclude_bands > 0) THEN
        !
        ALLOCATE(u_matrix_tmp(nbnd, num_wann, num_kpts))
        u_matrix_tmp = (0.d0, 0.d0)
        !
        ib_included = 0
        DO ib = 1, nbnd
          IF (ALL(exclude_bands /= ib)) THEN
              ! Band ib is not excluded.
              ib_included = ib_included + 1
              u_matrix_tmp(ib, :, :) = u_matrix(ib_included, :, :)
          ENDIF
        ENDDO
        !
        IF (ib_included /= num_bands) CALL errore("w90_read", &
          'exclude_bands and num_bands are not consistent', 1)
        !
        ! Copy u_matrix_tmp to u_matrix
        !
        DEALLOCATE(u_matrix)
        ALLOCATE(u_matrix(num_bands, num_wann, num_kpts))
        u_matrix = u_matrix_tmp
        DEALLOCATE(u_matrix_tmp)
        !
    ENDIF
    !
    ! Check u_matrix is a semi-unitary matrix
    ! u_matrix(:, :, ik)^H * u_matrix(:, :, ik) = I(num_wann)
    !
    ALLOCATE(tmp(num_wann, num_wann))
    DO ik = 1, num_kpts
        CALL ZGEMM('C', 'N', num_wann, num_wann, nbnd, &
          (1.d0, 0.d0), u_matrix(:, :, ik), nbnd, u_matrix(:, :, ik), nbnd, &
          (0.d0, 0.d0), tmp, num_wann)
        !
        DO ib = 1, num_wann
          tmp(ib, ib) = tmp(ib, ib) - 1.d0
        ENDDO
        !
        IF (SUM(ABS(tmp)) > 1.d-10) THEN
          WRITE(stdout, '(5x,a,I8,a,ES20.10)') "ik = ", ik, " Unitarity error = ", SUM(ABS(tmp))
        ENDIF
        !
    ENDDO
    DEALLOCATE(tmp)
    !
  END SUBROUTINE w90_read
  !
  SUBROUTINE w90_read_chkpt(seedname)
    !-------------------------------------------------------------------------
    !! Read unformatted checkpoint file written by Wannier90
    !-------------------------------------------------------------------------
    !
    USE kinds,      ONLY : DP
    USE io_global,  ONLY : stdout
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(IN) :: seedname
    !
    CHARACTER(LEN=33) :: header
    CHARACTER(LEN=20) :: checkpoint
    INTEGER :: chk_unit, i, j, k, l, nkp, ierr
    INTEGER :: nntot
    REAL(DP) :: omega_invariant
    REAL(DP) :: real_lattice(3, 3), recip_lattice(3, 3)
    COMPLEX(DP), ALLOCATABLE :: m_matrix(:,:,:,:)
    LOGICAL, ALLOCATABLE :: lwindow(:, :)
    INTEGER, ALLOCATABLE :: ndimwin(:)
    !
    WRITE (stdout, '(5x, 3a)') 'Opening file ', TRIM(seedname), '.chk'
    !
    OPEN(NEWUNIT=chk_unit, FILE=TRIM(seedname)//'.chk', status='old', form='unformatted', err=121)
    !
    ! Read comment line
    READ(chk_unit) header
    WRITE(stdout, '(5x, a)') TRIM(header)
    !
    ! num_bands
    READ(chk_unit) num_bands
    WRITE (stdout, '(5x, a, I8)') "Number of bands :", num_bands
    !
    ! num_exclude_bands
    READ(chk_unit) num_exclude_bands
    WRITE (stdout, '(5x, a, I8)') "Number of excluded bands :", num_exclude_bands
    IF (num_exclude_bands < 0) THEN
        CALL errore("w90_read_chkpt", 'Invalid value for num_exclude_bands', 1)
    ENDIF
    !
    ! exclude_bands
    ALLOCATE(exclude_bands(num_exclude_bands), STAT=ierr)
    IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating exclude_bands', 1)
    READ(chk_unit) (exclude_bands(i), i=1, num_exclude_bands) ! Excluded bands
    WRITE (stdout, '(5x, a)', ADVANCE='no') "Excluded bands :"
    IF (num_exclude_bands == 0) THEN
        WRITE (stdout, '(5x, a)') "none."
    ELSE
        DO i = 1, num_exclude_bands - 1
          WRITE (stdout, '(1x, I4, a)', ADVANCE='no') exclude_bands(i), ','
        ENDDO
        WRITE (stdout, '(1x, I4, a)') exclude_bands(num_exclude_bands), '.'
    ENDIF
    !
    ! real_lattice, recip_lattice
    READ (chk_unit) ((real_lattice(i, j), i=1, 3), j=1, 3)
    READ (chk_unit) ((recip_lattice(i, j), i=1, 3), j=1, 3)
    !
    ! k point mesh
    READ (chk_unit) num_kpts
    WRITE (stdout, '(5x, a, I8)') "Num kpts :", num_kpts
    READ (chk_unit) (mp_grid(i), i=1, 3)         ! M-P grid
    WRITE (stdout, '(5x, a, 3I8)') "mp_grid :", mp_grid
    !
    ALLOCATE (kpt_crys(3, num_kpts), STAT=ierr)
    IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating kpt_crys', 1)
    READ (chk_unit) ((kpt_crys(i, nkp), i=1, 3), nkp=1, num_kpts)
    !
    ! nntot
    READ (chk_unit) nntot
    WRITE (stdout, '(5x, a, I8)') "nntot :", nntot
    !
    ! num_wann
    READ (chk_unit) num_wann
    WRITE (stdout, '(5x, a, I8)') "num_wann :", num_wann
    !
    ! checkpoint
    READ (chk_unit) checkpoint
    checkpoint = ADJUSTL(TRIM(checkpoint))
    WRITE (stdout, '(5x, a)') "checkpoint: " // TRIM(checkpoint)
    !
    ! have_disentangled: whether a disentanglement has been performed
    READ (chk_unit) have_disentangled
    !
    IF (have_disentangled) THEN
        WRITE (stdout, '(5x, a)') "have_disentangled: TRUE"
        !
        READ (chk_unit) omega_invariant
        !
        ! lwindow
        ALLOCATE (lwindow(num_bands, num_kpts), STAT=ierr)
        IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating lwindow', 1)
        READ (chk_unit, err=122) ((lwindow(i, k), i=1, num_bands), k=1, num_kpts)
        WRITE (stdout, '(5x, a)') "lwindow: read."
        !
        ! ndimwin
        ALLOCATE (ndimwin(num_kpts), STAT=ierr)
        IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating ndimwin', 1)
        READ (chk_unit, err=123) (ndimwin(k), k=1, num_kpts)
        WRITE (stdout, '(5x, a)') "ndimwin: read."
        !
        ! U_matrix_opt
        ALLOCATE (u_matrix_opt(num_bands, num_wann, num_kpts), STAT=ierr)
        IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating u_matrix_opt', 1)
        READ (chk_unit, err=124) (((u_matrix_opt(i, j, k), i=1, num_bands), j=1, num_wann), k=1, num_kpts)
        WRITE (stdout, '(5x, a, 3I8)') "U_matrix_opt: read. shape = ", SHAPE(u_matrix_opt)
        !
    ELSE
        WRITE (stdout, '(5x, a)') "have_disentangled: FALSE"
    ENDIF
    !
    ! U_matrix
    !
    ALLOCATE (u_matrix(num_wann, num_wann, num_kpts), STAT=ierr)
    IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating u_matrix', 1)
    READ (chk_unit, err=125) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)
    WRITE (stdout, '(5x, a, 3I8)') "U_matrix: read. shape = ", SHAPE(u_matrix)
    !
    ! M_matrix
    !
    ALLOCATE (m_matrix(num_wann, num_wann, nntot, num_kpts), STAT=ierr)
    IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating m_matrix', 1)
    READ (chk_unit, err=126) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, nntot), l=1, num_kpts)
    WRITE (stdout, '(5x, a)') "M_matrix: read."
    !
    ! wann_centers
    !
    ALLOCATE (wann_centers(3, num_wann), STAT=ierr)
    IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating wann_centers', 1)
    READ (chk_unit, err=127) ((wann_centers(i, j), i=1, 3), j=1, num_wann)
    WRITE (stdout, '(5x, a)') "wannier_centres: read."
    !
    ! wannier spreads
    !
    ALLOCATE (wann_spreads(num_wann), STAT=ierr)
    IF (ierr /= 0) CALL errore("w90_read_chkpt", 'Error allocating wann_spreads', 1)
    READ (chk_unit, err=128) (wann_spreads(i), i=1, num_wann)
    WRITE (stdout, '(5x, a)') "wannier_spreads: read."
    !
    CLOSE(chk_unit)
    !
    WRITE(stdout, '(5x, a)') 'Completed reading Wannier90 checkpoint file'
    !
    RETURN
    !
121   CALL errore("w90_read_chkpt", "Error opening "//TRIM(seedname)//".chk", 1)
122   CALL errore("w90_read_chkpt", "Error reading lwindow from "//TRIM(seedname)//".chk", 1)
123   CALL errore("w90_read_chkpt", "Error reading ndimwin from "//TRIM(seedname)//".chk", 1)
124   CALL errore("w90_read_chkpt", "Error reading u_matrix_opt from "//TRIM(seedname)//".chk", 1)
125   CALL errore("w90_read_chkpt", "Error reading u_matrix from "//TRIM(seedname)//".chk", 1)
126   CALL errore("w90_read_chkpt", "Error reading m_matrix from "//TRIM(seedname)//".chk", 1)
127   CALL errore("w90_read_chkpt", "Error reading wannier_centres from "//TRIM(seedname)//".chk", 1)
128   CALL errore("w90_read_chkpt", "Error reading wannier_spreads from "//TRIM(seedname)//".chk", 1)
    !
  END SUBROUTINE w90_read_chkpt
  !
  !----------------------------------------------------------------------------
END MODULE w90_interface
