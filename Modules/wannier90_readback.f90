!
! Copyright (C) 2003-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE wannier90_readback
  !-----------------------------------------------------------------------
  !
  !! Readers for files produced by a wannier90.x run: the gauge (U)
  !! matrices written with write_u_matrices = .true. and the Wannier
  !! centres written with write_xyz = .true.
  !!
  !! Used by KCW, by pw2wannier90's decompose mode and by wannier2pw.
  !! Reading happens on ionode and the result is broadcast over the
  !! communicator supplied by the caller.
  !
  USE kinds,       ONLY : DP
  USE cell_base,   ONLY : bg
  USE io_global,   ONLY : ionode, ionode_id
  USE mp,          ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: w90rb_read_umat
  PUBLIC :: w90rb_read_centres_xyz
  !
CONTAINS
  !
  !----------------------------------------------------------------------
  SUBROUTINE w90rb_read_umat (filename, xk, comm, umat, nrows, ncols)
    !--------------------------------------------------------------------
    !
    !! Read one wannier90 U-matrix file (<seed>_u.mat, <seed>_u_dis.mat,
    !! <seed>_emp_u.mat, ...). File layout (wannier90 write_u_matrices):
    !! a header/date line; a line with num_kpts, num_wann, num_rows
    !! (num_rows = num_wann for u.mat, the number of bands for
    !! u_dis.mat); then, per k point, the k coordinates in crystal units
    !! followed by the matrix elements, row (first) index innermost, in
    !! format (f15.10,sp,f15.10).
    !!
    !! The k-point list of the file is checked against xk(3, nk); any
    !! mismatch is fatal.
    !! umat is allocated here to (nrows, ncols, nk) and broadcast.
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename
    !! File to read
    REAL(DP), INTENT(IN) :: xk(:, :)
    !! Expected k points, Cartesian coordinates in units of 2pi/alat (3, nk)
    INTEGER, INTENT(IN) :: comm
    !! Communicator for the broadcasts
    COMPLEX(DP), ALLOCATABLE, INTENT(OUT) :: umat(:, :, :)
    !! The matrices, (nrows, ncols, nk)
    INTEGER, INTENT(OUT) :: nrows
    !! Leading (band-like) dimension from the file header
    INTEGER, INTENT(OUT) :: ncols
    !! Second (Wannier-like) dimension from the file header
    !
    REAL(DP) :: xk_file(3)
    INTEGER :: nk, ik, i, j, nk_file, iun, ierr
    !
    nk = SIZE(xk, 2)
    !
    IF (ionode) THEN
      OPEN (NEWUNIT=iun, FILE=TRIM(filename), FORM='formatted', STATUS='old', IOSTAT=ierr)
      IF (ierr /= 0) CALL errore ('w90rb_read_umat', &
           'cannot open '//TRIM(filename)//': run wannier90 with write_u_matrices=.true.', ABS(ierr))
      READ (iun, *, IOSTAT=ierr)
      IF (ierr == 0) READ (iun, *, IOSTAT=ierr) nk_file, ncols, nrows
      IF (ierr /= 0) CALL errore ('w90rb_read_umat', 'error reading header of '//TRIM(filename), 1)
      IF (nk_file /= nk) CALL errore ('w90rb_read_umat', &
           'number of k points in '//TRIM(filename)//' does not match the calculation', nk)
    ENDIF
    CALL mp_bcast (nrows, ionode_id, comm)
    CALL mp_bcast (ncols, ionode_id, comm)
    !
    ALLOCATE (umat(nrows, ncols, nk))
    !
    IF (ionode) THEN
      DO ik = 1, nk
        READ (iun, *, IOSTAT=ierr) xk_file
        IF (ierr /= 0) CALL errore ('w90rb_read_umat', 'error reading k point in '//TRIM(filename), ik)
        CALL cryst_to_cart (1, xk_file, bg, 1)
        IF (MAXVAL(ABS(xk_file(:) - xk(:, ik))) > 1.0E-6_DP) &
             CALL errore ('w90rb_read_umat', &
             'k-point mismatch between '//TRIM(filename)//' and the calculation', ik)
        READ (iun, '(f15.10,sp,f15.10)', IOSTAT=ierr) ((umat(i, j, ik), i = 1, nrows), j = 1, ncols)
        IF (ierr /= 0) CALL errore ('w90rb_read_umat', 'error reading matrix in '//TRIM(filename), ik)
      ENDDO
      CLOSE (iun)
    ENDIF
    CALL mp_bcast (umat, ionode_id, comm)
    !
  END SUBROUTINE w90rb_read_umat
  !
  !----------------------------------------------------------------------
  SUBROUTINE w90rb_read_centres_xyz (filename, ncentres, centres, comm, found)
    !--------------------------------------------------------------------
    !
    !! Read the Wannier centres from a wannier90 .xyz file
    !! (<seed>_centres.xyz, written with write_xyz = .true.): first line
    !! total entry count, second line comment, then one 'X  x y z' line
    !! per Wannier centre (Cartesian, angstrom) followed by the atom
    !! lines. Exactly ncentres 'X' lines are expected: fewer, or an
    !! extra 'X' line after them, is fatal (a mismatch with num_wann).
    !!
    !! Centres are returned in Cartesian angstrom (no unit conversion)
    !! and broadcast over comm.
    !!
    !! If the optional argument `found` is present, a missing file sets
    !! found = .FALSE. and returns (centres untouched) instead of
    !! erroring.
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename
    !! File to read
    INTEGER, INTENT(IN) :: ncentres
    !! Number of Wannier-centre ('X') lines expected
    REAL(DP), INTENT(INOUT) :: centres(3, ncentres)
    !! The centres, Cartesian angstrom
    INTEGER, INTENT(IN) :: comm
    !! Communicator for the broadcasts
    LOGICAL, OPTIONAL, INTENT(OUT) :: found
    !! When present: .FALSE. + early return if the file does not exist
    !
    CHARACTER(LEN=256) :: line
    INTEGER :: n, iun, ierr
    LOGICAL :: exst
    !
    IF (PRESENT(found)) THEN
      IF (ionode) INQUIRE (FILE=TRIM(filename), EXIST=exst)
      CALL mp_bcast (exst, ionode_id, comm)
      found = exst
      IF (.NOT. exst) RETURN
    ENDIF
    !
    IF (ionode) THEN
      OPEN (NEWUNIT=iun, FILE=TRIM(filename), FORM='formatted', STATUS='old', IOSTAT=ierr)
      IF (ierr /= 0) CALL errore ('w90rb_read_centres_xyz', &
           'cannot open '//TRIM(filename)//': run wannier90 with write_xyz=.true.', ABS(ierr))
      READ (iun, *, IOSTAT=ierr)   ! entry count (unused)
      IF (ierr == 0) READ (iun, *, IOSTAT=ierr)   ! comment line
      IF (ierr /= 0) CALL errore ('w90rb_read_centres_xyz', 'error reading header of '//TRIM(filename), 1)
      DO n = 1, ncentres
        READ (iun, '(a)', IOSTAT=ierr) line
        IF (ierr /= 0) CALL errore ('w90rb_read_centres_xyz', &
             'unexpected end of '//TRIM(filename), n)
        IF (line(1:1) /= 'X') CALL errore ('w90rb_read_centres_xyz', &
             'X must precede each Wannier centre line in '//TRIM(filename), n)
        READ (line(2:), *, IOSTAT=ierr) centres(:, n)
        IF (ierr /= 0) CALL errore ('w90rb_read_centres_xyz', &
             'error parsing centre in '//TRIM(filename), n)
      ENDDO
      ! One more 'X' line here would mean ncentres does not match num_wann
      READ (iun, '(a)', IOSTAT=ierr) line
      IF (ierr == 0 .AND. line(1:1) == 'X') CALL errore ('w90rb_read_centres_xyz', &
           'more Wannier centres than expected in '//TRIM(filename)//': check num_wann', 1)
      CLOSE (iun)
    ENDIF
    CALL mp_bcast (centres, ionode_id, comm)
    !
  END SUBROUTINE w90rb_read_centres_xyz
  !
END MODULE wannier90_readback
