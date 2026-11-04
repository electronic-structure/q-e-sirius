!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_average_w()
  !-----------------------------------------------------------------------
  !! BZ average of every computed Coulomb matrix element (Wannier basis only),
  !! after the q loop and at the end of crpa_spectrum_only. Each element
  !! W_ij,kl(R1,R2)(w), (i j R1) a matrix-element pair and (k l R2) a
  !! perturbation pair, is averaged on its own as sum_q wq W_q / sum_q wq
  !! with the q weights of filU.q, and written with its frequencies to
  !! filU.dir-W/W_i_j_R1a_R1b_R1c_k_l_R2a_R2b_R2c.dat, the averaged bare V in
  !! the header. With crpa_reduce_q only the star-invariant combinations
  !! (trace, pair sum) of these averages are the BZ average.
  !!
  !! The matrices are read back from filU.bare#iq and filU.scrd#iq, since
  !! the q points can be on other images or skipped by crpa_coul_restart. A
  !! missing file (e.g. start_q/last_q) skips the average.
  !!
  !! With crpa_w_rspace the same pass also accumulates the Fourier transform
  !!   W_ij,kl(R1,R2)(R3) = sum_q wq exp(-iq.R3) W_q / sum_q wq
  !!     = int w*_i0(r) w_jR1(r) W(r,r') w_kR3(r') w*_l,R3+R2(r'),
  !! the third pair translated by R3 (drho_q of a pair translated by R3 is
  !! exp(-iq.R3) drho_q). It needs the full q grid with equal weights. The
  !! q=Gamma file already carries the q=0 divergence correction
  !! (crpa_q0correction, monopole x monopole elements only), so the
  !! transform adds the same wq(Gamma)*Delta/sum wq at every R3, as RESPACK
  !! does; it is not reapplied here. R3=0 is the average above. The values
  !! are periodic in R3 over the q-grid supercell; as R1 and R2 in
  !! crpa_find_close_WF_pairs, each R3 is written at the first image
  !! (w90_find_wigner_seitz) of shortest distance, here between the pair
  !! centres (w_i+w_j+R1)/2 and (w_k+w_l+R2)/2 + R3 (from filU.wcentres), to
  !! filU.dir-WR/W_i_j_R1a_R1b_R1c_k_l_R2a_R2b_R2c.dat.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout, meta_ionode
  USE mp_world,      ONLY : world_comm
  USE mp,            ONLY : mp_barrier
  USE clib_wrappers, ONLY : f_mkdir_safe
  USE constants,     ONLY : tpi
  USE cell_base,     ONLY : at
  USE crpacom,       ONLY : filU, crpa_reduce_q, crpa_w_rspace
  USE crpa_pert,     ONLY : pert_basis
  !
  IMPLICIT NONE
  !
  INTEGER :: nq, iq, nel, nf, iel, ifreq, iun, ios, k
  INTEGER :: nqg(3), nr, ir
  INTEGER, ALLOCATABLE :: idx0(:,:), idx(:,:), r3(:,:)
  REAL(DP) :: arg
  REAL(DP), ALLOCATABLE :: xq(:,:), wq(:), val(:,:), w0(:), eta(:)
  COMPLEX(DP), ALLOCATABLE :: vb(:), w(:,:), phase(:), vbr(:,:), wr(:,:,:)
  CHARACTER(LEN=256) :: dirname, filename
  CHARACTER(LEN=6) :: int_to_char
  LOGICAL :: ok, lrsp
  !
  IF (pert_basis /= 'wannier') RETURN
  !
  ! all the q points (on every image) must be written before reading them back
  !
  CALL mp_barrier(world_comm)
  IF (.NOT. meta_ionode) RETURN
  !
  filename = TRIM(filU) // '.q'
  OPEN(NEWUNIT=iun, FILE=TRIM(filename), STATUS='old', IOSTAT=ios)
  IF (ios /= 0) THEN
     WRITE(stdout, '(/5x,"W average skipped: ",a," not found")') TRIM(filename)
     RETURN
  ENDIF
  READ(iun, *) nqg(:)
  READ(iun, *) nq
  ALLOCATE(xq(3,nq), wq(nq))
  DO iq = 1, nq
     READ(iun, *) xq(:,iq), wq(iq)
  ENDDO
  CLOSE(iun)
  !
  ! R3 transform: the canonical R3 are the nq grid-cell vectors 0..nqg-1
  ! (crystal), the readin guard (crpa_reduce_q) makes the grid full; this
  ! check is for crpa_spectrum_only on files of another run
  !
  lrsp = crpa_w_rspace
  IF (lrsp) THEN
     nr = PRODUCT(nqg)
     IF (nq /= nr .OR. MAXVAL(ABS(wq - wq(1))) > 1.0d-10 * ABS(wq(1))) THEN
        WRITE(stdout, '(/5x,"WARNING: ",a," is not the full q grid with equal weights,", &
             &" W(R3) not computed")') TRIM(filename)
        lrsp = .FALSE.
     ENDIF
  ENDIF
  IF (lrsp) THEN
     ALLOCATE(r3(3,nr), phase(nr))
     DO k = 0, nr - 1
        r3(:,k+1) = (/ k / (nqg(2)*nqg(3)), MOD(k / nqg(3), nqg(2)), MOD(k, nqg(3)) /)
     ENDDO
  ENDIF
  !
  ! weighted sums; every q must have the same elements, in the same order
  !
  DO iq = 1, nq
     CALL read_rows(TRIM(filU) // '.bare' // TRIM(int_to_char(iq)), 2, ok)
     IF (.NOT. ok) RETURN
     IF (iq == 1) THEN
        nel = SIZE(idx, 2)
        idx0 = idx
        ALLOCATE(vb(nel))
        vb = (0.0_DP, 0.0_DP)
     ELSE IF (SIZE(idx, 2) /= nel) THEN
        CALL errore('crpa_average_w', 'different elements at q #' // int_to_char(iq), 1)
     ELSE IF (ANY(idx /= idx0)) THEN
        CALL errore('crpa_average_w', 'different elements at q #' // int_to_char(iq), 1)
     ENDIF
     vb = vb + wq(iq) * CMPLX(val(1,:), val(2,:), KIND=DP)
     !
     ! exp(-iq.R3), xq Cartesian in 2pi/alat, at in alat
     !
     IF (lrsp) THEN
        IF (iq == 1) THEN
           ALLOCATE(vbr(nel,nr))
           vbr = (0.0_DP, 0.0_DP)
        ENDIF
        DO ir = 1, nr
           arg = -tpi * SUM(MATMUL(xq(:,iq), at) * DBLE(r3(:,ir)))
           phase(ir) = wq(iq) * CMPLX(COS(arg), SIN(arg), KIND=DP)
           vbr(:,ir) = vbr(:,ir) + phase(ir) * CMPLX(val(1,:), val(2,:), KIND=DP)
        ENDDO
     ENDIF
     !
     ! screened: nf rows per element, the frequency innermost
     !
     CALL read_rows(TRIM(filU) // '.scrd' // TRIM(int_to_char(iq)), 4, ok)
     IF (.NOT. ok) RETURN
     IF (iq == 1) THEN
        nf = SIZE(idx, 2) / nel
        ALLOCATE(w(nf, nel), w0(nf), eta(nf))
        w = (0.0_DP, 0.0_DP)
        w0 = val(1, 1:nf)
        eta = val(2, 1:nf)
     ENDIF
     IF (SIZE(idx, 2) /= nel * nf) &
        CALL errore('crpa_average_w', 'different frequencies at q #' // int_to_char(iq), 1)
     DO iel = 1, nel
        IF (ANY(idx(:, (iel-1)*nf+1:iel*nf) /= SPREAD(idx0(:, iel), 2, nf))) &
           CALL errore('crpa_average_w', 'bare and screened elements differ at q #' &
                       // int_to_char(iq), 1)
     ENDDO
     w = w + wq(iq) * RESHAPE(CMPLX(val(3,:), val(4,:), KIND=DP), (/ nf, nel /))
     IF (lrsp) THEN
        IF (iq == 1) THEN
           WRITE(stdout, '(/5x,"W(R3): ",i6," R3 vectors, ",f10.1," MB")') &
              nr, 16.0_DP * nf * nel * nr / 1024.0_DP**2
           ALLOCATE(wr(nf,nel,nr))
           wr = (0.0_DP, 0.0_DP)
        ENDIF
        DO ir = 1, nr
           wr(:,:,ir) = wr(:,:,ir) + phase(ir) * &
              RESHAPE(CMPLX(val(3,:), val(4,:), KIND=DP), (/ nf, nel /))
        ENDDO
     ENDIF
  ENDDO
  vb = vb / SUM(wq)
  w = w / SUM(wq)
  IF (lrsp) THEN
     vbr = vbr / SUM(wq)
     wr = wr / SUM(wq)
  ENDIF
  !
  dirname = TRIM(filU) // '.dir-W'
  IF (f_mkdir_safe(TRIM(dirname)) > 0) &
     CALL errore('crpa_average_w', 'cannot create ' // TRIM(dirname), 1)
  !
  DO iel = 1, nel
     filename = TRIM(dirname) // '/W'
     DO k = 1, 10
        filename = TRIM(filename) // '_' // TRIM(int_to_char(idx0(k, iel)))
     ENDDO
     filename = TRIM(filename) // '.dat'
     OPEN(NEWUNIT=iun, FILE=TRIM(filename), STATUS='unknown', IOSTAT=ios)
     IF (ios /= 0) CALL errore('crpa_average_w', 'opening ' // TRIM(filename), ABS(ios))
     WRITE(iun, '("# q-averaged screened interaction W_ij,kl(R1,R2): sum_q w_q W_q / sum_q w_q", &
          &" over",i5," q points")') nq
     WRITE(iun, '("# i j R1 =",2i4,2x,3i4,"   k l R2 =",2i4,2x,3i4)') idx0(:, iel)
     WRITE(iun, '("# bare V (eV):  Re =",f20.12,"  Im =",f20.12)') vb(iel)
     IF (crpa_reduce_q) WRITE(iun, '("# WARNING: crpa_reduce_q, only star-invariant combinations", &
          &" (trace, pair sum) are the BZ average")')
     WRITE(iun, '("#   w0 (Ry)    eta (Ry)           Re W (eV)           Im W (eV)")')
     DO ifreq = 1, nf
        WRITE(iun, '(2F12.6,2F20.12)') w0(ifreq), eta(ifreq), w(ifreq, iel)
     ENDDO
     CLOSE(iun, STATUS='keep')
  ENDDO
  !
  WRITE(stdout, '(/5x,"q-averaged W of",i6," elements at",i5," frequencies over",i5, &
       &" q points written to ",a)') nel, nf, nq, TRIM(dirname)
  IF (crpa_reduce_q) WRITE(stdout, '(5x,"WARNING: crpa_reduce_q, only star-invariant", &
       &" combinations (trace, pair sum) of these elements are the BZ average")')
  !
  IF (lrsp) THEN
     CALL write_rspace()
     DEALLOCATE(r3, phase, vbr, wr)
  ENDIF
  !
  RETURN
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE write_rspace()
    !-----------------------------------------------------------------------
    !! Writes filU.dir-WR. With the convention of crpa_find_close_WF_pairs
    !! for R1 and R2, each canonical R3 is labelled R3+T, T the first
    !! q-grid supercell vector of w90_find_wigner_seitz minimising the
    !! distance |c_B + R3 + T - c_A| between the pair centres. Without
    !! filU.wcentres the centres are taken at 0.
    !
    USE cell_base,  ONLY : alat
    USE constants,  ONLY : BOHR_RADIUS_ANGS
    USE w90_wigner, ONLY : w90_wigner_allocate, w90_wigner_deallocate, &
                           w90_find_wigner_seitz, irvec
    !
    IMPLICIT NONE
    !
    INTEGER :: nwan, iw, jw, nrr
    REAL(DP) :: ca(3), cb(3), rcart(3), dmin
    REAL(DP), ALLOCATABLE :: tau(:,:)
    !
    ! Wannier centres (alat), written by crpa_setup_wannier
    !
    filename = TRIM(filU) // '.wcentres'
    OPEN(NEWUNIT=iun, FILE=TRIM(filename), STATUS='old', IOSTAT=ios)
    IF (ios == 0) THEN
       READ(iun, *)
       READ(iun, *) nwan
       ALLOCATE(tau(3, nwan))
       DO iw = 1, nwan
          READ(iun, *) jw, tau(:, iw)
       ENDDO
       CLOSE(iun)
       IF (MAXVAL(idx0(1:2, :)) > nwan .OR. MAXVAL(idx0(6:7, :)) > nwan) &
          CALL errore('crpa_average_w', 'Wannier index larger than num_wann of ' // TRIM(filename), 1)
    ELSE
       WRITE(stdout, '(5x,"WARNING: ",a," not found, R3 folded with the pair centres at 0")') &
          TRIM(filename)
       nwan = MAX(MAXVAL(idx0(1:2, :)), MAXVAL(idx0(6:7, :)))
       ALLOCATE(tau(3, nwan))
       tau = 0.0_DP
    ENDIF
    !
    dirname = TRIM(filU) // '.dir-WR'
    IF (f_mkdir_safe(TRIM(dirname)) > 0) &
       CALL errore('crpa_average_w', 'cannot create ' // TRIM(dirname), 1)
    !
    CALL w90_wigner_allocate(nqg(1), nqg(2), nqg(3))
    !
    DO iel = 1, nel
       !
       ! pair centres c_A = (t_i + t_j + R1)/2, c_B = (t_k + t_l + R2)/2
       !
       rcart = DBLE(idx0(3:5, iel))
       CALL cryst_to_cart(1, rcart, at, +1)
       ca = 0.5_DP * (tau(:, idx0(1, iel)) + tau(:, idx0(2, iel)) + rcart)
       rcart = DBLE(idx0(8:10, iel))
       CALL cryst_to_cart(1, rcart, at, +1)
       cb = 0.5_DP * (tau(:, idx0(6, iel)) + tau(:, idx0(7, iel)) + rcart)
       !
       filename = TRIM(dirname) // '/W'
       DO k = 1, 10
          filename = TRIM(filename) // '_' // TRIM(int_to_char(idx0(k, iel)))
       ENDDO
       filename = TRIM(filename) // '.dat'
       OPEN(NEWUNIT=iun, FILE=TRIM(filename), STATUS='unknown', IOSTAT=ios)
       IF (ios /= 0) CALL errore('crpa_average_w', 'opening ' // TRIM(filename), ABS(ios))
       WRITE(iun, '("# W_ij,kl(R1,R2)(R3) = sum_q w_q exp(-iq.R3) W_q / sum_q w_q over",i5, &
            &" q points")') nq
       WRITE(iun, '("#   = int w*_i(r) w_j(r-R1) W(r,r'') w_k(r''-R3) w*_l(r''-R3-R2)")')
       WRITE(iun, '("# i j R1 =",2i4,2x,3i4,"   k l R2 =",2i4,2x,3i4)') idx0(:, iel)
       WRITE(iun, '("# R3 (crystal) at the image of shortest pair-centre distance d")')
       WRITE(iun, '("# q=0 correction: the Gamma-point Delta (monopole x monopole elements),", &
            &" the same at every R3")')
       WRITE(iun, '("# R3a R3b R3c    d (A)     w0 (Ry)    eta (Ry)", &
            &"           Re V (eV)           Im V (eV)           Re W (eV)           Im W (eV)")')
       DO ir = 1, nr
          rcart = DBLE(r3(:, ir))
          CALL cryst_to_cart(1, rcart, at, +1)
          CALL w90_find_wigner_seitz(ca - cb - rcart, nrr, dmin)
          DO ifreq = 1, nf
             WRITE(iun, '(3i5,f12.6,2f12.6,4f20.12)') r3(:, ir) + irvec(:, 1), &
                dmin * alat * BOHR_RADIUS_ANGS, w0(ifreq), eta(ifreq), vbr(iel, ir), wr(ifreq, iel, ir)
          ENDDO
       ENDDO
       CLOSE(iun, STATUS='keep')
       !
    ENDDO
    !
    CALL w90_wigner_deallocate()
    DEALLOCATE(tau)
    !
    WRITE(stdout, '(5x,"W(R3) of",i6," elements at",i6," R3 written to ",a)') nel, nr, TRIM(dirname)
    WRITE(stdout, '(5x,"q=0 correction: the Gamma-point Delta enters every R3 as a constant", &
         &" (monopole x monopole elements only)")')
    !
  END SUBROUTINE write_rspace
  !
  !-----------------------------------------------------------------------
  SUBROUTINE read_rows(fname, nreal, ok)
    !-----------------------------------------------------------------------
    !! Reads the rows iw1 jw1 R1(3) iw2 jw2 R2(3) and nreal reals of the
    !! Coulomb-matrix file fname (after its '#' header) into idx and val.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: nreal
    LOGICAL, INTENT(OUT) :: ok
    !
    INTEGER :: iunf, nrow, irow
    CHARACTER(LEN=256) :: line
    !
    OPEN(NEWUNIT=iunf, FILE=fname, STATUS='old', IOSTAT=ios)
    ok = (ios == 0)
    IF (.NOT. ok) THEN
       WRITE(stdout, '(/5x,"W average skipped: ",a," not found")') fname
       RETURN
    ENDIF
    !
    nrow = 0
    DO
       READ(iunf, '(a)', IOSTAT=ios) line
       IF (ios /= 0) EXIT
       IF (line(1:1) /= '#') nrow = nrow + 1
    ENDDO
    IF (ALLOCATED(idx)) DEALLOCATE(idx, val)
    ALLOCATE(idx(10, nrow), val(nreal, nrow))
    !
    REWIND(iunf)
    READ(iunf, *)
    DO irow = 1, nrow
       READ(iunf, *, IOSTAT=ios) idx(:, irow), val(:, irow)
       IF (ios /= 0) CALL errore('crpa_average_w', 'unreadable line in ' // fname, ABS(ios))
    ENDDO
    CLOSE(iunf)
    !
  END SUBROUTINE read_rows
  !
END SUBROUTINE crpa_average_w
