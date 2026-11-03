!
! Copyright (C) 2003-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE wannier_decompose
  !-----------------------------------------------------------------------
  !
  !! Decomposition of Wannier-function densities about their centres onto
  !! Loewdin-orthonormalised Gaussian radial functions times real spherical
  !! harmonics (Himanen et al., Comput. Phys. Commun. 247, 106949 (2020)).
  !! Entry point: compute_decompose, called by pw2wannier90.x with
  !! wan_mode = 'decompose'.
  !!
  !! Radial part: phi_nl(r) = r^l exp(-alpha_nl r^2), with alpha_nl fixed
  !! by phi_nl(r_n) = 1e-3 for r_n equally spaced on [r_min, r_max], and
  !! g_nl = sum_n' beta_n'nl phi_n'l the orthonormalised combination.
  !! All lengths in Angstrom.
  !!
  !! Angular part: real spherical harmonics with the Condon-Shortley
  !! phase, evaluated at the inverted direction -r^, as in the descriptors
  !! of Schubert et al., npj Comput. Mater. 10, 299 (2024). They equal
  !! (-1)^(l+m) times the harmonics returned by ylmr2.
  !!
  !! Coefficients are ordered with n outermost, then l, then
  !! m = -l..l: index (n-1)*(l_max+1)^2 + l^2 + l + m + 1.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : fpi, sqrtpi
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: compute_decompose
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE compute_decompose
    !---------------------------------------------------------------------
    !
    !! Decompose each Wannier-function density |w_n(r)|^2 about its own
    !! centre onto the basis of this module. Called with
    !! wan_mode = 'decompose' after a wannier90.x run that wrote
    !! seedname_u.mat (and seedname_u_dis.mat) and seedname_centres.xyz.
    !!
    !! The density is periodic over the Born-von Karman supercell,
    !! rho_n(r) = sum_q e^{iqr} rho_q,n(r), with
    !! rho_q,n(r) = 1/Nk sum_k conj(u_k,n(r)) u_k+q,n(r) and u_k,n the
    !! periodic parts of the Bloch states in the Wannier gauge. Its Fourier
    !! coefficients rho_n(q+G) give the projection analytically:
    !! c_nlm = sum_{q,G} rho_n(q+G) e^{i(q+G).c} F_nlm(q+G).
    !!
    !! If decompose_centres_file is given, the sum of all Wannier-function
    !! densities is also decomposed about each centre listed there.
    !
    USE io_global,       ONLY : stdout, ionode, ionode_id
    USE mp,              ONLY : mp_bcast, mp_sum
    USE mp_pools,        ONLY : npool, intra_pool_comm
    USE mp_world,        ONLY : world_comm
    USE cell_base,       ONLY : at, bg, omega, tpiba
    USE constants,       ONLY : BOHR_RADIUS_ANGS
    USE klist,           ONLY : xk, ngk, igk_k
    USE wvfct,           ONLY : nbnd, npwx
    USE wavefunctions,   ONLY : evc
    USE io_files,        ONLY : nwordwfc, iunwfc
    USE buffers,         ONLY : get_buffer
    USE fft_base,        ONLY : dffts
    USE fft_interfaces,  ONLY : fwfft, invfft
    USE gvect,           ONLY : g, mill
    USE gvecs,           ONLY : ngms
    USE noncollin_module,ONLY : noncolin
    USE control_flags,   ONLY : gamma_only
    USE uspp,            ONLY : okvan
    USE wannier,         ONLY : seedname, iknum, ikstart, num_bands, excluded_band, &
                                n_wannier, decompose_n_max, decompose_l_max, &
                                decompose_r_min, decompose_r_max, decompose_centres_file
    USE wannier90_readback,   ONLY : w90rb_read_umat, w90rb_read_centres_xyz
    !
    IMPLICIT NONE
    !
    INTEGER :: num_wann
    !! Number of Wannier functions
    INTEGER :: n_ext
    !! Number of centres in decompose_centres_file
    COMPLEX(DP), ALLOCATABLE :: umat(:,:,:), udis(:,:,:)
    !! Gauge matrices U(k) and, with disentanglement, U_dis(k)
    REAL(DP), ALLOCATABLE :: centres(:,:)
    !! Wannier centres followed by the extra centres, angstrom (3, num_wann+n_ext)
    REAL(DP), ALLOCATABLE :: xk_cryst(:,:)
    !! k points of this spin, crystal coordinates
    COMPLEX(DP), ALLOCATABLE :: evc_w(:,:,:), evc_incl(:,:)
    !! Bloch states in the Wannier gauge (npwx, num_wann, iknum)
    COMPLEX(DP), ALLOCATABLE :: coeff(:,:)
    !! Coefficients of each Wannier-function density about its own centre
    COMPLEX(DP), ALLOCATABLE :: gcoeff(:,:)
    !! Coefficients of the summed density about each extra centre
    COMPLEX(DP), ALLOCATABLE :: psic_k(:), psic_kq(:), phase(:), rhoq(:,:), fk(:)
    REAL(DP), ALLOCATABLE :: alphas(:,:), betas(:,:,:), kq(:,:), kq2(:), ylm(:,:), power(:)
    !
    REAL(DP) :: xq(3), diff(3)
    INTEGER :: iq, ik, ikq, ibnd, m, n, e, ig, nhit
    INTEGER :: npw_k, npw_kq, n_coeff, lmax2, nr, nc, iun, ierr
    CHARACTER(LEN=256) :: line
    CHARACTER(LEN=5)   :: numstr
    CHARACTER(LEN=*), PARAMETER :: coeff_order = &
         'outer n (0..n_max-1), then l (0..l_max), then inner m (0..2l)'
    !
    CALL start_clock('compute_dcp')
    !
    IF (gamma_only) CALL errore('compute_decompose', 'gamma_only not implemented', 1)
    IF (noncolin)   CALL errore('compute_decompose', 'noncollinear not implemented', 1)
    IF (okvan)      CALL errore('compute_decompose', &
         'ultrasoft/PAW not implemented: use norm-conserving pseudopotentials', 1)
    IF (npool > 1)  CALL errore('compute_decompose', 'pool parallelism not implemented', npool)
    IF (decompose_n_max < 1 .OR. decompose_l_max < 0) &
         CALL errore('compute_decompose', 'invalid decompose_n_max/l_max', 1)
    IF (decompose_r_min <= 0.0_DP .OR. decompose_r_max < decompose_r_min) &
         CALL errore('compute_decompose', 'invalid decompose_r_min/r_max', 1)
    !
    ALLOCATE(xk_cryst(3, iknum))
    xk_cryst(:,:) = xk(:, ikstart:ikstart+iknum-1)
    CALL cryst_to_cart(iknum, xk_cryst, at, -1)
    !
    ! gauge matrices: wannier90 disentangles whenever num_bands > num_wann
    !
    CALL w90rb_read_umat(TRIM(seedname)//'_u.mat', xk(:, ikstart:ikstart+iknum-1), world_comm, &
                         umat, nr, num_wann)
    IF (n_wannier > 0 .AND. num_wann /= n_wannier) CALL errore('compute_decompose', &
         'num_wann in '//TRIM(seedname)//'_u.mat differs from the nnkp projections', 1)
    IF (num_bands < num_wann) CALL errore('compute_decompose', 'num_wann > num_bands', 1)
    IF (num_bands > num_wann) THEN
      CALL w90rb_read_umat(TRIM(seedname)//'_u_dis.mat', xk(:, ikstart:ikstart+iknum-1), &
                           world_comm, udis, nr, nc)
      IF (nr /= num_bands .OR. nc /= num_wann) CALL errore('compute_decompose', &
           'dimensions of '//TRIM(seedname)//'_u_dis.mat do not match num_bands, num_wann', 1)
    ENDIF
    !
    ! centres: those of the Wannier functions, then the optional extra ones
    ! (one per line, Cartesian angstrom, blank lines and '#' comments skipped)
    !
    n_ext = 0
    IF (ionode .AND. decompose_centres_file /= ' ') THEN
      OPEN(NEWUNIT=iun, FILE=TRIM(decompose_centres_file), FORM='formatted', &
           STATUS='old', IOSTAT=ierr)
      IF (ierr /= 0) CALL errore('compute_decompose', &
           'cannot open '//TRIM(decompose_centres_file), ABS(ierr))
      DO
        READ(iun, '(a)', IOSTAT=ierr) line
        IF (ierr /= 0) EXIT
        line = ADJUSTL(line)
        IF (line /= ' ' .AND. line(1:1) /= '#') n_ext = n_ext + 1
      ENDDO
      REWIND(iun)
    ENDIF
    CALL mp_bcast(n_ext, ionode_id, world_comm)
    !
    ALLOCATE(centres(3, num_wann + n_ext))
    CALL w90rb_read_centres_xyz(TRIM(seedname)//'_centres.xyz', num_wann, &
                                centres(:, 1:num_wann), world_comm)
    IF (ionode .AND. decompose_centres_file /= ' ') THEN
      e = 0
      DO WHILE (e < n_ext)
        READ(iun, '(a)') line
        line = ADJUSTL(line)
        IF (line == ' ' .OR. line(1:1) == '#') CYCLE
        e = e + 1
        READ(line, *, IOSTAT=ierr) centres(:, num_wann + e)
        IF (ierr /= 0) CALL errore('compute_decompose', &
             'error parsing '//TRIM(decompose_centres_file), e)
      ENDDO
      CLOSE(iun)
    ENDIF
    CALL mp_bcast(centres, ionode_id, world_comm)
    !
    WRITE(stdout, '(/,5x,a,i0,a,i0)') 'Wannier-density decomposition: n_max = ', &
         decompose_n_max, ', l_max = ', decompose_l_max
    WRITE(stdout, '(5x,a,f7.4,a,f7.4,a)') 'r_min =', decompose_r_min, &
         ', r_max =', decompose_r_max, ' angstrom'
    WRITE(stdout, '(5x,i0,a,i0,a,i0,a)') num_wann, ' Wannier functions from ', &
         num_bands, ' bands, ', n_ext, ' extra centres'
    !
    ! Bloch states in the Wannier gauge
    !
    ALLOCATE(evc_w(npwx, num_wann, iknum), evc_incl(npwx, num_bands))
    DO ik = 1, iknum
      CALL get_buffer(evc, nwordwfc, iunwfc, ikstart - 1 + ik)
      m = 0
      DO ibnd = 1, nbnd
        IF (excluded_band(ibnd)) CYCLE
        m = m + 1
        evc_incl(:, m) = evc(:, ibnd)
      ENDDO
      IF (num_bands > num_wann) THEN
        evc_w(:, :, ik) = MATMUL(evc_incl, MATMUL(udis(:, :, ik), umat(:, :, ik)))
      ELSE
        evc_w(:, :, ik) = MATMUL(evc_incl, umat(:, :, ik))
      ENDIF
    ENDDO
    DEALLOCATE(evc_incl, umat)
    IF (ALLOCATED(udis)) DEALLOCATE(udis)
    !
    lmax2 = (decompose_l_max + 1)**2
    n_coeff = decompose_n_max*lmax2
    ALLOCATE(alphas(decompose_n_max, 0:decompose_l_max))
    ALLOCATE(betas(decompose_n_max, decompose_n_max, 0:decompose_l_max))
    CALL wdcp_radial_params(decompose_n_max, decompose_l_max, decompose_r_min, &
                            decompose_r_max, alphas, betas)
    !
    ALLOCATE(coeff(n_coeff, num_wann), gcoeff(n_coeff, n_ext))
    ALLOCATE(psic_k(dffts%nnr), psic_kq(dffts%nnr), phase(dffts%nnr), rhoq(dffts%nnr, num_wann))
    ALLOCATE(kq(3, ngms), kq2(ngms), ylm(ngms, lmax2), fk(n_coeff))
    coeff = (0.0_DP, 0.0_DP)
    gcoeff = (0.0_DP, 0.0_DP)
    !
    DO iq = 1, iknum
      !
      xq(:) = xk_cryst(:, iq) - xk_cryst(:, 1)
      !
      ! rho_q,n(r) = 1/Nk sum_k conj(u_k,n(r)) u_k+q,n(r)
      !
      rhoq = (0.0_DP, 0.0_DP)
      DO ik = 1, iknum
        !
        ! k + q = k' + G
        DO ikq = 1, iknum
          diff(:) = xk_cryst(:, ik) + xq(:) - xk_cryst(:, ikq)
          IF (ALL(ABS(diff - NINT(diff)) < 1.0E-5_DP)) EXIT
        ENDDO
        IF (ikq > iknum) CALL errore('compute_decompose', &
             'k+q not found in the k list: the k mesh is not a uniform grid', ik)
        !
        ! phase = e^{-iGr}, so that u_k+q = e^{-iGr} u_k'
        phase = (0.0_DP, 0.0_DP)
        nhit = 0
        DO ig = 1, ngms
          IF (ALL(mill(:, ig) == -NINT(diff))) THEN
            phase(dffts%nl(ig)) = (1.0_DP, 0.0_DP)
            nhit = 1
          ENDIF
        ENDDO
        CALL mp_sum(nhit, intra_pool_comm)
        IF (nhit /= 1) CALL errore('compute_decompose', 'k+q folding vector not found', ik)
        CALL invfft('Rho', phase, dffts)
        !
        npw_k  = ngk(ikstart - 1 + ik)
        npw_kq = ngk(ikstart - 1 + ikq)
        DO n = 1, num_wann
          psic_k = (0.0_DP, 0.0_DP)
          psic_k(dffts%nl(igk_k(1:npw_k, ikstart - 1 + ik))) = evc_w(1:npw_k, n, ik)
          CALL invfft('Wave', psic_k, dffts)
          psic_kq = (0.0_DP, 0.0_DP)
          psic_kq(dffts%nl(igk_k(1:npw_kq, ikstart - 1 + ikq))) = evc_w(1:npw_kq, n, ikq)
          CALL invfft('Wave', psic_kq, dffts)
          rhoq(:, n) = rhoq(:, n) + CONJG(psic_k)*psic_kq*phase
        ENDDO
        !
      ENDDO
      DO n = 1, num_wann
        CALL fwfft('Rho', rhoq(:, n), dffts)
      ENDDO
      !
      ! c_nlm = sum_K rho_n(K) e^{iK.c} F_nlm(K), K = q + G in angstrom^-1
      !
      CALL cryst_to_cart(1, xq, bg, 1)
      DO ig = 1, ngms
        kq(:, ig) = (xq(:) + g(:, ig))*tpiba/BOHR_RADIUS_ANGS
        kq2(ig) = SUM(kq(:, ig)**2)
      ENDDO
      CALL ylmr2(lmax2, ngms, kq, kq2, ylm)
      DO ig = 1, ngms
        CALL wdcp_ft_kernel(SQRT(kq2(ig)), ylm(ig, :), decompose_n_max, decompose_l_max, &
                            alphas, betas, fk)
        DO n = 1, num_wann
          coeff(:, n) = coeff(:, n) + rhoq(dffts%nl(ig), n)*fk(:) &
               *EXP(CMPLX(0.0_DP, DOT_PRODUCT(kq(:, ig), centres(:, n)), KIND=DP))
        ENDDO
        DO e = 1, n_ext
          gcoeff(:, e) = gcoeff(:, e) + SUM(rhoq(dffts%nl(ig), :))*fk(:) &
               *EXP(CMPLX(0.0_DP, DOT_PRODUCT(kq(:, ig), centres(:, num_wann + e)), KIND=DP))
        ENDDO
      ENDDO
      !
    ENDDO
    !
    CALL mp_sum(coeff, intra_pool_comm)
    CALL mp_sum(gcoeff, intra_pool_comm)
    !
    ! rho_n(K) are Fourier coefficients over the Born-von Karman supercell of
    ! densities that integrate to one; convert them to angstrom^-3
    !
    coeff = coeff/(iknum**2*omega*BOHR_RADIUS_ANGS**3)
    gcoeff = gcoeff/(iknum**2*omega*BOHR_RADIUS_ANGS**3)
    IF (MAXVAL(ABS(AIMAG(coeff))) > 1.0E-6_DP .OR. MAXVAL(ABS(AIMAG(gcoeff))) > 1.0E-6_DP) &
         CALL errore('compute_decompose', 'coefficients are not real', 1)
    !
    IF (ionode) THEN
      ALLOCATE(power((decompose_l_max + 1)*decompose_n_max*(decompose_n_max + 1)/2))
      DO n = 1, num_wann
        CALL wdcp_power_orb(REAL(coeff(:, n), DP), decompose_n_max, decompose_l_max, power)
        WRITE(numstr, '(i5.5)') n
        CALL write_file(TRIM(seedname)//'_'//numstr//'.coeff', REAL(coeff(:, n), DP), &
             'coefficients', coeff_order, 'coefficients')
        CALL write_file(TRIM(seedname)//'_'//numstr//'.power', power, &
             'power spectrum (orbital-orbital)', 'outer n1 (0..n_max-1), then '// &
             'n2 (n1..n_max-1), then l (0..l_max); entry = sum_m c(n1,l,m) c(n2,l,m)', &
             'entries')
      ENDDO
      DO e = 1, n_ext
        WRITE(numstr, '(i5.5)') e
        CALL write_file(TRIM(seedname)//'_gc_'//numstr//'.coeff', REAL(gcoeff(:, e), DP), &
             'coefficients', coeff_order, 'coefficients')
      ENDDO
      DEALLOCATE(power)
    ENDIF
    !
    DEALLOCATE(xk_cryst, centres, evc_w, alphas, betas, coeff, gcoeff)
    DEALLOCATE(psic_k, psic_kq, phase, rhoq, kq, kq2, ylm, fk)
    !
    CALL stop_clock('compute_dcp')
    !
  CONTAINS
    !
    !----------------------------------------------------------------------
    SUBROUTINE write_file(fname, values, title, ordering, what)
      !--------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN) :: fname, title, ordering, what
      REAL(DP), INTENT(IN) :: values(:)
      INTEGER :: iun, ierr
      !
      OPEN(NEWUNIT=iun, FILE=fname, STATUS='replace', ACTION='write', IOSTAT=ierr)
      IF (ierr /= 0) CALL errore('compute_decompose', 'cannot open '//fname, ABS(ierr))
      WRITE(iun, '(a)') '# Wannier density decomposition '//title
      WRITE(iun, '(a,i0)') '# n_max = ', decompose_n_max
      WRITE(iun, '(a,i0)') '# l_max = ', decompose_l_max
      WRITE(iun, '(a,es24.15e3)') '# r_min = ', decompose_r_min
      WRITE(iun, '(a,es24.15e3)') '# r_max = ', decompose_r_max
      WRITE(iun, '(a)') '# ordering: '//ordering
      WRITE(iun, '(a,i0)') '# number of '//what//' = ', SIZE(values)
      WRITE(iun, '(es24.15e3)') values
      CLOSE(iun)
      !
    END SUBROUTINE write_file
    !
  END SUBROUTINE compute_decompose
  !
  !----------------------------------------------------------------------
  SUBROUTINE wdcp_radial_params (n_max, l_max, r_min, r_max, alphas, betas)
    !--------------------------------------------------------------------
    !
    !! Decay coefficients alpha_nl and Loewdin coefficients beta_nn'l.
    !! The overlap is analytic,
    !! S_nn'l = 0.5 * Gamma(l+3/2) * (alpha_nl + alpha_n'l)^-(l+3/2),
    !! and beta = S^-1/2 from its eigendecomposition.
    !
    INTEGER,  INTENT(IN)  :: n_max
    !! Number of radial functions
    INTEGER,  INTENT(IN)  :: l_max
    !! Maximum angular momentum
    REAL(DP), INTENT(IN)  :: r_min
    !! Smallest threshold radius (Angstrom)
    REAL(DP), INTENT(IN)  :: r_max
    !! Largest threshold radius (Angstrom)
    REAL(DP), INTENT(OUT) :: alphas(n_max, 0:l_max)
    !! Gaussian decay coefficients (Angstrom^-2)
    REAL(DP), INTENT(OUT) :: betas(n_max, n_max, 0:l_max)
    !! Loewdin orthonormalisation coefficients
    !
    REAL(DP), PARAMETER :: thr = 1.0E-3_DP
    !
    REAL(DP) :: r_thr(n_max)
    REAL(DP) :: evec(n_max, n_max), eval(n_max), work(3*n_max)
    INTEGER  :: n, np, l, info
    !
    r_thr(1) = r_min
    DO n = 2, n_max
      r_thr(n) = r_min + (r_max - r_min)*REAL(n - 1, DP)/REAL(n_max - 1, DP)
    ENDDO
    !
    DO l = 0, l_max
      alphas(:, l) = (l*LOG(r_thr(:)) - LOG(thr))/r_thr(:)**2
    ENDDO
    IF (ANY(alphas <= 0.0_DP)) CALL errore ('wdcp_radial_params', &
         'Gaussian exponent not positive: r_min is too small for this l_max', 1)
    !
    DO l = 0, l_max
      DO np = 1, n_max
        DO n = 1, n_max
          evec(n, np) = 0.5_DP*GAMMA(l + 1.5_DP)*(alphas(n, l) + alphas(np, l))**(-l - 1.5_DP)
        ENDDO
      ENDDO
      !
      CALL dsyev ('V', 'U', n_max, evec, n_max, eval, work, 3*n_max, info)
      IF (info /= 0) CALL errore ('wdcp_radial_params', 'dsyev failed', ABS(info))
      IF (ANY(eval <= 0.0_DP)) CALL errore ('wdcp_radial_params', &
           'overlap matrix not positive definite (basis too ill-conditioned)', 1)
      !
      DO np = 1, n_max
        DO n = 1, n_max
          betas(n, np, l) = SUM(evec(n, :)*evec(np, :)/SQRT(eval(:)))
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE wdcp_radial_params
  !
  !----------------------------------------------------------------------
  SUBROUTINE wdcp_ft_kernel (kmod, ylm, n_max, l_max, alphas, betas, fk)
    !--------------------------------------------------------------------
    !
    !! Fourier transform of every basis function at one wavevector K,
    !! F_nlm(K) = int f_nlm(r) e^{iK.r} dr. From the Rayleigh expansion,
    !!   F_nlm(K) = 4 pi (-i)^l (-1)^m ylm(K^) sum_n' beta_n'nl I_n'l(K),
    !!   I_nl(K) = int_0^inf r^(l+2) e^(-alpha_nl r^2) j_l(K r) dr
    !!           = sqrt(pi) K^l e^(-K^2/(4 alpha_nl)) / (2^(l+2) alpha_nl^(l+3/2)),
    !! with ylm the harmonics of ylmr2 (see the module header for the
    !! phase factors).
    !
    REAL(DP), INTENT(IN)  :: kmod
    !! |K| (Angstrom^-1)
    REAL(DP), INTENT(IN)  :: ylm(:)
    !! Real spherical harmonics of K^, as returned by ylmr2
    INTEGER,  INTENT(IN)  :: n_max
    !! Number of radial functions
    INTEGER,  INTENT(IN)  :: l_max
    !! Maximum angular momentum
    REAL(DP), INTENT(IN)  :: alphas(n_max, 0:l_max)
    !! Gaussian decay coefficients (Angstrom^-2)
    REAL(DP), INTENT(IN)  :: betas(n_max, n_max, 0:l_max)
    !! Loewdin coefficients
    COMPLEX(DP), INTENT(OUT) :: fk(n_max*(l_max + 1)**2)
    !! F_nlm(K)
    !
    REAL(DP) :: radint(n_max)
    COMPLEX(DP) :: pref
    INTEGER  :: n, l, m, lm
    !
    DO l = 0, l_max
      radint(:) = sqrtpi*kmod**l*EXP(-kmod**2/(4.0_DP*alphas(:, l))) &
                  /(2.0_DP**(l + 2)*alphas(:, l)**(l + 1.5_DP))
      DO n = 1, n_max
        pref = fpi*(0.0_DP, -1.0_DP)**l*DOT_PRODUCT(betas(:, n, l), radint)
        DO m = -l, l
          ! ylmr2 stores m = 0, then (+m, -m) pairs for m = 1..l
          lm = l*l + 2*ABS(m) + 1
          IF (m > 0) lm = lm - 1
          fk((n - 1)*(l_max + 1)**2 + l*l + l + m + 1) = pref*(-1)**ABS(m)*ylm(lm)
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE wdcp_ft_kernel
  !
  !----------------------------------------------------------------------
  SUBROUTINE wdcp_power_orb (coeff, n_max, l_max, power)
    !--------------------------------------------------------------------
    !
    !! Rotationally-invariant power spectrum
    !! p_{n1 n2 l} = sum_m c_{n1 l m} c_{n2 l m}, n1 <= n2.
    !! Ordering: outer n1, then n2 >= n1, then l.
    !
    INTEGER,  INTENT(IN)  :: n_max
    !! Number of radial functions
    INTEGER,  INTENT(IN)  :: l_max
    !! Maximum angular momentum
    REAL(DP), INTENT(IN)  :: coeff(n_max*(l_max + 1)**2)
    !! Expansion coefficients
    REAL(DP), INTENT(OUT) :: power((l_max + 1)*n_max*(n_max + 1)/2)
    !! Power-spectrum entries
    !
    INTEGER  :: n1, n2, l, ib1, ib2, ip
    !
    ip = 0
    DO n1 = 1, n_max
      DO n2 = n1, n_max
        DO l = 0, l_max
          ib1 = (n1 - 1)*(l_max + 1)**2 + l**2
          ib2 = (n2 - 1)*(l_max + 1)**2 + l**2
          ip = ip + 1
          power(ip) = SUM(coeff(ib1 + 1:ib1 + 2*l + 1)*coeff(ib2 + 1:ib2 + 2*l + 1))
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE wdcp_power_orb
  !
END MODULE wannier_decompose
