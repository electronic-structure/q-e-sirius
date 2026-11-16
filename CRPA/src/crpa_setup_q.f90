!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_setup_q()
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares several variables which are needed in the
  !  CRPA program at fixed q point:
  !  1) computes the total local potential (external+scf) on the smooth
  !     grid to be used in h_psi (and other places?)
  !  2) Set the nonlinear core correction variable
  !  3) Allocate the variable of the local magnetization
  !  4) Compute the derivative of the XC potential (dmuxc)
  !  5) Setup gradient correction stuff
  !  6) Compute the inverse of each matrix of the crystal symmetry group
  !  7) Computes the number of occupied bands for each k point
  !  8) Compute alpha_pv
  !  9) Set various symmetry-related variables
  !       time_reversal true if there is time-reversal symmetry
  !       gi            the G associated to each symmetry operation
  !       gimq          the G of the q -> -q+G symmetry
  !       nsymq         the order of the small group of q
  !       irotmq        the index of the q->-q+G symmetry
  !       minus_q       true if there is a symmetry sending q -> -q+G
  !       rtau          rtau = S\tau_a - \tau_b
  ! 10) Setup the parameters alpha_mix needed for the
  !     solution of the linear system
  ! 11) Initialize d1, d2, d3 to rotate the spherical harmonics
  !
  !  IMPORTANT NOTE ABOUT SYMMETRIES:
  !  nrot  is the number of sym.ops. of the Bravais lattice
  !        read from data file, only used in set_default_pw
  !  nsym  is the number of sym.ops. of the crystal symmetry group
  !        read from data file, should never be changed
  !  nsymq is the number of sym.ops. of the small group of q
  !        it is calculated in set_defaults_pw for each q
  !  The matrices "s" of sym.ops are ordered as follows:
  !   first the nsymq sym.ops. of the small group of q
  !   (the ordering is done in subroutine copy_sym in set_defaults_pw),
  !   followed by the remaining nsym-nsymq sym.ops. of the crystal group,
  !   followed by the remaining nrot-nsym sym.ops. of the Bravais  group
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,        ONLY : at, bg
  USE io_global,        ONLY : stdout
  USE lsda_mod,         ONLY : nspin, starting_magnetization
  USE scf,              ONLY : v, vrs, vltot, rho, kedtau
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE gvecs,            ONLY : doublegrid
  USE symm_base,        ONLY : nrot, nsym, s, ft, irt, time_reversal, &
                               inverse_s, d1, d2, d3
  USE uspp_param,       ONLY : upf
  USE uspp,             ONLY : nlcc_any, okvan, deeq_nc
  USE constants,        ONLY : degspin, pi, rytoev
  USE noncollin_module, ONLY : noncolin, domag, m_loc, nspin_mag, &
                               angle1, angle2, ux
  USE wvfct,            ONLY : nbnd, et
  USE control_flags,    ONLY : noinv
  USE eqv,              ONLY : dmuxc
  USE qpoint,           ONLY : xq
  USE control_lr,       ONLY : lgamma, flmixdpot, reduce_io, niter_ph, alpha_mix
  USE lr_symm_base,     ONLY : gi, gimq, irotmq, minus_q, invsymq, nsymq, rtau
  USE control_flags,    ONLY : modenum
  USE lr_nc_mag,        ONLY : deeq_nc_save
  USE dfunct,           ONLY : newd
  USE Coul_cut_2D,      ONLY : do_cutoff_2D
  USE Coul_cut_2D_ph,   ONLY : cutoff_fact_qg
  USE crpacom,          ONLY : crpa_nfreq, crpa_freq
  USE klist,            ONLY : ltetra
  !
  IMPLICIT NONE
  INTEGER :: ir, isym, ik, it, na
  LOGICAL :: sym(48), magnetic_sym
  !
  CALL start_clock ('crpa_setup_q')
  !
  ! 1) Compute the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! 2) Set the nonlinear core correction variable
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  ! 3) Allocate the variable of the local magnetization.
  !    This is needed in find_sym
  !
  IF (.NOT.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  !
  IF (noncolin.and.domag) THEN
     DO na = 1, nat
       !
       m_loc(1,na) = starting_magnetization(ityp(na)) * &
                     SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
       m_loc(2,na) = starting_magnetization(ityp(na)) * &
                     SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
       m_loc(3,na) = starting_magnetization(ityp(na)) * &
                     COS( angle1(ityp(na)) )
     END DO
     ux=0.0_DP
     !
     ! Change the sign of the magnetic field in the screened US coefficients
     ! and save also the coefficients computed with -B_xc.
     !
     IF (okvan) THEN
       deeq_nc_save(:,:,:,:,1)=deeq_nc(:,:,:,:)
       v%of_r(:,2:4)=-v%of_r(:,2:4)
       CALL newd()
       v%of_r(:,2:4)=-v%of_r(:,2:4)
       deeq_nc_save(:,:,:,:,2)=deeq_nc(:,:,:,:)
       deeq_nc(:,:,:,:)=deeq_nc_save(:,:,:,:,1)
     ENDIF
  ENDIF
  !
  ! 4) Compute the derivative of the XC potential (dmuxc)
  !
  CALL setup_dmuxc()
  !
  ! 5) Setup gradient correction stuff
  !
  CALL setup_dgc()
  !
  ! 6) Compute the inverse of each matrix of the crystal symmetry group
  !
  CALL inverse_s()
  !
  ! 7) Computes the number of occupied bands for each k point
  !
  CALL setup_nbnd_occ()
  !
  ! 8) Set various symmetry-related variables
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  ! The small group of q was already determined. At q\=0 it is calculated
  ! by set_nscf, at q=0 it coincides with the point group and we take nsymq=nsym
  !
  IF (lgamma) THEN
     !
     nsymq   = nsym
     !
     ! Note: condition on Re(w) is introduced to avoid forced real induced density
     ! symmetrization from symdvscf.f90, which is not true when Re(w)/= 0
     ! (for any of the frequencies)
     !
     IF ( time_reversal .AND. MAXVAL(ABS(DBLE(crpa_freq))) < 1.0d-10 ) THEN
         minus_q = .TRUE.
     ELSE
         minus_q = .FALSE.
     ENDIF
     !
  ENDIF
  !
  ! Same guard at q /= 0, where minus_q comes from the geometry (S q = -q+G)
  ! alone: symdvscf's minus_q branch averages drho with the conjugate of its
  ! rotated copy, i.e. drho(w) with drho(-w*) = conj(drho(w)), which keeps
  ! only Re(chi). Exact for real or purely imaginary w, wrong once Re(w)/=0.
  !
  IF (.NOT. lgamma .AND. MAXVAL(ABS(DBLE(crpa_freq))) > 1.0d-10) minus_q = .FALSE.
  !
  ! 9) Compute alpha_pv. With tetrahedra this also initializes the tetrahedra
  ! (dfpt_tetra_main -> opt_tetra_init), which needs nsymq, time_reversal and
  ! minus_q, so it must come after they are set - as in PHonon/PH/phq_setup.f90,
  ! which asserts nsymq/=0 on entry. No effect for lgauss or insulators.
  !
  CALL setup_alpha_pv()
  !
  ! Tetrahedra, static response at q = Gamma: make the Fermi-surface (Drude)
  ! weights of dfpt_tetra_mod independent of the basis inside degenerate
  ! multiplets. Every other case keeps the dfpt_tetra_mod weights unchanged.
  !
  IF (ltetra .AND. lgamma .AND. crpa_nfreq == 1 .AND. ABS(crpa_freq(1)) < 1.0d-10) &
     CALL crpa_tetra_degenerate_drude()
  !
  ! Calculate rtau (the Bravais lattice vector associated to a rotation)
  ! with the new symmetry order
  !
  CALL sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)
  !
  ! Calculate the vectors G associated to the symmetry Sq = q + G
  ! If minus_q=.true. calculate also irotmq and the G associated to Sq=-q+G
  !
  CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
  !
  ! 10) Setup the parameters alpha_mix
  !
  DO it = 2, niter_ph
     IF (alpha_mix(it).eq.0.d0) alpha_mix(it) = alpha_mix(it - 1)
  ENDDO
  !
  ! 11) Since the order of the S matrices is changed (for q\=0)
  !     we need to re-initialize d1, d2, d3 to rotate the spherical harmonics
  !
  CALL d_matrix( d1, d2, d3 )
  !
  ! Set flmixdpot
  !
  IF (reduce_io) THEN
     flmixdpot = ' '
  ELSE
     flmixdpot = 'mixd'
  ENDIF
  !
  ! For 2d calculations, we need to initialize the fact for the q+G
  ! component of the cutoff of the Coulomb interaction
  !
  IF (do_cutoff_2D) call cutoff_fact_qg()
  !
  CALL stop_clock ('crpa_setup_q')
  !
  RETURN
  !
END SUBROUTINE crpa_setup_q
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_tetra_degenerate_drude()
  !-----------------------------------------------------------------------
  !! Called by crpa_setup_q for tetrahedra with a static response (w_freq = 0)
  !! at q = Gamma only.
  !!
  !! dfpt_tetra_main adds the Fermi-surface (Drude) term -delta_i/2 (times
  !! alpha_pv) to the diagonal pair weight dfpt_tetra_beta(i,i) only, while the
  !! other pairs (i,j) of a degenerate multiplet at the same k keep the
  !! interband value. With a diagonal weight a different from the off-diagonal
  !! one b the response contains (a-b)*sum_i |<i|dV|i>|^2, which depends on the
  !! (arbitrary, layout-dependent) basis the eigensolver picked inside the
  !! multiplet. In the q -> 0 limit (f_i-f_j)/(e_i-e_j) -> -delta(e_F-e) for
  !! every pair of a degenerate block, which is what Gaussian smearing does
  !! (orthogonalize: theta = 1/2, w0gauss for all degenerate pairs). The same is
  !! done here: inside each degenerate block every pair (i,j) gets the block
  !! mean of the diagonal weights, and delta_i (used by the Fermi-energy shift
  !! and by the diagonal) its block mean. Non-degenerate bands are unchanged.
  !
  USE kinds,          ONLY : DP
  USE wvfct,          ONLY : nbnd, et
  USE qpoint,         ONLY : nksq, ikks
  USE dfpt_tetra_mod, ONLY : dfpt_tetra_beta, dfpt_tetra_delta
  USE io_global,      ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: thr = 1.0d-5   ! Ry, degeneracy threshold
  REAL(DP) :: bdiag, dbar
  INTEGER :: ik, ikk, i1, i2, n, ib
  !
  WRITE(stdout, '(/,5x,"WARNING: tetrahedra, static response at q = Gamma.",/, &
     & 5x,"The Fermi-surface (Drude) weights of degenerate multiplets are set to their",/, &
     & 5x,"block mean, as with Gaussian smearing, to remove the dependence on the basis",/, &
     & 5x,"inside the multiplets (MPI layout, diagonalizer). This treatment has not been",/, &
     & 5x,"cross-validated; Gaussian smearing is the safer choice for static U(q=Gamma).")')
  !
  DO ik = 1, nksq
     ikk = ikks(ik)
     i1 = 1
     DO WHILE (i1 <= nbnd)
        i2 = i1
        DO WHILE (i2 < nbnd)
           IF (ABS(et(i2+1,ikk) - et(i1,ikk)) > thr) EXIT
           i2 = i2 + 1
        ENDDO
        n = i2 - i1 + 1
        IF (n > 1) THEN
           bdiag = 0.0_DP
           DO ib = i1, i2
              bdiag = bdiag + dfpt_tetra_beta(ib, ib, ikk)
           ENDDO
           dfpt_tetra_beta(i1:i2, i1:i2, ikk) = bdiag / DBLE(n)
           dbar = SUM(dfpt_tetra_delta(i1:i2, ikk)) / DBLE(n)
           dfpt_tetra_delta(i1:i2, ikk) = dbar
        ENDIF
        i1 = i2 + 1
     ENDDO
  ENDDO
  !
END SUBROUTINE crpa_tetra_degenerate_drude
