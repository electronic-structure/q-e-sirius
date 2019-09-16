  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE elphon_shuffle_wrap
  !-----------------------------------------------------------------------
  !!
  !! Electron-phonon calculation with Wannier functions: load all phonon q's
  !!
  !! This subroutine is the main driver of the electron-phonon 
  !! calculation. It first calculates the electron-phonon matrix elements
  !! on the coarse mesh and then passes the data off to [[ephwann_shuffle]]
  !! to perform the interpolation.
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm, &
                            npool, inter_image_comm, world_comm  
  USE mp_images,     ONLY : my_image_id, nimage
  USE mp_world,      ONLY : mpime
  USE mp,            ONLY : mp_barrier, mp_bcast
  USE io_global,     ONLY : stdout, meta_ionode, meta_ionode_id
  USE us,            ONLY : nqxq, dq, qrad
  USE gvect,         ONLY : gcutm, ngm
  USE cellmd,        ONLY : cell_factor
  USE uspp_param,    ONLY : lmaxq, nbetam
  USE io_files,      ONLY : prefix, tmp_dir
  USE wavefunctions, ONLY : evc
  USE wvfct,         ONLY : npwx
  USE eqv,           ONLY : vlocq, dmuxc
  USE ions_base,     ONLY : nat, nsp, tau, ityp
  USE control_flags, ONLY : iverbosity
  USE io_epw,        ONLY : iuepb, iuqpeig
  USE pwcom,         ONLY : nks, nbnd, nkstot
  USE cell_base,     ONLY : at, bg
  USE symm_base,     ONLY : irt, s, nsym, ft, sname, invs, s_axis_to_cart, &
                            sr, nrot, copy_sym, set_sym_bl, find_sym, & 
                            inverse_s, remove_sym, allfrac
  USE start_k,       ONLY : nk1, nk2, nk3
  USE phcom,         ONLY : dpsi, dvpsi, evq, nq1, nq3, nq2 
  USE qpoint,        ONLY : igkq, xq, eigqts
  USE modes,         ONLY : nmodes, u, npert
  USE lr_symm_base,  ONLY : minus_q, rtau, gi, gimq, irotmq, nsymq, invsymq
  USE epwcom,        ONLY : epbread, epbwrite, epwread, lifc, etf_mem, vme, &
                            nbndsub, iswitch, kmaps, eig_read, dvscf_dir, lpolar
  USE elph2,         ONLY : epmatq, dynq, sumr, et_ks, &
                            zstar, epsi, cu, cuq, lwin, lwinq, bmat, igk_k_all, &
                            ngk_all, exband, wscache, umat, umat_all
  USE klist_epw,     ONLY : xk_all, et_loc, et_all
  USE constants_epw, ONLY : ryd2ev, zero, czero
  USE fft_base,      ONLY : dfftp
  USE control_ph,    ONLY : u_from_file
  USE noncollin_module, ONLY : m_loc, npol, noncolin
  USE iotk_module,   ONLY : iotk_open_read, iotk_scan_dat, iotk_free_unit, &
                            iotk_close_read
  USE division,      ONLY : fkbounds
  USE uspp,          ONLY : okvan
  USE spin_orb,      ONLY : lspinorb 
  USE lrus,          ONLY : becp1
  USE becmod,        ONLY : becp, deallocate_bec_type
  USE phus,          ONLY : int1, int1_nc, int2, int2_so, &
                            int4, int4_nc, int5, int5_so, alphap

#if defined(__NAG)
  USE f90_unix_io,   ONLY : flush
#endif
  !
  ! --------------------------------------------------------------
  !
  IMPLICIT NONE
  ! 
  INTEGER :: sym_smallq(48) 
  !! Set of all symmetries for the small group of one q.
  !! This is a subset of total crystal symmetries that remains
  !! after the q-point pertubation.
  INTEGER :: nqc_irr
  !! Number of qpoints in the irreducible wedge
  INTEGER :: nqc
  !! Number of qpoints on the uniform grid
  INTEGER :: maxvalue
  !! Temporary integer for max value
  INTEGER :: nqxq_tmp
  !! Maximum G+q length  
  INTEGER :: ibnd
  !! Band index
  INTEGER :: ik
  !! Total k-point index
  INTEGER :: ios
  !! Contains the state of the opened file 
  INTEGER :: dummy1
  !! Dummy variable
  INTEGER :: dummy2
  !! Dummy variable
  INTEGER :: ik_start
  !! Lower bound for the k-point of the coarse grid in parallel 
  INTEGER :: ik_stop
  !! Higher bound for the k-point of the coarse grid in parallel 
  INTEGER :: gmapsym(ngm,48)
  !! Correspondence G -> S(G)
  INTEGER :: nq
  !! Degeneracy of the star of q
  INTEGER :: isq(48)
  !! Index of q in the star of a given sym.op.
  INTEGER :: imq              
  !! Index of -q in the star of q (0 if not present)
  INTEGER :: sym_sgq(48)
  !! The symmetries giving the q point iq in the star
  INTEGER :: i
  !! Index for the star of q points
  INTEGER :: j
  !! Cartesian index
  INTEGER :: iq 
  !! q-point index
  INTEGER :: iq_irr
  !! Counter on irreducible q-points
  INTEGER :: isym
  !! Index of symmetry
  INTEGER :: iq_first
  !! First q in the star of q
  INTEGER :: jsym
  !! Symmetry index 
  INTEGER :: ism1
  !! Inverse of the symmetry
  INTEGER :: nsq
  !! The number of degeneracy of the small group for this iq in the star
  INTEGER :: ipol
  !! Polarization index
  INTEGER :: jpol
  !! Polarization index
  INTEGER :: ierr
  !! Error index when reading/writing a file
  INTEGER :: iunpun
  !! Unit of the file
  ! 
  REAL(kind=DP), ALLOCATABLE :: xqc_irr(:,:)
  !! The qpoints in the irr wedge
  REAL(kind=DP), ALLOCATABLE :: xqc(:,:)
  !! The qpoints in the uniform mesh
  REAL(kind=DP), ALLOCATABLE :: wqlist(:)
  !! The corresponding weigths
  REAL(kind=DP) :: sxq(3, 48)
  !! List of vectors in the star of q  
  REAL(kind=DP) :: et_tmp(nbnd, nkstot)
  !! Temporary array containing the eigenvalues (KS or GW) when read from files
  REAL(kind=DP) :: xq0(3) 
  !! Current coarse q-point coords.
  REAL(kind=DP) :: aq(3)
  !! Store the current q-point for symmetry multiplication
  REAL(kind=DP) :: saq(3)
  !! Rotated q-point
  REAL(kind=DP) :: raq(3)
  !! Rotate q-point in cartesian coordinate
  REAL(kind=DP) :: ft1
  !! Fractional translation x
  REAL(kind=DP) :: ft2
  !! Fractional translation y
  REAL(kind=DP) :: ft3
  !! Fractional translation z
  REAL(kind=DP) :: w_centers(3,nbndsub)
  !! Wannier centers
  REAL(kind=DP) :: qnorm_tmp
  !! Absolute value of xqc_irr
  !
  COMPLEX(kind=DP) :: eigv(ngm, 48)
  !! $e^{ iGv}$ for 1...nsym (v the fractional translation)
  COMPLEX(kind=DP) :: cz1(nmodes, nmodes)
  !! The eigenvectors for the first q in the star
  COMPLEX(kind=DP) :: cz2(nmodes, nmodes)
  !! The rotated eigenvectors, for the current q in the star
  !
  CHARACTER(len=256) :: tempfile
  !! Temporary .eig file
  CHARACTER(len=256) :: dirname
  !! Name of the directory
  CHARACTER(len=256) :: filename
  !! Name of the file
  CHARACTER(len=3) :: filelab
  !! Append the number of the core that works on that file
  CHARACTER(len=80)   :: line
  !! Use to read external eigenvalues
  CHARACTER(len=6), external :: int_to_char
  !! Transform an integer into a character
  !
  LOGICAL :: sym(48)
  !! Logical vectors that say which crystal symmetries exist in our system
  LOGICAL :: eqvect_strict
  !! This function tests if two tridimensional vectors are equal
  LOGICAL :: nog
  !! Find if G=0 or not in $$S(q_0)+G=q$$
  LOGICAL :: symmo
  !! Check whether the symmetry belongs to a symmorphic group
  LOGICAL :: exst
  !! Find if a file exists.
  !
  ! ---------------------------------------------------------------------
  !
  CALL start_clock( 'elphon_wrap' )
  !
  ! Read qpoint list from stdin
  !
  IF (meta_ionode) READ(5,*) nqc_irr
  CALL mp_bcast(nqc_irr, meta_ionode_id, world_comm)
  ALLOCATE (xqc_irr(3, nqc_irr))
  ALLOCATE (xqc(3, nq1 * nq2 * nq3))
  ALLOCATE (wqlist(nq1 * nq2 * nq3))
  xqc_irr(:, :) = zero
  xqc(:, :)     = zero
  wqlist(:)     = zero
  !  
  IF (meta_ionode) THEN
    DO iq_irr = 1, nqc_irr
      READ(5,*) xqc_irr(:,iq_irr)
    ENDDO
  ENDIF
  CALL mp_bcast(xqc_irr, meta_ionode_id, world_comm)
  !
  ! fix for uspp
  ! this is needed to get the correct size of the interpolation table 'qrad' 
  ! for the non-local part of the pseudopotential in PW/src/allocate_nlpot.f90
  !
  maxvalue = nqxq
  DO iq_irr = 1, nqc_irr
    qnorm_tmp = SQRT(xqc_irr(1, iq_irr)**2 + xqc_irr(2, iq_irr)**2 + &
                     xqc_irr(3, iq_irr)**2)
    nqxq_tmp = INT(((SQRT(gcutm) + qnorm_tmp) / dq + 4) * cell_factor)
    IF (nqxq_tmp > maxvalue)  maxvalue = nqxq_tmp
  ENDDO
  !
  IF (maxvalue > nqxq) THEN
    IF (ALLOCATED(qrad)) DEALLOCATE (qrad)
    ALLOCATE (qrad(maxvalue, nbetam * (nbetam + 1) / 2, lmaxq, nsp))
    qrad(:,:,:,:) = zero
    ! RM - need to call init_us_1 to re-calculate qrad 
    CALL init_us_1
  ENDIF
  ! 
  ! do not perform the check if restart
  IF (epwread .AND. .NOT. epbread) THEN
    CONTINUE
  ELSE
    IF (nkstot /= nk1 * nk2 * nk3) &
       CALL errore('elphon_shuffle_wrap','nscf run inconsistent with epw input',1)  
  ENDIF
  !
  ! Read in external electronic eigenvalues. e.g. GW 
  !
  ALLOCATE (et_ks(nbnd, nks))
  et_ks(:, :) = zero
  IF (eig_read) THEN
    IF (meta_ionode) THEN
      WRITE (stdout,'(5x,a,i5,a,i5,a)') "Reading external electronic eigenvalues (", &
           nbnd, ",", nkstot,")"
      tempfile = trim(prefix)//'.eig'
      OPEN(iuqpeig, FILE=tempfile, FORM='formatted', action='read', iostat=ios)
      IF (ios /= 0) CALL errore('elphon_shuffle_wrap','error opening' // tempfile, 1)
      READ(iuqpeig,'(a)') line
      DO ik=1, nkstot
        ! We do not save the k-point for the moment ==> should be read and
        ! tested against the current one  
        READ(iuqpeig,'(a)') line
        READ(iuqpeig,*) et_tmp(:,ik)
      ENDDO
      CLOSE(iuqpeig)
      ! from eV to Ryd
      et_tmp = et_tmp / ryd2ev
    ENDIF
    CALL mp_bcast(et_tmp, meta_ionode_id, world_comm)
    !
    CALL fkbounds(nkstot, ik_start, ik_stop)
    et_ks(:,:)  = et_loc(:,:)
    et_loc(:,:) = et_tmp(:,ik_start:ik_stop)
  ENDIF
  !
  ! Do not recompute dipole matrix elements
  IF ( epwread .AND. .NOT. epbread ) THEN 
    CONTINUE
  ELSE
    ! compute coarse grid dipole matrix elements.  Very fast 
    IF (.NOT. vme) CALL compute_pmn_para
  ENDIF
  !
  !  gather electronic eigenvalues for subsequent shuffle
  !  
  IF (eig_read) THEN
    et_all(:,:) = zero
    CALL poolgather(nbnd, nkstot, nks, et_loc(1:nbnd,1:nks), et_all)
  ENDIF
  !
  IF (.NOT. kmaps) THEN
    CALL start_clock('kmaps')
    CALL createkmap_pw2
    CALL stop_clock('kmaps')
    CALL print_clock('kmaps')
  ELSE
    ! 
    ! 26/06/2012 RM
    ! if we do not have epmatq already on file then epbread=.false.
    ! .kgmap is used from disk and .kmap is regenerated for each q-point 
    ! 
    WRITE(stdout,'(/5x,a)') 'Using kmap and kgmap from disk'
  ENDIF
  !
  ! Do not do symmetry stuff 
  IF (epwread .AND. .NOT. epbread) THEN
    CONTINUE
  ELSE
    !
    !  allocate dynamical matrix and ep matrix for all q's
    !
    ALLOCATE (dynq(nmodes, nmodes, nq1 * nq2 * nq3))
    ALLOCATE (epmatq(nbnd, nbnd, nks, nmodes, nq1 * nq2 * nq3))
    ALLOCATE (epsi(3, 3))
    ALLOCATE (zstar(3, 3, nat))
    ALLOCATE (bmat(nbnd, nbnd, nks, nq1 * nq2 * nq3))
    ALLOCATE (cu(nbnd, nbndsub, nks))
    ALLOCATE (cuq(nbnd, nbndsub, nks)) 
    ALLOCATE (lwin(nbnd, nks))
    ALLOCATE (lwinq(nbnd, nks))
    ALLOCATE (exband(nbnd))
    !
    dynq(:, :, :)         = czero
    epmatq(:, :, :, :, :) = czero
    epsi(:, :)            = zero
    zstar(:, :, :)        = zero
    bmat(:, :, :, :)      = czero
    cu(:, :, :)           = czero
    cuq(:, :, :)          = czero
    !
    ! read interatomic force constat matrix from q2r
    IF (lifc) CALL read_ifc
    !
    ! SP: The symmetries are now consistent with QE 5. This means that the order of the q in the star
    !     should be the same as in the .dyn files produced by QE 5.
    ! 
    !     First we start by setting up the lattice & crystal symm. as done in PHonon/PH/q2qstar.f90
    ! 
    ! ~~~~~~~~ setup Bravais lattice symmetry ~~~~~~~~
    CALL set_sym_bl( ) ! This should define the s matrix
    WRITE(stdout,'(5x,a,i3)') "Symmetries of Bravais lattice: ", nrot
    !
    ! ~~~~~~~~ setup crystal symmetry ~~~~~~~~ 
    CALL find_sym(nat, tau, ityp, .false., m_loc)
    IF ( .NOT. allfrac ) CALL remove_sym( dfftp%nr1, dfftp%nr2, dfftp%nr3 )
    WRITE(stdout,'(5x,a,i3)') "Symmetries of crystal:         ", nsym
    !   
    ! The following loop is required to propertly set up the symmetry matrix s. 
    ! We here copy the calls made in PHonon/PH/init_representations.f90 to have the same s as in QE 5.
    DO iq_irr=1, nqc_irr
      xq = xqc_irr(:, iq_irr)
      ! search for the small group of q
      CALL set_small_group_of_q(nsymq, invsymq, minus_q)
      ! calculate rtau with the new symmetry order
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      ! calculate the vectors G associated to the symmetry Sq = q + G
      ! if minus_q is true calculate also irotmq and the G associated to Sq=-g+G
      CALL set_giq(xq, s, nsymq, nsym, irotmq, minus_q, gi, gimq)
    ENDDO
  ENDIF ! epwread .and. .NOT. epbread
  ! 
  ! CV: if we read the .fmt files we don't need to read the .epb anymore
  !
  IF (.NOT. epbread .AND. .NOT. epwread) THEN
    ! 
    ALLOCATE (evq(npwx * npol, nbnd))
    IF (lifc) THEN
      ALLOCATE (wscache(-2*nq3:2*nq3, -2*nq2:2*nq2, -2*nq1:2*nq1, nat, nat))
      wscache(:,:,:,:,:) = zero      
    ENDIF
    ! 
    ! In the loop over irr q-point, we need to read the pattern that
    ! corresponds to the dvscf file computed with QE 5.
    !
    nqc = 0
    iq_first = 1
    DO iq_irr = 1, nqc_irr
      u_from_file = .TRUE.
      !tmp_dir_ph = './_ph0/'
      !
      !  read the displacement patterns
      !
      IF (u_from_file) THEN
         ierr = 0
         ! ... look for an empty unit (only ionode needs it)
         IF ( meta_ionode ) CALL iotk_free_unit(iunpun, ierr)
         dirname = TRIM(dvscf_dir) // TRIM(prefix) // '.phsave'
         filename = TRIM(dirname) // '/patterns.' // &
                    TRIM(int_to_char(iq_irr)) // '.xml'
         INQUIRE(FILE=TRIM(filename), EXIST=exst )
         IF ( .NOT. exst) CALL errore('elphon_shuffle_wrap', &
                   'cannot open file for reading or writing', ierr)
         CALL iotk_open_read(iunpun, file = TRIM(filename), &
                                     binary = .FALSE., ierr = ierr)
         CALL read_modes(iunpun, iq_irr, ierr)
         IF (ierr /= 0) CALL errore('elphon_shuffle_wrap', & 
                                    'problem with modes file', 1)
         IF (meta_ionode) CALL iotk_close_read(iunpun)
      ENDIF
      !  
      WRITE(stdout,'(//5x,a)') repeat('=',67) 
      WRITE(stdout,'(5x,"irreducible q point # ",i4)') iq_irr
      WRITE(stdout,'(5x,a/)') repeat('=',67) 
      FLUSH(stdout)
      !
      xq = xqc_irr(:,iq_irr)
      !
      ! SP : The following is largely inspiered by PHonon/PH/q2qstar.f90
      ! 
      ! ~~~~~~~~ setup small group of q symmetry ~~~~~~~~ 
      !
      minus_q = .true.
      sym = .false.
      sym(1:nsym) = .true.
      CALL smallg_q(xq, 0, at, bg, nsym, s, sym, minus_q) ! s is intent(in)
      !
      ! SP: Notice that the function copy_sym reshuffles the s matrix for each irr_q.  
      !     This is why we then need to call gmap_sym for each irr_q [see below]. 
      nsymq = copy_sym(nsym, sym)    
      !
      ! Recompute the inverses as the order of sym.ops. has changed
      CALL inverse_s( )       
      CALL s_axis_to_cart( )
      !
      ! This computes gi, gimq
      CALL set_giq(xq, s, nsymq, nsym, irotmq, minus_q, gi, gimq)
      WRITE(stdout,'(5x,a,i3)') "Symmetries of small group of q:", nsymq
      IF(minus_q) WRITE(stdout,'(10x,a)') "in addition sym. q -> -q+G:"
      ! 
      ! Finally this does some of the above again and also computes rtau...
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      ! ######################### star of q #########################
      ! 
      sym_smallq(:) = 0
      CALL star_q2(xq, at, bg, nsym, s, invs, nq, sxq, isq, imq, .true., sym_smallq)
      !
      ! The reason for xq instead of xq0 in the above is because xq is passed to QE through module  
      xq0 = xq
      !
      !  determine the G vector map S(G) -> G 
      !  SP: The mapping needs to be done for each irr_q because the QE 5 symmetry routine
      !      reshuffles the s matrix for each irr_q [putting the sym of the small group of q first].
      !
      !  [I checked that gmapsym(gmapsym(ig,isym),invs(isym)) = ig]
      CALL gmap_sym(nsym, s, ft, gmapsym, eigv, invs)
      !
      !  Re-set the variables needed for the pattern representation
      !  and the symmetries of the small group of irr-q
      !  (from phq_setup.f90)
      !
      DO isym = 1, nsym
        sym(isym) = .true.
      ENDDO
      !
      CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
      !
      IF ( .NOT. ALLOCATED(sumr) ) ALLOCATE ( sumr(2,3,nat,3) )
      IF (meta_ionode) THEN
        CALL readmat_shuffle2(iq_irr, nqc_irr, nq, iq_first, sxq, imq, isq, &
                              invs, s, irt, rtau)
      ENDIF
      CALL mp_bcast(zstar, meta_ionode_id, world_comm)
      CALL mp_bcast(epsi , meta_ionode_id, world_comm)
      CALL mp_bcast(dynq , meta_ionode_id, world_comm)
      CALL mp_bcast(sumr , meta_ionode_id, world_comm)
      !
      ! now dynq is the cartesian dyn mat (not divided by the masses)
      !
      minus_q = (iswitch > -3)  
      !
      !  loop over the q points of the star
      !
      DO iq = 1, nq
        ! SP: First the vlocq needs to be initialized properly with the first
        !     q in the star
        xq = xq0      
        CALL epw_init(.false.)
        !
        ! retrieve the q in the star
        xq = sxq(:,iq)                               
        !
        ! and populate the uniform grid
        nqc = nqc + 1
        xqc(:,nqc) = xq
        !
        IF (iq == 1) WRITE(stdout,*)
        WRITE(stdout,5) nqc, xq
        !
        !  prepare the gmap for the refolding
        !
        CALL createkmap( xq )                      
        !
        IF (iverbosity == 1) THEN
          !
          !   description of symmetries
          !
          WRITE(stdout, '(36x,"s",24x,"frac. trans.")')
          CALL s_axis_to_cart() ! give sr(:,:, isym)
          DO isym = 1, nsym
            WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
            IF ( ft(1,isym)**2 + ft(2,isym)**2 + ft(3,isym)**2 > 1.0d-8 ) THEN
                ft1 = at(1,1)*ft(1,isym) + at(1,2)*ft(2,isym) + at(1,3)*ft(3,isym)
                ft2 = at(2,1)*ft(1,isym) + at(2,2)*ft(2,isym) + at(2,3)*ft(3,isym)
                ft3 = at(3,1)*ft(1,isym) + at(3,2)*ft(2,isym) + at(3,3)*ft(3,isym)
                WRITE(stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (s(1,ipol,isym),ipol=1,3), ft(1,isym)
                WRITE(stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                            (s(2,ipol,isym),ipol=1,3), ft(2,isym)
                WRITE(stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                            (s(3,ipol,isym),ipol=1,3), ft(3,isym)
                WRITE(stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (sr(1,ipol,isym),ipol=1,3), ft1
                WRITE(stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                            (sr(2,ipol,isym),ipol=1,3), ft2
                WRITE(stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                            (sr(3,ipol,isym),ipol=1,3), ft3
            ELSE
                WRITE(stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                                       isym,  (s (1, ipol, isym) , ipol = 1,3)
                WRITE(stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
                WRITE(stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
                WRITE(stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                                                   isym, (sr(1,ipol,isym), ipol = 1, 3)
                WRITE(stdout, '(17x," (",3f11.7," )")')  (sr(2,ipol,isym), ipol = 1, 3)
                WRITE(stdout, '(17x," (",3f11.7," )"/)') (sr(3,ipol,isym), ipol = 1, 3)
            ENDIF
            ! 
          ENDDO
          !
        ENDIF
        !
        ! isq(isym)=iq means: when we apply symmetry isym to the originating q 
        ! of the star, we get the iq-th member of the star. There are as many 
        ! matches as the degeneracy of the star.
        !
        ! We now need to pick up the q in the small group of q* so that Sxq0+G=iq with G=0.
        ! If we choose another element in the small group
        ! the actual q-point may be Sq+G and we screw up the q-vector below to generate
        ! k+q from k and for the KB projectors
        !
        nsq = 0 ! nsq is the degeneracy of the small group for this iq in the star
        !
        DO jsym = 1, nsym
          IF ( isq(jsym) == iq ) THEN
             nsq = nsq + 1
             sym_sgq(nsq) = jsym
          ENDIF
        ENDDO
        IF ( nsq*nq .ne. nsym ) CALL errore('elphon_shuffle_wrap', 'wrong degeneracy', iq)
        ! 
        IF (iverbosity == 1) THEN
          !
          WRITE(stdout,*) 'iq, i, isym, nog, symmo'
          DO i = 1, nsq
            !
            isym = sym_sgq(i)
            ism1 = invs (isym)
            !
            !  check for G such that Sq = q* + G 
            ! 
            aq  = xq0
            saq = xq
            CALL cryst_to_cart(1, aq, at, -1)
            DO j = 1, 3
              raq(j) = s(j,1,ism1) * aq(1) &
                     + s(j,2,ism1) * aq(2) &
                     + s(j,3,ism1) * aq(3)
            ENDDO
            CALL cryst_to_cart(1, saq, at, -1)
            nog = eqvect_strict(raq, saq) 
            !
            !  check whether the symmetry belongs to a symmorphic group
            !
            symmo = ( ft(1,isym)**2 + ft(2,isym)**2 + ft(3,isym)**2 > 1.0d-8 )
            !
            WRITE(stdout,'(3i5,L3,L3)') iq, i, isym, nog, symmo
            !
          ENDDO  
          !
        ENDIF
        ! 
        ! SP: We now need to select one symmetry among the small group of q that has G=0 
        !     (i.e. that respects Sq0+G=q ). There should always be such symmetry. 
        !     We enforce this for later easiness. 
        ! 
        aq = xq0
        saq = xq
        CALL cryst_to_cart(1, aq, at, - 1)
        CALL cryst_to_cart(1, saq, at, -1)
        !write(*,*)'xq0 ',aq
        !write(*,*)'xq ',saq
        DO jsym = 1, nsq
          ism1 = invs(sym_sgq(jsym))
          raq = zero
          DO ipol = 1, 3
            DO jpol = 1, 3
              raq(ipol) = raq(ipol) + s(ipol,jpol,ism1) * aq(jpol)
            ENDDO
          ENDDO
          nog = eqvect_strict(raq,saq)
          IF (nog) THEN ! This is the symmetry such that Sq=q
            isym = sym_sgq(jsym)
            EXIT
          ENDIF
          ! If we enter into that loop it means that we have not found 
          ! such symmetry within the small group of q. 
          IF (jsym == nsq) THEN
            CALL errore( 'elphon_shuffle_wrap ', 'No sym. such that Sxq0=iq was found in the sgq !', 1 )
          ENDIF
        ENDDO
        !
        !
        CALL loadumat( nbnd, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq, exband, w_centers )
        !
        ! Calculate overlap U_k+q U_k^\dagger
        IF (lpolar) CALL compute_umn_c( nbnd, nbndsub, nks, cu, cuq, bmat(:,:,:,nqc) )
        !
        !   calculate the sandwiches
        !
        ! A more accurate way of doing this is to symmetrize the matrix element w.r.t.
        ! the small group of the given q in the star. I'm not doing this here.
        ! (but I checked that even without symm the result of full zone and irr zone
        ! are equal to 5+ digits).
        ! For any volunteers, please write to giustino@civet.berkeley.edu
        !
        CALL elphon_shuffle( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, xq0, .false. )
        !
        !  bring epmatq in the mode representation of iq_first, 
        !  and then in the cartesian representation of iq
        !
        CALL rotate_eigenm( iq_first, nqc, isym, s, invs, irt, rtau, xq, cz1, cz2 )
        !
        CALL rotate_epmat( cz1, cz2, xq, nqc, lwin, lwinq, exband )
  !DBSP
        !write(*,*)'epmatq(:,:,2,:,nqc)',SUM(epmatq(:,:,2,:,nqc))
        !write(*,*)'epmatq(:,:,2,:,nqc)**2',SUM((REAL(REAL(epmatq(:,:,2,:,nqc))))**2)+&
        !  SUM((REAL(AIMAG(epmatq(:,:,2,:,nqc))))**2)
        !print*,'dynq ', SUM(dynq(:,:,nqc))
        !print*,'et ', et_loc(:,2)
  !END
        ! SP: Now we treat separately the case imq == 0
        IF (imq == 0) THEN
          !
          ! SP: First the vlocq need to be initialized propertly with the first
          !     q in the star
          xq = -xq0
          CALL epw_init(.false.)
          !
          ! retrieve the q in the star
          xq = -sxq(:,iq)
          !
          ! and populate the uniform grid
          nqc = nqc + 1
          xqc(:,nqc) = xq
          !
          IF (iq == 1) write(stdout,*)
          WRITE(stdout,5) nqc, xq
          !
          !  prepare the gmap for the refolding
          !
          CALL createkmap( xq )
          !
          CALL loadumat( nbnd, nbndsub, nks, nkstot, xq, cu, cuq, lwin, lwinq, exband, w_centers )
          !
          ! Calculate overlap U_k+q U_k^\dagger
          IF (lpolar) CALL compute_umn_c( nbnd, nbndsub, nks, cu, cuq, bmat(:,:,:,nqc) )
          !
          xq0 = -xq0
          !
          CALL elphon_shuffle( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, xq0, .true. )
          !  bring epmatq in the mode representation of iq_first, 
          !  and then in the cartesian representation of iq
          !
          CALL rotate_eigenm( iq_first, nqc, isym, s, invs, irt, rtau, xq, cz1, cz2 )
          !
          CALL rotate_epmat( cz1, cz2, xq, nqc, lwin, lwinq, exband )
          !
          xq0 = -xq0
        ENDIF ! end imq == 0  
        !
      ENDDO
      !
      iq_first = iq_first + nq
      if (imq == 0) iq_first = iq_first + nq
      !
    ENDDO ! irr-q loop
    ! 
    IF (nqc.ne.nq1*nq2*nq3) &
       CALL errore('elphon_shuffle_wrap','nqc .ne. nq1*nq2*nq3',nqc)
    wqlist = dble(1) / dble(nqc)
    !
    IF (lifc) DEALLOCATE (wscache)
    DEALLOCATE (evc)
    DEALLOCATE (evq)
    DEALLOCATE (vlocq)
    DEALLOCATE (dmuxc)
    DEALLOCATE (eigqts)
    DEALLOCATE (rtau)
    DEALLOCATE (u)
    DEALLOCATE (npert)
    IF (okvan) THEN
      DEALLOCATE (int1)
      DEALLOCATE (int2)
      DEALLOCATE (int4)
      DEALLOCATE (int5)
      IF (noncolin) THEN 
        DEALLOCATE (int1_nc)
        DEALLOCATE (int4_nc)
        IF (lspinorb) THEN
          DEALLOCATE (int2_so)
          DEALLOCATE (int5_so)
        ENDIF
      ENDIF
    ENDIF
    DO ik = 1, nks
      DO ipol = 1, 3
        CALL deallocate_bec_type( alphap(ipol,ik) )
      ENDDO
    ENDDO
    DEALLOCATE (alphap)
    DO ik = 1, size(becp1)
      CALL deallocate_bec_type( becp1(ik) )
    ENDDO
    DEALLOCATE (becp1)
    CALL deallocate_bec_type ( becp )
  ENDIF ! IF (.NOT. epbread .AND. .NOT. epwread) THEN
  !
  IF (my_image_id == 0 ) THEN
    IF ( epbread .OR. epbwrite ) THEN
      !
      ! read/write the e-ph matrix elements and other info in the Bloch representation
      ! (coarse mesh) from/to .epb files (one for each pool)
      !
      tempfile = trim(tmp_dir) // trim(prefix) // '.epb' 
      CALL set_ndnmbr(0, my_pool_id+1, 1, npool, filelab)
      tempfile = trim(tmp_dir) // trim(prefix) // '.epb' // filelab
      !
      IF (epbread) THEN
         inquire(file = tempfile, exist=exst)
         IF ( .NOT.  exst) CALL errore( 'elphon_shuffle_wrap', 'epb files not found ', 1)
         OPEN(iuepb, file = tempfile, form = 'unformatted')
         WRITE(stdout,'(/5x,"Reading epmatq from .epb files"/)') 
         READ(iuepb) nqc, xqc, et_loc, dynq, epmatq, zstar, epsi
         CLOSE(iuepb)
         WRITE(stdout,'(/5x,"The .epb files have been correctly read"/)')
      ENDIF
      !
      IF (epbwrite) THEN
         OPEN(iuepb, file = tempfile, form = 'unformatted')
         WRITE(stdout,'(/5x,"Writing epmatq on .epb files"/)') 
         WRITE(iuepb) nqc, xqc, et_loc, dynq, epmatq, zstar, epsi
         CLOSE(iuepb)
         WRITE(stdout,'(/5x,"The .epb files have been correctly written"/)')
      ENDIF
    ENDIF
  ENDIF
  !
  ! In case of image parallelization we want to stop after writing the .epb file
  IF (nimage > 1) THEN
    WRITE(stdout,'(/5x,"Image parallelization. The code will stop now. "/)')
    WRITE(stdout,'(/5x,"You need to restart a calculation by reading the .epb "/)')
    WRITE(stdout,'(/5x,"                       with pool parallelization only. "/)')
    CALL stop_epw
  ENDIF
  !
  IF (  .NOT. epbread .AND. epwread ) THEN
  !  CV: need dummy nqc, xqc for the ephwann_shuffle call
     nqc = 1
     xqc = zero
     WRITE(stdout,'(/5x,"Do not need to read .epb files; read .fmt files"/)')
  !
  ENDIF
  !
  !CALL mp_barrier(inter_pool_comm)
  CALL mp_barrier(world_comm)
  !
  !   now dynq is the cartesian dyn mat ( not divided by the masses)
  !   and epmatq is the epmat in cartesian representation (rotation in elphon_shuffle)
  !
  ! free up some memory
  !
  DEALLOCATE (umat_all)
  DEALLOCATE (umat)
  DEALLOCATE (xqc_irr)
  DEALLOCATE (wqlist)
  ! 
  IF (maxvalue > nqxq) THEN
    DEALLOCATE (qrad)
  ENDIF
  !
  IF ( ASSOCIATED (igkq) )      NULLIFY    (igkq)
  IF ( ALLOCATED  (dvpsi))      DEALLOCATE (dvpsi)
  IF ( ALLOCATED  (dpsi) )      DEALLOCATE (dpsi)
  IF ( ALLOCATED  (sumr) )      DEALLOCATE (sumr)
  IF ( ALLOCATED  (cu) )        DEALLOCATE (cu)
  IF ( ALLOCATED  (cuq) )       DEALLOCATE (cuq)
  IF ( ALLOCATED  (lwin) )      DEALLOCATE (lwin)
  IF ( ALLOCATED  (lwinq) )     DEALLOCATE (lwinq)
  IF ( ALLOCATED  (bmat) )      DEALLOCATE (bmat)
  IF ( ALLOCATED  (igk_k_all) ) DEALLOCATE (igk_k_all)
  IF ( ALLOCATED  (ngk_all) )   DEALLOCATE (ngk_all)
  IF ( ALLOCATED  (exband) )    DEALLOCATE (exband)
  ! 
  CALL stop_clock ( 'elphon_wrap' )
!DBSP
!  DO iq = 1, nqc
!    write(*,*) iq, xqc(:,iq)
!    write(*,*)'epmatq(:,:,2,:,iq)',SUM(epmatq(:,:,2,:,iq))
!    write(*,*)'epmatq(:,:,2,:,iq)**2',SUM((REAL(REAL(epmatq(:,:,2,:,iq))))**2)+&
!               SUM((REAL(AIMAG(epmatq(:,:,2,:,iq))))**2)
!  ENDDO
!END
  !
  ! the electron-phonon wannier interpolation
  !
  IF(etf_mem == 0 .OR. etf_mem == 1 ) CALL ephwann_shuffle( nqc, xqc )
  IF(etf_mem == 2 ) THEN
#if defined(__MPI)         
    CALL ephwann_shuffle_mem( nqc, xqc )
#else
    WRITE(stdout,'(/5x,a)') 'WARNING: etf_mem==2 only works with MPI'
    WRITE(stdout,'(5x,a)')  '         Changing to etf_mem == 1 and continue ...'
    etf_mem = 1
    CALL ephwann_shuffle( nqc, xqc )
#endif
  ENDIF        
  DEALLOCATE (xqc)
  !
5 FORMAT (8x,"q(",i5," ) = (",3f12.7," )") 
  !
  RETURN
  END SUBROUTINE elphon_shuffle_wrap
  !
  !---------------------------------------------------------------------------
  SUBROUTINE irotate( x, s, sx)
  !---------------------------------------------------------------------------
  !!
  !! a simple symmetry operation in crystal coordinates ( s is INTEGER!)
  !!
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(kind=DP), INTENT(in) :: x(3)
  !! Input x
  INTEGER, INTENT(in) :: s(3,3)
  !! Symmetry matrix
  REAL(kind=DP), INTENT(out) :: sx(3)
  !! Output rotated x 
  !
  ! Local Variable
  INTEGER :: i
  !
  DO i = 1, 3
     sx(i) = dble( s(i,1) ) * x(1) &
           + dble( s(i,2) ) * x(2) &
           + dble( s(i,3) ) * x(3)
  ENDDO
  !
  RETURN
  ! 
  END SUBROUTINE irotate
  !---------------------------------------------------------------------------
  LOGICAL function eqvect_strict(x, y)
  !-----------------------------------------------------------------------
  !!
  !! This function test if two tridimensional vectors are equal
  !!
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(kind=DP), INTENT(in) :: x(3)
  !! input: input vector
  REAL(kind=DP), INTENT(in) :: y(3)
  !! input: second input vector
  REAL(kind=DP) :: accep
  !! acceptance parameter
  PARAMETER (accep = 1.0d-5)
  !
  eqvect_strict = abs( x(1)-y(1) ) < accep .AND. &
                  abs( x(2)-y(2) ) < accep .AND. &
                  abs( x(3)-y(3) ) < accep
  !
  END FUNCTION eqvect_strict
  !---------------------------------------------------------------------------
  SUBROUTINE read_modes(iunpun, current_iq, ierr)
  !---------------------------------------------------------------------------
  !!
  !! This routine reads the displacement patterns.
  !!
  USE modes,        ONLY : nirr, npert, u
  USE lr_symm_base, ONLY : minus_q, nsymq  
  USE iotk_module,  ONLY : iotk_index, iotk_scan_dat, iotk_scan_begin, &
                           iotk_scan_end
  USE io_global,    ONLY : meta_ionode, meta_ionode_id
  USE mp,           ONLY : mp_bcast
  USE mp_global,    ONLY : world_comm
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: current_iq
  !! Current q-point 
  INTEGER, INTENT(in) :: iunpun
  !! Current q-point 
  INTEGER, INTENT(out) :: ierr
  !! Error
  !
  ! Local variables
  INTEGER :: imode0, imode
  !! Counter on modes
  INTEGER :: irr
  !! Counter on irreducible representations
  INTEGER :: ipert
  !! Counter on perturbations at each irr
  INTEGER :: iq
  !! Current q-point 
  !
  ierr = 0
  IF (meta_ionode) THEN
     CALL iotk_scan_begin(iunpun, "IRREPS_INFO")
     !
     CALL iotk_scan_dat(iunpun, "QPOINT_NUMBER", iq)
  ENDIF
  CALL mp_bcast(iq,  meta_ionode_id, world_comm)
  IF (iq /= current_iq) CALL errore('read_modes', &
            'problems with current_iq', 1 )
  ! 
  IF (meta_ionode) THEN
     !
     CALL iotk_scan_dat(iunpun, "QPOINT_GROUP_RANK", nsymq)
     CALL iotk_scan_dat(iunpun, "MINUS_Q_SYM", minus_q)
     CALL iotk_scan_dat(iunpun, "NUMBER_IRR_REP", nirr)
     imode0 = 0
     DO irr = 1, nirr
        CALL iotk_scan_begin(iunpun, "REPRESENTION"// &
                                   TRIM( iotk_index(irr) ))
        CALL iotk_scan_dat(iunpun, "NUMBER_OF_PERTURBATIONS", npert(irr))
        DO ipert = 1, npert(irr)
           imode = imode0 + ipert
           CALL iotk_scan_begin(iunpun, "PERTURBATION"// &
                                   TRIM( iotk_index(ipert) ))
           CALL iotk_scan_dat(iunpun, "DISPLACEMENT_PATTERN", u(:,imode))
           CALL iotk_scan_end(iunpun, "PERTURBATION"// &
                                  TRIM( iotk_index(ipert) ))
        ENDDO
        imode0 = imode0 + npert(irr)
        CALL iotk_scan_end(iunpun, "REPRESENTION"// &
                                   TRIM( iotk_index(irr) ))
     ENDDO
     !
     CALL iotk_scan_end(iunpun, "IRREPS_INFO")
     !
  ENDIF
  !
  CALL mp_bcast(nirr   , meta_ionode_id, world_comm)
  CALL mp_bcast(npert  , meta_ionode_id, world_comm)
  CALL mp_bcast(nsymq  , meta_ionode_id, world_comm)
  CALL mp_bcast(minus_q, meta_ionode_id, world_comm)
  CALL mp_bcast(u      , meta_ionode_id, world_comm)
  !
  RETURN
  !
  END SUBROUTINE read_modes
