!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE elph_readin()
  !----------------------------------------------------------------------------
  !
  !  This routine reads the input parameters for ELPH and reads
  !  the data produced by PWscf.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : meta_ionode, meta_ionode_id, qestdin, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  ! USE mp_images,        ONLY : my_image_id
  ! USE check_stop,       ONLY : max_seconds
  USE read_namelists_module, ONLY : check_namelist_read
  USE io_files,         ONLY : tmp_dir, prefix, create_directory, check_tempdir
  USE open_close_input_file, ONLY : open_input_file, close_input_file
  ! USE paw_variables,    ONLY : okpaw
  USE control_flags,    ONLY : io_level
  USE ions_base,        ONLY : nat
  USE control_lr,       ONLY : reduce_io
  USE control_ph,       ONLY : tmp_dir_ph, tmp_dir_phq
  ! ethr_nscf, lrpa, reduce_io, niter_ph, maxter, alpha_mix, &
  !                              nmix_ph, tr2_ph, thresh_init, conv_thr_nscf
  ! USE qpoint,           ONLY : nq1, nq2, nq3, start_q, last_q
  USE elphcom,          ONLY : tmp_dir_save, wannier_seedname, qplot, folder_elph, &
                               folder_wan_Rr, folder_phonon, fildvscf, write_wan_Rr
  ! USE crpa_pert,        ONLY : pert_basis
  ! USE qpoint,           ONLY : x_q, nqs
  ! USE constrained_dfpt, ONLY : lcdfpt
  USE modes,            ONLY : nmodes
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst, parallelfs
  !! TODO
  CHARACTER(LEN=256) :: outdir
  ! CHARACTER(LEN=6) :: int_to_char
  ! CHARACTER(len=80) :: diagonalization
  !! Diagonalization method for NSCF calculation
  CHARACTER(LEN = 512) :: line
  !! Line in input file
  INTEGER :: ios
  !! integer variable for I/O control
  INTEGER :: ios2
  !! INTEGER variable for I/O control
  ! INTEGER :: iter
  ! !! counter on iterations
  ! !
  ! REAL(DP), ALLOCATABLE :: xqaux(:,:)
  ! INTEGER, ALLOCATABLE :: wqaux(:)
  ! INTEGER :: nqaux, iq, ipol
  ! !
  ! INTEGER :: nmix, niter_max
  ! REAL(DP) :: tr2
  !
  ! folder_elph : Output directory where the elph data will be stored
  ! folder_phonon : Folder where phonon data is stored.
  !                 (Default: outdir/_ph0)
  ! fildvscf : Filename for the phonon potential. Should be the same as in ph.x calculation.
  !            (Default: 'dvscf')
  !
  ! folder_wan_Rr : Folder where the collected real-space Wannier functions are saved.
  !                 If write_wan_Rr is .true. this folder is created and the Wannier functions
  !                 are written to it. If write_wan_Rr is .false. this folder should already
  !                 exist and contain the files prefix.wan_Rr.xxx and prefix.wan_R.dat.
  !                 (Default: tmp_dir_ph)
  ! write_wan_Rr : If .true. write the collected real-space Wannier functions to file.
  !                If .false. use the existing files in folder_wan_Rr. In this case,
  !                prefix.wan_Rr.xxx and prefix.wan_R.dat files must exist inside folder_wan_Rr.
  !                (Default: .true.)
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  NAMELIST / INPUTELPH / prefix, outdir, wannier_seedname, reduce_io, qplot, folder_elph, &
                         folder_wan_Rr, folder_phonon, fildvscf, write_wan_Rr
  ! , lrpa, niter_max, alpha_mix, nq1, nq2, nq3, &
                        !  nmix, start_q, last_q, thresh_init, tr2, active_space, &
                        !  active_bands_min, active_bands_max, diagonalization, iverbosity, &
                        !  conv_thr_nscf, reduce_io, pert_basis, filU, &
                        !  dist_thr, dist_thr_large, lcdfpt
! , skip_equivalence_q,   &
!                          conv_thr_chi, skip_atom, skip_type, equiv_type, iverbosity,  &
!                          background, find_atpert, max_seconds, rmax,     &
!                          , , compute_crpa, perturb_only_atom,   &
!                          , sum_pertq, num_neigh, lmin,      &
!                          determine_num_pert_only, disable_type_analysis, docc_thr,    &
!                          determine_q_mesh_only
  !
  ! Note: meta_ionode is a single processor that reads the input
  !       Data read from input is subsequently broadcast to all processors
  !       from meta_ionode_id (using the default communicator world_comm)
  !
  IF (meta_ionode) CALL input_from_file()
  !
  ! Set default values for variables in namelist
  !
  prefix             = 'pwscf'
!   conv_thr_chi       = 1.D-5
  ! thresh_init        = 1.D-14
  ! tr2                = 1.D-20
  ! diagonalization    = 'david'
  ! conv_thr_nscf      = 1.D-9
  wannier_seedname   = ''
  folder_wan_Rr      = ''
  folder_phonon      = ''
  folder_elph        = './'
  fildvscf           = 'dvscf'
!   docc_thr           = 5.D-5
!   rmax               = 100.D0
!   skip_atom(:)       = .FALSE.
!   skip_type(:)       = .FALSE.
!   perturb_only_atom(:)    = .FALSE.
!   skip_equivalence_q      = .FALSE.
!   determine_num_pert_only = .FALSE.
!   determine_q_mesh_only   = .FALSE.
!   disable_type_analysis   = .FALSE.
!   equiv_type(:)      = 0
!   find_atpert        = 1
!   background         = 'no'
!   compute_crpa         = .FALSE.
!   sum_pertq          = .FALSE.
!   num_neigh          = 6
!   lmin               = 2
!   nq1                = 1
!   nq2                = 1
!   nq3                = 1
!   start_q            = 1
!   last_q             = -1
!   iverbosity         = 1
!   niter_max          = 100
!   alpha_mix(:)       = 0.D0
!   alpha_mix(1)       = 0.7D0
!   nmix               = 12
!   lcdfpt             = .FALSE.
!   active_space       = ''
!   pert_basis         = 'bands'
!   active_bands_min   = 0
!   active_bands_max   = 0
  write_wan_Rr       = .TRUE.
  qplot              = .FALSE.
  reduce_io          = .TRUE.
!   filU               = 'crpa'
! !   max_seconds        = 1.E+7_DP
!   lrpa               = .FALSE.   ! Needed in dv_of_drho
!   dist_thr           = 6.D-4  ! same as eps_dist in PW/src/ldaU.f90
!   dist_thr_large     = 6.D-4  ! same as eps_dist in PW/src/ldaU.f90
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  ! Reading the namelist INPUTELPH and check
  !
  IF (meta_ionode) THEN
    !
    ! ... Input from file (ios=0) or standard input (ios=-1) on unit "qestdin"
    !
    ios = open_input_file (  )
    !
    READ(qestdin, INPUTELPH, IOSTAT = ios)
    !
  ENDIF ! meta_ionode
  !
  CALL check_namelist_read(ios, qestdin, "inputelph")
  CALL mp_bcast(ios, meta_ionode_id, world_comm)
  IF (ABS(ios) /= 0) CALL errore('elph_readin', 'reading inputelph namelist', ABS(ios))
  !
  ! Setup tmp_dir
  !
  tmp_dir = trimcheck (outdir)
  !
  ! Broadcast input parameters over all processors
  !
  CALL elph_bcast_input()
  !
  ! Read and set the q points
  !
  CALL elph_read_qpoints(qplot)
  !
  IF (meta_ionode) ios = close_input_file ()
  !
  ! IO option
  !
  IF (reduce_io) io_level = 0
  !
  ! Here we finished the reading of the input file.
  ! Now allocate space for pwscf variables, read and check them.
  !
  tmp_dir_save = tmp_dir
  !
  tmp_dir_ph = trimcheck(TRIM(tmp_dir) // '_elph0')
  CALL check_tempdir(tmp_dir_ph, exst, parallelfs)
  tmp_dir_phq = tmp_dir_ph
  !
  ! Initialize folders
  !
  IF (TRIM(folder_wan_Rr) == '') THEN
    ! If folder_wan_Rr is not set, use the default tmp_dir_ph
    folder_wan_Rr = tmp_dir_ph
  ELSE
    ! Otherwise, create the directory specified in folder_wan_Rr
    CALL create_directory(folder_wan_Rr)
  ENDIF
  folder_wan_Rr = trimcheck(folder_wan_Rr)
  !
  IF (TRIM(folder_phonon) == '') THEN
    folder_phonon = TRIM(tmp_dir) // '_ph0/'
  ENDIF
  folder_phonon = trimcheck(folder_phonon)
  !
  folder_elph = trimcheck(folder_elph)
  CALL create_directory(folder_elph)
  !
  ! Read various data produced by PWscf.
  ! In particular, read the unperturbed occupation matrices
  ! via calling the routine read_rho.
  ! read_file calls init_tab_atwfc which initializes tab_atwfc.
  !
  CALL read_file()
  !
  nmodes = 3 * nat
  !
  ! Make sure all the features used in the PWscf calculation
  ! are actually supported by CRPA.
  !
  CALL input_sanity()
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE input_sanity()
  !--------------------------------------------------------------------------
  !
  ! This subroutine aims to gather all of the input sanity checks
  ! (features enabled in PWscf which are unsupported in CRPA).
  !
!   USE klist,            ONLY : lgauss, ltetra, two_fermi_energies
!   USE control_flags,    ONLY : gamma_only, tqr
!   USE fixed_occ,        ONLY : tfixed_occ
!   USE cellmd,           ONLY : lmovecell
!   USE noncollin_module, ONLY : i_cons, noncolin
!   USE mp_bands,         ONLY : nbgrp
!   USE xc_lib,           ONLY : xclib_dft_is
! !   USE ldaU,             ONLY : lda_plus_u, Hubbard_projectors, lda_plus_u_kind, Hubbard_J0, &
! !                                Hubbard_V, is_hubbard_back
!   !
!   IMPLICIT NONE
!   !
!   IF (tr2_ph <= 0.D0) CALL errore ('elph_readin', ' Wrong tr2_ph', 1)
!   !
!   IF (thresh_init <= 0.D0) CALL errore ('elph_readin', ' Wrong thresh_init', 1)
!   !
!   IF (nq1 < 0 .OR. nq2 < 0 .OR. nq3 < 0) &
!        CALL errore('elph_readin','nq1, nq2, and nq3 must be greater than 0', 1)
!   !
!   IF (start_q <= 0 ) CALL errore('elph_readin', ' Wrong start_q ', 1)
!   !
!   IF (filU == '') CALL errore('elph_readin', ' Wrong filU', 1)
! !   !
! !   IF (compute_crpa .AND. ANY(perturb_only_atom(:))) &
! !      CALL errore ('elph_readin', 'compute_crpa and perturb_only_atom are not allowed to be true together', 1)
! !   !
! !   IF ( ANY(Hubbard_V(:,:,2).NE.0.d0) .OR. &
! !        ANY(Hubbard_V(:,:,3).NE.0.d0) .OR. &
! !        ANY(Hubbard_V(:,:,4).NE.0.d0) ) &
! !      CALL errore ('elph_readin', 'The CRPA code does not support DFT+U+V with the background', 1)
! !   !
! !   IF (ANY(is_hubbard_back(:))) CALL errore ("elph_readin", &
! !           &" Two (or more) Hubbard channels per atomic type is not implemented", 1)
! !   !
! !   IF (ANY(Hubbard_J0(:).NE.0.d0)) &
! !      CALL errore ('elph_readin', 'Hubbard_J0 /= 0 is not allowed.', 1)
! !   !
! !   IF (determine_q_mesh_only .AND. .NOT.ANY(perturb_only_atom(:))) &
! !      CALL errore ('elph_readin', 'determine_q_mesh_only can be set to .true. only if perturb_only_atom is .true. for some atom', 1)
! !   !
! !   IF (sum_pertq .AND. .NOT.ANY(perturb_only_atom(:))) &
! !      CALL errore ('elph_readin', 'sum_pertq can be set to .true. only if perturb_only_atom is .true. for some atom', 1)
! !   !
!   ! IF (niter_ph < 1 .OR. niter_ph > maxter) &
!   !   CALL errore ('elph_readin', ' Wrong niter_ph ', 1)
!   ! !
!   ! DO iter = 1, niter_ph
!   !    IF ( alpha_mix(iter) < 0.D0 .OR. alpha_mix(iter) > 1.D0 ) &
!   !       CALL errore ('elph_readin', ' Wrong alpha_mix ', iter)
!   ! ENDDO
! !   !
! !   IF (num_neigh.LT.1) CALL errore('elph_readin','Not allowed value of num_neigh',1)
! !   !
! !   IF (lmin.LT.0 .OR. lmin.GT.3) CALL errore('elph_readin','Not allowed value of lmin',1)
! !   !
! !   IF (nmix.LT.1) CALL errore ('elph_readin', ' Wrong nmix ', 1)
! !   !
! !   IF (ltetra) CALL errore ('elph_readin', 'CRPA with tetrahedra is not supported', 1)
! !   !
! !   IF (gamma_only) CALL errore('elph_readin',&
! !      & 'Cannot start from pw.x data file using Gamma-point tricks',1)
! !   !
! !   IF (.NOT.lda_plus_u) CALL errore('elph_readin',&
! !      & 'The CRPA code can be used only on top of DFT+Hubbard (i.e. when the HUBBARD card is used in pw.x)',1)
! !   !
! !   IF (lda_plus_u_kind.EQ.1) CALL errore("elph_readin", &
! !      & ' The CRPA code does not support the Liechtenstein formulation of DFT+U',1)
! !   !
! !   IF (Hubbard_projectors.NE."atomic" .AND. Hubbard_projectors.NE."ortho-atomic") &
! !      CALL errore("elph_readin", &
! !      " The CRPA code for this Hubbard_projectors type is not implemented",1)
! !   !
! !   IF (lmovecell) CALL errore('elph_readin','The CRPA code is not working after vc-relax',1)
! !   !
! !   IF (nbgrp > 1) CALL errore('elph_readin', &
! !      & 'band parallelization is not implemented in CRPA',1)
! !   !
! !   IF (i_cons /= 0) CALL errore('elph_readin',&
! !      & 'The CRPA code with constrained magnetization is not yet available',1)
! !   !
! !   IF (two_fermi_energies .AND. (ltetra .OR. lgauss)) CALL errore('elph_readin', &
! !      & 'The CRPA code with two Fermi energies is not available for metals',1)
! !   !
! !   IF (tqr) CALL errore('elph_readin',&
! !      & 'The CRPA code with Q in real space is not supported',1)
! !   !
! !   IF (tfixed_occ) CALL errore('elph_readin', &
! !      & 'The CRPA code with arbitrary occupations not tested',1)
! !   !
! !   IF ( xclib_dft_is('meta') ) CALL errore('elph_readin',&
! !      'The CRPA code with meta-GGA functionals is not yet available',1)
! !   !
! !   IF ( xclib_dft_is('hybrid') ) CALL errore('elph_readin',&
! !      'The CRPA code with hybrid functionals is not yet available',1)
!   !
!   ! Check sanity of active space definition
!   !
!   ! IF (pert_basis == '') CALL errore('elph_readin', 'pert_basis not defined', 1)
!   ! !
!   ! IF (lcdfpt) THEN
!   !   IF (active_space == '') CALL errore('elph_readin', 'active_space not defined', 1)
!   !   IF (active_space == 'bands') THEN
!   !     IF (active_bands_min < 1) CALL errore('elph_readin', 'active_bands_min is not set', 1)
!   !     IF (active_bands_max < 1) CALL errore('elph_readin', 'active_bands_max is not set', 1)
!   !     IF (active_bands_min > active_bands_max) CALL errore('elph_readin', &
!   !         'active_bands_min cannot be greater than active_bands_max', 1)
!   !   ELSEIF (active_space == 'wannier') THEN
!   !     IF (wannier_seedname == '') CALL errore('elph_readin', 'wannier_seedname not set', 1)
!   !   ELSE
!   !     CALL errore('elph_readin', 'Invalid active_space. Must be none or bands or wannier.', 1)
!   !   ENDIF
!   ! ENDIF
!   ! !
!   ! IF (pert_basis == 'bands') THEN
!   !   IF (active_bands_min < 1) CALL errore('elph_readin', 'active_bands_min is not set', 1)
!   !   IF (active_bands_max < 1) CALL errore('elph_readin', 'active_bands_max is not set', 1)
!   !   IF (active_bands_min > active_bands_max) CALL errore('elph_readin', &
!   !       'active_bands_min cannot be greater than active_bands_max', 1)
!   ! ELSEIF (pert_basis == 'wannier') THEN
!   !   IF (wannier_seedname == '') CALL errore('elph_readin', 'wannier_seedname not set', 1)
!   ! ELSE
!   !   CALL errore('elph_readin', 'Invalid pert_basis. Must be bands or wannier.', 1)
!   ! ENDIF
  !
  RETURN
  !
END SUBROUTINE input_sanity
  !
  !----------------------------------------------------------------------------
  SUBROUTINE elph_read_qpoints(qplot)
    !----------------------------------------------------------------------------
    !! Read q points from the input file. Use the same format as in ph.x.
    !! Assumes that the input file is already open on unit qestdin.
    !! The results are stored in module qpoint.
    !!
    !! Input:
    !! - If qplot = .TRUE. : Read a list of q points from the input file.
    !!
    !! Output (all in module qpoint):
    !! - nqs: Number of q points
    !! - x_q(3, nqs): Coordinates of q points in Cartesian coordinates
    !! - wq(nqs): Weights of q points
    !!
    !! TODO: Support options q2d, q_in_band_form, lqdir
    !----------------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE mp,               ONLY : mp_bcast
    USE mp_world,         ONLY : world_comm
    USE io_global,        ONLY : meta_ionode, meta_ionode_id, qestdin
    USE qpoint,           ONLY : x_q, nqs, wq
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: qplot
    !! If .TRUE., read the list of q points from the input file
    !
    INTEGER :: nqaux
    !! Number of q points read from the input file
    INTEGER :: ipol
    !! Index for directions
    INTEGER :: iq
    !! Index for q points
    INTEGER :: ios
    !! Integer variable for I/O control
    REAL(DP), ALLOCATABLE :: xqaux(:, :)
    !! Temporary array for q points read from the input file
    REAL(DP), ALLOCATABLE :: wqaux(:)
    !! Temporary array for weights of q points read from the input file
    !
    IF (.NOT. qplot) CALL errore('elph_read_qpoints', &
        'currently only qplot = .TRUE. is implemented', 1)
    !
    IF (qplot) THEN
      !
      ! Read a list of q points from the input file
      !
      IF (meta_ionode) THEN
        ios = 0
        READ(qestdin, *, IOSTAT=ios) nqaux
      ENDIF
      CALL mp_bcast(ios, meta_ionode_id, world_comm )
      CALL errore('elph_read_qpoints', 'reading nq', ABS (ios) )
      !
      CALL mp_bcast(nqaux, meta_ionode_id, world_comm )
      !
      ALLOCATE(xqaux(3, nqaux))
      ALLOCATE(wqaux(nqaux))
      !
      IF (meta_ionode) THEN
        DO iq = 1, nqaux
          READ (qestdin, *, IOSTAT=ios) (xqaux (ipol,iq), ipol=1,3), wqaux(iq)
        ENDDO
      ENDIF
      CALL mp_bcast(ios, meta_ionode_id, world_comm )
      CALL errore('elph_read_qpoints', 'reading xq', ABS (ios))
      !
      CALL mp_bcast(xqaux, meta_ionode_id, world_comm )
      CALL mp_bcast(wqaux, meta_ionode_id, world_comm )
      !
      ! IF (q2d) THEN
      !    nqs=wqaux(2)*wqaux(3)
      !    ALLOCATE(x_q(3,nqs))
      !    ALLOCATE(wq(nqs))
      !    CALL generate_k_in_plane(nqaux, xqaux, wqaux, x_q, wq, nqs)
      ! ELSEIF (q_in_band_form) THEN
      !    nqs=SUM(wqaux(1:nqaux-1))+1
      !    DO i=1,nqaux-1
      !       IF (wqaux(i)==0) nqs=nqs+1
      !    ENDDO
      !    ALLOCATE(x_q(3,nqs))
      !    ALLOCATE(wq(nqs))
      !    CALL generate_k_along_lines(nqaux, xqaux, wqaux, x_q, wq, nqs)
      ! ELSE
      !
      ! Use xqaux and wqaux directly as the q points and weights
      !
      nqs = nqaux
      ALLOCATE(x_q(3, nqs))
      ALLOCATE(wq(nqs))
      x_q(:, 1:nqs) = xqaux(:, 1:nqs)
      wq(:) = wqaux(:)
      CALL mp_bcast(x_q, meta_ionode_id, world_comm)
      ! ENDIF
      !
      DEALLOCATE(xqaux)
      DEALLOCATE(wqaux)
      !
    ENDIF ! qplot
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE elph_read_qpoints
  !----------------------------------------------------------------------------
  !
END SUBROUTINE elph_readin
