!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE crpa_readin()
  !----------------------------------------------------------------------------
  !
  !  This routine reads the input parameters for CRPA and reads
  !  the data produced by PWscf.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : meta_ionode, meta_ionode_id, qestdin, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE mp_images,        ONLY : my_image_id
  USE check_stop,       ONLY : max_seconds
  USE io_files,         ONLY : tmp_dir, prefix, create_directory
  USE open_close_input_file, ONLY : open_input_file, close_input_file
  USE paw_variables,    ONLY : okpaw
  USE control_flags,    ONLY : iverbosity, isolve, io_level
  USE control_lr,       ONLY : ethr_nscf, lrpa, reduce_io, niter_ph, maxter, alpha_mix, &
                               nmix_ph, tr2_ph, thresh_init, conv_thr_nscf
  USE qpoint,           ONLY : nq1, nq2, nq3, start_q, last_q
  USE crpacom,          ONLY : tmp_dir_crpa, tmp_dir_save, wannier_seedname, filU, &
                               active_bands_min, active_bands_max, dist_thr, &
                               dist_thr_large, folder_wan_Rr, write_wan_Rr, w_freq, &
                               q0div_treatment, crpa_mode
  USE crpa_pert,        ONLY : pert_basis
  USE crpa_qpoints,     ONLY : qplot
  USE qpoint,           ONLY : x_q, nqs
  USE constrained_dfpt, ONLY : cdfpt, cdfpt_active_space, cdfpt_bands_min, cdfpt_bands_max, &
                               cdfpt_validate_input
  ! conv_thr_chi, , find_atpert, skip_atom, skip_type, &
!                                equiv_type, background, compute_crpa,         &
!                                perturb_only_atom, sum_pertq, determine_num_pert_only,  &
!                                skip_equivalence_q, dist_thr,  &
!                                disable_type_analysis, docc_thr, num_neigh, lmin, rmax, &
!                                determine_q_mesh_only
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: outdir
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(len=80) :: diagonalization
  !! Diagonalization method for NSCF calculation
  CHARACTER(LEN = 512) :: line
  !! Line in input file
  INTEGER :: ios
  !! integer variable for I/O control
  INTEGER :: ios2
  !! INTEGER variable for I/O control
  INTEGER :: iter
  !! counter on iterations
  !
  REAL(DP), ALLOCATABLE :: xqaux(:,:)
  INTEGER, ALLOCATABLE :: wqaux(:)
  INTEGER :: nqaux, iq, ipol
  !
  INTEGER :: nmix, niter_max
  REAL(DP) :: tr2
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  ! dist_thr: Wannier function distance threshold for computing self-consistent response.
  !           The number of DFPT perturbations is determined by this parameter.
  !           (Default: 6.D-4 bohr)
  ! dist_thr_large: Wannier function distance threshold for computing matrix elements.
  !                 This parameter can be much larger than dist_thr, still does not affect
  !                 computational cost much because they only affect the matrix element
  !                 calculation after the self-consistent response calculation.
  !                 (Default: same as dist_thr)
  ! folder_wan_Rr : Folder where the collected real-space Wannier functions are saved.
  !                 If write_wan_Rr is .true. this folder is created and the Wannier functions
  !                 are written to it. If write_wan_Rr is .false. this folder should already
  !                 exist and contain the files prefix.wan_Rr.xxx and prefix.wan_R.dat.
  !                 (Default: tmp_dir_crpa)
  ! write_wan_Rr : If .true. write the collected real-space Wannier functions to file.
  !                If .false. use the existing files in folder_wan_Rr. In this case,
  !                prefix.wan_Rr.xxx and prefix.wan_R.dat files must exist inside folder_wan_Rr.
  !                (Default: .true.)
  ! w_freq : Complex frequency value for finite-frequency DFPT. (Default: 0 (static DFPT))
  ! q0div_treatment : Correction scheme for the q->0 divergence of the Coulomb kernel.
  !                   (Default: 'HL' - Hybertsen Louie correction)
  ! crpa_mode : Calculation mode of the Uijkl(R) elements.
  !             (Default: 'full' (all possible elements within the given distance thresholds)
  !
  NAMELIST / INPUTCRPA / prefix, outdir, lrpa, niter_max, alpha_mix, nq1, nq2, nq3, &
                         nmix, start_q, last_q, thresh_init, tr2, &
                         active_bands_min, active_bands_max, diagonalization, iverbosity, &
                         conv_thr_nscf, reduce_io, wannier_seedname, pert_basis, filU, &
                         qplot, dist_thr, dist_thr_large, cdfpt, cdfpt_active_space, &
                         cdfpt_bands_min, cdfpt_bands_max, folder_wan_Rr, &
                         write_wan_Rr, w_freq, q0div_treatment, crpa_mode
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
  IF (meta_ionode) CALL input_from_file ()
  !
  ! Set default values for variables in namelist
  !
  prefix             = 'pwscf'
!   conv_thr_chi       = 1.D-5
  thresh_init        = 1.D-14
  tr2                = 1.D-20
  diagonalization    = 'david'
  conv_thr_nscf      = 1.D-9
  wannier_seedname   = ''
  folder_wan_Rr      = ''
  write_wan_Rr       = .TRUE.
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
  nq1                = 1
  nq2                = 1
  nq3                = 1
  start_q            = 1
  last_q             = -1
  iverbosity         = 1
  niter_max          = 100
  alpha_mix(:)       = 0.D0
  alpha_mix(1)       = 0.7D0
  nmix               = 12
  cdfpt              = .FALSE.
  cdfpt_active_space = ''
  cdfpt_bands_min    = 0
  cdfpt_bands_max    = 0
  pert_basis         = 'bands'
  active_bands_min   = 0
  active_bands_max   = 0
  reduce_io          = .FALSE.
  filU               = 'crpa'
!   max_seconds        = 1.E+7_DP
  lrpa               = .FALSE.   ! Needed in dv_of_drho
  dist_thr           = 6.D-4
  dist_thr_large     = -999.0
  w_freq             = (0.d0, 0.d0)
  q0div_treatment    = 'HL'
  crpa_mode          = 'full'
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  ! Reading the namelist inputcrpa and check
  !
  IF (meta_ionode) THEN
    !
    ! ... Input from file (ios=0) or standard input (ios=-1) on unit "qestdin"
    !
    ios = open_input_file (  )
    !
    READ(qestdin, inputcrpa, IOSTAT = ios)
    ios2 = 0
    IF (ios /= 0) THEN
      BACKSPACE(qestdin)
      READ(qestdin, '(A512)', IOSTAT = ios2) line
    ENDIF
    IF (ios2 /= 0) CALL errore('crpa_readin', 'Could not find namelist &inputcrpa', 2)
    IF (ios /= 0) THEN
      CALL errore('crpa_readin', 'Bad line in namelist &inputcrpa: "' &
                 & //TRIM(line)//'" (error could be in the previous line)', 1)
    ENDIF
  ENDIF ! meta_ionode
  !
  ! diagonalization option
  !
  SELECT CASE(TRIM(diagonalization))
  CASE ('david','davidson')
    isolve = 0
  CASE ('cg')
    isolve = 1
  CASE DEFAULT
    CALL errore('crpa_readin','diagonalization '//trim(diagonalization)//' not implemented',1)
  END SELECT
  !
  ! Setup tmp_dir
  !
  tmp_dir = trimcheck (outdir)
  !
  ! The name of the lrcom variables have _ph suffix, so we copy to them and use.
  ! TODO: Change PHonon to do this conversion
  !
  tr2_ph = tr2
  nmix_ph = nmix
  niter_ph = niter_max
  !
  IF (dist_thr_large < 0.d0) dist_thr_large = dist_thr
  !
  ! Broadcast input parameters over all processors
  !
  CALL crpa_bcast_input()
  !
  ! If qplot, read broadcast q points
  !
  IF (qplot) THEN
    IF (meta_ionode) THEN
      ios = 0
      READ (qestdin, *, iostat = ios) nqaux
    ENDIF
    !
    CALL mp_bcast(nqaux, meta_ionode_id, world_comm )
    CALL mp_bcast(ios, meta_ionode_id, world_comm )
    CALL errore ('crpa_readin', 'reading nq', ABS (ios) )
    !
    ALLOCATE(xqaux(3,nqaux))
    ALLOCATE(wqaux(nqaux))
    IF (meta_ionode) THEN
      DO iq=1, nqaux
        READ (qestdin, *, iostat = ios) (xqaux (ipol,iq), ipol=1,3), wqaux(iq)
      ENDDO
    ENDIF
    CALL mp_bcast(ios, meta_ionode_id, world_comm )
    CALL errore ('crpa_readin', 'reading xq', ABS (ios) )
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
    nqs = nqaux
    ALLOCATE(x_q(3,nqs))
    x_q(:,1:nqs)=xqaux(:,1:nqs)
    CALL mp_bcast(x_q, meta_ionode_id, world_comm)
    ! ENDIF
    DEALLOCATE(xqaux)
    DEALLOCATE(wqaux)
    !
  ENDIF ! qplot
  !
  IF (meta_ionode) ios = close_input_file ()
  !
  ! IO option
  !
  IF (reduce_io) io_level = 0
  !
  ! Initialize folders
  !
  tmp_dir_save = tmp_dir
  tmp_dir_crpa = trimcheck( TRIM (tmp_dir) // '_crpa' // int_to_char(my_image_id) )
  CALL create_directory(tmp_dir_crpa)
  !
  IF (TRIM(folder_wan_Rr) == '') THEN
    ! If folder_wan_Rr is not set, use the default tmp_dir_crpa
    folder_wan_Rr = tmp_dir_crpa
  ELSE
    ! Otherwise, create the directory specified in folder_wan_Rr
    ! FIXME: When using image parallelization, create_directory checks file writable on
    !        ionode (not meta_ionode), so the root of every image does the check. This
    !        creates conflict and IO error. This should be fine if folder_wan_Rr is
    !        different for each image, but if it is the same, it will fail.
    !        As a workaround, we call create_directory only on the root image.
    IF (my_image_id == 0) CALL create_directory(folder_wan_Rr)
  ENDIF
  folder_wan_Rr = trimcheck(folder_wan_Rr)
  !
  ! Here we finished the reading of the input file.
  ! Now allocate space for pwscf variables, read and check them.
  !
  ! Read various data produced by PWscf.
  ! In particular, read the unperturbed occupation matrices
  ! via calling the routine read_rho.
  ! read_file calls init_tab_atwfc which initializes tab_atwfc.
  !
  CALL read_file()
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
  USE klist,            ONLY : lgauss, ltetra, two_fermi_energies
  USE control_flags,    ONLY : gamma_only, tqr
  USE fixed_occ,        ONLY : tfixed_occ
  USE cellmd,           ONLY : lmovecell
  USE noncollin_module, ONLY : i_cons, noncolin
  USE mp_bands,         ONLY : nbgrp
  USE xc_lib,           ONLY : xclib_dft_is
!   USE ldaU,             ONLY : lda_plus_u, Hubbard_projectors, lda_plus_u_kind, Hubbard_J0, &
!                                Hubbard_V, is_hubbard_back
  !
  IMPLICIT NONE
  !
  IF (tr2_ph <= 0.D0) CALL errore ('crpa_readin', ' Wrong tr2_ph', 1)
  !
  IF (thresh_init <= 0.D0) CALL errore ('crpa_readin', ' Wrong thresh_init', 1)
  !
  IF (nq1 < 0 .OR. nq2 < 0 .OR. nq3 < 0) &
       CALL errore('crpa_readin','nq1, nq2, and nq3 must be greater than 0', 1)
  !
  IF (start_q <= 0 ) CALL errore('crpa_readin', ' Wrong start_q ', 1)
  !
  IF (filU == '') CALL errore('crpa_readin', ' Wrong filU', 1)
!   !
!   IF (compute_crpa .AND. ANY(perturb_only_atom(:))) &
!      CALL errore ('crpa_readin', 'compute_crpa and perturb_only_atom are not allowed to be true together', 1)
!   !
!   IF ( ANY(Hubbard_V(:,:,2).NE.0.d0) .OR. &
!        ANY(Hubbard_V(:,:,3).NE.0.d0) .OR. &
!        ANY(Hubbard_V(:,:,4).NE.0.d0) ) &
!      CALL errore ('crpa_readin', 'The CRPA code does not support DFT+U+V with the background', 1)
!   !
!   IF (ANY(is_hubbard_back(:))) CALL errore ("crpa_readin", &
!           &" Two (or more) Hubbard channels per atomic type is not implemented", 1)
!   !
!   IF (ANY(Hubbard_J0(:).NE.0.d0)) &
!      CALL errore ('crpa_readin', 'Hubbard_J0 /= 0 is not allowed.', 1)
!   !
!   IF (determine_q_mesh_only .AND. .NOT.ANY(perturb_only_atom(:))) &
!      CALL errore ('crpa_readin', 'determine_q_mesh_only can be set to .true. only if perturb_only_atom is .true. for some atom', 1)
!   !
!   IF (sum_pertq .AND. .NOT.ANY(perturb_only_atom(:))) &
!      CALL errore ('crpa_readin', 'sum_pertq can be set to .true. only if perturb_only_atom is .true. for some atom', 1)
!   !
  IF (niter_ph < 1 .OR. niter_ph > maxter) &
    CALL errore ('crpa_readin', ' Wrong niter_ph ', 1)
  !
  DO iter = 1, niter_ph
     IF ( alpha_mix(iter) < 0.D0 .OR. alpha_mix(iter) > 1.D0 ) &
        CALL errore ('crpa_readin', ' Wrong alpha_mix ', iter)
  ENDDO
!   !
!   IF (num_neigh.LT.1) CALL errore('crpa_readin','Not allowed value of num_neigh',1)
!   !
!   IF (lmin.LT.0 .OR. lmin.GT.3) CALL errore('crpa_readin','Not allowed value of lmin',1)
!   !
!   IF (nmix.LT.1) CALL errore ('crpa_readin', ' Wrong nmix ', 1)
!   !
!   IF (ltetra) CALL errore ('crpa_readin', 'CRPA with tetrahedra is not supported', 1)
!   !
!   IF (gamma_only) CALL errore('crpa_readin',&
!      & 'Cannot start from pw.x data file using Gamma-point tricks',1)
!   !
!   IF (.NOT.lda_plus_u) CALL errore('crpa_readin',&
!      & 'The CRPA code can be used only on top of DFT+Hubbard (i.e. when the HUBBARD card is used in pw.x)',1)
!   !
!   IF (lda_plus_u_kind.EQ.1) CALL errore("crpa_readin", &
!      & ' The CRPA code does not support the Liechtenstein formulation of DFT+U',1)
!   !
!   IF (Hubbard_projectors.NE."atomic" .AND. Hubbard_projectors.NE."ortho-atomic") &
!      CALL errore("crpa_readin", &
!      " The CRPA code for this Hubbard_projectors type is not implemented",1)
!   !
!   IF (lmovecell) CALL errore('crpa_readin','The CRPA code is not working after vc-relax',1)
!   !
!   IF (nbgrp > 1) CALL errore('crpa_readin', &
!      & 'band parallelization is not implemented in CRPA',1)
!   !
!   IF (i_cons /= 0) CALL errore('crpa_readin',&
!      & 'The CRPA code with constrained magnetization is not yet available',1)
!   !
!   IF (two_fermi_energies .AND. (ltetra .OR. lgauss)) CALL errore('crpa_readin', &
!      & 'The CRPA code with two Fermi energies is not available for metals',1)
!   !
!   IF (tqr) CALL errore('crpa_readin',&
!      & 'The CRPA code with Q in real space is not supported',1)
!   !
!   IF (tfixed_occ) CALL errore('crpa_readin', &
!      & 'The CRPA code with arbitrary occupations not tested',1)
!   !
!   IF ( xclib_dft_is('meta') ) CALL errore('crpa_readin',&
!      'The CRPA code with meta-GGA functionals is not yet available',1)
!   !
!   IF ( xclib_dft_is('hybrid') ) CALL errore('crpa_readin',&
!      'The CRPA code with hybrid functionals is not yet available',1)
  !
  ! Check sanity of perturbation basis definition
  !
  IF (pert_basis == '') CALL errore('crpa_readin', 'pert_basis not defined', 1)
  !
  IF (pert_basis == 'bands') THEN
    IF (active_bands_min < 1) CALL errore('crpa_readin', 'active_bands_min is not set', 1)
    IF (active_bands_max < 1) CALL errore('crpa_readin', 'active_bands_max is not set', 1)
    IF (active_bands_min > active_bands_max) CALL errore('crpa_readin', &
        'active_bands_min cannot be greater than active_bands_max', 1)
  ELSEIF (pert_basis == 'wannier') THEN
    IF (wannier_seedname == '') CALL errore('crpa_readin', 'wannier_seedname not set', 1)
  ELSE
    CALL errore('crpa_readin', 'Invalid pert_basis. Must be bands or wannier.', 1)
  ENDIF
  !
  ! Validate constrained DFPT parameters
  !
  IF (cdfpt) CALL cdfpt_validate_input()
  !
  ! Validate q->0 correction scheme
  !
  IF (TRIM(q0div_treatment) /= 'HL' .AND. TRIM(q0div_treatment) /= 'GB') THEN
     CALL errore('crpa_readin', 'Invalid q0div_treatment. Only HL and GB are allowed.', 1)
  END IF
  !
  ! Validate matrix elements calculation mode
  !
  IF (TRIM(crpa_mode) /= 'full' .AND. TRIM(crpa_mode) /= 'dHP') THEN
     CALL errore('crpa_readin', 'Invalid crpa_mode. Only full and dHP are currently allowed.', 1)
  END IF
  !
  RETURN
  !
END SUBROUTINE input_sanity
  !
END SUBROUTINE crpa_readin
