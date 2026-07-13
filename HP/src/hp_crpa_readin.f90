!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE hp_crpa_readin()
  !----------------------------------------------------------------------------
  !
  !  This routine reads the CRPA-specific input (the &CRPA namelist and,
  !  if qplot, the q-point card). It is called from hp_readin when
  !  method = 'CRPA', after the &INPUTHP namelist has been read and
  !  broadcast, and before read_file().
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : meta_ionode, meta_ionode_id
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE mp_images,        ONLY : my_image_id
  USE io_files,         ONLY : tmp_dir, create_directory
  USE control_flags,    ONLY : isolve, io_level
  USE control_lr,       ONLY : lrpa, reduce_io, niter_ph, maxter, alpha_mix, &
                               nmix_ph, tr2_ph, thresh_init, conv_thr_nscf
  USE qpoint,           ONLY : nq1, nq2, nq3, start_q, last_q, x_q, nqs
  USE ldaU_hp,          ONLY : nq1_hp => nq1, nq2_hp => nq2, nq3_hp => nq3,     &
                               start_q_hp => start_q, last_q_hp => last_q,      &
                               alpha_mix_hp => alpha_mix, nmix_hp => nmix,      &
                               niter_max_hp => niter_max,                       &
                               thresh_init_hp => thresh_init,                   &
                               tmp_dir_hp, tmp_dir_save
  USE crpacom,          ONLY : wannier_seedname, active_bands_min, active_bands_max, &
                               wan_dist_thr, wan_dist_thr_mels, folder_wan_Rr,  &
                               write_wan_Rr, w_freq
  USE crpa_pert,        ONLY : pert_basis
  USE crpa_qpoints,     ONLY : qplot
  USE constrained_dfpt, ONLY : cdfpt, cdfpt_active_space, cdfpt_bands_min, cdfpt_bands_max
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=6) :: int_to_char
  CHARACTER(len=80) :: diagonalization
  !! Diagonalization method for NSCF calculation
  CHARACTER(LEN = 512) :: line
  !! Line in input file
  INTEGER :: ios
  !! integer variable for I/O control
  INTEGER :: ios2
  !! INTEGER variable for I/O control
  !
  REAL(DP), ALLOCATABLE :: xqaux(:,:)
  INTEGER, ALLOCATABLE :: wqaux(:)
  INTEGER :: nqaux, iq, ipol
  !
  REAL(DP) :: tr2
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  ! wan_dist_thr: Wannier function distance threshold for computing self-consistent response.
  !           The number of DFPT perturbations is determined by this parameter.
  !           (Default: 6.D-4 bohr)
  ! wan_dist_thr_mels: Wannier function distance threshold for computing matrix elements.
  !                 This parameter can be much larger than wan_dist_thr, still does not affect
  !                 computational cost much because they only affect the matrix element
  !                 calculation after the self-consistent response calculation.
  !                 (Default: same as wan_dist_thr)
  ! folder_wan_Rr : Folder where the collected real-space Wannier functions are saved.
  !                 If write_wan_Rr is .true. this folder is created and the Wannier functions
  !                 are written to it. If write_wan_Rr is .false. this folder should already
  !                 exist and contain the files prefix.wan_Rr.xxx and prefix.wan_R.dat.
  !                 (Default: tmp_dir_hp)
  ! write_wan_Rr : If .true. write the collected real-space Wannier functions to file.
  !                If .false. use the existing files in folder_wan_Rr. In this case,
  !                prefix.wan_Rr.xxx and prefix.wan_R.dat files must exist inside folder_wan_Rr.
  !                (Default: .true.)
  ! w_freq : Complex frequency value for finite-frequency DFPT. (Default: 0 (static DFPT))
  !
  NAMELIST / CRPA / pert_basis, active_bands_min, active_bands_max,          &
                    wannier_seedname, folder_wan_Rr, write_wan_Rr,           &
                    wan_dist_thr, wan_dist_thr_mels, tr2, lrpa, w_freq,      &
                    diagonalization, conv_thr_nscf, reduce_io,               &
                    cdfpt, cdfpt_active_space, cdfpt_bands_min, cdfpt_bands_max, &
                    qplot
  !
  ! Set default values for variables in namelist
  !
  tr2                = 1.D-20
  diagonalization    = 'david'
  conv_thr_nscf      = 1.D-9
  wannier_seedname   = ''
  folder_wan_Rr      = ''
  write_wan_Rr       = .TRUE.
  cdfpt              = .FALSE.
  cdfpt_active_space = ''
  cdfpt_bands_min    = 0
  cdfpt_bands_max    = 0
  pert_basis         = 'bands'
  active_bands_min   = 0
  active_bands_max   = 0
  reduce_io          = .FALSE.
  lrpa               = .FALSE.   ! Needed in dv_of_drho
  wan_dist_thr       = 6.D-4
  wan_dist_thr_mels  = -999.0
  w_freq             = (0.d0, 0.d0)
  qplot              = .FALSE.
  !
  ! Reading the namelist &CRPA from the input unit already opened by
  ! hp_readin (input_from_file connects it to unit 5), positioned right
  ! after the &INPUTHP namelist.
  !
  IF (meta_ionode) THEN
    !
    READ(5, crpa, IOSTAT = ios)
    ios2 = 0
    IF (ios /= 0) THEN
      BACKSPACE(5)
      READ(5, '(A512)', IOSTAT = ios2) line
    ENDIF
    IF (ios2 /= 0) CALL errore('hp_crpa_readin', 'Could not find namelist &crpa', 2)
    IF (ios /= 0) THEN
      CALL errore('hp_crpa_readin', 'Bad line in namelist &crpa: "' &
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
    CALL errore('hp_crpa_readin','diagonalization '//trim(diagonalization)//' not implemented',1)
  END SELECT
  !
  ! The names of the lrcom variables have the _ph suffix, so we copy to them and use.
  ! The DFPT-iteration knobs (niter_max, nmix, alpha_mix, thresh_init) are read from
  ! &INPUTHP into ldaU_hp; the CRPA code paths use the control_lr copies.
  !
  tr2_ph = tr2
  nmix_ph = nmix_hp
  niter_ph = niter_max_hp
  alpha_mix(1:maxter) = alpha_mix_hp(1:maxter)
  thresh_init = thresh_init_hp
  !
  ! Copy the q-grid variables from ldaU_hp (read from &INPUTHP) to the
  ! qpoint module used by the CRPA code paths.
  !
  nq1 = nq1_hp
  nq2 = nq2_hp
  nq3 = nq3_hp
  start_q = start_q_hp
  last_q = last_q_hp
  !
  ! NOTE: Both conv_thr_nscf and ethr_nscf exist in control_lr. The CRPA NSCF
  ! calculation uses conv_thr_nscf (read from &CRPA), while the LRT (HP) NSCF
  ! uses ethr_nscf (read from &INPUTHP).
  !
  IF (wan_dist_thr_mels < 0.d0) wan_dist_thr_mels = wan_dist_thr
  !
  ! Broadcast input parameters over all processors
  !
  CALL hp_crpa_bcast_input()
  !
  ! If qplot, read and broadcast q points
  !
  IF (qplot) THEN
    IF (meta_ionode) THEN
      ios = 0
      READ (5, *, iostat = ios) nqaux
    ENDIF
    !
    CALL mp_bcast(nqaux, meta_ionode_id, world_comm )
    CALL mp_bcast(ios, meta_ionode_id, world_comm )
    CALL errore ('hp_crpa_readin', 'reading nq', ABS (ios) )
    !
    ALLOCATE(xqaux(3,nqaux))
    ALLOCATE(wqaux(nqaux))
    IF (meta_ionode) THEN
      DO iq=1, nqaux
        READ (5, *, iostat = ios) (xqaux (ipol,iq), ipol=1,3), wqaux(iq)
      ENDDO
    ENDIF
    CALL mp_bcast(ios, meta_ionode_id, world_comm )
    CALL errore ('hp_crpa_readin', 'reading xq', ABS (ios) )
    CALL mp_bcast(xqaux, meta_ionode_id, world_comm )
    CALL mp_bcast(wqaux, meta_ionode_id, world_comm )
    !
    nqs = nqaux
    ALLOCATE(x_q(3,nqs))
    x_q(:,1:nqs)=xqaux(:,1:nqs)
    CALL mp_bcast(x_q, meta_ionode_id, world_comm)
    !
    DEALLOCATE(xqaux)
    DEALLOCATE(wqaux)
    !
  ENDIF ! qplot
  !
  ! IO option
  !
  IF (reduce_io) io_level = 0
  !
  ! Initialize folders (image-aware HP-style scratch directory)
  !
  tmp_dir_save = tmp_dir
  tmp_dir_hp = trimcheck( TRIM (tmp_dir) // 'HP' // int_to_char(my_image_id) )
  CALL create_directory(tmp_dir_hp)
  !
  IF (TRIM(folder_wan_Rr) == '') THEN
    ! If folder_wan_Rr is not set, use the default tmp_dir_hp
    folder_wan_Rr = tmp_dir_hp
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
  RETURN
  !
END SUBROUTINE hp_crpa_readin
