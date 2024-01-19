!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kcw_run_nscf (do_band)
  !-----------------------------------------------------------------------
  !
  ! This is the main driver of the PWscf program called from the KCW code.
  !
  USE control_flags,   ONLY : conv_ions, restart, iverbosity, isolve
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : tmp_dir, wfc_dir
  USE fft_types,       ONLY : fft_type_allocate
  USE fft_base,        ONLY : dffts, dfftp
  USE cell_base,       ONLY : at, bg
  USE gvect,           ONLY : gcutm
  USE gvecs,           ONLY : gcutms
  USE mp_bands,        ONLY : intra_bgrp_comm, nyfft
  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : ethr_nscf
  USE control_kcw,     ONLY : tmp_dir_kcwq
  USE klist,           ONLY : nelec
  USE paw_variables,   ONLY : okpaw
  USE mp_pools,        ONLY : inter_pool_comm, intra_pool_comm, npool
  USE mp_images,       ONLY : nproc_image, intra_image_comm
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  INTEGER :: verbosity_save
  INTEGER :: ik
  !
  CALL start_clock( 'kcw_run_nscf' )
  !
  CALL clean_pw (.FALSE.)
  !
  CALL close_files (.TRUE.)
  !
  ! From now on, work only on the KC directory
  !
  wfc_dir = tmp_dir_kcwq
  tmp_dir = tmp_dir_kcwq
  !
  ! Setting the values for the NSCF run
  !
  startingconfig = 'input'
  starting_pot   = 'file'
  starting_wfc   = 'atomic'
  restart        = .FALSE.
  conv_ions      = .TRUE.
  ethr_nscf      = 1.0D-9 / nelec 
  isolve         = 0
  !
  ! iverbosity is used by the PWscf routines
  IF (iverbosity.LE.2) THEN
     ! temporarily change the value of iverbosity
     ! in order to have less output from the PWscf routines
     verbosity_save = iverbosity
     iverbosity = 0
  ENDIF
  !
  CALL fft_type_allocate ( dfftp, at, bg, gcutm,  intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft)
  !
  CALL setup_nscf ( .FALSE., xq, .TRUE. )
  !
  CALL init_run()
  !
  IF (do_band) THEN
<<<<<<< HEAD
!#if defined(__SIRIUS)
!     CALL setup_sirius()
!     !
!     ! create k-point set
!     ! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice and
!     !          without x2 multiplication for the lsda case
!     CALL sirius_create_kset(sctx, num_kpoints, kpoints, wkpoints, .FALSE., ks_handler1)    
!     CALL sirius_initialize_kset(ks_handler1)
!     !
!     ! create ground-state class    
!     CALL sirius_create_ground_state(ks_handler1, gs_handler1)
!     CALL put_density_to_sirius(gs_handler1)
!     IF (okpaw) THEN
!       CALL put_density_matrix_to_sirius(gs_handler1)
!       CALL sirius_generate_density(gs_handler1, paw_only=.TRUE.)
!     ENDIF
!     CALL sirius_generate_effective_potential(gs_handler1)
!     CALL sirius_initialize_subspace(gs_handler1, ks_handler1)
!     CALL sirius_find_eigen_states(gs_handler1, ks_handler1, iter_solver_tol=1.d-13)!, iter_solver_steps=100)
!     !save wfs
!     CALL get_wave_functions_from_sirius(ks_handler1)
!#else 
=======
#if defined(__SIRIUS)
     !WARNING: This is a copy of the last part of setup sirius. Factorize?
     !WRITE(*,*) "Gonna create kset"
     WRITE(*,*) "nkpt", num_kpoints
     WRITE(*,*) kpoints(:,:)
     WRITE(*,*) wkpoints(:)

     CALL setup_sirius()
    !
    ! create context of simulation
     !WRITE(*,*) "Gonna create ctx"
     !CALL sirius_initialize_context(sctx)

     !
     ! create k-point set
     ! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice and
     !          without x2 multiplication for the lsda case
     WRITE(*,*) "Gonna create kset"
     CALL sirius_create_kset(sctx, num_kpoints, kpoints, wkpoints, .FALSE., ks_handler1)    !
     CALL sirius_initialize_kset(ks_handler1)
     !
     ! create ground-state class    
     WRITE(*,*) "Gonna create ground state"
     CALL sirius_create_ground_state(ks_handler1, gs_handler1)
     WRITE(*,*) "Gonna put density to sirius"
     CALL put_density_to_sirius(gs_handler1)
     IF (okpaw) THEN
       CALL put_density_matrix_to_sirius(gs_handler1)
       CALL sirius_generate_density(gs_handler1, paw_only=.TRUE.)
     ENDIF
     WRITE(*,*) "Gonna generate effective potential"
     CALL sirius_generate_effective_potential(gs_handler1)
     WRITE(*,*) "Gonna initialize subspace"
     CALL sirius_initialize_subspace(gs_handler1, ks_handler1)
     WRITE(*,*) "Gonna find eigenstates"
     CALL sirius_find_eigen_states(gs_handler1, ks_handler1, iter_solver_tol=1.d-13)!, iter_solver_steps=100)
#else 
>>>>>>> c065fddbb542b256b77efd7b9228d0b5243f7f8b
     CALL non_scf()
!#endif 
     CALL punch( 'all' )
  ENDIF
  !
  IF (iverbosity.EQ.0) iverbosity = verbosity_save 
  !
  CALL close_files(.TRUE.)
  !
  CALL stop_clock( 'kcw_run_nscf' )
  !
  RETURN
  ! 
END SUBROUTINE kcw_run_nscf
