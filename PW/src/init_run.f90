!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE init_run()
  !----------------------------------------------------------------------------
  !
  USE klist,              ONLY : nkstot
  USE symme,              ONLY : sym_rho_init
  USE wvfct,              ONLY : nbnd, et, wg, btype
  USE control_flags,      ONLY : lmd, gamma_only, smallmem, ts_vdw, mbd_vdw, io_level
  USE gvect,              ONLY : g, gg, mill, gcutm, ig_l2g, ngm, ngm_g, &
                                 gshells, gstart ! to be comunicated to the Solvers if gamma_only
  USE gvecs,              ONLY : gcutms, ngms
  USE cell_base,          ONLY : at, bg, set_h_ainv
  USE cellmd,             ONLY : lmovecell
  USE dynamics_module,    ONLY : allocate_dyn_vars
  USE paw_variables,      ONLY : okpaw
  USE paw_init,           ONLY : paw_init_onecenter, allocate_paw_internals
#if defined(__MPI)
  USE paw_init,           ONLY : paw_post_init
#endif
  USE bp,                 ONLY : allocate_bp_efield, bp_global_map
  USE fft_base,           ONLY : dfftp, dffts
  USE recvec_subs,        ONLY : ggen, ggens
  USE wannier_new,        ONLY : use_wannier    
  USE dfunct,             ONLY : newd
  USE esm,                ONLY : do_comp_esm, esm_init
  USE tsvdw_module,       ONLY : tsvdw_initialize
  USE Coul_cut_2D,        ONLY : do_cutoff_2D, cutoff_fact 
  USE lsda_mod,           ONLY : nspin
  USE xc_lib,             ONLY : xclib_dft_is_libxc, xclib_init_libxc
  USE buffers,            ONLY : open_buffer
  USE io_files,           ONLY : iunwfc, nwordwfc
  USE mod_sirius
  !
  USE control_flags,      ONLY : use_gpu
  USE dfunct_gpum,        ONLY : newd_gpu
  USE wvfct_gpum,         ONLY : using_et, using_wg, using_wg_d
  USE gvect_gpum,         ONLY : using_g, using_gg, using_g_d, using_gg_d, &
                                 using_mill, using_mill_d
  !
  IMPLICIT NONE
  LOGICAL exst_file,exst_mem
  !
  CALL start_clock( 'init_run' )
  !
  ! ... calculate limits of some indices, used in subsequent allocations
  !
  CALL pre_init()
  !
  ! ... determine the data structure for fft arrays
  !
  CALL data_structure( gamma_only )
  !
  ! ... print a summary and a memory estimate before starting allocating
  !
  CALL summary()
  CALL memory_report()
  !
  ! ... allocate memory for G- and R-space fft arrays
  !
  CALL allocate_fft()
  !
  ! ... generate reciprocal-lattice vectors and fft indices
  !
  IF( smallmem ) THEN
     CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
          g, gg, mill, ig_l2g, gstart, no_global_sort = .TRUE. )
  ELSE
     CALL ggen( dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
       g, gg, mill, ig_l2g, gstart )
  END IF
  CALL ggens( dffts, gamma_only, at, g, gg, mill, gcutms, ngms )
  if (gamma_only) THEN
     ! ... Solvers need to know gstart
     call export_gstart_2_solvers(gstart)
  END IF

#if defined(__CUDA)
  ! All these variables are actually set by ggen which has intent out
  CALL using_mill(2); CALL using_mill_d(0); ! updates mill indices,
  CALL using_g(2);    CALL using_g_d(0);    ! g and gg that are used almost only after
  CALL using_gg(2);   CALL using_gg_d(0)    ! a single initialization .
                                            ! This is a trick to avoid checking for sync everywhere.
#endif
  !
  IF (do_comp_esm) CALL esm_init()
  !
  ! ... setup the 2D cutoff factor
  !
  IF (do_cutoff_2D) CALL cutoff_fact()
  !
  CALL gshells ( lmovecell )
  !
  ! ... variable initialization for parallel symmetrization
  !
  CALL sym_rho_init (gamma_only )
  !
  ! ... allocate memory for all other arrays (potentials, wavefunctions etc)
  !
  CALL allocate_nlpot()
  IF (okpaw) THEN
     CALL allocate_paw_internals()
     CALL paw_init_onecenter()
  ENDIF
  CALL allocate_locpot()
  CALL allocate_bp_efield()
  CALL bp_global_map()
  !
  call plugin_initbase()
  !
  ALLOCATE( et( nbnd, nkstot ) , wg( nbnd, nkstot ), btype( nbnd, nkstot ) )
  !
  et(:,:) = 0.D0
  CALL using_et(2)
  !
  wg(:,:) = 0.D0
  CALL using_wg(2)
#if defined(__CUDA)
  ! Sync here. Shouldn't be done and will be removed ASAP.
  CALL using_wg_d(0)
#endif
  !
  btype(:,:) = 1
  !
  IF (ts_vdw .or. mbd_vdw) THEN
     CALL tsvdw_initialize()
     CALL set_h_ainv()
  END IF
  !
  CALL allocate_wfc_k()
  CALL openfil()
  !
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_init_libxc( nspin )
  !
  CALL hinit0()
  !
  CALL potinit()
  !
  IF (use_sirius_scf.OR.always_setup_sirius) THEN
    CALL clear_sirius
    CALL setup_sirius
  ENDIF
  !
  IF (use_sirius_scf.OR.always_setup_sirius) THEN
    CALL sirius_initialize_kset(ks_handler)
  ENDIF
  !
  IF ( use_gpu ) THEN
    !
    CALL newd_gpu()
    !
    CALL wfcinit_gpu()
    !
  ELSE
    !
    CALL newd()
    !
    IF (use_sirius_scf.OR.use_sirius_nlcg) THEN
      CALL sirius_initialize_subspace(gs_handler, ks_handler)
      CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
    ELSE
    CALL wfcinit()
    ENDIF
    !
  END IF
  !
  IF(use_wannier) CALL wannier_init()
  !
#if defined(__MPI)
  ! Cleanup PAW arrays that are only used for init
  !IF (.NOT.(use_sirius_scf.OR.use_sirius_nlcg)) THEN
  IF (okpaw) CALL paw_post_init() ! only parallel!
  !ENDIF
#endif
  !
  IF ( lmd ) CALL allocate_dyn_vars()
  !
  CALL stop_clock( 'init_run' )
  !
  RETURN
  !
END SUBROUTINE init_run
  !
!----------------------------------------------------------------------------
SUBROUTINE pre_init()
  !----------------------------------------------------------------------------
  !
  USE ions_base,        ONLY : nat, nsp, ityp
  USE uspp_param,       ONLY : upf, lmaxkb, nh, nhm, nbetam
  USE uspp,             ONLY : nkb, nkbus
  IMPLICIT NONE
  INTEGER :: na, nt, nb
  !
  !     calculate the number of beta functions for each atomic type
  !
  lmaxkb = - 1
  DO nt = 1, nsp
     !
     nh (nt) = 0
     !
     ! do not add any beta projector if pseudo in 1/r fmt (AF)
     IF ( upf(nt)%tcoulombp ) CYCLE 
     !
     DO nb = 1, upf(nt)%nbeta
        nh (nt) = nh (nt) + 2 * upf(nt)%lll(nb) + 1
        lmaxkb = MAX (lmaxkb, upf(nt)%lll(nb) )
     ENDDO
     !
  ENDDO
  !
  ! calculate the maximum number of beta functions
  !
  nhm = MAXVAL (nh (1:nsp))
  nbetam = MAXVAL (upf(:)%nbeta)
  !
  ! calculate the number of beta functions of the solid
  !
  nkb = 0
  nkbus = 0
  do na = 1, nat
     nt = ityp(na)
     nkb = nkb + nh (nt)
     if (upf(nt)%tvanp) nkbus = nkbus + nh (nt)
  enddo


END SUBROUTINE pre_init
