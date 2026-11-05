!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE crpa_run_nscf (do_band)
  !-----------------------------------------------------------------------
  !
  ! This is the main driver of the PWscf program called from the CRPA code.
  !
  USE mp,              ONLY : mp_barrier
  USE mp_bands,        ONLY : intra_bgrp_comm, nyfft
  USE mp_images,       ONLY : intra_image_comm
  USE fft_types,       ONLY : fft_type_allocate
  USE fft_base,        ONLY : dffts, dfftp
  USE control_flags,   ONLY : conv_ions, restart, iverbosity, io_level
  USE starting_scf,    ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir, wfc_dir
  USE lsda_mod,        ONLY : nspin
  USE check_stop,      ONLY : check_stop_now
  USE cell_base,       ONLY : at, bg
  USE gvect,           ONLY : gcutm
  USE gvecs,           ONLY : gcutms
  USE io_global,       ONLY : stdout
  USE scf,             ONLY : vrs
  USE ldaU,            ONLY : lda_plus_u_kind
  USE klist,           ONLY : nelec
  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : lgamma, ethr_nscf, conv_thr_nscf
  USE lr_symm_base,    ONLY : nsymq, invsymq
  USE crpacom,         ONLY : tmp_dir_save, tmp_dir_crpa
  USE crpa_pert,       ONLY : pert_basis
  USE start_k,         ONLY : nks_start, wk_start
  USE noncollin_module,ONLY : noncolin, domag
  USE ldaU,            ONLY : lda_plus_u_kind
  USE rism_module,     ONLY : lrism, rism_set_restart
  USE kinds,           ONLY : DP
  USE random_numbers,  ONLY : randy
  USE noncollin_module,ONLY : noncolin, domag
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  LOGICAL :: exst
  REAL(DP) :: arg
  INTEGER :: verbosity_save
  !
  CALL start_clock( 'crpa_run_nscf' )
  !
  CALL clean_pw (.FALSE.)
  !
  CALL close_files (.TRUE.)
  !
  ! From now on, work only on the CRPA directory
  !
  wfc_dir = tmp_dir_crpa
  tmp_dir = tmp_dir_crpa
  !
  ! Setting the values for the NSCF run
  !
  startingconfig = 'input'
  starting_pot   = 'file'
  starting_wfc   = 'atomic'
  restart        = .FALSE.
  conv_ions      = .TRUE.
  ethr_nscf      = MAX(1.D-13, 0.1D0 * MIN( 1.D-2, conv_thr_nscf / nelec ))
  !
  ! Reseed the random-number generator, as PHonon does at every q (set_irr:
  ! arg = randy(0)). When the pseudopotentials carry no atomic wavefunctions
  ! the NSCF starts from random vectors (wfcinit); reseeding makes them
  ! independent of whatever drew random numbers earlier in the run (e.g.
  ! another q point), so a q point gives the same result in any run.
  arg = randy(0)
  !
  IF (lrism) CALL rism_set_restart()
  !
  ! iverbosity is used by the PWscf routines
  IF (iverbosity <= 2) THEN
     ! temporarily change the value of iverbosity
     ! in order to have less output from the PWscf routines
     verbosity_save = iverbosity
     iverbosity = 0
  ENDIF
  !
  IF (lgamma) THEN
     WRITE(stdout,'(/5x,"Performing NSCF calculation at all points k...")')
  ELSE
     WRITE(stdout,'(/5x,"Performing NSCF calculation at all points k and k+q...")')
  ENDIF
  !
  CALL fft_type_allocate ( dfftp, at, bg, gcutm,  intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft)
  !
  ! Wannier basis: keep exactly the k points of the SCF (the Wannier mesh) plus
  ! k+q. Without dont_use_ibz, irreducible_BZ treats them as the IBZ of the full
  ! lattice group and adds rotated, off-mesh k with split weights whenever the
  ! mesh is not invariant under the lattice rotations (e.g. 3x1x1): there the
  ! Wannier interpolation is not exact and the pair density, which weights every
  ! k uniformly, is wrong. irreducible_BZ also normalizes the weights to 1, so
  ! skipping it needs the input weights (e.g. 1.0 per k of an explicit list)
  ! normalized here.
  !
  IF (pert_basis == 'wannier' .AND. nks_start > 0) &
     wk_start(1:nks_start) = wk_start(1:nks_start) / SUM(wk_start(1:nks_start))
  CALL setup_nscf ( .FALSE., xq, pert_basis == 'wannier' )
  !
  CALL init_run()
  !
  IF (do_band) THEN
     IF (noncolin.AND.domag) THEN
        ! this subroutine calls c_bands_ph, which applies
        ! the time-reversal operator
        CALL non_scf_ph()
     ELSE
        CALL non_scf()
     ENDIF
     CALL punch( 'all' )
  ENDIF
  !
  IF (iverbosity == 0) iverbosity = verbosity_save
  !
  IF (io_level > 0) THEN
     CALL close_files(.TRUE.)
  ELSE
     CALL mp_barrier(intra_image_comm)
  ENDIF
  !
  WRITE(stdout,'(5x,"Done!")')
  !
  CALL stop_clock( 'crpa_run_nscf' )
  !
  RETURN
  !
END SUBROUTINE crpa_run_nscf
