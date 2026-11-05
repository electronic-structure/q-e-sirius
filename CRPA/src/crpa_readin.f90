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
  USE ramanm,           ONLY : eth_rps
  USE output,           ONLY : fildrho
  USE control_ph,       ONLY : tmp_dir_ph, ext_recover
  USE io_global,        ONLY : meta_ionode, meta_ionode_id, qestdin, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE mp_images,        ONLY : my_image_id
  USE check_stop,       ONLY : max_seconds
  USE io_files,         ONLY : tmp_dir, prefix, create_directory
  USE open_close_input_file, ONLY : open_input_file, close_input_file
  USE paw_variables,    ONLY : okpaw
  USE control_flags,    ONLY : iverbosity, isolve, io_level
  USE control_lr,       ONLY : ethr_nscf, lrpa, lnoloc, reduce_io, niter_ph, maxter, alpha_mix, &
                               nmix_ph, tr2_ph, thresh_init, conv_thr_nscf
  USE qpoint,           ONLY : nq1, nq2, nq3, start_q, last_q
  USE crpacom,          ONLY : tmp_dir_crpa, tmp_dir_save, wannier_seedname, filU, &
                               active_bands_min, active_bands_max, dist_thr, &
                               dist_thr_large, folder_wan_Rr, write_wan_Rr, w_freq, &
                               q0div_treatment, crpa_mode, crpa_epsil, crpa_done_epsil, crpa_epsil_only, &
                               crpa_diel_llanczos, crpa_diel_lanczos_steps, crpa_diel_lanczos_restart_step, &
                               crpa_diel_lanczos_restart, crpa_diel_lanczos_save_final, &
                               crpa_freq_start, crpa_freq_end, crpa_nfreq, crpa_freq_broadening, &
                               crpa_diel_extrapolation, crpa_diel_lanczos_steps_ext, &
                               crpa_llanczos, crpa_coul_lanczos_steps, &
                               crpa_coul_extrapolation, crpa_coul_lanczos_steps_ext, &
                               crpa_coul_restart, crpa_coul_lanczos_save_final, &
                               crpa_coul_lanczos_restart_step, &
                               crpa_spectrum_only, crpa_reduce_q, crpa_w_rspace, crpa_diel_intraband, &
                               crpa_diel_flat, crpa_diel_deg_thr, &
                               crpa_xcmode, crpa_freq, current_w
  USE crpa_freq_grid,   ONLY : crpa_build_freq_grid
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
  ! w_freq : Complex frequency (Ry) of the calculation, used when no frequency grid
  !             is set (crpa_nfreq = 0). (Default: 0 (static DFPT))
  !             The frequencies are the same for every path (CG or Lanczos, Coulomb
  !             and dielectric): either w_freq, or the crpa_freq_* grid below.
  ! q0div_treatment : Correction scheme for the q->0 divergence of the Coulomb kernel.
  !                   (Default: 'HL' - Hybertsen Louie correction)
  ! crpa_mode : Calculation mode of the Uijkl(R) elements.
  !             (Default: 'full' (all possible elements within the given distance thresholds)
  ! crpa_epsil : Logical for the evaluation of the complex macroscopic dielectric tensor eps_ij(w)
  !             (Default: 'false')
  ! crpa_done_epsil : Logical that checks if the dielectric tensor calculation is completed
  !             (Default: 'false')
  ! crpa_diel_intraband : For metals (lgauss .or. ltetra in the underlying pw.x run) only:
  !             whether the intraband/Drude contribution -omega_p^2/(w*(w+i*eta)) is added
  !             to the head of epsilon_w(_diel), on top of the interband CG/Lanczos result
  !             (see crpa_plasma_frequency.f90). Set to .false. to isolate the interband-only
  !             result for comparison. No effect for insulators. (Default: .true.)
  ! crpa_diel_flat : For metals (Gaussian smearing or tetrahedra) only: build the electric-field
  !             auxiliary function with the occupation-free projector (single metallic
  !             weight in eps, instead of W^3) and compute the plasma frequency as the trace
  !             over degenerate blocks. No effect for insulators. (Default: .false.)
  ! crpa_diel_deg_thr : Degeneracy threshold (Ry) used by crpa_diel_flat: pairs closer
  !             than this are moved from the interband part to the plasma frequency.
  !             (Default: 1.0d-5)
  ! crpa_diel_llanczos : If .true. epsilon_w_diel(3,3,crpa_nfreq) is computed by
  !             building a non-Hermitian Lanczos chain instead of solving the Sternheimer
  !             equation self-consistently at every frequency. Evaluates the complex
  !             dielectric tensor at all the frequencies (w_freq or the crpa_freq_*
  !             grid) from a single chain per polarization. (Default: 'false')
  ! crpa_diel_lanczos_steps : Number of steps of the Lanczos chain. The
  !             non-Hermitian chain converges more slowly than the earlier
  !             pseudo-Hermitian shortcut (self-consistent interaction is fed
  !             back every step through alternating q/p components, rather
  !             than every other step through a single vector), so this
  !             typically needs several hundred steps. (Default: 800)
  !
  !             IMPORTANT: unlike a plane-wave cutoff or k-point grid, there is
  !             no universally "safe" default here: convergence depends on the
  !             system and, especially, on how small crpa_freq_broadening is
  !             (resolving the response very close to the real axis requires
  !             finely resolving near-degenerate low-energy transitions, which
  !             can need MANY more steps). ALWAYS check convergence by
  !             increasing crpa_diel_lanczos_steps and confirming epsilon_w_diel
  !             stops changing; a chain that looks perfectly well-behaved
  !             (bounded, smooth beta/gamma) can still be far from converged.
  !             The practical way to afford deep convergence is
  !             crpa_diel_extrapolation: run a moderate real chain, then
  !             extrapolate it much further (see below) rather than paying for
  !             thousands of real (expensive) Lanczos steps.
  ! crpa_diel_lanczos_restart_step : Interval (in Lanczos steps) for writing
  !             intermediate restart points of the dielectric chain. 0 means
  !             none. (Default: 0)
  ! crpa_diel_lanczos_save_final : If .true., a checkpoint is written at the end
  !             of the dielectric chain, so that a finished chain can later be
  !             extended. (Default: .true.)
  ! crpa_diel_lanczos_restart : If .true., resume the dielectric chain from its
  !             checkpoint: an interrupted chain continues, a finished one is
  !             extended or truncated to crpa_diel_lanczos_steps. Needs the same
  !             input and number of MPI processes and pools. (Default: .false.)
  ! crpa_nfreq : Number of points of the frequency grid. 0 means no grid: the single
  !             frequency w_freq is used, and setting both is an error. With
  !             crpa_nfreq >= 1, w_i = crpa_freq_start + (i-1)*step + i*crpa_freq_broadening,
  !             step = (crpa_freq_end - crpa_freq_start)/(crpa_nfreq - 1). (Default: 0)
  ! crpa_freq_start, crpa_freq_end : Start and end (Ry) of the real-frequency grid.
  !             (Default: 0.d0, 0.d0)
  ! crpa_freq_broadening : Imaginary broadening eta (Ry) added to every grid point.
  !             The Lanczos chains need Im(w) > 0. (Default: 0.01d0)
  ! crpa_diel_extrapolation : 'no' (default) uses the beta/gamma Lanczos chain
  !             as computed. 'constant' or 'oscillate' extrapolate beta/gamma
  !             past crpa_diel_lanczos_steps out to crpa_diel_lanczos_steps_ext
  !             (assuming they have settled into a stable, possibly two-periodic
  !             pattern), reducing the number of expensive chain steps needed
  !             for a converged continued fraction. See calc_chi.f90/extrapolate()
  !             in thermo_pw for the same technique.
  !             Validated on a small test case: extrapolating an under-converged
  !             800-step chain out to 2000 steps reproduced (to ~1e-4) the result
  !             of an actual 3000-step chain, at a small fraction of the cost.
  !             Only extrapolates a pattern already present in the real chain,
  !             though: it cannot substitute for a chain that is still actively
  !             evolving (not yet showing a stable pattern) at
  !             crpa_diel_lanczos_steps.
  ! crpa_diel_lanczos_steps_ext : Total chain length used in the continued-fraction
  !             solve when crpa_diel_extrapolation /= 'no'. Must be >=
  !             crpa_diel_lanczos_steps. (Default: crpa_diel_lanczos_steps)
  ! crpa_llanczos : If .true. the screened Coulomb matrix elements are computed by
  !             building one non-Hermitian Lanczos chain per bare perturbation
  !             (crpa_do_lanczos_coulomb.f90) instead of solving the Sternheimer
  !             equation self-consistently for every perturbation
  !             (crpa_solve_linear_system.f90). Evaluates the screened Coulomb
  !             matrix at all the frequencies (w_freq or the crpa_freq_* grid),
  !             from a single chain per bare perturbation. Currently only implemented for
  !             pert_basis = 'wannier'. (Default: 'false')
  ! crpa_coul_lanczos_steps : Number of steps of the Lanczos chain used to compute
  !             the screened Coulomb matrix elements when crpa_llanczos = .true.
  !             (Default: 800)
  ! crpa_coul_extrapolation : Same scheme as crpa_diel_extrapolation ('no',
  !             'constant', 'oscillate'), applied independently to each of the
  !             npert_tot screened Coulomb Lanczos chains. (Default: 'no')
  ! crpa_coul_lanczos_steps_ext : Total chain length used in the continued-fraction
  !             solve when crpa_coul_extrapolation /= 'no'. Must be >=
  !             crpa_coul_lanczos_steps. (Default: crpa_coul_lanczos_steps)
  ! crpa_coul_restart : If .true., reuse the screened-Coulomb results of a
  !             previous run with the same outdir/filU. CG: a completed q point
  !             (marker filU.scrd#iq.done) is skipped. Lanczos: a q point whose
  !             filU.coul_lanczos_chain#iq holds enough steps is rebuilt from it;
  !             otherwise completed chains are reused and chains with a vector
  !             checkpoint (crpa_coul_lanczos_save_final/_restart_step) resumed or
  !             extended. Needs the same input and, for the vector checkpoints,
  !             number of MPI processes and pools. (Default: .false.)
  ! crpa_coul_lanczos_save_final : If .true., the Lanczos vectors are written at the
  !             end of every screened-Coulomb chain, so that the chains can later be
  !             extended. Costs disk: 4 wavefunction-sized arrays per perturbation
  !             and q point. (Default: .false.)
  ! crpa_coul_lanczos_restart_step : Interval (in Lanczos steps) for writing
  !             intermediate vector checkpoints of the screened-Coulomb chains.
  !             0 means none. (Default: 0)
  ! crpa_spectrum_only : If .true., skip SCF/NSCF/Wannier setup and the q-point
  !             loop entirely, and instead read back previously saved raw
  !             Lanczos chain coefficients (diel_lanczos_chain.dat if
  !             crpa_epsil, filU.coul_lanczos_chain#iq if crpa_llanczos) and
  !             just redo the frequency-dependent continued-fraction solve, at
  !             whatever frequencies (w_freq or crpa_freq_*) and
  !             extrapolation scheme are given in the CURRENT input (which may
  !             differ from the run that produced the saved chains). Requires
  !             crpa_epsil and/or crpa_llanczos. (Default: 'false')
  ! crpa_w_rspace : If .true., also write W_ij,kl(R1,R2)(R3) = sum_q w_q exp(-iq.R3) W_q / sum_q w_q,
  !             the element of w_i(0) w_j(R1) with w_k(R3) w_l(R3+R2), for every R3 of
  !             the Wigner-Seitz supercell of the q grid, to filU.dir-WR (crpa_average_w.f90).
  !             Needs the full q grid: not compatible with crpa_reduce_q or qplot.
  !             Only pert_basis = 'wannier'. (Default: 'false')
  !
  NAMELIST / INPUTCRPA / prefix, outdir, lrpa, lnoloc, niter_max, alpha_mix, nq1, nq2, nq3, &
                         nmix, start_q, last_q, thresh_init, tr2, &
                         active_bands_min, active_bands_max, diagonalization, iverbosity, &
                         conv_thr_nscf, reduce_io, wannier_seedname, pert_basis, filU, &
                         qplot, dist_thr, dist_thr_large, cdfpt, cdfpt_active_space, &
                         cdfpt_bands_min, cdfpt_bands_max, folder_wan_Rr, &
                         write_wan_Rr, w_freq, q0div_treatment, crpa_mode, crpa_epsil, &
                         crpa_done_epsil, crpa_epsil_only, crpa_diel_llanczos, &
                         crpa_diel_lanczos_steps, crpa_diel_lanczos_restart_step, &
                         crpa_diel_lanczos_restart, crpa_diel_lanczos_save_final, &
                         crpa_freq_start, crpa_freq_end, crpa_nfreq, &
                         crpa_freq_broadening, crpa_diel_extrapolation, crpa_diel_lanczos_steps_ext, &
                         crpa_llanczos, crpa_coul_lanczos_steps, &
                         crpa_coul_extrapolation, crpa_coul_lanczos_steps_ext, &
                         crpa_coul_restart, crpa_coul_lanczos_save_final, &
                         crpa_coul_lanczos_restart_step, &
                         crpa_spectrum_only, crpa_reduce_q, crpa_w_rspace, crpa_diel_intraband, &
                         crpa_diel_flat, crpa_diel_deg_thr, &
                         crpa_xcmode
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
  crpa_reduce_q      = .FALSE.
  crpa_w_rspace      = .FALSE.
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
  lnoloc             = .FALSE.
  dist_thr           = 6.D-4
  dist_thr_large     = -999.0
  w_freq             = (0.d0, 0.d0)
  q0div_treatment    = 'HL'
  crpa_mode          = 'full'
  crpa_epsil         = .false.
  crpa_done_epsil    = .false.
  crpa_epsil_only    = .false.
  crpa_diel_intraband = .true.
  crpa_diel_flat      = .false.
  crpa_diel_deg_thr   = 1.0d-5
  crpa_diel_llanczos          = .false.
  crpa_diel_lanczos_steps     = 800
  crpa_diel_lanczos_restart_step = 0
  crpa_diel_lanczos_restart      = .false.
  crpa_diel_lanczos_save_final   = .true.
  crpa_freq_start             = 0.d0
  crpa_freq_end               = 0.d0
  crpa_nfreq                  = 0
  crpa_freq_broadening        = 0.01d0
  crpa_diel_extrapolation     = 'no'
  ! -1 = "not set by the user": resolved to crpa_diel_lanczos_steps (whatever
  ! value that ends up with, AFTER the namelist below is read) further down,
  ! since a simple "= crpa_diel_lanczos_steps" assignment here would only
  ! ever capture today's default (800), not a value the user sets below.
  crpa_diel_lanczos_steps_ext = -1
  crpa_llanczos               = .false.
  crpa_coul_lanczos_steps     = 800
  crpa_coul_extrapolation     = 'no'
  ! -1 = "not set by the user", resolved the same way as
  ! crpa_diel_lanczos_steps_ext above.
  crpa_coul_lanczos_steps_ext = -1
  crpa_coul_restart              = .false.
  crpa_coul_lanczos_save_final   = .false.
  crpa_coul_lanczos_restart_step = 0
  crpa_spectrum_only          = .false.
  crpa_xcmode                 = 'TCTC'
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
    !
    ! crpa_xcmode: kernel on the external legs of the screened interaction W
    ! (see crpacom.f90). 'TETE' puts f_xc there and only makes sense when the
    ! screening also carries f_xc, i.e. lrpa = .false.
    !
    IF (TRIM(crpa_xcmode) /= 'TCTC' .AND. TRIM(crpa_xcmode) /= 'TETE') &
       CALL errore('crpa_readin', "crpa_xcmode must be 'TCTC' or 'TETE'", 1)
    IF (TRIM(crpa_xcmode) == 'TETE' .AND. lrpa) &
       CALL errore('crpa_readin', 'crpa_xcmode = TETE requires lrpa = .false.', 1)
    !
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
  eth_rps = 1.0D-9
  fildrho = ' '
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
  tmp_dir_ph = trimcheck( TRIM (tmp_dir) // '_epsw' // int_to_char(my_image_id) )
  CALL create_directory(tmp_dir_ph)
  ext_recover = .FALSE.
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
  ! Check sanity of perturbation basis definition. Not needed in
  ! crpa_epsil_only mode: that mode never builds a screened-interaction
  ! perturbation basis (bands- or Wannier-projected), so it should not
  ! require active_bands_min/max or wannier_seedname to be set.
  !
  IF (.NOT. crpa_epsil_only) THEN
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
  ENDIF
  !
  ! Validate epsilon-only mode
  !
  IF (crpa_epsil_only .AND. .NOT. crpa_epsil) &
     CALL errore('crpa_readin', &
        'crpa_epsil_only = .true. requires crpa_epsil = .true. in the input', 1)
  !
  ! crpa_diel_llanczos and crpa_llanczos are tied together at a general
  ! level, outside crpa_epsil_only mode: crpa_llanczos = .true. always
  ! means crpa_diel_llanczos = .true. too (crpa_q0correction.f90 needs to
  ! be able to call the dielectric-function Lanczos chain and share its
  ! frequency grid), and crpa_llanczos = .false. always means
  ! crpa_diel_llanczos = .false. too. They are only genuinely independent
  ! in crpa_epsil_only mode (which has no crpa_llanczos concept at all:
  ! the q-loop that would use it is never reached).
  !
  IF (.NOT. crpa_epsil_only) crpa_diel_llanczos = crpa_llanczos
  !
  ! Validate the Lanczos-chain dielectric-function calculation
  !
  IF (crpa_diel_llanczos) THEN
     IF (TRIM(crpa_diel_extrapolation) /= 'no' .AND. &
         TRIM(crpa_diel_extrapolation) /= 'constant' .AND. &
         TRIM(crpa_diel_extrapolation) /= 'oscillate') CALL errore('crpa_readin', &
        'crpa_diel_extrapolation must be one of: no, constant, oscillate', 1)
     ! Resolve the "not set by the user" sentinel now that
     ! crpa_diel_lanczos_steps has its final (namelist-read) value.
     IF (crpa_diel_lanczos_steps_ext < 0) &
        crpa_diel_lanczos_steps_ext = crpa_diel_lanczos_steps
     IF (crpa_diel_lanczos_steps_ext < crpa_diel_lanczos_steps) CALL errore('crpa_readin', &
        'crpa_diel_lanczos_steps_ext must be >= crpa_diel_lanczos_steps', 1)
     ! crpa_extrapolate_diel only fills diel_beta/gamma_store_ext beyond
     ! crpa_diel_lanczos_steps when extrapolation is requested: silently
     ! force steps_ext = steps otherwise, to avoid leaving the tail of
     ! diel_beta/gamma_store_ext uninitialized.
     IF (TRIM(crpa_diel_extrapolation) == 'no') &
        crpa_diel_lanczos_steps_ext = crpa_diel_lanczos_steps
  ENDIF
  !
  ! Validate the Lanczos-chain screened-Coulomb-matrix-elements calculation.
  ! Since crpa_diel_llanczos = crpa_llanczos is enforced above, the
  ! crpa_diel_llanczos block already resolved its own steps_ext/extrapolation - only
  ! the genuinely Coulomb-specific settings need checking here.
  !
  IF (crpa_llanczos) THEN
     IF (pert_basis /= 'wannier') CALL errore('crpa_readin', &
        'crpa_llanczos = .true. is currently only implemented for pert_basis = wannier', 1)
     IF (crpa_coul_lanczos_steps < 1) CALL errore('crpa_readin', &
        'crpa_coul_lanczos_steps must be >= 1', 1)
     IF (TRIM(crpa_coul_extrapolation) /= 'no' .AND. &
         TRIM(crpa_coul_extrapolation) /= 'constant' .AND. &
         TRIM(crpa_coul_extrapolation) /= 'oscillate') CALL errore('crpa_readin', &
        'crpa_coul_extrapolation must be one of: no, constant, oscillate', 1)
     ! Resolve the "not set by the user" sentinel now that
     ! crpa_coul_lanczos_steps has its final (namelist-read) value.
     IF (crpa_coul_lanczos_steps_ext < 0) &
        crpa_coul_lanczos_steps_ext = crpa_coul_lanczos_steps
     IF (crpa_coul_lanczos_steps_ext < crpa_coul_lanczos_steps) CALL errore('crpa_readin', &
        'crpa_coul_lanczos_steps_ext must be >= crpa_coul_lanczos_steps', 1)
     ! crpa_extrapolate_coul only fills coul_beta/gamma_store_ext beyond
     ! crpa_coul_lanczos_steps when extrapolation is requested: silently
     ! force steps_ext = steps otherwise, to avoid leaving the tail of
     ! coul_beta/gamma_store_ext uninitialized.
     IF (TRIM(crpa_coul_extrapolation) == 'no') &
        crpa_coul_lanczos_steps_ext = crpa_coul_lanczos_steps
  ENDIF
  !
  ! Frequencies of every path (CG or Lanczos, Coulomb and dielectric):
  ! w_freq, or the crpa_freq_* grid if crpa_nfreq >= 1
  !
  IF (crpa_nfreq < 0) CALL errore('crpa_readin', 'crpa_nfreq must be >= 0', 1)
  IF (crpa_nfreq > 0) THEN
     IF (ABS(w_freq) > 1.0d-10) CALL errore('crpa_readin', &
        'set either w_freq or the crpa_freq_* grid (crpa_nfreq >= 1), not both', 1)
     IF (crpa_nfreq > 1 .AND. crpa_freq_end <= crpa_freq_start) CALL errore('crpa_readin', &
        'crpa_freq_end must be > crpa_freq_start when crpa_nfreq > 1', 1)
     IF (crpa_freq_broadening < 0.d0) CALL errore('crpa_readin', &
        'crpa_freq_broadening must be >= 0', 1)
     IF (crpa_nfreq > 1 .AND. .NOT. crpa_epsil_only .AND. pert_basis /= 'wannier') &
        CALL errore('crpa_readin', 'a frequency grid is only implemented for pert_basis = wannier', 1)
     CALL crpa_build_freq_grid(crpa_nfreq, crpa_freq_start, crpa_freq_end, &
                               crpa_freq_broadening, crpa_freq)
  ELSE
     crpa_nfreq = 1
     CALL crpa_build_freq_grid(1, DBLE(w_freq), DBLE(w_freq), AIMAG(w_freq), crpa_freq)
  ENDIF
  IF ((crpa_llanczos .OR. crpa_diel_llanczos) .AND. ANY(AIMAG(crpa_freq) <= 0.d0)) &
     CALL errore('crpa_readin', 'the Lanczos chains need a broadening > 0: Im(w_freq) &
     &or crpa_freq_broadening', 1)
  current_w = crpa_freq(1)
  !
  ! Validate crpa_spectrum_only: needs at least one saved chain to read back
  ! (crpa_epsil for diel_lanczos_chain.dat, crpa_llanczos for
  ! filU.coul_lanczos_chain#iq). Compatible with crpa_epsil_only: crpa_main.f90
  ! branches to crpa_run_spectrum_only() before crpa_setup_wannier() is ever
  ! called, so crpa_epsil_only's value has no bearing on this mode - the
  ! combination just means "reconstruct epsilon(w) only, from a chain that
  ! was itself produced by a crpa_epsil_only run" (crpa_llanczos, and hence
  ! any Coulomb/U(w) reconstruction, stays under the user's control via the
  ! pert_basis = 'wannier' requirement below, exactly as outside
  ! crpa_spectrum_only).
  !
  IF (crpa_spectrum_only) THEN
     IF (.NOT. (crpa_epsil .OR. crpa_llanczos)) CALL errore('crpa_readin', &
        'crpa_spectrum_only = .true. requires crpa_epsil and/or crpa_llanczos', 1)
     IF (crpa_diel_lanczos_restart) CALL errore('crpa_readin', &
        'crpa_diel_lanczos_restart continues a chain: not compatible with crpa_spectrum_only', 1)
     IF (crpa_coul_restart) CALL errore('crpa_readin', &
        'crpa_coul_restart: not compatible with crpa_spectrum_only', 1)
  ENDIF
  IF (crpa_diel_lanczos_restart_step < 0) CALL errore('crpa_readin', &
     'crpa_diel_lanczos_restart_step must be >= 0', 1)
  !
  ! The R3 transform needs W_q on every point of the q grid
  !
  IF (crpa_w_rspace) THEN
     IF (crpa_reduce_q) CALL errore('crpa_readin', 'crpa_w_rspace needs the full q grid: &
        &set crpa_reduce_q = .false. (SCF and NSCF with nosym and noinv)', 1)
     IF (qplot) CALL errore('crpa_readin', 'crpa_w_rspace needs the full q grid: not compatible with qplot', 1)
     IF (pert_basis /= 'wannier') CALL errore('crpa_readin', &
        'crpa_w_rspace is only implemented for pert_basis = wannier', 1)
  ENDIF
  IF (crpa_coul_lanczos_restart_step < 0) CALL errore('crpa_readin', &
     'crpa_coul_lanczos_restart_step must be >= 0', 1)
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
  IF (TRIM(crpa_mode) /= 'full' .AND. TRIM(crpa_mode) /= 'dHP' .AND. TRIM(crpa_mode) /= 'debug') THEN
     CALL errore('crpa_readin', 'Invalid crpa_mode. Only full and dHP are currently allowed.', 1)
  END IF
  !
  RETURN
  !
END SUBROUTINE input_sanity
  !
END SUBROUTINE crpa_readin
