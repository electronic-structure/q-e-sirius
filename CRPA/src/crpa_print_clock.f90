!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_print_clock
  !-----------------------------------------------------------------------

  USE io_global,  ONLY : stdout
  USE uspp,       ONLY : okvan
  USE fft_base,   ONLY : dffts
  !
  IMPLICIT NONE
  !
  WRITE( stdout, * )
  CALL print_clock ('init_vloc')
  CALL print_clock ('init_us_1')
  CALL print_clock ('newd')
  CALL print_clock ('add_vuspsi')
  !
  WRITE( stdout, * )
  WRITE( stdout, * )  '    PRINTING TIMING FROM CRPA ROUTINES: '
  WRITE( stdout, * )
  !
  CALL print_clock ('crpa_nscf_rotate')
  CALL print_clock ('crpa_setup_q')
  CALL print_clock ('crpa_init_q')
  CALL print_clock ('crpa_solve_linear_system')
  CALL print_clock ('crpa_dvpsi_pert')
  CALL print_clock ('crpa_dnsq')
  ! CALL print_clock ('crpa_symdnsq')
  ! CALL print_clock ('crpa_dnstot_sum_q')
  ! CALL print_clock ('crpa_rotate_dnsq')
  ! CALL print_clock ('crpa_calc_chi')
  ! CALL print_clock ('crpa_vpsifft')
  CALL print_clock ('ef_shift')
  CALL print_clock ('crpa_run_nscf')
  CALL print_clock ('crpa_active')
  CALL print_clock ('crpa_postproc')
  !
#if defined (__MPI)
  CALL print_clock ('crpa_psymdvscf')
#else
  CALL print_clock ('crpa_symdvscf')
#endif
  !
  WRITE( stdout, * )
  WRITE( stdout, * )  '    PRINTING TIMING FROM Wannier90 INTERFACE ROUTINES: '
  WRITE( stdout, * )
  CALL print_clock("w90_rotate")
  CALL print_clock("w90_interpolate")
  !
  WRITE( stdout, * )
  CALL print_clock("w90_wigner")
  CALL print_clock("w90_readwrite")
  CALL print_clock("w90_mult_iqr")
  CALL print_clock("w90_shift_wfc")
  !
  WRITE( stdout, * )
  WRITE( stdout, * )  '    PRINTING TIMING FROM LR MODULE: '
  WRITE( stdout, * )
  !
  CALL print_clock ('sth_kernel')
  CALL print_clock ('apply_dpot_b')
  CALL print_clock ('ortho')
  CALL print_clock ('cgsolve')
  CALL print_clock ('ch_psi')
  CALL print_clock ('incdrhoscf')
  CALL print_clock ('localdos')
  CALL print_clock ('dv_of_drho')
  CALL print_clock ('mix_pot')
  CALL print_clock ('setup_dgc')
  CALL print_clock ('setup_dmuxc')
  CALL print_clock ('setup_nbnd_occ')
  CALL print_clock ('lr_orthoUwfc')
  CALL print_clock ('lr_matrix_elements')
  IF (dffts%has_task_groups) THEN
     CALL print_clock ('cft_wave_tg')
  ELSE
     CALL print_clock ('cft_wave')
  ENDIF
  !
  IF (okvan) THEN
     WRITE( stdout, * )
     WRITE( stdout, * )  '    USPP ROUTINES: '
     WRITE( stdout, * )
     CALL print_clock ('newdq')
     CALL print_clock ('adddvscf')
     CALL print_clock ('addusdbec')
     CALL print_clock ('addusldos')
     CALL print_clock ('crpa_addusddens')
  ENDIF
  !
  RETURN
  !
END SUBROUTINE crpa_print_clock
