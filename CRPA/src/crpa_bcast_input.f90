!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_bcast_input ( )
  !-----------------------------------------------------------------------
  !
  ! In this routine the first processor sends the input parameters to all
  ! the other processors
  !
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE io_files,         ONLY : tmp_dir, prefix
  USE control_flags,    ONLY : iverbosity, isolve
  USE check_stop,       ONLY : max_seconds
  USE io_global,        ONLY : meta_ionode_id
  USE control_lr,       ONLY : lrpa, ethr_nscf, reduce_io, alpha_mix, nmix_ph, niter_ph,     &
                               tr2_ph, conv_thr_nscf, thresh_init
  USE qpoint,           ONLY : nq1, nq2, nq3, start_q, last_q
  USE constrained_dfpt, ONLY : cdfpt, cdfpt_active_space, cdfpt_bands_min, cdfpt_bands_max
  USE crpacom,          ONLY : filU, wannier_seedname, dist_thr, dist_thr_large, &
                               active_bands_min, active_bands_max, &
                               folder_wan_Rr, write_wan_Rr, w_freq, q0div_treatment, crpa_mode
  USE crpa_pert,        ONLY : pert_basis
  USE crpa_qpoints,     ONLY : qplot
  !
  IMPLICIT NONE
  !
  ! Logicals
  !
  CALL mp_bcast (lrpa, meta_ionode_id, world_comm)
  CALL mp_bcast (reduce_io, meta_ionode_id, world_comm)
  CALL mp_bcast (qplot, meta_ionode_id, world_comm)
  CALL mp_bcast (write_wan_Rr, meta_ionode_id, world_comm)
  !
  ! Integers
  !
  CALL mp_bcast (nq1, meta_ionode_id, world_comm)
  CALL mp_bcast (nq2, meta_ionode_id, world_comm)
  CALL mp_bcast (nq3, meta_ionode_id, world_comm)
  CALL mp_bcast (start_q, meta_ionode_id, world_comm)
  CALL mp_bcast (last_q, meta_ionode_id, world_comm)
  CALL mp_bcast (iverbosity, meta_ionode_id, world_comm)
  CALL mp_bcast (isolve, meta_ionode_id, world_comm)
  CALL mp_bcast (niter_ph, meta_ionode_id, world_comm)
  CALL mp_bcast (nmix_ph, meta_ionode_id, world_comm)
  !
  ! Real*8
  !
  CALL mp_bcast (thresh_init, meta_ionode_id, world_comm)
  CALL mp_bcast (conv_thr_nscf, meta_ionode_id, world_comm)
  CALL mp_bcast (alpha_mix, meta_ionode_id, world_comm)
  CALL mp_bcast (tr2_ph, meta_ionode_id, world_comm)
  CALL mp_bcast (dist_thr, meta_ionode_id, world_comm)
  CALL mp_bcast (dist_thr_large, meta_ionode_id, world_comm)
  !
  ! Complex
  !
  CALL mp_bcast (w_freq, meta_ionode_id, world_comm)
  !
  ! Characters
  !
  CALL mp_bcast (prefix, meta_ionode_id, world_comm)
  CALL mp_bcast (tmp_dir, meta_ionode_id, world_comm)
  CALL mp_bcast (filU, meta_ionode_id, world_comm)
  CALL mp_bcast (folder_wan_Rr, meta_ionode_id, world_comm)
  CALL mp_bcast (wannier_seedname, meta_ionode_id, world_comm)
  CALL mp_bcast (q0div_treatment, meta_ionode_id, world_comm)
  CALL mp_bcast (crpa_mode, meta_ionode_id, world_comm)
  !
  ! Constrained DFPT parameters (from LR_Modules)
  CALL mp_bcast (cdfpt,              meta_ionode_id, world_comm)
  CALL mp_bcast (cdfpt_active_space, meta_ionode_id, world_comm)
  CALL mp_bcast (cdfpt_bands_min,    meta_ionode_id, world_comm)
  CALL mp_bcast (cdfpt_bands_max,    meta_ionode_id, world_comm)
  !
  ! CRPA perturbation basis parameters (from crpacom)
  CALL mp_bcast (pert_basis,       meta_ionode_id, world_comm)
  CALL mp_bcast (active_bands_min, meta_ionode_id, world_comm)
  CALL mp_bcast (active_bands_max, meta_ionode_id, world_comm)
  !
  RETURN
  !
END SUBROUTINE crpa_bcast_input
