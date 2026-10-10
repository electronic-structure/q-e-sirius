!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE elph_bcast_input ( )
  !-----------------------------------------------------------------------
  !
  ! In this routine the first processor sends the input parameters to all
  ! the other processors
  !
  USE mp,               ONLY : mp_bcast
  USE mp_world,         ONLY : world_comm
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : meta_ionode_id
  USE control_lr,       ONLY : reduce_io
  USE elphcom,          ONLY : wannier_seedname, qplot, folder_elph, folder_wan_Rr, &
                               folder_phonon, fildvscf, write_wan_Rr
  !
  IMPLICIT NONE
  !
  CALL mp_bcast(prefix, meta_ionode_id, world_comm)
  CALL mp_bcast(tmp_dir, meta_ionode_id, world_comm)
  !
  CALL mp_bcast(reduce_io, meta_ionode_id, world_comm)
  CALL mp_bcast(qplot, meta_ionode_id, world_comm)
  CALL mp_bcast(write_wan_Rr, meta_ionode_id, world_comm)
  !
  CALL mp_bcast(folder_wan_Rr, meta_ionode_id, world_comm)
  CALL mp_bcast(folder_elph, meta_ionode_id, world_comm)
  CALL mp_bcast(wannier_seedname, meta_ionode_id, world_comm)
  CALL mp_bcast(folder_phonon, meta_ionode_id, world_comm)
  CALL mp_bcast(fildvscf, meta_ionode_id, world_comm)
  !
END SUBROUTINE elph_bcast_input
