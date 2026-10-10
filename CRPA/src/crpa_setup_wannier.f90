!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
SUBROUTINE crpa_setup_wannier()
   !----------------------------------------------------------------------------
   !
   ! This routine prepares the perturbations.
   ! The preparation specific to each q point is done in
   ! crpa_prepare_perturbations_q.
   !
   !----------------------------------------------------------------------------
   !
   USE io_global,         ONLY : stdout
   USE io_files,          ONLY : tmp_dir
   USE w90_interface,     ONLY : w90_read
   USE w90_real_space,    ONLY : w90_write_real_space
   USE w90_wan_Rr_buffer, ONLY : w90_wan_Rr_save_buffer
   USE crpacom,           ONLY : wannier_seedname, tmp_dir_save, folder_wan_Rr, write_wan_Rr
   USE crpa_pert,         ONLY : pert_basis
   !
   IMPLICIT NONE
   !
   IF (pert_basis == 'wannier') THEN
      !
      WRITE (stdout, '(/5x, a)') "========================================"
      !
      WRITE (stdout, '(5x, a)') 'Reading Wannier90 checkpoint file'
      CALL w90_read(wannier_seedname)
      !
      IF (write_wan_Rr) THEN
         WRITE (stdout, '(5x, a)') 'Writing real-space Wannier functions in collected form'
         CALL w90_write_real_space(tmp_dir_save, folder_wan_Rr)
      ELSE
         WRITE (stdout, '(5x, a)') 'Using collected real-space Wannier functions in ' // TRIM(folder_wan_Rr)
      ENDIF
      !
      WRITE (stdout, '(5x, a)') 'Loading real-space Wannier functions to buffer in distributed form'
      CALL w90_wan_Rr_save_buffer(folder_wan_Rr, tmp_dir)
      !
      WRITE (stdout, '(5x, a)') "========================================"
      !
   ENDIF
   !
   !----------------------------------------------------------------------------
END SUBROUTINE crpa_setup_wannier
!------------------------------------------------------------------------------
