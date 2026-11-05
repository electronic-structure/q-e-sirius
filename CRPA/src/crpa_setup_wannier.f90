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
   USE io_global,         ONLY : stdout, meta_ionode
   USE io_files,          ONLY : tmp_dir
   USE w90_interface,     ONLY : w90_read, num_wann, wann_centers
   USE w90_real_space,    ONLY : w90_write_real_space
   USE w90_wan_Rr_buffer, ONLY : w90_wan_Rr_save_buffer
   USE crpacom,           ONLY : wannier_seedname, tmp_dir_save, folder_wan_Rr, write_wan_Rr, filU
   USE crpa_pert,         ONLY : pert_basis
   !
   IMPLICIT NONE
   !
   INTEGER :: iun, ios, iw
   !
   IF (pert_basis == 'wannier') THEN
      !
      WRITE (stdout, '(/5x, a)') "========================================"
      !
      WRITE (stdout, '(5x, a)') 'Reading Wannier90 checkpoint file'
      CALL w90_read(wannier_seedname)
      !
      ! Wannier centres for crpa_average_w (WS folding of R3), which also
      ! runs in crpa_spectrum_only, where this routine is not called
      !
      IF (meta_ionode) THEN
         OPEN(NEWUNIT=iun, FILE=TRIM(filU)//'.wcentres', STATUS='unknown', IOSTAT=ios)
         IF (ios /= 0) CALL errore('crpa_setup_wannier', 'opening '//TRIM(filU)//'.wcentres', ABS(ios))
         WRITE(iun, '("# Wannier centres, Cartesian, alat units")')
         WRITE(iun, '(i8)') num_wann
         DO iw = 1, num_wann
            WRITE(iun, '(i8,3e24.15)') iw, wann_centers(:, iw)
         ENDDO
         CLOSE(iun, STATUS='keep')
      ENDIF
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
