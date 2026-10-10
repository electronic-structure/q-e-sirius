!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE crpa_openfil_q()
  !--------------------------------------------------------------------------
  !
  ! This subroutine opens all necessary files necessary.
  !
  USE io_files,         ONLY : prefix, tmp_dir, iunhub, nwordwfcU
  USE control_flags,    ONLY : io_level
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : npol, noncolin, domag
  USE buffers,          ONLY : open_buffer
  USE qpoint,           ONLY : nksq
  USE control_lr,       ONLY : lgamma
  USE units_lr,         ONLY : iuwfc, lrwfc, iuatswfc, iudwf, lrdwf, iudvwfc, lrdvwfc
  USE ldaU,             ONLY : nwfcU, lda_plus_u
  USE crpacom,          ONLY : recalc_sym, tmp_dir_save, tmp_dir_crpa, iuwan, lrwan
  USE crpa_pert,        ONLY : pert_basis
  USE w90_interface,    ONLY : num_wann
  !
  IMPLICIT NONE
  LOGICAL :: exst, exst_mem
  ! logical variable to check file exists
  ! logical variable to check file exists in memory
  !
  IF (LEN_TRIM(prefix) == 0) CALL errore ('crpa_openfil_q', 'wrong prefix', 1)
  !
  IF (lgamma .AND. .NOT.recalc_sym.and. .not.(noncolin.and.domag)) THEN
     tmp_dir = tmp_dir_save
  ELSE
     tmp_dir = tmp_dir_crpa
  ENDIF
  !
  ! Open a file to read the unperturbed KS wavefunctions
  !
  iuwfc = 20
  lrwfc = nbnd * npwx * npol
  IF (io_level > 0) THEN
     CALL open_buffer (iuwfc, 'wfc', lrwfc, io_level, exst_mem, exst, tmp_dir)
     IF (.NOT.exst .AND. .NOT.exst_mem) &
        CALL errore ('crpa_openfil_q', 'file '//trim(prefix)//'.wfc not found', 1)
  ELSE
     iuwfc = 10
  ENDIF
  !
  ! From now on all files are written to tmp_dit_crpa
  !
  tmp_dir = tmp_dir_crpa
  !
  ! Open a file to write/read deltaV_{bare} * psi
  ! (i.e. perturbing potential times the unperterbed wfct)
  !
  iudvwfc = 21
  lrdvwfc = nbnd * npwx * npol
  CALL open_buffer (iudvwfc, 'dvwfc', lrdvwfc, io_level, exst_mem, exst, tmp_dir)
  !
  ! Open a file to write/read a solution of the linear system (dpsi)
  !
  iudwf = 22
  lrdwf = nbnd * npwx * npol
  CALL open_buffer (iudwf, 'dwfc', lrdwf, io_level, exst_mem, exst, tmp_dir)
  !
  IF (lda_plus_u) THEN
     !
     ! Open a file to write/read S*phi at k and k+q (atomic wfct's)
     !
     iuatswfc  = 23
     nwordwfcU = npwx * nwfcU * npol
     CALL open_buffer (iuatswfc, 'satwfc', nwordwfcU, io_level, exst_mem, exst, tmp_dir)
     !
  ENDIF
  !
  IF (lgamma) THEN
     !
     ! If q = Gamma, open unit iunhub which contain S*phi at k.
     ! Unit iunhub is used in commutator_Vhubx_psi.f90.
     !
     CALL open_buffer(iunhub, 'hub', nwordwfcU, io_level, exst_mem, exst, tmp_dir)
     !
  ENDIF
  !
  IF (pert_basis == 'wannier') THEN
     iuwan = 24
     lrwan = num_wann * npwx * npol
     CALL open_buffer(iuwan, 'wan', lrwan, io_level, exst_mem, exst, tmp_dir)
  ENDIF
  !
END SUBROUTINE crpa_openfil_q
