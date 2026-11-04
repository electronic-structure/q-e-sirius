!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine crpa_init()
  !----------------------------------------------------------------------
  !
  ! Setup various variables for the CRPA calculation.
  !
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE io_global,     ONLY : stdout
  USE lsda_mod,      ONLY : nspin
  USE ldaU,          ONLY : Hubbard_lmax, is_hubbard
  USE qpoint,        ONLY : nq1, nq2, nq3
  USE crpacom,       ONLY : chi0, chi, nqsh, code
                           !  determine_num_pert_only
                           !  perturbed_atom, todo_atom
  USE ph_restart,    ONLY : ph_readfile
  USE control_lr,    ONLY : rec_code_read
  !
  IMPLICIT NONE
  INTEGER :: na, nt, ierr
  !
!   ALLOCATE (todo_atom(nat))
!   ALLOCATE (perturbed_atom(nat))
  !
  ! Determine the number of q-points in the q-mesh without symmetry
  !
  nqsh = nq1*nq2*nq3
  !
  ! Try to read a previous save file (written by crpa_phescf_w after the
  ! dielectric calculation completes). If found, rec_code_read is updated
  ! and crpa_phescf_w will skip the E-field calculation on restart.
  !
  CALL ph_readfile('status_ph', 0, 0, ierr)
  IF (ierr == 0) THEN
     WRITE(stdout, '(/,5x,"Restart file found. Resuming from rec_code = ", i5)') rec_code_read
  ENDIF
!   !
!   ! Determine the total number of Hubbard atoms in the virtual supercell
!   !
!   nath_sc = nath * nqsh
!   !
!   ! Find inequivalent sites
!   !
!   CALL crpa_find_inequiv_sites()
!   !
!   IF (.NOT.determine_num_pert_only) THEN
!      !
!      ! The first dimension of chi0 and chi runs over all possible
!      ! real+virtual atoms (nath_sc), whereas the second dimension runs over
!      ! real atoms in the primitive cell which can be perturbed (nat).
!      !
!      ALLOCATE (chi0(nath_sc, nat))
!      ALLOCATE (chi(nath_sc, nat))
!      chi0(:,:) = 0.0d0
!      chi(:,:) = 0.0d0
!      !
!   ENDIF
  !
  RETURN
  !
end subroutine crpa_init
