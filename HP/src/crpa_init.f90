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
  ! Setup various variables for the cRPA calculation.
  !
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE io_global,     ONLY : stdout
  USE lsda_mod,      ONLY : nspin
  USE ldaU,          ONLY : Hubbard_lmax, is_hubbard, lda_plus_u
  USE qpoint,        ONLY : nq1, nq2, nq3
  USE crpacom,       ONLY : chi0, chi, nqsh
  USE ldaU_hp,       ONLY : code
  USE uspp,          ONLY : okvan
  USE paw_variables, ONLY : okpaw
  USE symm_base,     ONLY : nosym
  USE control_flags, ONLY : noinv
  USE noncollin_module, ONLY : noncolin
                           !  determine_num_pert_only
                           !  perturbed_atom, todo_atom
  !
  IMPLICIT NONE
  INTEGER :: na, nt
  !
  ! cRPA is implemented only for norm-conserving pseudopotentials.
  !
  IF (okvan) CALL errore('crpa_init', &
     'cRPA is not implemented for ultrasoft pseudopotentials', 1)
  IF (okpaw) CALL errore('crpa_init', &
     'cRPA is not implemented for PAW pseudopotentials', 1)
  IF (lda_plus_u) CALL errore('crpa_init', &
     'cRPA is not implemented for DFT+Hubbard', 1)
  IF (noncolin) CALL errore('crpa_init', &
     'cRPA is not implemented for noncollinear or spin-orbit calculations', 1)
  IF (nspin == 2) CALL errore('crpa_init', &
     'cRPA is not implemented for spin-polarized (LSDA) calculations', 1)
  !
  ! The response is symmetrized assuming a totally symmetric perturbation, which a
  ! Wannier/band pair density is not. Require the pw.x run to disable symmetry.
  !
  IF (.NOT. nosym) CALL errore('crpa_init', &
     'cRPA requires nosym = .true. in the pw.x calculation', 1)
  IF (.NOT. noinv) CALL errore('crpa_init', &
     'cRPA requires noinv = .true. in the pw.x calculation', 1)
  !
!   ALLOCATE (todo_atom(nat))
!   ALLOCATE (perturbed_atom(nat))
  !
  ! Determine the number of q-points in the q-mesh without symmetry
  !
  nqsh = nq1*nq2*nq3
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
