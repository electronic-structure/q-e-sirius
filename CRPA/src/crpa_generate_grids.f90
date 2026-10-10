!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_generate_grids()
   !-----------------------------------------------------------------------
   !
   ! This routine generates the q-points and R-points grids.
   ! comp_iq = .TRUE. if this q point is calculated in this run, .FALSE. otherwise
   !
   USE mp_images,       ONLY : inter_image_comm
   USE io_files,        ONLY : tmp_dir
   USE ions_base,       ONLY : nat
   USE lsda_mod,        ONLY : nspin
   USE scf,             ONLY : rho
   USE io_rho_xml,      ONLY : write_scf
   USE qpoint,          ONLY : start_q, last_q, nqs, comp_iq
   USE lr_symm_base,    ONLY : rtau
   USE ldaU,            ONLY : Hubbard_lmax
   USE ldaU_lr,         ONLY : dnsscf
   USE crpacom,         ONLY : dns0, dnsscf_tot, dns0_tot, nqsh, tmp_dir_crpa
   !
   IMPLICIT NONE
   !
   INTEGER :: iq
   !! counter over q points
   INTEGER :: start_q_loc, last_q_loc
   !! First and last q point to compute in this image
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   !
   tmp_dir = tmp_dir_crpa
   !
   ! Write the ground-state density
   !
   CALL write_scf( rho, nspin )
   !
   ! Copy the k point grid to xk_orig
   !
   CALL crpa_copy_original_xk()
   !
   ! Generate the q-points grid
   !
   CALL crpa_q_points()
   !
   IF ( last_q < 1 .OR. last_q > nqs ) last_q = nqs
   !
   ! Distribute the q points
   !
   CALL divide(inter_image_comm, last_q - start_q + 1, start_q_loc, last_q_loc)
   start_q_loc = start_q_loc + start_q - 1
   last_q_loc = last_q_loc + start_q - 1
   !
   ! Set the q points for which the calculation must be performed.
   !
   ALLOCATE(comp_iq(nqs))
   comp_iq(:) = .FALSE.
   !
   DO iq = start_q_loc, last_q_loc
      comp_iq(iq) = .TRUE.
   ENDDO
   !
   ! Generate the R-points grid
   !
   CALL crpa_R_points()
   !
   ! Allocate some arrays
   !
   ALLOCATE (dnsscf(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, nqs))
   ALLOCATE (dns0  (2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, nqs))
   dnsscf = (0.0d0, 0.0d0)
   dns0   = (0.0d0, 0.0d0)
   !
   ALLOCATE (dnsscf_tot(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, nqsh))
   ALLOCATE (dns0_tot  (2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, nqsh))
   dnsscf_tot = (0.0d0, 0.0d0)
   dns0_tot   = (0.0d0, 0.0d0)
   !
   ALLOCATE (rtau(3,48,nat))
   !
   RETURN
   !
END SUBROUTINE crpa_generate_grids
