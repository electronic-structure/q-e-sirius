!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_summary
  !-----------------------------------------------------------------------
  !
  ! This routine writes on output a summary of the variables which do not
  ! depend on the q point.
  !
  USE constants,        ONLY : rytoev
  USE ions_base,        ONLY : nat, ityp, atm, tau, ntyp => nsp, amass
  USE io_global,        ONLY : stdout
  USE cell_base,        ONLY : at, bg, ibrav, alat, omega, celldm
  USE gvecs,            ONLY : dual
  USE gvecw,            ONLY : ecutwfc
  USE funct,            ONLY : write_dft_name
  USE control_lr,       ONLY : ethr_nscf
  USE ldaU,             ONLY : is_hubbard, Hubbard_U, lda_plus_u_kind, Hubbard_V, num_uc
  USE crpacom,          ONLY : conv_thr_chi, active_space, active_bands_min, active_bands_max
  USE constrained_dfpt, ONLY : lcdfpt

  IMPLICIT NONE
  !
  INTEGER :: i, ipol, apol, na, nb, nt, isymq, isym, ik, nsymtot, dimn
  ! generic counter
  ! counter on polarizations
  ! counter on polarizations
  ! counter on atoms
  ! counter on atomic types
  ! counter on symmetries
  ! counter on symmetries
  ! counter on k points
  !
  WRITE( stdout, * )
  !
  WRITE( stdout, 100) ibrav, alat, omega, nat, ntyp, &
                      ecutwfc, ecutwfc * dual, ethr_nscf, conv_thr_chi
100 FORMAT (/5x, 'bravais-lattice index     =  ',i12,/,5x, &
               & 'lattice parameter (alat)  =  ',f12.4,' (a.u.)',/,5x, &
               & 'unit-cell volume          =  ',f12.4,' (a.u.)^3',/,5x, &
               & 'number of atoms/cell      =  ',i12,/,5x, &
               & 'number of atomic types    =  ',i12,/,5x, &
               & 'kinetic-energy cut-off    =  ',f12.2,' (Ry)',/,5x, &
               & 'charge density cut-off    =  ',f12.2,' (Ry)',/,5x, &
               & 'conv. thresh. for NSCF    =  ',3x,1pe9.1,/,5x, &
               & 'conv. thresh. for chi     =  ',3x,1pe9.1)
  !
  ! Info about Hubbard parameters
  !
  WRITE (stdout,'(5x,a)') 'Input Hubbard parameters (in eV):'
  IF (lda_plus_u_kind.EQ.0) THEN
     DO nt = 1, ntyp
        IF (is_hubbard(nt)) THEN
           WRITE (stdout,'(7x,a,i2,a,21x,a,1x,1pe12.5)') &
                  & 'U (',nt,')','= ', Hubbard_U(nt)*rytoev
        ENDIF
     ENDDO
  ELSEIF (lda_plus_u_kind.EQ.2) THEN
     ! Number of atoms in the 3x3x3 supercell
     dimn = num_uc * nat
     DO na = 1, nat
        DO nb = 1, dimn
           IF ( ABS(Hubbard_V(na,nb,1)).GE.1.d-15 .OR. nb==na) THEN
               WRITE (stdout,'(7x,a,i3,a,i4,a,2x,a,1x,f7.4)') &
                  & 'V (',na,',',nb,')','= ', Hubbard_V(na,nb,1)*rytoev
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  ! Description of the unit cell
  !
  WRITE( stdout, '(/,2(3x,3(2x,"celldm(",i1,") =",f9.5),/))') (i, celldm(i), i=1,6)
  !
  ! Description of the direct and reciprocal lattice vectors
  !
  WRITE( stdout, '(5x, "crystal axes: (cart. coord. in units of alat)",/, &
                  & 3(15x,"a(",i1,") = (",3f8.4," )  ",/ ) )')  &
                  & (apol, (at(ipol,apol),ipol=1,3), apol=1,3)
  WRITE( stdout, '(5x, "reciprocal axes: (cart. coord. in units 2 pi/alat)",/, &
                  & 3(15x,"b(",i1,") = (",3f8.4," )  ",/ ) )')  &
                  & (apol, (bg(ipol,apol),ipol=1,3), apol=1,3)
  !
  ! Description of the atoms inside the unit cell
  !
  WRITE( stdout, '(5x,"Atoms inside the unit cell (Cartesian axes):")')
  WRITE( stdout, '(5x,"site n.  atom      mass ", "          positions (alat units)")')
  WRITE( stdout, '(6x,i3,3x,a6,3x,f8.4,"   tau(",i3, ") = (",3f9.5,"  )")')  &
         & (na, atm(ityp(na)), amass(ityp(na)), na, (tau(ipol,na), ipol=1,3), na=1,nat)
  !
  ! Print info about the active space
  !
  IF (lcdfpt) THEN
     WRITE(stdout, '()')
     WRITE(stdout, '(5x,"Active space for constrained DFPT")')
     IF (TRIM(active_space) == 'none') THEN
        WRITE(stdout, '(5x,5x,"Active space is not used. Perform ordinary DFPT.")')
     ELSEIF (TRIM(active_space) == 'bands') THEN
        WRITE(stdout, '(5x,5x,"Active space defined by the band index.")')
        WRITE(stdout, '(5x,5x,a,i5,a,i5)') 'Active bands: ', active_bands_min, ' - ', active_bands_max
     ENDIF
     WRITE(stdout, '()')
  ENDIF
  !
  RETURN
  !
END SUBROUTINE crpa_summary
