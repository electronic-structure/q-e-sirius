!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE crpa_postprocess_coulomb(iq)
!------------------------------------------------------------
   !
   ! Postprocess and write Coulomb matrix of the current q to file.
   !
   USE kinds,         ONLY : DP
   USE io_global,     ONLY : ionode, stdout
   USE qpoint,        ONLY : nksq, nksqtot
   USE lsda_mod,      ONLY : nspin
   USE io_files,      ONLY : prefix, tmp_dir
   USE crpacom,       ONLY : v_coul_bare, v_coul_scrd, filU
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iq
   ! index of the q point
   !
   CHARACTER(LEN=256) :: filename
   !! Filename to write the U matrix
   INTEGER :: iun
   !! Unit number to write to
   INTEGER :: ios
   !! I/O status
   CHARACTER(LEN=6) :: int_to_char
   !
   CALL start_clock('crpa_postproc')
   !
   IF (ionode) THEN
      !
      WRITE(stdout, '(/5x, a)') "Writing Coulomb matrix elements to file"
      !
      filename = TRIM(filU) // ".bare" // TRIM(int_to_char(iq))
      OPEN (NEWUNIT = iun, FILE = TRIM(filename), ERR = 100, IOSTAT = ios)
      CALL write_coulomb(iun, v_coul_bare)
      CLOSE(iun, STATUS = "keep")
      !
      filename = TRIM(filU) // ".scrd" // TRIM(int_to_char(iq))
      OPEN (NEWUNIT = iun, FILE = TRIM(filename), ERR = 100, IOSTAT = ios)
      CALL write_coulomb(iun, v_coul_scrd)
      CLOSE(iun, STATUS = "keep")
      !
      WRITE(stdout, '(5x, a/)') "Done writing Coulomb matrix elements to file"
      !
   ENDIF
   !
   CALL stop_clock('crpa_postproc')
   !
   RETURN
100  CALL errore('crpa_postprocess_coulomb', 'opening file' // TRIM(filename), ABS(ios))
   !
CONTAINS
!
SUBROUTINE write_coulomb(iun, v)
   !
   USE crpa_pert,     ONLY : npert_tot, nmels_tot
   USE crpa_pert,     ONLY : pert_basis
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iun
   !! Unit number to write to
   COMPLEX(DP), INTENT(IN) :: v(nmels_tot, npert_tot)
   !
   IF (pert_basis == 'bands') THEN
      CALL write_coulomb_bands(iun, v)
   ELSEIF (pert_basis == 'wannier') THEN
      CALL write_coulomb_wannier(iun, v)
   ENDIF
   !
END SUBROUTINE write_coulomb
!
SUBROUTINE write_coulomb_bands(iun, v)
   !
   USE constants,     ONLY : RYTOEV
   USE crpa_pert,     ONLY : npert_tot, nmels_tot, pert_nbnd
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iun
   !! Unit number to write to
   COMPLEX(DP), INTENT(IN) :: v(nmels_tot, npert_tot)
   !
   INTEGER :: ib1, ib2, jb1, jb2, ik1, ik2, ipert1, ipert2
   !
   WRITE(iun, '(a)') "# ib1 jb1 ik1 ib2 jb2 ik2      Re V (eV)      Im V (eV)"
   !
   DO ipert2 = 1, npert_tot
      !
      ib2 = MOD(ipert2 - 1, pert_nbnd) + 1
      jb2 = MOD((ipert2 - 1) / pert_nbnd, pert_nbnd) + 1
      ik2 = (ipert2 - 1) / (pert_nbnd * pert_nbnd) + 1
      !
      DO ipert1 = 1, npert_tot
         !
         ib1 = MOD(ipert1 - 1, pert_nbnd) + 1
         jb1 = MOD((ipert1 - 1) / pert_nbnd, pert_nbnd) + 1
         ik1 = (ipert1 - 1) / (pert_nbnd * pert_nbnd) + 1
         !
         WRITE(iun, '(6I4, 2F20.12)') ib1, jb1, ik1, ib2, jb2, ik2, v(ipert1, ipert2) * RYTOEV
      ENDDO
   ENDDO
   !
END SUBROUTINE write_coulomb_bands
!
SUBROUTINE write_coulomb_wannier(iun, v)
   !
   USE constants,     ONLY : RYTOEV
   USE crpa_pert,     ONLY : npert_tot, pert_iwlist, pert_jwlist, pert_Rlist, &
                             nmels_tot, mels_iwlist, mels_jwlist, mels_Rlist
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iun
   !! Unit number to write to
   COMPLEX(DP), INTENT(IN) :: v(nmels_tot, npert_tot)
   !
   INTEGER :: iw1, iw2, jw1, jw2, R1(3), R2(3), ipert1, ipert2
   !
   WRITE(iun, '(a)') "# iw1 jw1 R1a R1b R1c iw2 jw2 R2a R2b R2c      Re V (eV)      Im V (eV)"
   !
   DO ipert2 = 1, npert_tot
      !
      iw2 = pert_iwlist(ipert2)
      jw2 = pert_jwlist(ipert2)
      R2  = pert_Rlist(:, ipert2)
      !
      DO ipert1 = 1, nmels_tot
         !
         iw1 = mels_iwlist(ipert1)
         jw1 = mels_jwlist(ipert1)
         R1  = mels_Rlist(:, ipert1)
         !
         WRITE(iun, '(10I4, 2F20.12)') iw1, jw1, R1, iw2, jw2, R2, v(ipert1, ipert2) * RYTOEV
      ENDDO
   ENDDO
   !
END SUBROUTINE write_coulomb_wannier
  !
END SUBROUTINE crpa_postprocess_coulomb
