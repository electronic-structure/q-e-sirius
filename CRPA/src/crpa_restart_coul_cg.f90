!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! q-level restart of the CG (Sternheimer) screened-Coulomb path
! (crpa_solve_linear_system.f90), with crpa_coul_restart = .true.
!
! All the results of a q point are in filU.bare#iq and filU.scrd#iq, written
! by crpa_postprocess_coulomb, and nothing after the q loop uses them in
! memory: a finished q point is skipped, with no NSCF and no Sternheimer
! solves. After both files are written the ionode creates the marker
! filU.scrd#iq.done, so a run interrupted while writing them leaves no
! marker. As for TDDFPT, the restart needs the same input.
!
!-----------------------------------------------------------------------
MODULE crpa_coul_cg_restart
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: crpa_coul_cg_q_done, crpa_coul_cg_write_done
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  FUNCTION result_file(iq, ext)
    !-----------------------------------------------------------------------
    USE crpacom,  ONLY : filU
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iq
    CHARACTER(LEN=*), INTENT(IN) :: ext
    CHARACTER(LEN=320) :: result_file
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    result_file = TRIM(filU) // ext // TRIM(int_to_char(iq))
  END FUNCTION result_file
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_coul_cg_write_done(iq)
    !-----------------------------------------------------------------------
    !! Creates the marker of q point iq, after crpa_postprocess_coulomb(iq).
    !
    USE io_global,  ONLY : ionode, ionode_id
    USE mp,         ONLY : mp_bcast
    USE mp_images,  ONLY : intra_image_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iq
    !
    INTEGER :: iun, ios
    !
    ios = 0
    IF (ionode) THEN
       OPEN(NEWUNIT=iun, FILE=TRIM(result_file(iq, '.scrd'))//'.done', STATUS='replace', IOSTAT=ios)
       IF (ios == 0) WRITE(iun, '(a)') '# CRPA CG screened Coulomb: q point completed'
       IF (ios == 0) CLOSE(iun, STATUS='keep')
    ENDIF
    CALL mp_bcast(ios, ionode_id, intra_image_comm)
    IF (ios /= 0) CALL errore('crpa_coul_cg_write_done', 'cannot write the marker of q point', iq)
    !
  END SUBROUTINE crpa_coul_cg_write_done
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_coul_cg_q_done(iq, done)
    !-----------------------------------------------------------------------
    !! done = .true. if q point iq was completed by a previous run: its marker
    !! and result files exist, and the caller skips it.
    !
    USE io_global,  ONLY : ionode, ionode_id, stdout
    USE mp,         ONLY : mp_bcast
    USE mp_images,  ONLY : intra_image_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iq
    LOGICAL, INTENT(OUT) :: done
    !
    done = .FALSE.
    IF (ionode) THEN
       INQUIRE(FILE=TRIM(result_file(iq, '.scrd'))//'.done', EXIST=done)
       IF (done) INQUIRE(FILE=TRIM(result_file(iq, '.bare')), EXIST=done)
       IF (done) INQUIRE(FILE=TRIM(result_file(iq, '.scrd')), EXIST=done)
    ENDIF
    CALL mp_bcast(done, ionode_id, intra_image_comm)
    !
    IF (done) WRITE(stdout, '(/,5x,"q point #",i4,": screened Coulomb taken from ",a, &
         &" (crpa_coul_restart), no NSCF and no Sternheimer solves")') &
         iq, TRIM(result_file(iq, '.scrd'))
    !
  END SUBROUTINE crpa_coul_cg_q_done
  !
END MODULE crpa_coul_cg_restart
