!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE elph_prepare_q (iq, do_iq, setup_pw)
  !-----------------------------------------------------------------------
  !
  !  This routine prepares a few variables that are needed to control
  !  the ELPH run after the q point has been decided, but before
  !  doing the NSCF calculation. In particular it sets:
  !  xq       : the coordinate of the current q point
  !  do_iq    : if .true. q point has to be calculated
  !  lgamma   : if this is a gamma point calculation (i.e. q=0)
  !  setup_pw : if .true. the NSCF calculation will be performed
  !
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : prefix, create_directory, tmp_dir
  USE klist,            ONLY : ltetra
  USE dfpt_tetra_mod,   ONLY : dfpt_tetra_linit
  USE qpoint,           ONLY : xq, x_q
  USE control_lr,       ONLY : lgamma
  USE noncollin_module, ONLY : noncolin, domag
  USE control_ph,       ONLY : tmp_dir_ph, tmp_dir_phq, current_iq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: iq
  !! Index of the q point to be prepared
  LOGICAL, INTENT(OUT) :: do_iq
  !! If .TRUE. the q point has to be calculated
  LOGICAL, INTENT(OUT) :: setup_pw
  !! TODO
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  do_iq = .TRUE.
  !
  WRITE( stdout, '(//,5X,"=-------------------------------------------------------------=")')
  WRITE( stdout, '(/,5X,"Calculation for q #",i4," = (", 3F12.7, " )")') iq, x_q(:,iq)
  WRITE( stdout, '(/,5X, "=-------------------------------------------------------------=")')
  !
  ! Set the coordinates of the current q point
  !
  xq(1:3) = x_q(1:3,iq)
  !
  ! Check if it is the Gamma point (0,0,0)
  !
  lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
  !
  ! Check if the NSCF calculation must be performed
  !
  setup_pw = .NOT.lgamma  &            ! perform NCSF calculation for all q\=0
            ! .OR. (lgamma .AND. recalc_sym) &   ! perform NCSF calculation also for q=0
             .OR.(noncolin.AND.domag)      ! only if the number of symmetries was reduced
                                           ! (which may lead to the increase in the number of k points)
  !
  ! If setup_pw=.TRUE. then tetra is initialized through setup_nscf.
  ! If setup_pw=.FALSE. then we need to set here dfpt_tetra_linit=.TRUE.
  ! which is needed in dfpt_tetra_setup in order to initizalize tetra.
  !
  dfpt_tetra_linit = .FALSE.
  IF ( (.NOT.setup_pw) .AND. ltetra) dfpt_tetra_linit = .TRUE.
!   !
!   IF (lgamma .AND. recalc_sym) &
!      WRITE( stdout, '(/,5X,"Do NSCF calculation at q=0 because the number of symmetries was reduced")')
  !
  current_iq = iq
  !
  tmp_dir_phq = trimcheck(TRIM(tmp_dir_ph) // TRIM(prefix) // '.q_' // int_to_char(current_iq))
  CALL create_directory(tmp_dir_phq)
  tmp_dir = tmp_dir_phq
  !
END SUBROUTINE elph_prepare_q
