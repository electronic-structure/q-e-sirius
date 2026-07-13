!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE crpa_init_q()
  !----------------------------------------------------------------------------
  !
  ! This subroutine prepares several variables which are needed in the
  ! CRPA program at fixed q point:
  !
  ! 1) Compute the phase factor (for the USPP case)
  ! 2) Compute becp = <vkb|evc> for every k point
  ! 3) Calculate and write the S\phi for k and k+q
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, tau
  USE becmod,               ONLY : calbec
  USE constants,            ONLY : eps8, tpi
  USE klist,                ONLY : xk, ngk, igk_k, nks
  USE io_global,            ONLY : stdout
  USE buffers,              ONLY : get_buffer, save_buffer
  USE wavefunctions,        ONLY : evc
  USE uspp,                 ONLY : vkb, okvan
  USE eqv,                  ONLY : evq
  USE lrus,                 ONLY : becp1, becpt
  USE control_lr,           ONLY : lgamma
  USE units_lr,             ONLY : lrwfc, iuwfc
  USE qpoint,               ONLY : xq, nksq, eigqts, ikks, ikqs, ikmks
  USE noncollin_module,     ONLY : npol, noncolin, domag
  USE uspp_init,            ONLY : init_us_2
  USE wvfct,                ONLY : nbnd, npwx
  USE ldaU,                 ONLY : lda_plus_u
  USE crpacom,              ONLY : lrwan, iuwan
  USE crpa_pert,            ONLY : pert_basis
  USE w90_interface,        ONLY : num_wann
  USE w90_interpolate,      ONLY : w90_interpolate_wfc
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: ik, ikk, ikq, na
    ! generic counter
    ! counter on k points
    ! counter on k+q points
    ! counter on atoms
  INTEGER :: npw
    ! number of plane waves at k
  REAL(DP) :: arg
    ! the argument of the phase
  COMPLEX(DP), ALLOCATABLE :: tevc(:,:)
  COMPLEX(DP), ALLOCATABLE :: wan(:, :)
  !
  CALL start_clock( 'crpa_init_q' )
  !
  ! 1) USPP: Compute the phase factor exp(-i q*\tau)
  !
  IF (okvan) THEN
     DO na = 1, nat
        arg = ( xq(1) * tau(1,na) + &
                xq(2) * tau(2,na) + &
                xq(3) * tau(3,na) ) * tpi
        eigqts(na) = CMPLX( COS( arg ), - SIN( arg ) ,kind=DP)
     ENDDO
  ENDIF
  !
  IF (noncolin.AND.domag) ALLOCATE(tevc(npwx*npol,nbnd))
  !
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     !
     ! Read the wavefunctions evc (at k) and evq (at k+q).
     ! Note: this is important because if nksq=1 then evc and evq are read only
     ! once (here) and then used throughout the code. This may happen e.g.
     ! when the ratio of the total number of k points (without k+q) and k pools
     ! is not an integer number (as a consequence some k pools will have nksq=1).
     !
     CALL get_buffer (evc, lrwfc, iuwfc, ikk)
     !
     IF (noncolin .and. domag) &
         CALL get_buffer (tevc, lrwfc, iuwfc, ikmks(ik) )
     !
     IF (.NOT.lgamma .AND. nksq.EQ.1) CALL get_buffer (evq, lrwfc, iuwfc, ikq)
     !
     ! 2) USPP: Compute the becp terms which are used in the rest of the code
     !
     IF (okvan) THEN
        !
        ! Compute the beta function vkb(k+G)
        !
        CALL init_us_2 (npw, igk_k(1,ikk), xk(1,ikk), vkb, .true.)
        !
        !$acc update host(vkb)
        !
        ! becp1 = <vkb|evc>
        !
        CALL calbec (npw, vkb, evc, becp1(ik))
        IF (noncolin .and. domag) &
           CALL calbec (npw, vkb, tevc, becpt(ik) )
        !
     ENDIF
     !
  ENDDO
  !
  IF (noncolin.AND.domag) DEALLOCATE( tevc )
  !
  IF (lda_plus_u) THEN
     !
     ! 3) Calculate and write to file S\phi for k and k+q
     !
     CALL lr_orthoUwfc (.FALSE.)
  ENDIF
  !
  ! 4) Calculate and write to file Wannier functions for k and k+q
  !
  IF (pert_basis == 'wannier') THEN
     !
     WRITE( stdout, '(/5x, A, I0, A)') "Computing Wannier functions at k and k+q, ", nks, " points"
     !
     ALLOCATE(wan(npwx*npol, num_wann))
     !
     ! Loop over both k and k+q at once
     !
     DO ik = 1, nks
        WRITE(stdout, '(I8)', ADVANCE='NO') ik
        IF (MOD(ik, 10) == 0) WRITE(stdout, '()')
        CALL w90_interpolate_wfc(ik, wan)
        CALL save_buffer(wan, lrwan, iuwan, ik)
     ENDDO ! ik
     !
     DEALLOCATE(wan)
     !
     IF (MOD(nks, 10) /= 0) WRITE(stdout, '()')
     WRITE(stdout, '(5x, A/)') "Done"
     !
  ENDIF
  !
  CALL stop_clock ( 'crpa_init_q' )
  !
  RETURN
  !
END SUBROUTINE crpa_init_q
