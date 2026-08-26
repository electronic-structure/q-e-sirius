!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------
SUBROUTINE elph_q_points ( )
!------------------------------------------------
   !
   ! Generate the q points
   !
   USE kinds,            ONLY : DP
   USE io_global,        ONLY : stdout, meta_ionode
   USE symm_base,        ONLY : nsym, s, time_reversal, t_rev, invs
   USE cell_base,        ONLY : at, bg
   USE qpoint,           ONLY : nq1, nq2, nq3, x_q, wq, nqs, lgamma_iq

   implicit none

   integer :: i, iq, ios, iun
   logical :: exist_gamma, check
   logical, external :: check_q_points_sym
   real(DP), allocatable :: xq(:,:)
   INTEGER :: nqmax
   !
   ALLOCATE(lgamma_iq(nqs))
   DO iq = 1, nqs
      lgamma_iq(iq)= ( ABS(x_q(1,iq)) < 1.0e-10_dp ) .AND. &
                     ( ABS(x_q(2,iq)) < 1.0e-10_dp ) .AND. &
                     ( ABS(x_q(3,iq)) < 1.0e-10_dp )
   ENDDO
   !
   !
   ! Write the q points in the output
   !
   write(stdout, '(//5x,"The grid of q-points (", 2(i2,","),i2,")",2x, &
   & "(",i3," q-points ) :")') nq1, nq2, nq3, nqs
   write(stdout, '(5x,"  N       xq(1)         xq(2)         xq(3)       wq" )')
   do iq = 1, nqs
      write(stdout, '(5x,i3, 4f14.9)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq), wq(iq)
   enddo
!    !
!    !  Check that the q point grid is compatible with the symmetry.
!    !
!    IF (nsym>1) THEN
!       !
!       check = check_q_points_sym(nqs, x_q, at, bg, nsym, s, invs, nq1, nq2, nq3)
!       !
!       IF (.NOT.check) THEN
!          WRITE(stdout, '(/,5x,"This q-mesh breaks symmetry!")')
!          WRITE(stdout, '(/,5x,"Try to disable the symmetry (nosym=.true. and noinv=.true. in PWscf)!")')
!          WRITE(stdout, '(5x,"Or try to choose different nq1, nq2, nq3")')
!          CALL errore('elph_q_points', 'q-mesh breaks symmetry', 1)
!       ENDIF
!       !
!    ENDIF
!    !
!    ! Write the q points to filU.q file
!    !
!    IF (meta_ionode) THEN
!       OPEN(NEWUNIT=iun, file=TRIM(filU)//'.q', status='unknown', ERR=100, IOSTAT=ios)
!       WRITE(iun, '(3i8)' ) nq1, nq2, nq3
!       WRITE(iun, '( i8)' ) nqs
!       DO iq = 1, nqs
!          WRITE(iun, '(3e24.15)') x_q(1, iq), x_q(2, iq), x_q(3, iq)
!       ENDDO
!       CLOSE(iun, STATUS = 'keep')
!    ENDIF
!    !
!    deallocate(wq)
!    !
!    RETURN
!    !
! 100 CALL errore('elph_q_points', 'opening file' // TRIM(filU) // ".q", ABS(ios))
   !
END SUBROUTINE elph_q_points
