!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------
SUBROUTINE crpa_q_points ( )
!------------------------------------------------
   !
   ! Generate the q points
   !
   USE kinds,            ONLY : DP
   USE io_global,        ONLY : stdout, meta_ionode
   USE symm_base,        ONLY : nsym, s, time_reversal, t_rev, invs
   USE cell_base,        ONLY : at, bg
   USE qpoint,           ONLY : nq1, nq2, nq3, x_q, nqs
   USE crpacom,          ONLY : skip_equivalence_q, filU
   USE crpacom,          ONLY : crpa_reduce_q, crpa_mode
   USE crpa_qpoints,     ONLY : qplot

   implicit none

   integer :: i, iq, ios, iun
   logical :: exist_gamma, check
   logical, external :: check_q_points_sym
   real(DP), allocatable :: xq(:,:), wq(:)
   INTEGER :: nqmax
   !
   ! Crystal point group used to classify the q grid when crpa_reduce_q is set.
   ! It is recovered from the structure (crpa_crystal_sym) WITHOUT touching
   ! symm_base, so the DFPT still runs unsymmetrized. Used only in the
   ! crpa_reduce_q branches below; when not reducing, the punch-file symmetry
   ! (symm_base) is used directly, exactly as before.
   !
   INTEGER :: nsym_q, s_q(3,3,48), t_rev_q(48), invs_q(48)
   !
   IF (crpa_reduce_q) THEN
      CALL crpa_crystal_sym(nsym_q, s_q, t_rev_q, invs_q)
      WRITE(stdout, '(/5x,a,i3,a)') &
         "crpa_reduce_q: using the crystal point group (", nsym_q, &
         " operations) to select the irreducible q points"
      IF (TRIM(crpa_mode) == 'debug') WRITE(stdout, '(5x,a,/5x,a,/5x,a)') &
         "WARNING: crpa_mode='debug' computes a single Wannier matrix element,", &
         "which is NOT invariant under the star permutation. A weighted sum over", &
         "the irreducible q will NOT reproduce the full-grid average. Use 'dHP'."
   ENDIF
   !
   IF (qplot) THEN
      !
      ! q points in x_q are set in crpa_readin
      !
      ! We temporarily allocate wq. It is used only for printing.
      !
      ALLOCATE(wq(nqs))
      wq = 1.d0
      !
      ! FIXME: Unify interface
      !
   ELSE  ! .NOT. qplot
      !
      if ( nq1 <= 0 .or. nq2 <= 0 .or. nq3 <= 0 ) &
         call errore('crpa_q_points','nq1 or nq2 or nq3 <= 0',1)
      !
      nqmax = nq1 * nq2 * nq3
      !
      allocate (wq(nqmax))
      allocate (xq(3,nqmax))
      !
      IF (crpa_reduce_q) THEN
         ! Reduce to the irreducible q using the crystal group, with time
         ! reversal (q <-> -q) enabled and reduction forced on.
         CALL kpoint_grid( nsym_q, .TRUE., .FALSE., s_q, t_rev_q, bg, nqmax,&
            0,0,0, nq1,nq2,nq3, nqs, xq, wq )
      ELSE
         CALL kpoint_grid( nsym, time_reversal, skip_equivalence_q, s, t_rev, bg, nqmax,&
            0,0,0, nq1,nq2,nq3, nqs, xq, wq )
      ENDIF
      !
      allocate(x_q(3,nqs))
      x_q(:,:) = xq(:,1:nqs)
      deallocate (xq)
      !
   ENDIF
   !
   ! Check if the Gamma point is one of the points and put
   ! it in the first position (it should already be the first)
   !
   exist_gamma = .false.
   !
   do iq = 1, nqs
      if ( abs(x_q(1,iq)) .lt. 1.0e-10_dp .and. &
            abs(x_q(2,iq)) .lt. 1.0e-10_dp .and. &
            abs(x_q(3,iq)) .lt. 1.0e-10_dp ) then
         exist_gamma = .true.
         if (iq .ne. 1) then
            do i = 1, 3
               x_q(i,iq) = x_q(i,1)
               x_q(i,1) = 0.0_dp
            end do
         end if
      end if
   end do
   !
   ! Write the q points in the output
   !
   write(stdout, '(//5x,"The grid of q-points (", 2(i2,","),i2,")",2x, &
   & "(",i3," q-points ) :")') nq1, nq2, nq3, nqs
   write(stdout, '(5x,"  N       xq(1)         xq(2)         xq(3)       wq" )')
   do iq = 1, nqs
      write(stdout, '(5x,i3, 4f14.9)') iq, x_q(1,iq), x_q(2,iq), x_q(3,iq), wq(iq)
   enddo
   !
   !  Check that the q point grid is compatible with the symmetry.
   !
   check = .TRUE.
   IF (crpa_reduce_q) THEN
      IF (nsym_q>1) check = check_q_points_sym(nqs, x_q, at, bg, nsym_q, s_q, invs_q, nq1, nq2, nq3)
   ELSE
      IF (nsym>1)   check = check_q_points_sym(nqs, x_q, at, bg, nsym, s, invs, nq1, nq2, nq3)
   ENDIF
   !
   IF (.NOT.check) THEN
      WRITE(stdout, '(/,5x,"This q-mesh breaks symmetry!")')
      WRITE(stdout, '(/,5x,"Try to disable the symmetry (nosym=.true. and noinv=.true. in PWscf)!")')
      WRITE(stdout, '(5x,"Or try to choose different nq1, nq2, nq3")')
      CALL errore('crpa_q_points', 'q-mesh breaks symmetry', 1)
   ENDIF
   !
   ! Write the q points to filU.q file
   !
   IF (meta_ionode) THEN
      OPEN(NEWUNIT=iun, file=TRIM(filU)//'.q', status='unknown', ERR=100, IOSTAT=ios)
      WRITE(iun, '(3i8)' ) nq1, nq2, nq3
      WRITE(iun, '( i8)' ) nqs
      DO iq = 1, nqs
         WRITE(iun, '(4e24.15)') x_q(1, iq), x_q(2, iq), x_q(3, iq), wq(iq)
      ENDDO
      CLOSE(iun, STATUS = 'keep')
   ENDIF
   !
   deallocate(wq)
   !
   RETURN
   !
100 CALL errore('crpa_q_points', 'opening file ' // TRIM(filU) // ".q", ABS(ios))
   !
END SUBROUTINE crpa_q_points
!
!------------------------------------------------------------------------------
SUBROUTINE crpa_crystal_sym(nsym_c, s_c, t_rev_c, invs_c)
   !---------------------------------------------------------------------------
   !! Return the crystal point group (rotations s_c, time-reversal flags t_rev_c,
   !! inverse map invs_c, order nsym_c) of the current structure.
   !!
   !! The SCF/NSCF are run with nosym (CRPA forces it, see setup.f90 "IF (lecrpa)
   !! nosym=.TRUE."), so the punch file only records nsym=1: it keeps all nrot
   !! LATTICE operations but tags only the identity as a crystal symmetry
   !! (qexsd_init.f90, "IF (i<=nsym) crystal_symmetry"), and the crystal-group
   !! count is not stored. It therefore cannot be read back and must be
   !! recomputed from the lattice + basis with the standard PW routines
   !! set_sym_bl + find_sym (exactly what setup.f90 does at start-up).
   !!
   !! Those routines write their result into the global symm_base module. Since
   !! the rest of the CRPA run (in particular the DFPT solve) must keep seeing
   !! the unsymmetrized nsym=1 state, this routine SAVES the symm_base variables
   !! they touch, recomputes the symmetry, copies out what the caller needs, and
   !! RESTORES symm_base to its exact previous contents.
   !---------------------------------------------------------------------------
   USE kinds,            ONLY : DP
   USE ions_base,        ONLY : nat, tau, ityp
   USE symm_base,        ONLY : nrot, nsym, nsym_ns, nsym_na, s, sr, ft, sname, &
                                t_rev, invs, irt, invsym, time_reversal, &
                                set_sym_bl, find_sym, inverse_s
   USE noncollin_module, ONLY : noncolin, domag
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(OUT) :: nsym_c, s_c(3,3,48), t_rev_c(48), invs_c(48)
   !
   ! saved symm_base state
   INTEGER :: nrot_s, nsym_s, nsym_ns_s, nsym_na_s
   INTEGER :: s_s(3,3,48), t_rev_s(48), invs_s(48)
   REAL(DP) :: sr_s(3,3,48), ft_s(3,48)
   CHARACTER(LEN=45) :: sname_s(48)
   LOGICAL :: invsym_s, time_reversal_s, magnetic_sym
   INTEGER, ALLOCATABLE :: irt_s(:,:)
   REAL(DP), ALLOCATABLE :: m_loc0(:,:)
   !
   ! Save everything set_sym_bl / find_sym / inverse_s may overwrite.
   !
   nrot_s = nrot;  nsym_s = nsym;  nsym_ns_s = nsym_ns;  nsym_na_s = nsym_na
   s_s = s;  sr_s = sr;  ft_s = ft;  sname_s = sname
   t_rev_s = t_rev;  invs_s = invs
   invsym_s = invsym;  time_reversal_s = time_reversal
   IF (ALLOCATED(irt)) THEN
      ALLOCATE(irt_s(SIZE(irt,1), SIZE(irt,2)))
      irt_s = irt
   ENDIF
   !
   ! Recompute the crystal symmetry from the structure. m_loc drives magnetic
   ! symmetry analysis only; pass zeros (Wannier-basis CRPA is nonmagnetic).
   !
   magnetic_sym = noncolin .AND. domag
   ALLOCATE(m_loc0(3,nat));  m_loc0 = 0.0_DP
   !
   CALL set_sym_bl()
   CALL find_sym(nat, tau, ityp, magnetic_sym, m_loc0)
   CALL inverse_s()
   !
   nsym_c  = nsym
   s_c     = s
   t_rev_c = t_rev
   invs_c  = invs
   !
   DEALLOCATE(m_loc0)
   !
   ! Restore symm_base exactly as it was.
   !
   nrot = nrot_s;  nsym = nsym_s;  nsym_ns = nsym_ns_s;  nsym_na = nsym_na_s
   s = s_s;  sr = sr_s;  ft = ft_s;  sname = sname_s
   t_rev = t_rev_s;  invs = invs_s
   invsym = invsym_s;  time_reversal = time_reversal_s
   IF (ALLOCATED(irt_s)) THEN
      irt = irt_s
      DEALLOCATE(irt_s)
   ENDIF
   !
END SUBROUTINE crpa_crystal_sym
