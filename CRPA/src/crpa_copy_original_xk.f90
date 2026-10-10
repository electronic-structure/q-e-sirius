!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_copy_original_xk()
   !-----------------------------------------------------------------------
   !
   ! Copy the list of k points (directly read from the previous SCF calculation) to
   ! nks_orig and xk_orig. These are the full collected k points (not the distributed ones).
   !
   ! This step is necessary as DFPT calculation at q /= 0 overwrites xk and changes the
   ! order of k vectors.
   !
   USE kinds,           ONLY : DP
   USE klist,           ONLY : nks, nkstot, xk
   USE crpacom,         ONLY : nks_orig, xk_orig
   !
   IMPLICIT NONE
   !
   REAL(DP), ALLOCATABLE :: xk_orig_loc(:, :)
   !! Local list of k points
   !
   ALLOCATE(xk_orig_loc(3, nks))
   ALLOCATE(xk_orig(3, nkstot))
   !
   xk_orig_loc = xk(:, 1:nks)
   !
   CALL poolgather(3, nkstot, nks, xk_orig_loc, xk_orig)
   nks_orig = nkstot
   !
   DEALLOCATE(xk_orig_loc)
   !
END SUBROUTINE crpa_copy_original_xk


!--------------------------------------------------------------------
SUBROUTINE poolgather(nsize, nkstot, nks, f_in, f_out)
   !--------------------------------------------------------------------
   !!
   !! Gather the kpoints and the electronic eigenvalues
   !! across the pools
   !! doesn't work with the double grid (k and k+q)
   !!
   USE kinds,     only : DP
   USE mp_global, ONLY : my_pool_id, inter_pool_comm, kunit,npool, my_pool_id
   USE mp,        ONLY : mp_barrier, mp_bcast,mp_sum
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(in) :: nsize
   !! first dimension of vectors f_in and f_out
   INTEGER, INTENT(in) :: nks
   !! number of k-points per pool
   INTEGER, INTENT(in) :: nkstot
   !! total number of k-points
   REAL(KIND = DP), INTENT(in) :: f_in(nsize, nks)
   !! input ( only for k-points of mypool )
   REAL(KIND = DP), INTENT(out) :: f_out(nsize, nkstot)
   !! output  ( contains values for all k-point )
   !
   ! Local variables
#if defined(__MPI)
   INTEGER :: rest
   !! the rest of the INTEGER division nkstot / npo
   INTEGER :: nbase
   ! the position in the original list
   !
   rest = nkstot / kunit - (nkstot / kunit / npool) * npool
   !
   nbase = nks * my_pool_id
   !
   IF ((my_pool_id + 1) > rest) nbase = nbase + rest * kunit
   f_out = 0.d0
   f_out(:, (nbase + 1):(nbase + nks)) = f_in(:, 1:nks)
   !
   ! Reduce across the pools
   CALL mp_sum(f_out, inter_pool_comm)
   !
#else
   f_out(:, :) = f_in(:, :)
   !
#endif
   !
   !--------------------------------------------------------------------
END SUBROUTINE poolgather
!--------------------------------------------------------------------
