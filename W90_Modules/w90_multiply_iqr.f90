!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE w90_multiply_iqr(dfft, xq, func)
!----------------------------------------------------------------------------
!!
!! Multiply exp(i*q*r) to func.
!! Real space indexing is adapted from Modules/compute_dipole.f90
!!
!----------------------------------------------------------------------------
  USE kinds,               ONLY : DP
  USE constants,           ONLY : tpi
  USE cell_base,           ONLY : at
  USE fft_types,           ONLY : fft_type_descriptor, fft_index_to_3d
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfft
  ! fft_type_descriptor. dffts or dfftp
  REAL(DP), INTENT(in) :: xq(3)
  ! q vector in cartesian coordinate
  COMPLEX(DP), INTENT(inout) :: func(dfft%nnr)
  ! input, output: real-space function
  !
  LOGICAL :: offrange
  INTEGER :: ir, i, j, k, ir_end
  REAL(DP) :: arg, xq_crys(3)
  COMPLEX(DP) :: phase
  !
  CALL start_clock('w90_mult_iqr')
  !
  xq_crys = xq
  CALL cryst_to_cart(1, xq_crys, at, -1)
  !
#if defined (__MPI)
  ir_end = MIN(dfft%nnr, dfft%nr1x*dfft%my_nr2p*dfft%my_nr3p)
#else
  ir_end = dfft%nnr
#endif
  !
  DO ir = 1, ir_end
    !
    CALL fft_index_to_3d(ir, dfft, i, j, k, offrange)
    IF ( offrange ) CYCLE
    !
    ! (i,j,k) is the zero-based coordinate of the real-space grid
    arg = tpi * (   xq_crys(1) * REAL(i, DP) / REAL(dfft%nr1, DP) &
                  + xq_crys(2) * REAL(j, DP) / REAL(dfft%nr2, DP) &
                  + xq_crys(3) * REAL(k, DP) / REAL(dfft%nr3, DP)  )
    phase = CMPLX( COS(arg), SIN(arg), kind=DP )
    !
    func(ir) = func(ir) * phase
    !
  ENDDO ! ir
  !
  CALL stop_clock('w90_mult_iqr')
  !
!----------------------------------------------------------------------------
END SUBROUTINE w90_multiply_iqr
!----------------------------------------------------------------------------
