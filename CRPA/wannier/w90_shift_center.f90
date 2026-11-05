!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE w90_shift_wfc(ik, shift, wfc_r)
  !----------------------------------------------------------------------------
  !! Shift evc in real space by shift.
  !----------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE fft_interfaces,   ONLY : invfft, fwfft
  USE fft_base,         ONLY : dffts
  USE klist,            ONLY : igk_k, ngk
  USE wvfct,            ONLY : npwx
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k point index
  REAL(DP), INTENT(IN) :: shift(3)
  !! Distance to shift in real space (in Cartesian coordinates)
  COMPLEX(DP), INTENT(INOUT) :: wfc_r(dffts%nnr, npol)
  !! Real space wavefunction to shift
  !
  INTEGER :: ipol
  !! Spin index
  INTEGER :: npw
  !! Number of plane waves
  COMPLEX(DP), ALLOCATABLE :: wfc_G(:)
  !! Buffer for wavefunction in reciprocal space
  !
  CALL start_clock('w90_shift_wfc')
  !
  ALLOCATE(wfc_G(npwx*npol))
  !
  npw = ngk(ik)
  !
  wfc_G = (0.d0, 0.d0)
  DO ipol = 1, npol
    CALL fwfft('Wave', wfc_r(:, ipol), dffts)
    wfc_G((ipol-1)*npwx+1 : (ipol-1)*ipol+npw) = wfc_r(dffts%nl(igk_k(1:npw, ik)), ipol)
  ENDDO
  !
  CALL w90_shift_wfc_G(ik, shift, wfc_G)
  !
  wfc_r = (0.d0, 0.d0)
  DO ipol = 1, npol
    wfc_r(dffts%nl(igk_k(1:npw, ik)), ipol) = wfc_G((ipol-1)*npwx+1 : (ipol-1)*ipol+npw)
    CALL invfft('Wave', wfc_r(:, ipol), dffts)
  ENDDO
  !
  DEALLOCATE(wfc_G)
  !
  CALL stop_clock('w90_shift_wfc')
  !
  !----------------------------------------------------------------------------
END SUBROUTINE w90_shift_wfc
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE w90_shift_wfc_G(ik, shift, wfc_G)
  !----------------------------------------------------------------------------
  !! Shift evc in real space by shift.
  !----------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE constants,        ONLY : tpi
  USE klist,            ONLY : xk, igk_k, ngk
  USE gvect,            ONLY : g
  USE wvfct,            ONLY : npwx
  USE noncollin_module, ONLY : noncolin, npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k point index
  REAL(DP), INTENT(IN) :: shift(3)
  !! Distance to shift in real space (in Cartesian coordinates)
  COMPLEX(DP), INTENT(INOUT) :: wfc_G(npwx * npol)
  !! Wavefunction to shift in reciprocal space
  !
  INTEGER :: ig
  !! G vector index
  INTEGER :: npw
  !! Number of plane waves
  REAL(DP) :: arg
  !! Argument of the exponent
  REAL(DP) :: gk(3)
  !! k + G vector
  COMPLEX(DP) :: phase
  !! exp(-i(k+G) * shift) factor
  !
  CALL start_clock('w90_shift_wfc_G')
  !
  npw = ngk(ik)
  !
  ! wfc_G(G) *= exp(-i * (k + G) * shift)
  !
  DO ig = 1, npw
    !
    gk = xk(:, ik) + g(:, igk_k(ig, ik))
    arg = tpi * SUM(shift * gk)
    phase = CMPLX(COS(arg), -SIN(arg), KIND=DP)
    !
    wfc_G(ig) = wfc_G(ig) * phase
    !
    IF (noncolin) wfc_G(ig+npwx) = wfc_G(ig+npwx) * phase
    !
  ENDDO
  !
  CALL stop_clock('w90_shift_wfc_G')
  !
  !----------------------------------------------------------------------------
END SUBROUTINE w90_shift_wfc_G
!----------------------------------------------------------------------------
