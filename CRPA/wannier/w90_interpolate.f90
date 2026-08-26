!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE w90_interpolate
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE w90_interpolate_wfc(ik, evc_kp_G)
    !----------------------------------------------------------------------------
    !! Using real-space Wannier gauge wavefunctions on the coarse grid, interpolate
    !! the wavefunctions to an arbitrary k point.
    !----------------------------------------------------------------------------
    USE kinds,            ONLY : DP
    USE constants,        ONLY : tpi
    USE fft_interfaces,   ONLY : fwfft
    USE fft_base,         ONLY : dffts
    USE klist,            ONLY : xk, igk_k, ngk
    USE wvfct,            ONLY : npwx
    USE cell_base,        ONLY : at
    USE noncollin_module, ONLY : npol
    USE w90_interface,    ONLY : num_wann, wann_centers
    USE w90_real_space,   ONLY : nR_ws, iRlist_ws
    USE w90_wan_Rr_buffer,ONLY : w90_wan_Rr_get_buffer
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! k point index to interpolate to
    COMPLEX(DP), INTENT(out) :: evc_kp_G(npwx*npol, num_wann)
    !! Interpolated wavefunction at kp
    !
    INTEGER :: iw
    !! Wannier function index
    INTEGER :: ipol
    !! Spin index
    INTEGER :: iR
    !! Real-space lattice vector index
    INTEGER :: npw
    !! Number of plane waves at xkp
    REAL(DP) :: kdotR
    !! Dot product of k and R
    REAL(DP) :: xkp(3)
    !! k point to interpolate the wavefunctions to (in Cartesian coordinates).
    REAL(DP) :: R(3)
    !! Real-space lattice vector
    COMPLEX(DP) :: fac
    !! Factor for Fourier transform
    COMPLEX(DP), ALLOCATABLE :: evc_kp_r(:, :, :)
    !! Interpolated wavefunction at kp in real space
    COMPLEX(DP), ALLOCATABLE :: wan_Rr(:, :, :)
    !! Wavefunction in the Wannier gauge in real space
    !
    CALL start_clock("w90_interpolate")
    !
    xkp = xk(:, ik)
    npw = ngk(ik)
    !
    ALLOCATE(evc_kp_r(dffts%nnr, npol, num_wann))
    ALLOCATE(wan_Rr(dffts%nnr, npol, num_wann))
    !
    evc_kp_r = (0.d0, 0.d0)
    !
    ! Read wan_Rr wavefunction from file
    !
    DO iR = 1, nR_ws
      R(:) = REAL(iRlist_ws(:, iR), KIND=DP)
      CALL cryst_to_cart(1, R, at, +1)
      !
      CALL w90_wan_Rr_get_buffer(iR, wan_Rr)
      !
      kdotR = tpi * SUM(xkp * R)
      fac = CMPLX(COS(kdotR), SIN(kdotR))
      !
      evc_kp_r = evc_kp_r + fac * wan_Rr
      !
    ENDDO ! iR
    !
    ! Convert from psi_k(r) to u_k(r) = psi_k(r) * exp(-i k r)
    !
    DO iw = 1, num_wann
      DO ipol = 1, npol
        CALL w90_multiply_iqr(dffts, -xkp, evc_kp_r(:, ipol, iw))
      ENDDO
    ENDDO
    !
    evc_kp_G = (0.d0, 0.d0)
    !
    DO iw = 1, num_wann
      DO ipol = 1, npol
        CALL fwfft('Wave', evc_kp_r(:, ipol, iw), dffts)
        evc_kp_G((ipol-1)*npwx+1 : (ipol-1)*npwx+npw, iw) = evc_kp_r(dffts%nl(igk_k(1:npw, ik)), ipol, iw)
      ENDDO
    ENDDO
    !
    ! Shift the wavefunction back to its original center
    ! (see the comment in w90_real_space.f90)
    !
    DO iw = 1, num_wann
      CALL w90_shift_wfc_G(ik, wann_centers(:, iw), evc_kp_G(:, iw))
    ENDDO
    !
    CALL stop_clock("w90_interpolate")
    !
  !----------------------------------------------------------------------------
  END SUBROUTINE w90_interpolate_wfc
  !----------------------------------------------------------------------------
!
END MODULE w90_interpolate
