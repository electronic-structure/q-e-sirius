!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE crpa_freq_grid
  !----------------------------------------------------------------------------
  !! Wraps crpa_build_freq_grid in a module (rather than a bare external
  !! subroutine) so that its ALLOCATABLE dummy argument gets an explicit
  !! interface automatically at every call site (USE crpa_freq_grid): an
  !! external subroutine call with no explicit interface is non-conforming
  !! for an allocatable dummy argument and corrupts the call.
  !
  IMPLICIT NONE
  !
  CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE crpa_build_freq_grid(nfreq, freq_start, freq_end, broadening, freq)
    !-----------------------------------------------------------------------
    !! Builds a real-frequency grid (with a fixed imaginary broadening added
    !! to every point):
    !! freq(i) = freq_start + (i-1)*step + i*broadening (i = imaginary unit),
    !! step = (freq_end - freq_start) / (nfreq - 1) (0 if nfreq = 1).
    !! Used by crpa_readin to build crpacom%crpa_freq, the frequencies of
    !! every path (CG or Lanczos, Coulomb and dielectric).
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfreq
    !! Number of points of the frequency grid
    REAL(DP), INTENT(IN) :: freq_start
    !! Start (Ry) of the real-frequency grid
    REAL(DP), INTENT(IN) :: freq_end
    !! End (Ry) of the real-frequency grid (ignored if nfreq = 1)
    REAL(DP), INTENT(IN) :: broadening
    !! Imaginary broadening eta (Ry) added to every frequency
    COMPLEX(DP), ALLOCATABLE, INTENT(INOUT) :: freq(:)
    !! Complex frequency grid built by this routine. Dimension: (nfreq)
    !
    INTEGER :: ifreq
    REAL(DP) :: dfreq
    !
    IF (ALLOCATED(freq)) DEALLOCATE(freq)
    ALLOCATE(freq(nfreq))
    !
    dfreq = 0.0_DP
    IF (nfreq > 1) dfreq = (freq_end - freq_start) / (nfreq - 1)
    !
    DO ifreq = 1, nfreq
       freq(ifreq) = CMPLX(freq_start + (ifreq-1)*dfreq, broadening, KIND=DP)
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE crpa_build_freq_grid
  !
END MODULE crpa_freq_grid
