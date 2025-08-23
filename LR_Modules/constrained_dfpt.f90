!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------------------
MODULE constrained_dfpt
   !----------------------------------------------------------------------------------------
   !! This module contains routines for the constrained DFPT method.
   !! The constrained DFPT method is used to compute the response of a system where an
   !! "active-space" susceptibility is removed from the full susceptibility.
   !!
   !! Before calling dfpt_kernel, the high-level routine must allocate and set chi_active
   !! to be the active-space susceptibility in the band basis.
   !! - If the active space is defined by band indices, one can use cdfpt_chi_active_by_bands.
   !! - If the active space is defined otherwise (e.g. by Wannier functions), calculation
   !!   of chi_active is left to the user.
   !!
   !! Subroutine cdfpt_subtract_active_wfc uses chi_active and subtracts the active-space
   !! contribution from the wavefunction response. It is called after cgsolve_all
   !! (the Sternheimer solver) inside subroutine sternheimer_kernel.
   !----------------------------------------------------------------------------------------
   !
   USE kinds,                ONLY : DP
   !
   SAVE
   !
   LOGICAL :: lcdfpt = .FALSE.
   !! Logical flag to enable constrained DFPT
   !
   COMPLEX(DP), ALLOCATABLE :: chi_active(:, :, :)
   !! Bare active-space susceptibility, to be subtracted from the full susceptibility.
   !
   CONTAINS
   !
   !----------------------------------------------------------------------------------------
   SUBROUTINE cdfpt_subtract_active_wfc(ik, dvpsi, dpsi)
   !----------------------------------------------------------------------------------------
   !!
   !! This routine updates the wavefunction response by subtracting the
   !! active-space contribution as needed in the constrained DFPT method.
   !!
   !! - dvpsi : (Input) dV * psi, before orthogonalization.
   !! - dvscf : (In/Out) In : Wavefunction response from DFPT.
   !!                    Out: Wavefunction response minus the active-space contribution.
   !!
   !----------------------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE mp,                   ONLY : mp_sum
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE noncollin_module,     ONLY : npol
      USE wvfct,                ONLY : nbnd, npwx, et
      USE klist,                ONLY : lgauss, degauss
      USE qpoint,               ONLY : ikks, ikqs
      USE control_lr,           ONLY : nbnd_occ
      USE eqv,                  ONLY : evq
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ik
      !! k point index
      COMPLEX(DP), INTENT(IN) :: dvpsi(npwx * npol, nbnd)
      !! (dV_bare + dV_induced) * psi, before orthogonaliztaion.
      COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx * npol, nbnd)
      !! Wavefunction renormalization
      !
      INTEGER :: ikk, ikq, ibnd, jbnd
      !! indices
      REAL(DP) :: deltae
      !! e_mkq - e_nk
      REAL(DP) :: theta
      !! f(e_mkq - e_nk)
      COMPLEX(DP), ALLOCATABLE :: mel(:, :)
      !! Matrix element of the total perturbation potential
      COMPLEX(DP), ALLOCATABLE :: ps(:, :)
      !! Matrix element to be used for subtracting the active space contribution
      REAL(DP), EXTERNAL :: wgauss
      !! Occupation function
      !
      CALL start_clock('cdfpt_subtract')
      !
      ikk = ikks(ik)
      ikq = ikqs(ik)
      !
      ALLOCATE(mel(nbnd, nbnd))
      ALLOCATE(ps(nbnd, nbnd_occ(ikk)))
      !
      ! Compute the full matrix elements (bare + induced)
      ! mel(j, i) = <psi_{k+q,j} | dV_tot * psi_{k,i}>
      !
      CALL ZGEMM('C', 'N', nbnd, nbnd, npwx*npol, (1.d0, 0.d0), &
                evq, npwx*npol, dvpsi, npwx*npol, (0.d0, 0.d0), &
                mel, nbnd)
      CALL mp_sum(mel, intra_bgrp_comm)
      !
      ! Set ps to be the coefficients of the active-space wavefunction response
      !
      ps = (0.d0, 0.d0)
      IF (lgauss) THEN
         DO ibnd = 1, nbnd_occ(ikk)
            DO jbnd = 1, nbnd
               deltae = et(jbnd, ikq) - et(ibnd, ikk)
               theta = wgauss(deltae / degauss, 0)
               ps(jbnd, ibnd) = mel(jbnd, ibnd) * chi_active(jbnd, ibnd, ik) * theta
            ENDDO
         ENDDO
      ELSE
         DO ibnd = 1, nbnd_occ(ikk)
            DO jbnd = nbnd_occ(ikq) + 1, nbnd
               ps(jbnd, ibnd) = mel(jbnd, ibnd) * chi_active(jbnd, ibnd, ik)
            ENDDO
         ENDDO
      ENDIF ! lgauss
      !
      ! dpsi(:, i) = dpsi(:, i) - evq(:, j) * ps(j, i)
      !
      CALL ZGEMM('N', 'N', npwx*npol, nbnd_occ(ikk), nbnd, &
         (-1.d0, 0.d0), evq, npwx*npol, ps, nbnd, &
         ( 1.d0, 0.d0), dpsi, npwx*npol)
      !
      DEALLOCATE(mel)
      DEALLOCATE(ps)
      !
      CALL stop_clock('cdfpt_subtract')
      !
   END SUBROUTINE cdfpt_subtract_active_wfc
   !----------------------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------------------
   SUBROUTINE cdfpt_chi_active_by_bands(active_bands_min, active_bands_max)
   !----------------------------------------------------------------------------------------
   !! Compute the active-space susceptibility for active space defined by band indices
   !! from active_bands_min to active_bands_max.
   !----------------------------------------------------------------------------------------
      USE kinds,                ONLY : DP
      USE wvfct,                ONLY : nbnd
      USE qpoint,               ONLY : ikks, ikqs, nksq
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: active_bands_min
      !! Minimum band index in the active space
      INTEGER, INTENT(IN) :: active_bands_max
      !! Maximum band index in the active space
      !
      INTEGER :: ik, ikk, ikq
      !! k and k+q point indices
      INTEGER :: ibnd, jbnd
      !! Band indices
      REAL(DP) :: factor
      !! Susceptibility factor
      !
      IF (.NOT. ALLOCATED(chi_active)) THEN
         CALL errore('cdfpt_chi_active_by_bands', 'chi_active not allocated', 1)
      ENDIF
      IF (ANY(SHAPE(chi_active) /= (/ nbnd, nbnd, nksq /))) THEN
         CALL errore('cdfpt_chi_active_by_bands', 'chi_active has wrong size', 1)
      ENDIF
      !
      chi_active = (0.d0, 0.d0)
      !
      DO ik = 1, nksq
         !
         ikk = ikks(ik)
         ikq = ikqs(ik)
         !
         DO ibnd = active_bands_min, active_bands_max
            DO jbnd = active_bands_min, active_bands_max
               factor = susceptibility_factor(ibnd, jbnd, ikk, ikq)
               chi_active(jbnd, ibnd, ik) = CMPLX(factor, 0.d0, KIND = DP)
            ENDDO
         ENDDO
         !
      ENDDO
      !
   !----------------------------------------------------------------------------------------
   END SUBROUTINE cdfpt_chi_active_by_bands
   !----------------------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------------------
   SUBROUTINE cdfpt_allocate()
      USE wvfct,                ONLY : nbnd
      USE qpoint,               ONLY : nksq
      IMPLICIT NONE
      ALLOCATE(chi_active(nbnd, nbnd, nksq))
   END SUBROUTINE cdfpt_allocate
   !----------------------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------------------
   SUBROUTINE cdfpt_deallocate()
      IMPLICIT NONE
      DEALLOCATE(chi_active)
   END SUBROUTINE cdfpt_deallocate
   !----------------------------------------------------------------------------------------
   !
   !----------------------------------------------------------------------------------------
   FUNCTION susceptibility_factor(ibnd, jbnd, ikk, ikq) RESULT(factor)
   !----------------------------------------------------------------------------------------
   !! Compute the susceptibility factor (f_mk+q - f_nk) / (e_mk+q - e_nk)
   !! Treat the metallic and insulating cases seprately.
   !----------------------------------------------------------------------------------------
      !
      USE kinds,                ONLY : DP
      USE ener,                 ONLY : ef
      USE klist,                ONLY : lgauss, degauss, ngauss
      USE wvfct,                ONLY : et
      USE control_lr,           ONLY : nbnd_occ
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: ibnd, jbnd
      !! Band indices
      INTEGER, INTENT(IN) :: ikk, ikq
      !! k and k+q point indices
      !
      REAL(DP) :: factor
      !! Output: Susceptibility factor
      !
      REAL(DP) :: ek, ekq
      !! Energy of the k and k+q points
      REAL(DP) :: fk, fkq
      !! Occupation function at k and k+q
      REAL(DP), EXTERNAL :: wgauss
      !! Occupation function
      REAL(DP), EXTERNAL :: w0gauss
      !! Derivative of the occupation function
      !
      ek  = et(ibnd, ikk)
      ekq = et(jbnd, ikq)
      !
      IF (lgauss) THEN
         !
         ! Metallic case
         !
         fk  = wgauss((ef - ek ) / degauss, ngauss)
         fkq = wgauss((ef - ekq) / degauss, ngauss)
         !
         IF (ABS(ekq - ek) > 1.d-5) THEN
            factor = (fkq - fk) / (ekq - ek)
         ELSE
            !
            ! This is a special case, where the energy difference is very small
            ! Take the limit of 0/0, which gives the derivative of the occupation function
            !factor = -w0gauss((ef - (ek + ekq) / 2) / degauss, ngauss) / degauss
            factor = -w0gauss((ef - ek) / degauss, ngauss) / degauss
            !
         ENDIF
      ELSE
         !
         ! Insulating case
         !
         IF ((jbnd > nbnd_occ(ikq)) .AND. (ibnd <= nbnd_occ(ikk))) THEN
            ! fkq = 0, fk = 1
            factor = -1.d0 / (ekq - ek)
         ELSEIF ((jbnd <= nbnd_occ(ikq)) .AND. (ibnd > nbnd_occ(ikk))) THEN
            ! fkq = 1, fk = 0
            factor = 1.d0 / (ekq - ek)
         ELSE
            factor = 0.d0
         ENDIF
         !
      ENDIF ! lgauss
      !
   END FUNCTION susceptibility_factor
   !
END MODULE constrained_dfpt
