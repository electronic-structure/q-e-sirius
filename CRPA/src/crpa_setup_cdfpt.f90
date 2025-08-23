!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_setup_cdfpt()
   !-----------------------------------------------------------------------
   !! Setup for constrained DFPT.
   !! See LR_Modules/constrained_dfpt.f90 for details.
   !
   USE kinds,                ONLY : DP
   USE crpacom,              ONLY : active_space, active_bands_min, active_bands_max
   USE constrained_dfpt,     ONLY : lcdfpt, chi_active, cdfpt_allocate, &
                                    cdfpt_chi_active_by_bands
   !
   IMPLICIT NONE
   !
   CALL cdfpt_allocate()
   !
   IF (active_space == 'bands') THEN
      CALL cdfpt_chi_active_by_bands(active_bands_min, active_bands_max)
   ELSEIF (active_space == 'wannier') THEN
      ! User must set chi_active
      CALL errore('crpa_setup_cdfpt', 'wannier not yet implemented', 1)
   ELSE
      CALL errore('crpa_setup_cdfpt', 'Invalid active_space', 1)
   ENDIF
   !
END SUBROUTINE crpa_setup_cdfpt
