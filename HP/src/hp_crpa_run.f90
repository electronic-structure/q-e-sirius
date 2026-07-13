!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hp_crpa_run()
  !-----------------------------------------------------------------------
  !
  ! Driver of the CRPA calculation (hp.x with method = 'CRPA').
  ! STUB: the solver is ported in a later merge stage; for now this
  ! routine only echoes the parsed &CRPA input and returns.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE control_flags,    ONLY : isolve
  USE control_lr,       ONLY : lrpa, tr2_ph, conv_thr_nscf, niter_ph, nmix_ph
  USE qpoint,           ONLY : nq1, nq2, nq3, nqs, x_q
  USE crpacom,          ONLY : wannier_seedname, active_bands_min, active_bands_max, &
                               wan_dist_thr, wan_dist_thr_mels, folder_wan_Rr,       &
                               write_wan_Rr, w_freq
  USE crpa_pert,        ONLY : pert_basis
  USE crpa_qpoints,     ONLY : qplot
  USE constrained_dfpt, ONLY : cdfpt, cdfpt_active_space, cdfpt_bands_min, cdfpt_bands_max
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, ipol
  !
  WRITE(stdout, '(/,5x,69("="))')
  WRITE(stdout, '(/,25x,"CRPA calculation",/)')
  WRITE(stdout, '(5x,69("="))')
  !
  ! Echo the parsed &CRPA input
  !
  WRITE(stdout, '(/,5x,"Perturbation basis: ",a)') TRIM(pert_basis)
  IF (TRIM(pert_basis) == 'bands') THEN
     WRITE(stdout, '(5x,5x,"Perturbation bands: ",i5," - ",i5)') active_bands_min, active_bands_max
  ELSEIF (TRIM(pert_basis) == 'wannier') THEN
     WRITE(stdout, '(5x,5x,"Wannier seedname:   ",a)') TRIM(wannier_seedname)
     WRITE(stdout, '(5x,5x,"folder_wan_Rr:      ",a)') TRIM(folder_wan_Rr)
     WRITE(stdout, '(5x,5x,"write_wan_Rr:       ",l1)') write_wan_Rr
     WRITE(stdout, '(5x,5x,"wan_dist_thr:       ",1pe10.2," bohr")') wan_dist_thr
     WRITE(stdout, '(5x,5x,"wan_dist_thr_mels:  ",1pe10.2," bohr")') wan_dist_thr_mels
  ENDIF
  !
  WRITE(stdout, '(/,5x,"lrpa:              ",l1)') lrpa
  WRITE(stdout, '(  5x,"tr2:               ",1pe10.2)') tr2_ph
  WRITE(stdout, '(  5x,"conv_thr_nscf:     ",1pe10.2)') conv_thr_nscf
  WRITE(stdout, '(  5x,"niter_max:         ",i5)') niter_ph
  WRITE(stdout, '(  5x,"nmix:              ",i5)') nmix_ph
  WRITE(stdout, '(  5x,"diagonalization:   ",a)') MERGE('david', 'cg   ', isolve == 0)
  WRITE(stdout, '(  5x,"w_freq:          (",1pe10.2,",",1pe10.2," ) Ry")') w_freq
  !
  IF (cdfpt) THEN
     WRITE(stdout, '(/,5x,"Constrained DFPT active space: ",a)') TRIM(cdfpt_active_space)
     IF (TRIM(cdfpt_active_space) == 'bands') THEN
        WRITE(stdout, '(5x,5x,"cDFPT active bands: ",i5," - ",i5)') cdfpt_bands_min, cdfpt_bands_max
     ENDIF
  ENDIF
  !
  IF (qplot) THEN
     WRITE(stdout, '(/,5x,"q points from input (qplot), nqs = ",i5)') nqs
     DO iq = 1, nqs
        WRITE(stdout, '(5x,5x,"q(",i4,") = (",3f12.7," )")') iq, (x_q(ipol,iq), ipol=1,3)
     ENDDO
  ELSE
     WRITE(stdout, '(/,5x,"q mesh: ",i4," x ",i4," x ",i4)') nq1, nq2, nq3
  ENDIF
  !
  WRITE(stdout, '(/,5x,"hp.x method=''CRPA'': solver not yet implemented (merge stage 4 stub)",/)')
  !
  RETURN
  !
END SUBROUTINE hp_crpa_run
