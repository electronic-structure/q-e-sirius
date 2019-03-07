!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine stres_har (sigmahar)
  !----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : e2, fpi
  USE cell_base, ONLY: omega, tpiba2
  USE ener,      ONLY: ehart
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft
  USE gvect,     ONLY: ngm, gstart, g, gg
  USE scf,       ONLY: rho
  USE control_flags,        ONLY: gamma_only
  USE wavefunctions, ONLY : psic
  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum
  USE Coul_cut_2D,  ONLY: do_cutoff_2D, cutoff_stres_sigmahar

  implicit none
  !
  real(DP) :: sigmahar (3, 3), shart, g2
  real(DP), parameter :: eps = 1.d-8
  integer :: ig, l, m

  sigmahar(:,:) = 0.d0
  psic (:) = CMPLX (rho%of_r(:,1), KIND=dp)

  CALL fwfft ('Rho', psic, dfftp)
  ! psic contains now the charge density in G space
  ! the  G=0 component is not computed
  IF (do_cutoff_2D) THEN  
    call cutoff_stres_sigmahar(psic, sigmahar)
  ELSE
  do ig = gstart, ngm
     g2 = gg (ig) * tpiba2
     shart = psic (dfftp%nl (ig) ) * CONJG(psic (dfftp%nl (ig) ) ) / g2
     do l = 1, 3
        do m = 1, l
           sigmahar (l, m) = sigmahar (l, m) + shart * tpiba2 * 2 * &
                g (l, ig) * g (m, ig) / g2
        enddo
     enddo
  enddo
  ENDIF 
  !
  call mp_sum(  sigmahar, intra_bgrp_comm )
  !
  if (gamma_only) then
     sigmahar(:,:) =       fpi * e2 * sigmahar(:,:)
  else
     sigmahar(:,:) = 0.5d0 * fpi * e2 * sigmahar(:,:)
  end if
  do l = 1, 3
     sigmahar (l, l) = sigmahar (l, l) - ehart / omega
  enddo
  do l = 1, 3
     do m = 1, l - 1
        sigmahar (m, l) = sigmahar (l, m)
     enddo
  enddo

  sigmahar(:,:) = -sigmahar(:,:)

  return
end subroutine stres_har

