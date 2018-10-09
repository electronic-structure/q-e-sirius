!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine force_corr (forcescc)
  !-----------------------------------------------------------------------
  !   This routine calculates the force term vanishing at full
  !     self-consistency. It follows the suggestion of Chan-Bohnen-Ho
  !     (PRB 47, 4771 (1993)). The true charge density is approximated
  !     by means of a free atom superposition.
  !     (alessio f.)
  ! Uses superposition of atomic charges contained in the array rho_at
  ! and read from pseudopotential files
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : msh, rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, ngl, gl, igtongl
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vnew
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  use mod_sirius
  !
  implicit none
  !
  real(DP) :: forcescc (3, nat)
  !
  real(DP), allocatable :: rhocgnt (:), aux (:)
  ! work space
  real(DP) ::  gx, arg, fact
  ! temp factors
  integer :: ir, isup, isdw, ig, nt, na, ipol, ndm
  ! counters
  !
  ! vnew is V_out - V_in, psic is the temp space
  !
  if (use_sirius.and.use_sirius_stress) then
    call sirius_get_forces(gs_handler, string("scf_corr"), forcescc(1, 1))
    forcescc = forcescc * 2 ! convert to Ry
    return
  endif
  if (nspin == 1 .or. nspin == 4) then
     psic(:) = vnew%of_r (:, 1)
  else
     isup = 1
     isdw = 2
     psic(:) = (vnew%of_r (:, isup) + vnew%of_r (:, isdw)) * 0.5d0
  end if
  !
  ndm = MAXVAL ( msh(1:ntyp) )
  allocate ( aux(ndm), rhocgnt(ngl) )

  forcescc(:,:) = 0.d0

  CALL fwfft ('Rho', psic, dfftp)

  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if

  do nt = 1, ntyp
     !
     ! Here we compute the G.ne.0 term
     !
     do ig = gstart, ngl
        gx = sqrt (gl (ig) ) * tpiba
        do ir = 1, msh (nt)
           if (rgrid(nt)%r(ir) .lt.1.0d-8) then
              aux (ir) = upf(nt)%rho_at (ir)
           else
              aux (ir) = upf(nt)%rho_at (ir) * &
                         sin(gx*rgrid(nt)%r(ir)) / (rgrid(nt)%r(ir)*gx)
           endif
        enddo
        call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgnt (ig) )
     enddo
     do na = 1, nat
        if (nt.eq.ityp (na) ) then
           do ig = gstart, ngm
              arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
                   + g (3, ig) * tau (3, na) ) * tpi
              do ipol = 1, 3
                 forcescc (ipol, na) = forcescc (ipol, na) + fact * &
                      rhocgnt (igtongl(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g(ipol,ig) * tpiba * CONJG(psic(dfftp%nl(ig)))
              enddo
           enddo
        endif
     enddo
  enddo
  !
  call mp_sum(  forcescc, intra_bgrp_comm )
  !
  deallocate ( aux, rhocgnt )

  return
end subroutine force_corr

