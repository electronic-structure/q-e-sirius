!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine force_cc (forcecc)
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : rgrid, msh
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, gg, ngl, gl, igtongl, mill
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  use mod_sirius
  !
  implicit none
  !
  !   first the dummy variable
  !
  real(DP) :: forcecc (3, nat)
  ! output: the local forces on atoms

  integer :: ipol, ig, ir, nt, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms


  real(DP), allocatable :: vxc (:,:), rhocg (:), rho_core_g(:)
  ! exchange-correlation potential
  ! radial fourier transform of rho core
  real(DP)  ::  arg, fact
  
  if (use_sirius.and.use_sirius_forces) then
    call sirius_get_forces(gs_handler, string("core"), forcecc(1, 1))
    forcecc = forcecc * 2 ! convert to Ry
    return
  endif
  !
  forcecc(:,:) = 0.d0
  if ( ANY ( upf(1:ntyp)%nlcc ) ) go to 15
  return
  !
15 continue
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  !
  ! recalculate the exchange-correlation potential
  !
  allocate ( vxc(dfftp%nnr,nspin) )
  !
  call v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
  !
  psic=(0.0_DP,0.0_DP)
  if (nspin == 1 .or. nspin == 4) then
     do ir = 1, dfftp%nnr
        psic (ir) = vxc (ir, 1)
     enddo
  else
     do ir = 1, dfftp%nnr
        psic (ir) = 0.5d0 * (vxc (ir, 1) + vxc (ir, 2) )
     enddo
  endif
  deallocate (vxc)
  CALL fwfft ('Rho', psic, dfftp)
  !
  ! psic contains now Vxc(G)
  !
  allocate ( rhocg(ngl) )
  allocate(rho_core_g(ngm))
  !
  ! core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then

        if (use_sirius.and.use_sirius_rho_core) then
          call sirius_get_pw_coeffs_real(sctx, atom_type(nt)%label, string("rhoc"), rho_core_g(1),&
                                        &ngm, mill(1, 1), intra_bgrp_comm)
        else
          call drhoc (ngl, gl, omega, tpiba2, msh(nt), rgrid(nt)%r,&
                      rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
          rho_core_g(:) = rhocg(igtongl(:))
        endif
        do na = 1, nat
           if (nt.eq.ityp (na) ) then
              do ig = gstart, ngm
                 arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
                      + g (3, ig) * tau (3, na) ) * tpi
                 do ipol = 1, 3
                    forcecc (ipol, na) = forcecc (ipol, na) + tpiba * omega * &
                         rho_core_g(ig) * CONJG(psic (dfftp%nl (ig) ) ) * &
                         CMPLX( sin (arg), cos (arg) ,kind=DP) * g (ipol, ig) * fact
                 enddo
              enddo
           endif
        enddo
     endif
  enddo
  !
  call mp_sum(  forcecc, intra_bgrp_comm )
  !
  deallocate (rhocg)
  deallocate(rho_core_g)
  !
  return
end subroutine force_cc
