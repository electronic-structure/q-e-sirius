!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_vloc()
  !----------------------------------------------------------------------
  !
  !    This routine computes the fourier coefficient of the local
  !    potential vloc(ig,it) for each type of atom
  !
  USE atom,       ONLY : msh, rgrid
  USE m_gth,      ONLY : vloc_gth
  USE kinds,      ONLY : dp
  USE uspp_param, ONLY : upf
  USE ions_base,  ONLY : ntyp => nsp
  USE cell_base,  ONLY : omega, tpiba2
  USE vlocal,     ONLY : vloc
  USE gvect,      ONLY : ngl, gl, ngm, mill, igtongl
  USE Coul_cut_2D, ONLY : do_cutoff_2D, cutoff_lr_Vloc
  use mp_bands, only : intra_bgrp_comm
  use mod_sirius
  !
  implicit none
  !
  integer :: nt, i
  real(8), allocatable :: tmp(:)
  ! counter on atomic types
  !
  call start_clock ('init_vloc')
  call sirius_start_timer(string("qe|init_vloc"))
  vloc(:,:) = 0._dp
  do nt = 1, ntyp
     !
     ! compute V_loc(G) for a given type of atom
     !
     IF ( .NOT. ASSOCIATED ( upf(nt)%vloc ) ) THEN
        !
        IF ( upf(nt)%is_gth ) THEN
           !
           ! special case: GTH pseudopotential
           !
           call vloc_gth (nt, upf(nt)%zp, tpiba2, ngl, gl, omega, vloc (1, nt) )
           !
        ELSE
           !
           ! special case: pseudopotential is coulomb 1/r potential
           !
           call vloc_coul (upf(nt)%zp, tpiba2, ngl, gl, omega, vloc (1, nt) )
           !
        END IF
        !
     ELSE
        !
        ! normal case
        !
        if (use_sirius.and.use_sirius_vloc) then
          allocate(tmp(ngm))
          call sirius_get_pw_coeffs_real(sctx, atom_type(nt)%label, string("vloc"), tmp(1), ngm, mill(1, 1), intra_bgrp_comm)
          do i = 1, ngm
            vloc(igtongl(i), nt) = tmp(i) * 2 ! convert to Ry
          enddo
          deallocate(tmp)
        else
          call vloc_of_g (rgrid(nt)%mesh, msh (nt), rgrid(nt)%rab, rgrid(nt)%r, &
              upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngl, gl, omega, vloc (1, nt) )
        endif
        !
     END IF
  enddo
  ! in 2D calculations the long range part of vloc(g) (erf/r part)
  ! was not re-added in g-space because everything is caclulated in
  ! radial coordinates, which is not compatible with 2D cutoff. 
  ! It will be re-added each time vloc(g) is used in the code. 
  ! Here, this cutoff long-range part of vloc(g) is computed only once
  ! by the routine below and stored
  IF (do_cutoff_2D) call cutoff_lr_Vloc() 
  call sirius_stop_timer(string("qe|init_vloc"))
  call stop_clock ('init_vloc')
  return
end subroutine init_vloc

