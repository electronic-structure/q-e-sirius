subroutine get_band_energies_from_sirius
  !
  use wvfct,    only : nbnd, et
  use klist,    only : kset_id, nkstot, nks
  use lsda_mod, only : nspin
  use sirius
  !
  implicit none
  !
  integer, external :: global_kpoint_index
  !
  real(8), allocatable :: band_e(:,:), tmp(:)
  integer :: ik, nk, nb, nfv
  
  allocate(band_e(nbnd, nkstot))
  
  ! get band energies
  if (nspin.ne.2) then
    ! non-magnetic or non-collinear case
    do ik = 1, nkstot
      call sirius_get_band_energies(kset_id, ik, band_e(1, ik), nbnd)
    end do
  else
    ! collinear magnetic case
    nb = nbnd * 2
    nk = nkstot / 2
    allocate(tmp(nb))
    ! get band energies
    do ik = 1, nk
      call sirius_get_band_energies(kset_id, ik, tmp(1), nb)
      band_e(1 : nbnd, ik) = tmp(1 : nbnd)
      band_e(1 : nbnd, nk + ik) = tmp(nbnd + 1 : 2 * nbnd)
    end do
  
    deallocate(tmp)
  
  endif
  
  ! convert to Ry
  do ik = 1, nks
    et(:, ik) = 2.d0 * band_e(:, global_kpoint_index(nkstot, ik))
  enddo
  
  deallocate(band_e)

end subroutine get_band_energies_from_sirius

subroutine put_band_occupancies_to_sirius
  !
  use wvfct,    only : nbnd, wg
  use klist,    only : kset_id, nkstot, nks, wk
  use lsda_mod, only : nspin
  use mp_pools, only : inter_pool_comm
  use parallel_include
  use sirius
  !
  implicit none
  !
  integer, external :: global_kpoint_index
  !
  real(8), allocatable :: bnd_occ(:, :), tmp(:)
  real(8) :: maxocc
  integer :: ik, ierr, nk, nb
  
  ! compute occupancies
  allocate(bnd_occ(nbnd, nkstot))
  bnd_occ = 0.d0
  ! define a maximum band occupancy (2 in case of spin-unpolarized, 1 in case of spin-polarized)
  maxocc = 2.d0
  if (nspin.gt.1) then
    maxocc = 1.d0
  endif
  do ik = 1, nks
    bnd_occ(:, global_kpoint_index(nkstot, ik)) = maxocc * wg(:, ik) / wk(ik)
  enddo
  call mpi_allreduce(MPI_IN_PLACE, bnd_occ(1, 1), nbnd * nkstot, MPI_DOUBLE, MPI_SUM, inter_pool_comm, ierr)
  
  if (nspin.ne.2) then
    ! set band occupancies
    do ik = 1, nkstot
      call sirius_set_band_occupancies(kset_id, ik, bnd_occ(1, ik), nbnd)
    enddo
  else 
    nk = nkstot / 2
    nb = nbnd * 2
    allocate(tmp(nb))
    do ik = 1, nk
      tmp(1 : nbnd) = bnd_occ(1 : nbnd, ik)
      tmp(nbnd + 1 : nb) = bnd_occ(1 : nbnd, ik + nk)
      call sirius_set_band_occupancies(kset_id, ik, tmp(1), nb)
    enddo
    deallocate(tmp)
  endif
  
  deallocate(bnd_occ)

end subroutine put_band_occupancies_to_sirius

subroutine get_density_from_sirius
  !
  use scf,        only : rho
  use gvect,      only : mill, ngm
  use mp_bands,   only : intra_bgrp_comm
  use lsda_mod,   only : nspin
  use ions_base,  only : nat, nsp, ityp
  use uspp_param, only : nhm, nh
  use sirius
  !
  implicit none
  !
  complex(8), allocatable :: dens_mtrx(:,:,:)
  integer iat, ig, ih, jh, ijh, na, ispn
  complex(8) z1, z2
  
  ! get rho(G)
  call sirius_get_pw_coeffs(c_str("rho"), rho%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
  
  if (nspin.eq.2) then
    call sirius_get_pw_coeffs(c_str("magz"), rho%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
    ! convert to rho_{up}, rho_{dn}
    do ig = 1, ngm
      z1 = rho%of_g(ig, 1)
      z2 = rho%of_g(ig, 2)
      rho%of_g(ig, 1) = 0.5 * (z1 + z2)
      rho%of_g(ig, 2) = 0.5 * (z1 - z2)
    enddo
  endif
  if (nspin.eq.4) then
    call sirius_get_pw_coeffs(c_str("magx"), rho%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_get_pw_coeffs(c_str("magy"), rho%of_g(1, 3), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_get_pw_coeffs(c_str("magz"), rho%of_g(1, 4), ngm, mill(1, 1), intra_bgrp_comm)
  endif
  ! get density matrix
  ! complex density matrix in SIRIUS has at maximum three components
  allocate(dens_mtrx(nhm, nhm, 3))
  do iat = 1, nsp
    do na = 1, nat
      if (ityp(na).eq.iat.and.allocated(rho%bec)) then
        rho%bec(:, na, :) = 0.d0
        call sirius_get_density_matrix(na, dens_mtrx(1, 1, 1), nhm)

        ijh = 0
        do ih = 1, nh(iat)
          do jh = ih, nh(iat)
            ijh = ijh + 1
            if (nspin.le.2) then
              do ispn = 1, nspin
                rho%bec(ijh, na, ispn) = dreal(dens_mtrx(ih, jh, ispn))
              enddo
            endif
            if (nspin.eq.4) then
              rho%bec(ijh, na, 1) = dreal(dens_mtrx(ih, jh, 1) + dens_mtrx(ih, jh, 2))
              rho%bec(ijh, na, 4) = dreal(dens_mtrx(ih, jh, 1) - dens_mtrx(ih, jh, 2))
              rho%bec(ijh, na, 2) = 2.d0 * dreal(dens_mtrx(ih, jh, 3))
              rho%bec(ijh, na, 3) = -2.d0 * dimag(dens_mtrx(ih, jh, 3))
            endif
            ! off-diagonal elements have a weight of 2
            if (ih.ne.jh) then
              do ispn = 1, nspin
                rho%bec(ijh, na, ispn) = rho%bec(ijh, na, ispn) * 2.d0
              enddo
            endif
          enddo
        enddo
      endif
    enddo
  enddo
  deallocate(dens_mtrx)
end subroutine get_density_from_sirius

!subroutine put_density_to_sirius
!  !
!  use scf,        only : rho
!  use gvect,      only : mill, ngm
!  use mp_bands,   only : intra_bgrp_comm
!  use lsda_mod,   only : nspin
!  use ions_base,  only : nat, nsp, ityp
!  use uspp_param, only : nhm, nh
!  use sirius
!  implicit none
!  !
!  complex(8), allocatable :: rho_tot(:), mag(:)
!  complex(8), allocatable :: dens_mtrx(:,:,:)
!  integer iat, ig, ih, jh, ijh, na, ispn
!  real(8) :: fact
!  !
!  if (nspin.eq.1.or.nspin.eq.4) then
!    call sirius_set_pw_coeffs(c_str("rho"), rho%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
!  endif
!
!  if (nspin.eq.2) then
!    allocate(rho_tot(ngm))
!    allocate(mag(ngm))
!    do ig = 1, ngm
!      rho_tot(ig) = rho%of_g(ig, 1) + rho%of_g(ig, 2)
!      mag(ig) = rho%of_g(ig, 1) - rho%of_g(ig, 2)
!    enddo
!    call sirius_set_pw_coeffs(c_str("rho"), rho_tot(1), ngm, mill(1, 1), intra_bgrp_comm)
!    call sirius_set_pw_coeffs(c_str("magz"), mag(1), ngm, mill(1, 1), intra_bgrp_comm)
!    deallocate(rho_tot)
!    deallocate(mag)
!  endif
!
!  if (nspin.eq.4) then
!    call sirius_set_pw_coeffs(c_str("magx"), rho%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
!    call sirius_set_pw_coeffs(c_str("magy"), rho%of_g(1, 3), ngm, mill(1, 1), intra_bgrp_comm)
!    call sirius_set_pw_coeffs(c_str("magz"), rho%of_g(1, 4), ngm, mill(1, 1), intra_bgrp_comm)
!  endif
!
!  !!== ! set density matrix
!  !!== ! complex density matrix in SIRIUS has at maximum three components
!  !!== allocate(dens_mtrx(nhm, nhm, 3))
!  !!== do iat = 1, nsp
!  !!==   do na = 1, nat
!  !!==     if (ityp(na).eq.iat.and.allocated(rho%bec)) then
!  !!==       dens_mtrx = (0.d0, 0.d0)
!  !!==       ijh = 0
!  !!==       do ih = 1, nh(iat)
!  !!==         do jh = ih, nh(iat)
!  !!==           ijh = ijh + 1
!  !!==           ! off-diagonal elements have a weight of 2
!  !!==           if (ih.ne.jh) then
!  !!==             fact = 0.5d0
!  !!==           else
!  !!==             fact = 1.d0
!  !!==           endif
!  !!==           if (nspin.le.2) then
!  !!==             do ispn = 1, nspin
!  !!==               dens_mtrx(ih, jh, ispn) = fact * rho%bec(ijh, na, ispn)
!  !!==               dens_mtrx(jh, ih, ispn) = fact * rho%bec(ijh, na, ispn)
!  !!==             enddo
!  !!==           endif
!  !!==           if (nspin.eq.4) then
!  !!==             ! 0.5 * (rho + mz)
!  !!==             dens_mtrx(ih, jh, 1) = fact * 0.5 * (rho%bec(ijh, na, 1) + rho%bec(ijh, na, 4))
!  !!==             dens_mtrx(jh, ih, 1) = fact * 0.5 * (rho%bec(ijh, na, 1) + rho%bec(ijh, na, 4))
!  !!==             ! 0.5 * (rho - mz)
!  !!==             dens_mtrx(ih, jh, 2) = fact * 0.5 * (rho%bec(ijh, na, 1) - rho%bec(ijh, na, 4))
!  !!==             dens_mtrx(jh, ih, 2) = fact * 0.5 * (rho%bec(ijh, na, 1) - rho%bec(ijh, na, 4))
!  !!==             ! 0.5 * (mx - I * my)
!  !!==             dens_mtrx(ih, jh, 3) = fact * 0.5 * dcmplx(rho%bec(ijh, na, 2), -rho%bec(ijh, na, 3))
!  !!==             dens_mtrx(jh, ih, 3) = fact * 0.5 * dcmplx(rho%bec(ijh, na, 2), -rho%bec(ijh, na, 3))
!  !!==           endif
!  !!==         enddo
!  !!==       enddo
!  !!==       call sirius_set_density_matrix(na, dens_mtrx(1, 1, 1), nhm)
!  !!==     endif
!  !!==   enddo
!  !!== enddo
!  !!== deallocate(dens_mtrx)
!
!end subroutine put_density_to_sirius

subroutine put_potential_to_sirius
  use scf,                  only : v, vltot
  use gvect,                only : mill, ngm
  use mp_bands,             only : intra_bgrp_comm
  use lsda_mod,             only : nspin
  use noncollin_module,     only : nspin_mag
  USE wavefunctions_module,  ONLY: psic

  use fft_base,             only : dfftp
  USE fft_interfaces,       ONLY : fwfft
  use uspp,                 only : deeq
  use uspp_param,           only : nhm
  use paw_variables,        only : okpaw
  use ions_base,            only : nat
  use sirius
  !
  implicit none
  !
  complex(8), allocatable :: vxcg(:)
  complex(8) :: z1, z2
  integer ig, is, ir, i, ia, j
  character(10) label
  real(8), allocatable :: deeq_tmp(:,:)
  real(8) :: d1,d2
  !
  if (nspin.eq.1.or.nspin.eq.4) then
    ! add local part of the potential and transform to PW domain
    psic(:) = v%of_r(:, 1) + vltot(:)
    call fwfft('Rho', psic, dfftp)
    ! convert to Hartree
    do ig = 1, ngm
      v%of_g(ig, 1) = psic(dfftp%nl(ig)) * 0.5d0
    enddo
    ! set effective potential
    call sirius_set_pw_coeffs(c_str("veff"), v%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
  endif

  if (nspin.eq.2) then
    do is = 1, 2
      ! add local part of the potential and transform to PW domain
      psic(:) = v%of_r(:, is) + vltot(:)
      call fwfft('Rho', psic, dfftp)
      ! convert to Hartree
      do ig = 1, ngm
         v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
      enddo
    enddo

    do ig = 1, ngm
      z1 = v%of_g(ig, 1) 
      z2 = v%of_g(ig, 2)
      v%of_g(ig, 1) = 0.5 * (z1 + z2)
      v%of_g(ig, 2) = 0.5 * (z1 - z2)
    enddo
    ! set effective potential and magnetization
    call sirius_set_pw_coeffs(c_str("veff"),v%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_set_pw_coeffs(c_str("bz"),v%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
  endif

  if (nspin.eq.4) then
    do is = 2, nspin_mag
      psic(:) = v%of_r(:, is)
      call fwfft('Rho', psic, dfftp)
      ! convert to Hartree
      do ig = 1, ngm
        v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
      enddo
      if (is.eq.2) label="bx"
      if (is.eq.3) label="by"
      if (is.eq.4) label="bz"

      call sirius_set_pw_coeffs(c_str(label),v%of_g(1, is), ngm, mill(1, 1), intra_bgrp_comm)
    enddo
  endif

!  ! convert Vxc to plane-wave domain
!  if (nspin.eq.1.or.nspin.eq.4) then
!    do ir = 1, dfftp%nnr
!      psic(ir) = vxc(ir, 1)
!    enddo
!  else
!    do ir = 1, dfftp%nnr
!      psic(ir) = 0.5d0 * (vxc(ir, 1) + vxc(ir, 2))
!    enddo
!  endif
!  call fwfft('Rho', psic, dfftp)
!  allocate(vxcg(ngm))
!  ! convert to Hartree
!  do ig = 1, ngm
!     vxcg(ig) = psic(dfftp%nl(ig)) * 0.5d0
!  end do
!  ! set XC potential
!  call sirius_set_pw_coeffs(c_str("vxc"), vxcg(1), ngm, mill(1, 1), intra_bgrp_comm)
!  deallocate(vxcg)
!
  ! update D-operator matrix
  !call sirius_generate_d_operator_matrix()
  if (okpaw) then
    allocate(deeq_tmp(nhm, nhm))
    !! get D-operator matrix
    !do ia = 1, nat
    !  do is = 1, nspin
    !    call sirius_get_d_operator_matrix(ia, is, deeq(1, 1, ia, is), nhm)
    !  enddo
    !  if (nspin.eq.2) then
    !    do i = 1, nhm
    !      do j = 1, nhm
    !        d1 = deeq(i, j, ia, 1)
    !        d2 = deeq(i, j, ia, 2)
    !        deeq(i, j, ia, 1) = d1 + d2
    !        deeq(i, j, ia, 2) = d1 - d2
    !      enddo
    !    enddo
    !  endif
    !  ! convert to Ry
    !  deeq(:, :, ia, :) = deeq(:, :, ia, :) * 2
    !enddo
    !call add_paw_to_deeq(deeq)
    do ia = 1, nat
      do is = 1, nspin
        if (nspin.eq.2.and.is.eq.1) then
          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) + deeq(:, :, ia, 2)) / 2 ! convert to Ha
        endif
        if (nspin.eq.2.and.is.eq.2) then
          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) - deeq(:, :, ia, 2)) / 2 ! convert to Ha
        endif
        if (nspin.eq.1.or.nspin.eq.4) then
          deeq_tmp(:, :) = deeq(:, :, ia, is) / 2 ! convert to Ha
        endif
        call sirius_set_d_operator_matrix(ia, is, deeq_tmp(1, 1), nhm)
      enddo
    enddo
    deallocate(deeq_tmp)
  endif

end subroutine put_potential_to_sirius

!subroutine put_vltot_to_sirius
!  use scf,       only : vltot
!  use gvect, only : mill, ngm
!  use wavefunctions_module, only : psic
!  use fft_interfaces,       only : fwfft, invfft
!  use fft_base,             only : dfftp
!  use mp_bands, only : intra_bgrp_comm
!  use sirius
!  !
!  implicit none
!  !
!  complex(8), allocatable :: vg(:)
!  integer ig
!  !
!  allocate(vg(ngm))
!  psic(:) = vltot(:)
!  call fwfft('Rho', psic, dfftp)
!  ! convert to Hartree
!  do ig = 1, ngm
!    vg(ig) = psic(dfftp%nl(ig)) * 0.5d0 ! convert to Ha
!  enddo
!  ! set local potential
!  call sirius_set_pw_coeffs(c_str("vloc"), vg(1), ngm, mill(1, 1), intra_bgrp_comm)
!  deallocate(vg)
!end subroutine put_vltot_to_sirius

!subroutine get_rhoc_from_sirius
!  use uspp_param,only : upf
!  use ener,      only : etxcc
!  use scf,       only : rho_core, rhog_core
!  use control_flags, only : gamma_only
!  use wavefunctions_module, only : psic
!  use gvect, only : mill, ngm
!  use scf, only : rho_core, rhog_core
!  use mp_bands, only : intra_bgrp_comm
!  use ions_base, only : ntyp => nsp
!  use fft_interfaces,only : invfft
!  use fft_base,  only : dfftp
!  use sirius
!  !
!  implicit none
!  
!  etxcc = 0.0d0
!  if (any(upf(1:ntyp)%nlcc)) then
!    call sirius_get_pw_coeffs(c_str("rhoc"), rhog_core(1), ngm, mill(1, 1), intra_bgrp_comm)
!    psic(:) = (0.d0, 0.d0)
!    psic(dfftp%nl(:)) = rhog_core(:)
!    if (gamma_only) psic(dfftp%nlm(:)) = conjg(rhog_core(:))
!    call invfft ('Rho', psic, dfftp)
!    rho_core(:) = psic(:)
!  else 
!    rhog_core(:) = 0.0d0
!    rho_core(:)  = 0.0d0
!  endif
!
!end subroutine get_rhoc_from_sirius

!subroutine set_vloc_sirius
!use sirius
!use gvect, only : ngm, mill, igtongl, ngl
!use vlocal, only : vloc
!use mp_bands, only : intra_bgrp_comm
!use ions_base, only : atm
!integer nt,i
!real(8), allocatable :: tmp(:)
!
!allocate(tmp(ngm))
!vloc(:,:) = 0.d0
!do nt = 1, ntyp
!  call sirius_get_pw_coeffs_real(c_str(atm(nt)), c_str("vloc"), tmp(1), ngm, mill(1, 1), intra_bgrp_comm)
!  do i = 1, ngm
!    vloc(igtongl(i), nt) = tmp(i) * 2 ! convert to Ry
!  enddo
!enddo
!
!deallocate(tmp, tmp1)
!
!call set_vloc_sirius
!CALL setlocal()
!end subroutine

!subroutine get_vloc_from_sirius
!  use wavefunctions_module, only : psic
!  use gvect, only : mill, ngm, gg
!  use scf, only: vltot, v_of_0
!  use fft_interfaces, only : fwfft, invfft
!  use fft_base, only : dfftp
!  use constants, only : eps8
!  use control_flags, only : gamma_only
!  use mp_bands, only : intra_bgrp_comm
!  use mp, only : mp_bcast, mp_sum
!  use sirius
!  !
!  implicit none
!  !
!  complex(8), allocatable :: vpw(:)
!  allocate(vpw(ngm))
!  call sirius_get_pw_coeffs(c_str("vloc"), vpw(1), ngm, mill(1, 1), intra_bgrp_comm)
!  psic(:) = 0.d0
!  psic(dfftp%nl(:)) = vpw(:)
!  if (gamma_only) psic(dfftp%nlm(:)) = conjg(vpw(:))
!  call invfft('Rho', psic, dfftp)
!  vltot(:) = dble(psic(:)) * 2 ! convert to Ry
!  v_of_0=0.d0
!  IF (gg(1) < eps8) v_of_0 = dble(vpw(1))
!  !
!  call mp_sum(v_of_0, intra_bgrp_comm)
!  deallocate(vpw)
!
!end subroutine get_vloc_from_sirius

