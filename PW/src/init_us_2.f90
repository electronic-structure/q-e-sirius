!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_us_2 (npw_, igk_, q_, vkb_)
  !----------------------------------------------------------------------
  !
  !   Calculates beta functions (Kleinman-Bylander projectors), with
  !   structure factor, for all atoms, in reciprocal space. On input:
  !      npw_       : number of PWs 
  !      igk_(npw_) : indices of G in the list of q+G vectors
  !      q_(3)      : q vector (2pi/a units)
  !  On output:
  !      vkb_(npwx,nkb) : beta functions (npw_ <= npwx)
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,  ONLY : tpiba
  USE constants,  ONLY : tpi
  USE gvect,      ONLY : eigts1, eigts2, eigts3, mill, g
  USE wvfct,      ONLY : npwx
  USE us,         ONLY : nqx, dq, tab, tab_d2y, spline_ps
  USE m_gth,      ONLY : mk_ffnl_gth
  USE splinelib
  USE uspp,       ONLY : nkb, nhtol, nhtolm, indv
  USE uspp_param, ONLY : upf, lmaxkb, nhm, nh
  USE constants,    ONLY : fpi
  USE cell_base,    ONLY : omega, bg
  USE mod_sirius
  !
  implicit none
  !
  INTEGER, INTENT (IN) :: npw_, igk_ (npw_)
  REAL(dp), INTENT(IN) :: q_(3)
  COMPLEX(dp), INTENT(OUT) :: vkb_ (npwx, nkb)
  !
  !     Local variables
  !
  integer :: i0,i1,i2,i3, ig, ig_orig, lm, na, nt, nb, ih, jkb

  real(DP) :: px, ux, vx, wx, arg
  real(DP), allocatable :: gk (:,:), qg (:), vq (:), ylm (:,:), vkb1(:,:)

  complex(DP) :: phase, pref
  complex(DP), allocatable :: sk(:)

  real(DP), allocatable :: xdata(:)
  integer :: iq
  integer, allocatable :: gvl(:,:)
  real(8) vkl(3),t1

  ! cache blocking parameters
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: iblock, numblock, realblocksize
  !
  if (lmaxkb.lt.0) return
  if (use_sirius.and.use_sirius_beta_projectors) then
    call invert_mtrx(bg, bg_inv)
    vkl = matmul(bg_inv, q_)
    allocate(gvl(3, npw_))
    do ig = 1, npw_
      gvl(:,ig) = mill(:, igk_(ig))
    enddo
    stop 'not implemented'
    ! this has to be done corretly for the case when k-point on QE side and k-point
    ! on SIIRUS side are not located on the same MPI rank
    !call sirius_get_beta_projectors_by_kp(kset_id, vkl(1), npw_, gvl(1, 1), vkb_(1, 1), npwx, nkb)
    deallocate(gvl)
    return
  endif

  call start_clock ('init_us_2')

!   write(*,'(3i4,i5,3f10.5)') size(tab,1), size(tab,2), size(tab,3), size(vq), q_

  ! setting cache blocking size
  numblock  = (npw_+blocksize-1)/blocksize

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif

!$omp parallel private(vkb1, sk, qg, vq, ylm, gk, ig_orig, &
!$omp                  realblocksize, jkb, px, ux, vx, wx, &
!$omp                  i0, i1, i2, i3, lm, arg, phase, pref)
  !
  allocate (vkb1(blocksize,nhm))
  allocate (  sk(blocksize))
  allocate (  qg(blocksize))
  allocate (  vq(blocksize))
  allocate ( ylm(blocksize, (lmaxkb + 1) **2))
  allocate (  gk(3, blocksize))
  !
!$omp do
  DO iblock = 1, numblock
     !
     realblocksize = MIN(npw_-(iblock-1)*blocksize,blocksize)
     !
     do ig = 1, realblocksize
        ig_orig = (iblock-1)*blocksize+ig
        gk (1,ig) = q_(1) + g(1, igk_(ig_orig) )
        gk (2,ig) = q_(2) + g(2, igk_(ig_orig) )
        gk (3,ig) = q_(3) + g(3, igk_(ig_orig) )
        qg (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
     enddo
     !
     call ylmr2 ((lmaxkb+1)**2, realblocksize, gk, qg, ylm(1:realblocksize,:))
     !
     ! set now qg=|q+G| in atomic units
     !
     do ig = 1, realblocksize
        qg(ig) = sqrt(qg(ig))*tpiba
     enddo

     ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
     jkb = 0
     do nt = 1, ntyp
        ! calculate beta in G-space using an interpolation table f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
        do nb = 1, upf(nt)%nbeta
           if ( upf(nt)%is_gth ) then
              call mk_ffnl_gth( nt, nb, realblocksize, qg, vq )
           else
              do ig = 1, realblocksize
                 if (spline_ps) then
                    vq(ig) = splint(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), qg(ig))
                 else
                    px = qg (ig) / dq - int (qg (ig) / dq)
                    ux = 1.d0 - px
                    vx = 2.d0 - px
                    wx = 3.d0 - px
                    i0 = INT( qg (ig) / dq ) + 1
                    i1 = i0 + 1
                    i2 = i0 + 2
                    i3 = i0 + 3
                    if (use_sirius.and.use_sirius_radial_integrals_beta.and.sirius_context_initialized(sctx)) then
                      vq(ig) = sirius_get_radial_integral(sctx, atom_type(nt)%label, string("beta"), qg(ig), nb)
                      vq(ig) = vq(ig) * fpi / sqrt(omega)
                    else
                      vq (ig) = tab (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                                tab (i1, nb, nt) * px * vx * wx / 2.d0 - &
                                tab (i2, nb, nt) * px * ux * wx / 2.d0 + &
                                tab (i3, nb, nt) * px * ux * vx / 6.d0
                    endif
                 endif
              enddo
           endif
           ! add spherical harmonic part  (Y_lm(q)*f_l(q)) 
           do ih = 1, nh (nt)
              if (nb.eq.indv (ih, nt) ) then
                 !l = nhtol (ih, nt)
                 lm =nhtolm (ih, nt)
                 do ig = 1, realblocksize
                    vkb1 (ig,ih) = ylm (ig, lm) * vq (ig)
                 enddo
              endif
           enddo
        enddo
        !
        ! vkb1 contains all betas including angular part for type nt
        ! now add the structure factor and factor (-i)^l
        !
        do na = 1, nat
           ! ordering: first all betas for atoms of type 1
           !           then  all betas for atoms of type 2  and so on
           if (ityp (na) .eq.nt) then
              arg = (q_(1) * tau (1, na) + &
                     q_(2) * tau (2, na) + &
                     q_(3) * tau (3, na) ) * tpi
              phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
              do ig = 1, realblocksize
                 ig_orig = (iblock-1)*blocksize+ig
                 sk (ig) = eigts1 (mill(1,igk_(ig_orig)), na) * &
                           eigts2 (mill(2,igk_(ig_orig)), na) * &
                           eigts3 (mill(3,igk_(ig_orig)), na)
              enddo
              do ih = 1, nh (nt)
                 jkb = jkb + 1
                 pref = (0.d0, -1.d0) **nhtol (ih, nt) * phase
                 do ig = 1, realblocksize
                    vkb_((iblock-1)*blocksize+ig, jkb) = vkb1 (ig,ih) * sk (ig) * pref
                 enddo
                 ! clean up garbage in the last block
                 if (iblock.eq.numblock) then
                    do ig = npw_+1, npwx
                       vkb_(ig, jkb) = (0.0_dp, 0.0_dp)
                    enddo
                 endif
              enddo
           endif
        enddo
     enddo
  enddo
!$omp end do nowait
  deallocate (gk)
  deallocate (ylm)
  deallocate (vq)
  deallocate (qg)
  deallocate (sk)
  deallocate (vkb1)
!$omp end parallel
  !
  if (spline_ps) deallocate(xdata)
  !
  call stop_clock ('init_us_2')
  return
end subroutine init_us_2

