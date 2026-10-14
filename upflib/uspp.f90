!
! Copyright (C) 2004-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined(__CUDA)
#define PINMEM ,PINNED 
#else
#define PINMEM
#endif
!
!
MODULE uspp
  !
  !! Variables for ultrasoft PPs:
  !! 1) Clebsch-Gordan coefficients "ap", auxiliary variables "lpx", "lpl"
  !! 2) Atomic D_lm, Q_lm functions and related indices
  !! FIXME: maybe 1) and 2) should better stay in module uspp_param
  !! 3) Beta_l for the solid in G-space (vkb) and related indices
  !!    (also: beta and dbeta for the Car-Parrinello code)
  !! 4) self-consistent variables becsum, ebecsum, deeq
  !! FIXME: Variables in 3) and 4) should better be defined and allocated
  !!        outside upflib/
  !!
  USE upf_kinds,   ONLY: DP
  USE upf_params,  ONLY: lmaxx, lqmax
  USE upf_spinorb, ONLY: is_spinorbit, fcoef
  IMPLICIT NONE
  PRIVATE
  SAVE
  !
  PUBLIC :: nlx, lpx, lpl, ap, aainit, indv, nhtol, nhtolm, ofsbeta, &
            nkb, nkbus, vkb, dvan, deeq, qq_at, qq_nt, nhtoj, ijtoh, beta, &
            becsum, ebecsum
  PUBLIC :: okvan, nlcc_any
  PUBLIC :: qq_so,   dvan_so,   deeq_nc,   fcoef 
  PUBLIC :: dbeta
  !
  PUBLIC :: allocate_uspp, deallocate_uspp
  !
  ! Vars
  !
  INTEGER, PARAMETER :: &
       nlx  = (lmaxx+1)**2, &! maximum number of combined angular momentum
       mx   = 2*lqmax-1      ! maximum magnetic angular momentum of Q
  !
  INTEGER ::             &! for each pair of combined momenta lm(1),lm(2): 
       lpx(nlx,nlx),     &! maximum combined angular momentum LM
       lpl(nlx,nlx,mx)    ! list of combined angular momenta  LM
  REAL(DP) :: ap(lqmax*lqmax,nlx,nlx)
                          ! Clebsch-Gordan coefficients for spherical harmonics
  !
  INTEGER :: nkb,        &! total number of beta functions, with struct.fact.
             nkbus        ! as above, for US-PP only
  !
  INTEGER, ALLOCATABLE PINMEM ::&
       indv(:,:),        &! index linking  atomic beta's to beta's in the solid
       nhtol(:,:),       &! correspondence n <-> angular momentum l
       nhtolm(:,:),      &! correspondence n <-> combined lm index for (l,m)
       ijtoh(:,:,:),     &! correspondence beta indexes ih,jh -> composite index ijh
       ofsbeta(:)      ! first beta (index in the solid) for each atom 
  !
  LOGICAL :: &
       okvan = .FALSE.,&  ! if .TRUE. at least one pseudo is Vanderbilt
       nlcc_any=.FALSE.   ! if .TRUE. at least one pseudo has core corrections
  ! 
  !!FIXME: vkb should be created and then computed and used on device only
  !!FIXME: this should be done (if it works) with "$acc declare create(vkb)"
  !
  COMPLEX(DP), ALLOCATABLE, TARGET PINMEM :: &
       vkb(:,:)                ! all beta functions in reciprocal space
  REAL(DP), ALLOCATABLE :: &
       becsum(:,:,:)           ! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
  REAL(DP), ALLOCATABLE :: &
       ebecsum(:,:,:)          ! \sum_i f(i) et(i) <psi(i)|beta_l><beta_m|psi(i)>
  REAL(DP), ALLOCATABLE PINMEM :: &
       dvan(:,:,:),           &! the D functions of the solid
       deeq(:,:,:,:),         &! the integral of V_eff and Q_{nm} 
       qq_nt(:,:,:),          &! the integral of q functions in the solid (ONE PER NTYP) used to be the qq array
       qq_at(:,:,:),          &! the integral of q functions in the solid (ONE PER ATOM !!!!) 
       nhtoj(:,:)              ! correspondence n <-> total angular momentum
  !
  COMPLEX(DP), ALLOCATABLE :: & ! variables for spin-orbit/noncolinear case:
       qq_so(:,:,:,:),           &! Q_{nm}
       dvan_so(:,:,:,:),         &! D_{nm}
       deeq_nc(:,:,:,:)           ! \int V_{eff}(r) Q_{nm}(r) dr 
  !
  ! spin-orbit coupling: qq and dvan are complex, qq has additional spin index
  ! noncolinear magnetism: deeq is complex (even in absence of spin-orbit)
  !
  REAL(DP), ALLOCATABLE PINMEM :: &
       beta(:,:,:)           ! beta functions for CP (without struct.factor)
  REAL(DP), ALLOCATABLE PINMEM :: &
       dbeta(:,:,:,:,:)      ! derivative of beta functions w.r.t. cell for CP (without struct.factor)
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  subroutine aainit(lli)
    !-----------------------------------------------------------------------
    ! Product coefficients of real spherical harmonics, in ylmr2 order:
    ! Y_li(r) Y_lj(r) = sum_l ap(l,li,lj) Y_l(r).
    ! Evaluate the Gaunt integral analytically instead of inverting a matrix
    ! of harmonics at random points. See DLMF 34.3.22 and 34.2.4 (3-j sum).
    implicit none
    integer :: lli
    integer :: llx, li, lj, lm, l1, l2, l3, i
    real(dp) :: factorial(0:4*lmaxx+1)
    integer :: angular(lqmax*lqmax), magnetic(2,lqmax*lqmax)
    complex(dp) :: transform(2,lqmax*lqmax)

    if (lli < 0) call upf_error('aainit','lli not allowed',lli)
    if (lli*lli > nlx) call upf_error('aainit','nlx is too small ',lli*lli)
    llx = (2*lli-1)**2
    if (2*lli-1 > lqmax) &
         call upf_error('aainit','ap leading dimension is too small',llx)

    ap = 0.0_dp
    lpx = 0
    lpl = 0
    if (lli == 0) return
    factorial(0) = 1.0_dp
    do i = 1, ubound(factorial,1)
       factorial(i) = real(i,dp)*factorial(i-1)
    end do
    do lm = 1, llx
       call real_harmonic_components(lm,angular(lm),magnetic(:,lm),transform(:,lm))
    end do
    do li = 1, lli*lli
       l2 = angular(li)
       do lj = 1, li
          l3 = angular(lj)
          do lm = 1, llx
             l1 = angular(lm)
             if (l1 < abs(l2-l3) .or. l1 > l2+l3) cycle
             if (mod(l1+l2+l3,2) /= 0) cycle
             ap(lm,li,lj) = real_gaunt(l1,l2,l3,magnetic(:,lm),magnetic(:,li),magnetic(:,lj), &
                                      transform(:,lm),transform(:,li),transform(:,lj),factorial)
             ap(lm,lj,li) = ap(lm,li,lj)
             if (abs(ap(lm,li,lj)) > 1.d-3) then
                lpx(li,lj) = lpx(li,lj) + 1
                if (lpx(li,lj) > mx) &
                     call upf_error('aainit','mx dimension too small',lpx(li,lj))
                lpl(li,lj,lpx(li,lj)) = lm
             end if
          end do
          lpx(lj,li) = lpx(li,lj)
          lpl(lj,li,:) = lpl(li,lj,:)
       end do
    end do
  end subroutine aainit

  subroutine real_harmonic_components(lm,l,m,c)
    ! ylmr2 orders each l as m=0, cos(phi), sin(phi), cos(2phi), ...
    ! With Condon-Shortley complex harmonics:
    ! R_cos = (Y_m + (-1)^m Y_-m)/sqrt(2)
    ! R_sin = (Y_m - (-1)^m Y_-m)/(i sqrt(2)).
    implicit none
    integer, intent(in) :: lm
    integer, intent(out) :: l, m(2)
    complex(dp), intent(out) :: c(2)
    integer :: k, phase
    real(dp) :: q
    l = int(sqrt(real(lm-1,dp)))
    k = lm-l*l-1
    m = 0
    c = (0.0_dp,0.0_dp)
    if (k == 0) then
       c(1) = (1.0_dp,0.0_dp)
    else
       m(1) = (k+1)/2
       m(2) = -m(1)
       phase = 1-2*mod(m(1),2)
       q = 1.0_dp/sqrt(2.0_dp)
       if (mod(k,2) == 1) then
          c(1) = cmplx(q,0.0_dp,dp)
          c(2) = cmplx(phase*q,0.0_dp,dp)
       else
          c(1) = cmplx(0.0_dp,-q,dp)
          c(2) = cmplx(0.0_dp,phase*q,dp)
       end if
    end if
  end subroutine real_harmonic_components

  real(dp) function real_gaunt(l1,l2,l3,m1,m2,m3,c1,c2,c3,f) result(value)
    use upf_const, only: fpi
    implicit none
    integer, intent(in) :: l1,l2,l3,m1(2),m2(2),m3(2)
    complex(dp), intent(in) :: c1(2),c2(2),c3(2)
    real(dp), intent(in) :: f(0:)
    integer :: i,j,k
    real(dp) :: prefactor
    complex(dp) :: total
    prefactor = sqrt(real((2*l1+1)*(2*l2+1)*(2*l3+1),dp)/fpi) * &
                wigner_3j(l1,l2,l3,0,0,0,f)
    total = (0.0_dp,0.0_dp)
    do i=1,2
       do j=1,2
          do k=1,2
             if (m1(i)+m2(j)+m3(k) /= 0) cycle
             total = total+c1(i)*c2(j)*c3(k)*wigner_3j(l1,l2,l3,m1(i),m2(j),m3(k),f)
          end do
       end do
    end do
    value = prefactor*real(total,dp)
  end function real_gaunt

  real(dp) function wigner_3j(j1,j2,j3,m1,m2,m3,f) result(value)
    ! Finite Racah sum for integer angular momenta. Only the small bounded
    ! angular momenta supported by upflib are used (j1+j2+j3 <= 4*lmaxx).
    implicit none
    integer, intent(in) :: j1,j2,j3,m1,m2,m3
    real(dp), intent(in) :: f(0:)
    integer :: z,zmin,zmax,phase
    real(dp) :: prefactor,series
    value = 0.0_dp
    if (m1+m2+m3 /= 0) return
    if (abs(m1)>j1 .or. abs(m2)>j2 .or. abs(m3)>j3) return
    if (j3<abs(j1-j2) .or. j3>j1+j2) return
    phase = 1-2*mod(abs(j1-j2-m3),2)
    prefactor = real(phase,dp)*sqrt( &
         f(j1+j2-j3)*f(j1-j2+j3)*f(-j1+j2+j3)/f(j1+j2+j3+1) * &
         f(j1+m1)*f(j1-m1)*f(j2+m2)*f(j2-m2)*f(j3+m3)*f(j3-m3))
    zmin = max(0,j2-j3-m1,j1-j3+m2)
    zmax = min(j1+j2-j3,j1-m1,j2+m2)
    series = 0.0_dp
    do z=zmin,zmax
       series = series+real(1-2*mod(z,2),dp)/(f(z)*f(j1+j2-j3-z)* &
                f(j1-m1-z)*f(j2+m2-z)*f(j3-j2+m1+z)*f(j3-j1-m2+z))
    end do
    value = prefactor*series
  end function wigner_3j
  !
  !-----------------------------------------------------------------------
  subroutine allocate_uspp(use_gpu,noncolin,lspinorb,tqr,nhm,nsp,nat,nspin)
    !-----------------------------------------------------------------------
    implicit none
    logical, intent(in) :: use_gpu
    logical, intent(in) :: noncolin,lspinorb,tqr
    integer, intent(in) :: nhm,nsp,nat,nspin
    !
    allocate( nhtolm(nhm,nsp) )
    allocate( indv(nhm,nsp)   )
    !$acc enter data create(indv)
    allocate( nhtol(nhm,nsp)  )
    !$acc enter data create(nhtol)
    allocate( nhtoj(nhm,nsp)  )
    !$acc enter data create(nhtoj)
    allocate( ijtoh(nhm,nhm,nsp) )
    !$acc enter data create(ijtoh)
    allocate( deeq(nhm,nhm,nat,nspin) )
    !$acc enter data create(deeq)
    if ( noncolin ) then
       allocate( deeq_nc(nhm,nhm,nat,nspin) )
       !$acc enter data create(deeq_nc)
    endif
    allocate( qq_at(nhm,nhm,nat) )
    !$acc enter data create(qq_at)
    allocate( qq_nt(nhm,nhm,nsp) )
    ! set the internal spin-orbit flag
    is_spinorbit = lspinorb
    if ( lspinorb ) then
       allocate( qq_so(nhm,nhm,4,nsp) )
       !$acc enter data create(qq_so)
       allocate( dvan_so(nhm,nhm,nspin,nsp) )
       !$acc enter data create(dvan_so)
       allocate( fcoef(nhm,nhm,2,2,nsp) )
       !$acc enter data create(fcoef)
    else
       allocate( dvan(nhm,nhm,nsp) )
       !$acc enter data create(dvan)
    endif
    allocate(becsum( nhm*(nhm+1)/2, nat, nspin))
    !$acc enter data create(becsum)
    if (tqr) then
       allocate(ebecsum( nhm*(nhm+1)/2, nat, nspin))
       !$acc enter data create(ebecsum)
    endif
    allocate( ofsbeta(nat) )
    !
  end subroutine allocate_uspp
  !
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_uspp()
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !
    IF( ALLOCATED( ofsbeta ) )    DEALLOCATE( ofsbeta )
    IF( ALLOCATED( nhtolm ) )      DEALLOCATE( nhtolm )
    IF( ALLOCATED( nhtol ) ) THEN
        !$acc exit data delete( nhtol ) 
        DEALLOCATE( nhtol )
    END IF
    IF( ALLOCATED( indv ) ) THEN
        !$acc exit data delete( indv ) 
        DEALLOCATE( indv )
    END IF
    IF( ALLOCATED( ijtoh ) ) THEN
        !$acc exit data delete( ijtoh ) 
        DEALLOCATE( ijtoh )
    END IF
    IF( ALLOCATED( nhtoj ) ) THEN
        !$acc exit data delete( nhtoj ) 
        DEALLOCATE( nhtoj )
    END IF
    IF( ALLOCATED( vkb ) ) THEN
        !$acc exit data delete(vkb ) 
        DEALLOCATE( vkb )
    END IF 
    !$acc exit data delete( becsum )
    IF( ALLOCATED( becsum ) )     DEALLOCATE( becsum )
    IF( ALLOCATED( ebecsum ) ) THEN
       !$acc exit data delete( ebecsum )
       DEALLOCATE( ebecsum )
    END IF
    IF( ALLOCATED( qq_at ) ) THEN
      !$acc exit data delete( qq_at )
      DEALLOCATE( qq_at )
    ENDIF
    IF( ALLOCATED( qq_nt ) )      DEALLOCATE( qq_nt )
    IF( ALLOCATED( dvan ) ) THEN
      !$acc exit data delete( dvan )
      DEALLOCATE( dvan )
    END IF
    IF( ALLOCATED( deeq ) ) THEN
      !$acc exit data delete( deeq )
      DEALLOCATE( deeq )
    ENDIF
    IF( ALLOCATED( qq_so ) ) THEN
      !$acc exit data delete( qq_so )
      DEALLOCATE( qq_so )
    ENDIF
    IF( ALLOCATED( dvan_so ) ) THEN
      !$acc exit data delete( dvan_so )
      DEALLOCATE( dvan_so )
    END IF
    IF( ALLOCATED( deeq_nc ) ) THEN
      !$acc exit data delete( deeq_nc )
      DEALLOCATE( deeq_nc )
    ENDIF
    IF( ALLOCATED( fcoef ) ) THEN
      !$acc exit data delete( fcoef )
      DEALLOCATE( fcoef )
    ENDIF
    IF( ALLOCATED( beta ) )       DEALLOCATE( beta )
    IF( ALLOCATED( dbeta ) )      DEALLOCATE( dbeta )
    !
  END SUBROUTINE deallocate_uspp
  !
END MODULE uspp

