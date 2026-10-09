! Copyright (C) 2026 Quantum ESPRESSO group
! This file is distributed under the terms of the GNU General Public License.
! See the file License in the root directory of the present distribution.
program test_aainit
  use upf_kinds, only: dp
  use upf_const, only: pi, fpi
  use upf_params, only: lmaxx, lqmax
  use upf_parallel_include
  use uspp, only: aainit, ap, lpx, lpl, nlx
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: nq=lqmax*lqmax, nd=37
  real(dp) :: ref(nq,nlx,nlx), refined(nq,nlx,nlx), saved(nq,nlx,nlx)
  real(dp) :: xyz(3,nd), rr(nd), y(nd,nq), residual, sparse(nd), before, after
  integer :: i,j,k,l,n,lli,rep,failures,rank,info,all_failures
  character(32) :: argument
  failures=0
  rank=0
#if defined(__MPI)
  call MPI_INIT(info)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,info)
#endif
  call get_command_argument(1,argument)
  if (trim(argument)=='negative' .or. trim(argument)=='oversized') then
     if (trim(argument)=='negative') call aainit(-1)
     if (trim(argument)=='oversized') call aainit(lmaxx+2)
     ! The error harness must reject incorrectly accepted input.
#if defined(__MPI)
     call MPI_FINALIZE(info)
#endif
     stop 0
  end if
  if (len_trim(argument)>0) stop 3

  call reference(24,64,ref)
  call reference(48,128,refined)
  call check('reference refinement',maxval(abs(ref-refined)),1.e-13_dp)
  do i=1,nd
     xyz(3,i)=1.0_dp-2.0_dp*(real(i,dp)-0.5_dp)/nd
     xyz(1,i)=sqrt(1.0_dp-xyz(3,i)**2)*cos(real(i,dp)*2.399963229728653_dp)
     xyz(2,i)=sqrt(1.0_dp-xyz(3,i)**2)*sin(real(i,dp)*2.399963229728653_dp)
  end do
  rr=sum(xyz*xyz,dim=1)
  call ylmr2(nq,nd,xyz,rr,y)
  call cpu_time(before)
  do rep=1,16
     do lli=1,lmaxx+1
        n=(2*lli-1)**2
        call aainit(lli)
        if (.not.all(ieee_is_finite(ap))) failures=failures+1
        call check('quadrature coefficients',maxval(abs(ap(:n,:lli**2,:lli**2)-refined(:n,:lli**2,:lli**2))),1.e-12_dp)
        residual=0.0_dp
        do j=1,lli**2
           do i=1,lli**2
              call check('exchange symmetry',maxval(abs(ap(:n,i,j)-ap(:n,j,i))),1.e-14_dp)
              residual=max(residual,maxval(abs(matmul(y(:,:n),ap(:n,i,j))-y(:,i)*y(:,j))))
              if (lpx(i,j)<0 .or. lpx(i,j)>size(lpl,3)) then
                 failures=failures+1
                 cycle
              end if
              if (lpx(i,j)/=count(abs(ap(:n,i,j))>1.e-3_dp)) failures=failures+1
              sparse=0.0_dp
              do k=1,lpx(i,j)
                 l=lpl(i,j,k)
                 if (l<1 .or. l>n) then
                    failures=failures+1
                    cycle
                 end if
                 if (k>1) then
                    if (l<=lpl(i,j,k-1)) failures=failures+1
                 end if
                 sparse=sparse+ap(l,i,j)*y(:,l)
              end do
              call check('sparse reconstruction',maxval(abs(sparse-y(:,i)*y(:,j))),1.e-11_dp)
           end do
        end do
        call check('dense reconstruction',residual,1.e-11_dp)
        if (lli==lmaxx+1) then
           if (rep==1) saved=ap
           call check('repeat initialization',maxval(abs(ap-saved)),1.e-14_dp)
        end if
     end do
  end do
  call cpu_time(after)
  if(rank==0) write(*,'(a,es14.6)') '80 initialization/check cycles CPU seconds: ',after-before
  do lli=lmaxx+1,1,-1
     call aainit(lli)
     n=(2*lli-1)**2
     call check('reverse initialization order', &
          maxval(abs(ap(:n,:lli**2,:lli**2)-refined(:n,:lli**2,:lli**2))),1.e-12_dp)
  end do
  call aainit(lmaxx+1)
  call check('Y00 product',abs(ap(1,1,1)-1.0_dp/sqrt(fpi)),1.e-14_dp)
  call check('pz squared',abs(ap(5,2,2)-1.0_dp/sqrt(5.0_dp*pi)),1.e-14_dp)
  call check('px py product',abs(ap(9,3,4)-sqrt(3.0_dp/(20.0_dp*pi))),1.e-14_dp)
  do j=1,nlx
     do i=1,nlx
        do l=1,nq
           block
             integer :: ll,li,lj
             ll=int(sqrt(real(l-1,dp)))
             li=int(sqrt(real(i-1,dp)))
             lj=int(sqrt(real(j-1,dp)))
             if (ll<abs(li-lj) .or. ll>li+lj .or. mod(ll+li+lj,2)/=0) &
                 call check('selection rule',abs(ap(l,i,j)),1.e-14_dp)
           end block
        end do
     end do
  end do
#if defined(__MPI)
  saved=ap
  call MPI_BCAST(saved,size(saved),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
  call check('MPI coefficient agreement',maxval(abs(ap-saved)),0.0_dp)
#endif
  call aainit(0)
  call check('zero angular range',maxval(abs(ap)),0.0_dp)
  all_failures=failures
#if defined(__MPI)
  call MPI_ALLREDUCE(failures,all_failures,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_FINALIZE(info)
#endif
  if (all_failures/=0) then
     if(rank==0) write(*,'(a,i0)') 'FAILURES: ',all_failures
     stop 1
  end if
  if(rank==0) write(*,'(a)') 'PASS: aainit coefficients, reconstruction, sparsity, state and symmetry'
contains
  subroutine check(label,error,tolerance)
    character(*),intent(in)::label
    real(dp),intent(in)::error,tolerance
    if (.not.ieee_is_finite(error) .or. error>tolerance) then
       failures=failures+1
       if(rank==0 .and. failures<=20) write(*,'(a,1x,es14.6,1x,es14.6)') label,error,tolerance
    end if
  end subroutine

  subroutine reference(nt,np,coefficients)
    ! Independent integration of three real harmonics. Gauss-Legendre in
    ! cos(theta), trapezoidal in phi; no 3-j identities are used here.
    integer,intent(in)::nt,np
    real(dp),intent(out)::coefficients(nq,nlx,nlx)
    real(dp)::r(3,nt*np),norms(nt*np),weights(nt*np),values(nt*np,nq)
    real(dp)::x(nt),w(nt),phi,s
    integer::it,ip,g,a,b,c
    call legendre_rule(nt,x,w)
    g=0
    do it=1,nt
       s=sqrt(1.0_dp-x(it)**2)
       do ip=1,np
          g=g+1
          phi=(real(ip-1,dp)+0.37_dp)*2.0_dp*pi/np
          r(:,g)=[s*cos(phi),s*sin(phi),x(it)]
          weights(g)=w(it)*2.0_dp*pi/np
       end do
    end do
    norms=sum(r*r,dim=1)
    call ylmr2(nq,nt*np,r,norms,values)
    do b=1,nlx
       do a=1,nlx
          do c=1,nq
             coefficients(c,a,b)=sum(weights*values(:,c)*values(:,a)*values(:,b))
          end do
       end do
    end do
  end subroutine

  subroutine legendre_rule(n,x,w)
    integer,intent(in)::n
    real(dp),intent(out)::x(n),w(n)
    real(dp)::z,previous,p0,p1,p2,derivative
    integer::i,j,iteration
    do i=1,(n+1)/2
       z=cos(pi*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
       do iteration=1,50
          p0=1.0_dp
          p1=z
          do j=2,n
             p2=((2*j-1)*z*p1-(j-1)*p0)/j
             p0=p1
             p1=p2
          end do
          derivative=n*(z*p1-p0)/(z*z-1.0_dp)
          previous=z
          z=previous-p1/derivative
          if(abs(z-previous)<2.0_dp*epsilon(z)) exit
       end do
       if(iteration>50) stop 4
       x(i)=-z
       x(n+1-i)=z
       w(i)=2.0_dp/((1.0_dp-z*z)*derivative*derivative)
       w(n+1-i)=w(i)
    end do
  end subroutine
end program test_aainit
