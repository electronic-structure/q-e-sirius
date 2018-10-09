module mod_spline

logical, parameter :: use_spline = .true.
integer, parameter :: spline_integration_method = -1 ! -1: spline, -2: simpson

contains

subroutine spline_v3(n, r, f, cf)
implicit none
integer, intent(in)  :: n
real(8), intent(in)  :: r(n)
real(8), intent(in)  :: f(n)
real(8), intent(out) :: cf(3, n)
!
real(8) dl(n - 1)
real(8) du(n - 1)
real(8) d(n)
real(8) x(n)
real(8) dy(n - 1)
real(8) dr(n - 1)
integer i
real(8) h0, h1

do i = 1, n - 1
  dr(i) = r(i + 1) - r(i)
  dy(i) = (f(i + 1) - f(i)) / dr(i)
enddo
do i = 1, n - 2
  x(i + 1) = 6.d0 * (dy(i + 1) - dy(i))
enddo
x(1) = x(2)
x(n) = x(n-1)
do i = 1, n - 2
  d(i + 1) = 2.d0 * (r(i+2) - r(i))
enddo
do i = 1, n - 1
  du(i) = dr(i)
  dl(i) = dr(i)
enddo
h0 = dr(1)
h1 = dr(2)
d(1) = h0 - (h1 / h0) * h1
du(1) = h1 * ((h1 / h0) + 1) + 2 * (h0 + h1)

h0 = dr(n-1)
h1 = dr(n-2)
d(n) = h0 - (h1 / h0) * h1;
dl(n-1) = h1 * ((h1 / h0) + 1) + 2 * (h0 + h1);

call dgtsv(n, 1, dl, d, du, x, n, i)
do i = 1, n - 1
  cf(2, i) = x(i) / 2.d0
  h0 = (x(i + 1) - x(i)) / 6.d0
  cf(1, i) = dy(i) - (cf(2, i) + h0) * dr(i)
  cf(3, i) = h0 / dr(i)
enddo
cf(:,n) = 0.d0
end subroutine


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine integrate(n,f,x,res)
! !INPUT/OUTPUT PARAMETERS:
!   n : number of points (in,integer)
!   x : abscissa array (in,real(n))
!   f : function array (in,real(n))
!   res : resulting integral
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n),f(n)
real(8), intent(out) :: res
! local variables
integer i
real(8) x0,x1,x2,dx,third
! automatic arrays
real(8) cf(3,n)
real(8) cf1(3,n)
real(8) g(n)
if (n.le.0) then
  write(*,*)
  write(*,'("Error(integrate): invalid number of points : ",I8)') n
  write(*,*)
  stop
end if
select case(spline_integration_method)
case(-3)
! low accuracy trapezoidal integration
  g(1)=0.d0
  do i=1,n-1
    g(i+1)=g(i)+0.5d0*(x(i+1)-x(i))*(f(i+1)+f(i))
  end do
case(-2)
! medium accuracy Simpson integration
  g(1)=0.d0
  do i=1,n-2
    x0=x(i)
    x1=x(i+1)
    x2=x(i+2)
    g(i+1)=g(i)+(x0-x1)*(f(i+2)*(x0-x1)**2+f(i+1)*(x2-x0)*(x0+2.d0*x1-3.d0*x2) &
     +f(i)*(x2-x1)*(2.d0*x0+x1-3.d0*x2))/(6.d0*(x0-x2)*(x1-x2))
  end do
  x0=x(n)
  x1=x(n-1)
  x2=x(n-2)
  g(n)=g(n-1)+(x1-x0)*(f(n-2)*(x1-x0)**2+f(n)*(x1-x2)*(3.d0*x2-x1-2.d0*x0) &
   +f(n-1)*(x0-x2)*(3.d0*x2-2.d0*x1-x0))/(6.d0*(x2-x1)*(x2-x0))
case(-1)
  third = 1.d0 / 3
  call spline_v3(n,x,f,cf)
  !call sirius_spline(n, x, f, cf1)
  g(1)=0.d0
  do i=1,n-1
    dx=x(i+1)-x(i)
    g(i+1)=g(i)+(((0.25d0*cf(3,i)*dx + third*cf(2,i))*dx + 0.5d0*cf(1,i))*dx+f(i))*dx
  end do
case default
  write(*,*)
  write(*,'("Error(integrate): wrong integration type: ",I8)') spline_integration_method
  write(*,*)
  stop
end select
res = g(n)
return
end subroutine

end module
