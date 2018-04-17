module mod_spline

logical, parameter :: use_spline = .true.
integer, parameter :: spline_integration_method = -1 ! -1: spline, -2: simpson

contains

Subroutine spline_v2 (n, x, f, cf)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   ld : leading dimension (in,integer)
!   f  : input data array (in,real(ld,n))
!   cf : cubic spline coefficients (out,real(3,n))
! !DESCRIPTION:
!   Calculates the coefficients of a cubic spline fitted to input data. In other
!   words, given a set of data points $f_i$ defined at $x_i$, where
!   $i=1\ldots n$, the coefficients $c_j^i$ are determined such that
!   $$ y_i(x)=f_i+c_1^i(x-x_i)+c_2^i(x-x_i)^2+c_3^i(x-x_i)^3, $$
!   is the interpolating function for $x\in[x_i,x_{i+1})$. This is done by
!   determining the end-point coefficients $c_2^1$ and $c_2^n$ from the first
!   and last three points, and then solving the tridiagonal system
!   $$ d_{i-1}c_2^{i-1}+2(d_{i-1}+d_i)c_2^i+d_ic_2^{i+1}
!    =3\left(\frac{f_{i+1}-f_i}{d_i}-\frac{f_i-f_{i-1}}{d_{i-1}}\right), $$
!   where $d_i=x_{i+1}-x_i$, for the intermediate coefficients.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!   Improved speed and accuracy, April 2006 (JKD)
!   Optimisations and improved end-point coefficients, February 2008 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: n
      Real (8), Intent (In) :: x (n)
      Real (8), Intent (In) :: f (n)
      Real (8), Intent (Out) :: cf (3, n)
! local variables
      Integer :: i
      Real (8) :: t1, t2, t3, t4
! automatic arrays
      Real (8) :: w (n)
      Do i = 1, n - 1
         w (i) = 1.d0 / (x(i+1)-x(i))
         cf (1, i) = w (i) * (f(i+1)-f(i))
      End Do
      cf (2, 1) = 1.d0
! estimate second derivative at the first point
      cf (3, 1) = (cf(1, 2)-cf(1, 1)) / (x(3)-x(1))
! use Gaussian elimination to solve tridiagonal system
      t1 = (x(2)-x(1)) * w (2)
      t2 = t1 * cf (2, 1)
      t3 = 1.d0 / (2.d0*(t1+1.d0))
      cf (2, 2) = t3
      t4 = 3.d0 * (cf(1, 2)-cf(1, 1)) * w (2) - t2 * cf (3, 1)
      cf (3, 2) = t4
      Do i = 3, n - 1
         t1 = (x(i)-x(i-1)) * w (i)
         t2 = t1 * t3
         t3 = 1.d0 / (2.d0*t1+2.d0-t2)
         cf (2, i) = t3
         t4 = 3.d0 * (cf(1, i)-cf(1, i-1)) * w (i) - t2 * t4
         cf (3, i) = t4
      End Do
      cf (2, n) = 1.d0
! estimate second derivative at the last point
      t1 = (cf(1, n-1)-cf(1, n-2)) / (x(n)-x(n-2))
      cf (3, n) = t1
      Do i = n - 1, 2, - 1
         t1 = cf (3, i) - t1 * cf (2, i+1)
         cf (3, i) = t1
      End Do
! compute coefficients
      Do i = 1, n
         cf (2, i) = cf (2, i) * cf (3, i)
      End Do
      Do i = 1, n - 1
         t1 = 0.3333333333333333333d0 * (cf(2, i+1)-cf(2, i)) * w (i)
         cf (3, i) = t1
         t2 = x (i+1) - x (i)
         cf (1, i) = cf (1, i) - (cf(2, i)+t1*t2) * t2
      End Do
! determine end-point coefficients
      t1 = x (n) - x (n-1)
      cf (1, n) = cf (1, n-1) + (2.d0*cf(2, n-1)+3.d0*cf(3, n-1)*t1) * &
     & t1
      cf (3, n) = cf (3, n-1)
      Return
End Subroutine
!EOC

! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine cubic_spline(n,x,f,cf)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   f  : input data array (in,real(n))
!   cf : cubic spline coefficients (out,real(3,n))
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n),f(n)
real(8), intent(out) :: cf(3,n)
! local variables
integer i
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
if (n.le.0) then
  write(*,*)
  write(*,'("Error(spline): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
if (n.eq.1) then
  cf(:,1)=0.d0
  return
end if
if (n.eq.2) then
  cf(1,1)=(f(2)-f(1))/(x(2)-x(1))
  cf(2:3,1)=0.d0
  cf(1,2)=cf(1,1)
  cf(2:3,2)=0.d0
  return
end if
if (n.eq.3) then
  x0=x(1)
  x1=x(2)-x0
  x2=x(3)-x0
  y0=f(1)
  y1=f(2)-y0
  y2=f(3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t1=x1*y2
  t2=x2*y1
  c1=t0*(x2*t2-x1*t1)
  c2=t0*(t1-t2)
  cf(1,1)=c1
  cf(2,1)=c2
  cf(3,1)=0.d0
  t3=2.d0*c2
  cf(1,2)=c1+t3*x1
  cf(2,2)=c2
  cf(3,2)=0.d0
  cf(1,3)=c1+t3*x2
  cf(2,3)=c2
  cf(3,3)=0.d0
  return
end if
y0=f(1)
y1=f(2)-y0
y2=f(3)-y0
y3=f(4)-y0
x0=x(1)
x1=x(2)-x0
x2=x(3)-x0
x3=x(4)-x0
t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
t1=x1*x2*y3
t2=x2*x3*y1
t3=x3*x1*y2
t4=x1**2
t5=x2**2
t6=x3**2
y1=t3*t6-t1*t5
y3=t2*t5-t3*t4
y2=t1*t4-t2*t6
c1=t0*(x1*y1+x2*y2+x3*y3)
c2=-t0*(y1+y2+y3)
c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
cf(1,1)=c1
cf(2,1)=c2
cf(3,1)=c3
cf(1,2)=c1+2.d0*c2*x1+3.d0*c3*t4
cf(2,2)=c2+3.d0*c3*x1
cf(3,2)=c3
if (n.eq.4) then
  cf(1,3)=c1+2.d0*c2*x2+3.d0*c3*t5
  cf(2,3)=c2+3.d0*c3*x2
  cf(3,3)=c3
  cf(1,4)=c1+2.d0*c2*x3+3.d0*c3*t6
  cf(2,4)=c2+3.d0*c3*x3
  cf(3,4)=c3
  return
end if
do i=3,n-2
  y0=f(i)
  y1=f(i-1)-y0
  y2=f(i+1)-y0
  y3=f(i+2)-y0
  x0=x(i)
  x1=x(i-1)-x0
  x2=x(i+1)-x0
  x3=x(i+2)-x0
  t1=x1*x2*y3
  t2=x2*x3*y1
  t3=x3*x1*y2
  t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
  c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
  t4=x1**2
  t5=x2**2
  t6=x3**2
  y1=t3*t6-t1*t5
  y2=t1*t4-t2*t6
  y3=t2*t5-t3*t4
  cf(1,i)=t0*(x1*y1+x2*y2+x3*y3)
  cf(2,i)=-t0*(y1+y2+y3)
  cf(3,i)=c3
end do
c1=cf(1,n-2)
c2=cf(2,n-2)
c3=cf(3,n-2)
cf(1,n-1)=c1+2.d0*c2*x2+3.d0*c3*t5
cf(2,n-1)=c2+3.d0*c3*x2
cf(3,n-1)=c3
cf(1,n)=c1+2.d0*c2*x3+3.d0*c3*t6
cf(2,n)=c2+3.d0*c3*x3
cf(3,n)=c3
return
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
real(8) x0,x1,x2,dx
! automatic arrays
real(8) cf(3,n)
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
  !call cubic_spline(n,x,f,cf)
  call spline_v2(n,x,f,cf)
  g(1)=0.d0
  do i=1,n-1
    dx=x(i+1)-x(i)
    g(i+1)=g(i)+(((0.25d0*cf(3,i)*dx+0.3333333333333333333d0*cf(2,i))*dx &
     +0.5d0*cf(1,i))*dx+f(i))*dx
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
