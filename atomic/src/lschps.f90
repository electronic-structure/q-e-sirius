!
! Copyright (C) 2004-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------
SUBROUTINE lschps_meta (mode, z, eps, grid, nin, n, l, e, v, vtau, &
                        u, nstop)
  !----------------------------------------------------------------
  !
  ! Meta-GGA version of lschps
  ! vtau is the meta-GGA potential
  !
  USE kinds, ONLY : DP
  USE radial_grids, ONLY: radial_grid_type
  USE ld1inc, ONLY : cau_fact
  USE atomic_radial, ONLY : aei, aeo, aii, aio, derV
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER, INTENT (in) :: mode, n, l
  real(DP), INTENT(in) :: z, eps
  TYPE (radial_grid_type), INTENT(in) :: grid
  real(DP), INTENT(in) :: v(grid%mesh), vtau(grid%mesh)
  INTEGER, INTENT(inout) :: nin
  real(DP), INTENT(inout) :: e
  INTEGER, INTENT(out) :: nstop
  real (DP), INTENT(out) :: u(grid%mesh)
  !
  ! local variables
  !
  INTEGER, PARAMETER :: maxter=60
  real(DP), EXTERNAL:: d2u, int_0_inf_dr
  ! arrays used as work space
  real(DP), ALLOCATABLE :: up(:),upp(:),cf(:),dv(:),fr(:),frp(:),dvtau(:)
  real(DP):: al, als, cn
  real(DP):: de, emax, emin
  real(DP):: fss, gamma, ro, sc
  real(DP):: sls, sn, uld, uout,  upin, upout
  real(DP):: xkap, mm, mvt, fac1, fact2, alpha
  INTEGER :: i, it, mmax, n_it, node, mch, ierr
  LOGICAL ::  rel, find_e, fixed_e
  !
  !
  IF (mode<1.or.mode>5) STOP 'lschps: wrong mode'
  !
  nstop=0
  al   = grid%dx
  mmax = grid%mesh
  !
  ALLOCATE(up(mmax), stat=ierr)
  ALLOCATE(upp(mmax), stat=ierr)
  ALLOCATE(cf(mmax), stat=ierr)
  ALLOCATE(dv(mmax), stat=ierr)
  ALLOCATE(fr(mmax), stat=ierr)
  ALLOCATE(frp(mmax), stat=ierr)
  ALLOCATE(dvtau(mmax), stat=ierr)
  !
  uld=0.0_dp
  !
  rel    = ( mode == 1 .or. mode == 3 )
  find_e = ( mode == 1 .or. mode == 2 )
  fixed_e= ( mode == 3 .or. mode == 5 )
  !
  ! fss  : square of the fine structure constant
  ! gamma: u(r)=r^gamma is the asymptotic behavior of u(r) for r->0
  !
  IF (rel) THEN
     ! relativistic calculation
     !    fss=(1.0_dp/137.036_dp)**2
     fss=(1.0_dp/cau_fact)**2
     IF(l == 0) THEN
        gamma=sqrt(1.0_dp-fss*z**2)
     ELSE
        gamma=(l*sqrt(l**2-fss*z**2) + &
              (l+1)*sqrt((l+1)**2-fss*z**2))/(2*l+1)
     ENDIF
  ELSE
     ! nonrelativistic calculation
     fss=0.0_dp
     gamma=l+1
  ENDIF
  sls=l*(l+1)
  !
  ! emin, emax = estimated bounds for e
  !
  IF(find_e) THEN
     emax=v(mmax)+sls/grid%r(mmax)**2
     emin=0.0_dp
     DO i=1,mmax
        emin=min(emin,v(i)+sls/grid%r(i)**2)
     ENDDO
     IF(e > emax) e=1.25_dp*emax
     IF(e < emin) e=0.75_dp*emin
     IF(e > emax) e=0.5_dp*(emax+emin)
  ELSEIF(mode == 4) THEN
     emax=e + 10.0_dp
     emin=e - 10.0_dp
  ENDIF
  !
  DO i=1,mmax
     u(i)=0.0_dp
     up(i)=0.0_dp
     upp(i)=0.0_dp
  ENDDO
  als=al**2
  !
  ! calculate dv = dv/dr for darwin correction
  ! calculate dvtau = d tau / dr
  !
  CALL derV(mmax,al,grid%r,v,dv)
  CALL derV(mmax,al,grid%r,vtau,dvtau)
  !
  !     starting of loop on energy for bound state
  !
  DO n_it = 1, maxter
     !
     !     find classical turning point for matching
     !
     DO i=1,mmax
        cf(i)=(sls + (v(i)-e)*grid%r(i)**2)
     ENDDO
     IF(find_e) THEN
        DO i=mmax,2,-1
           IF(cf(i-1) <= 0.0_dp .and. cf(i) > 0.0_dp) THEN
              mch=i
              GOTO 40
           ENDIF
        ENDDO
        !PRINT '('' warning: wfc '',2i2,'' no turning point'')', n, l
        e=0.0
        DO i=1,mmax
           u (i)=0.0
           ! up(i)=0.0
        ENDDO
        nstop=1
        GOTO 999
     ELSE
        mch=nin
     ENDIF
40   CONTINUE
     !
     ! coefficient array for u in differential eq.
     !
     DO i=1,mmax
        mm=1.0_dp+0.25_dp*fss*(e-v(i))
        mvt=1.0_dp/(1.0_dp+0.5_dp*vtau(i))
        fac1=grid%r(i)*fss*0.25_dp*dv(i)*mvt/mm
        !
        cf(i)=als*(sls+mvt*(mm*(v(i)-e)*grid%r(i)**2                     &
             &           +grid%r(i)*0.5_dp*dvtau(i)))
        !
        !            cf(i)=als*(sls+mvt*(mm*(v(i)-e)*grid%r(i)**2))
        !
        !     relativistic coefficient arrays for u (fr) and up (frp).
        !
        fr(i)=als*fac1
        frp(i)=-al*(0.5*mvt*grid%r(i)*dvtau(i)+fac1)
        !            frp(i)=-al*fac1
        !
     ENDDO
     !
     !     start wavefunction with series
     !
     fac1=1.0_dp+0.5_dp*vtau(1)
     alpha=-z/(l+1.0_dp)/fac1*grid%r(1)
     fact2=-l/(2.0_dp*(l+1.0_dp))
     DO i=1,4
        u(i) = grid%r(i)**gamma*exp(alpha)*(1+0.5_dp*vtau(i))**fact2
        up(i)=al*(gamma+(fact2-z/(l+1.0_dp))*grid%r(i)/fac1)*u(i)
        upp(i)=d2u(al,cf(i),fr(i),frp(i),u(i),up(i))
        fac1=(1.0_dp+0.5*vtau(i+1))
        alpha = alpha - z/(l+1.0_dp)/fac1*(grid%r(i+1)-grid%r(i))
     ENDDO
     !
     !     outward integration using predictor once, corrector twice
     !
     node=0
     DO i=4,mch-1
        u (i+1)=u (i)+aeo(up,i)
        up(i+1)=up(i)+aeo(upp,i)
        DO it=1,2
           upp(i+1)=d2u(al,cf(i+1),fr(i+1),frp(i+1),u(i+1),up(i+1))
           up (i+1)=up(i)+aio(upp,i)
           u  (i+1)=u (i)+aio(up,i)
        ENDDO
        IF(u(i+1)*u(i) <= 0.0_dp) node=node+1
     ENDDO
     !
     uout=u(mch)
     upout=up(mch)
     !
     IF(node-n+l+1 == 0 .or. fixed_e) THEN
        !
        IF (find_e) THEN
           !
           ! good number of nodes: start inward integration
           ! at 10*classical turning point with simple exponential
           !
           nin=mch+2.3_dp/al
           IF(nin+4 > mmax ) nin=mmax-4
           !
           mm=1.0_dp+0.25_dp*fss*(e-v(nin))
           fac1=0.5_dp*dvtau(nin)/(1.0_dp+mm*0.5_dp*vtau(nin))
           fact2=0.5_dp*dvtau(nin)/grid%r(nin)
           alpha=1.0_dp+mm*0.5_dp*vtau(nin)
           xkap=sqrt(fac1**2.0_dp + 4.0_dp*(sls/grid%r(nin)**2 + &
                   (v(nin)-e+fact2)/alpha))
           !
           xkap=(-fac1+xkap)/2.0_dp
           DO i=nin,nin+4
              u (i) =exp(-xkap*(grid%r(i)-grid%r(nin)))
              up(i) =-grid%r(i)*al*xkap*u(i)
              upp(i)=d2u(al,cf(i),fr(i),frp(i),u(i),up(i))
           ENDDO
           !
           ! integrate inward
           !
           DO i=nin,mch+1,-1
              u (i-1)=u (i)+aei(up,i)
              up(i-1)=up(i)+aei(upp,i)
              DO it=1,2
                 upp(i-1)=d2u(al,cf(i-1),fr(i-1),frp(i-1),u(i-1),up(i-1))
                 up(i-1)=up(i)+aii(upp,i)
                 u (i-1)=u (i)+aii(up,i)
              ENDDO
           ENDDO
           !
           ! scale outside wf for continuity
           !
           sc=uout/u(mch)
           DO i=mch,nin
              up(i)=sc*up(i)
              u (i)=sc*u (i)
           ENDDO
           !
           upin=up(mch)
           !
        ELSE
           !
           ! this is used only in fixed-logarithmic-derivative calculation
           !
           upin=uld*uout
           !
        ENDIF
        !
        ! normalization: upp is used to store u^2
        !
        DO i=1,nin
           upp(i)=u(i)**2
        ENDDO
        !
        ! integral over dr (includes series expansion for r->0) ...
        ! assumes u(r)=r^(l+1) behaviour at the origin
        !
        sn=int_0_inf_dr ( upp, grid, nin-3, 2*l+2 )
        !
        ! ...extrapolation for the asymptotic behaviour (maybe useless...)
        !
        sn=sn + al*(23.0_dp*grid%r(nin-2)*upp(nin-2)     &
                  + 28.0_dp*grid%r(nin-1)*upp(nin-1)     &
                  +  9.0_dp*grid%r(nin  )*upp(nin  ))/24.0_dp
        !
        ! ...normalize u
        !
        cn=1.0_dp/sqrt(sn)
        uout=cn*uout
        upout=cn*upout
        upin=cn*upin
        !
        DO i=1,nin
           up(i)=cn*up(i)
           u(i)=cn*u(i)
        ENDDO
        DO i=nin+1,mmax
           u(i) =0.0_dp
           up(i)=0.0_dp
        ENDDO
        !
        !
        ! exit for fixed-energy calculation
        !
        IF (fixed_e)  GOTO 999
        !
        ! perturbation theory for energy shift
        de=uout*(upout-upin)/(al*grid%r(mch))
        !
        ! convergence test and possible exit
        !
        IF ( abs(de) < max(abs(e),0.2_dp)*eps) GOTO 999
        !
        IF(de > 0.0_dp) THEN
           emin=e
        ELSE
           emax=e
        ENDIF
        e=e+de
        IF(e > emax .or. e < emin) e=0.5_dp*(emax+emin)
        !
     ELSEIF(node-n+l+1 < 0) THEN
        !
        ! too few nodes
        !
        emin=e
        e=0.5_dp*(emin+emax)
        !
     ELSE
        !
        ! too many nodes
        !
        emax=e
        e=0.5_dp*(emin+emax)
        !
     ENDIF
     !
     ! loop back to converge e
     !
  ENDDO
  !
  !PRINT '('' warning: wfc '',2i2,'' not converged'')', n, l
  nstop=1
  u=0.0_dp
  nstop=1
  !
  ! deallocate arrays and exit
  !
999 CONTINUE

  DEALLOCATE(frp)
  DEALLOCATE(fr)
  DEALLOCATE(dv)
  DEALLOCATE(cf)
  DEALLOCATE(upp)
  DEALLOCATE(up)
  !do i=1,mmax
  !   up(i)=up(i)/grid%r(i)/al
  !enddo
  RETURN
  !
END SUBROUTINE lschps_meta
!
!---------------------------------------------------------------
FUNCTION d2u(al,cf,fr,frp,u,up)
  !---------------------------------------------------------------
  ! second derivative from radial KS equation
  USE kinds, ONLY : DP
  IMPLICIT NONE
  real(dp):: d2u
  real(dp), INTENT(in):: al,cf,fr,frp,u,up
  !
  d2u = (al+frp)*up + (cf+fr)*u
  RETURN
END FUNCTION d2u

!-----------------------------------------------------------
FUNCTION estimatealpha(mmax,u,up,al,r)
  !-----------------------------------------------------------
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: mmax
  real(dp) :: estimatealpha
  real(dp), INTENT(in) :: u(mmax),up(mmax),al,r(mmax)
  INTEGER i,istart,iend
  estimatealpha = 0.0_dp
  istart=5
  iend=100
  DO i=istart,iend
     IF(u(i) > 1.0d-8) THEN
        estimatealpha=estimatealpha+(1.0_dp-up(i)/u(i)/al)/r(i)
     ENDIF
  ENDDO
  estimatealpha=estimatealpha/(iend-istart+1)
  RETURN
END FUNCTION estimatealpha
