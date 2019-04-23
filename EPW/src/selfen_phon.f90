  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_phon_q(iqq, iq, totq)
  !-----------------------------------------------------------------------
  !!
  !!  compute the imaginary part of the phonon self energy due to electron-
  !!  phonon interaction in the Migdal approximation. This corresponds to 
  !!  the phonon linewidth (half width). The phonon frequency is taken into
  !!  account in the energy selection rule.
  !!
  !!  Use matrix elements, electronic eigenvalues and phonon frequencies
  !!  from ep-wannier interpolation.  This routine is similar to the one above
  !!  but it is ONLY called from within ephwann_shuffle and calculates 
  !!  the selfenergy for one phonon at a time.  Much smaller footprint on
  !!  the disk
  !!
  !!  RM 24/02/2014
  !!  redefined the size of coskkq, vkk, vkq within the fermi windwow
  !!  cleaned up the subroutine
  !!
  !-----------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  use phcom,      ONLY : nmodes
  use epwcom,     ONLY : nbndsub, fsthick, efermi_read, fermi_energy,  &
                         eptemp, ngaussw, degaussw, shortrange,        &
                         nsmear, delta_smear, eps_acustic, specfun_ph, &
                         delta_approx, vme
  use pwcom,      ONLY : nelec, ef
  USE klist_epw,  ONLY : isk_dummy
  use elph2,      ONLY : epf17, ibndmax, ibndmin, etf, wkf, xqf, wqf, nkqf, &
                         nkf, wf, nkqtotf, xqf, lambda_all, lambda_v_all,   &
                         dmef, vmef, gamma_all,gamma_v_all, efnew
  USE constants_epw, ONLY : ryd2mev, ryd2ev, two, zero, pi, eps4, eps6, eps8
  use mp,         ONLY : mp_barrier, mp_sum
  use mp_global,  ONLY : inter_pool_comm
  !
  implicit none
  !
  INTEGER, INTENT (in) :: iqq
  !! Current q-point index from the selecq
  INTEGER, INTENT (in) :: iq
  !! Current q-point index
  INTEGER, INTENT (in) :: totq
  !! Total number of q-points in selecq.fmt
  ! 
  ! Local variables 
  !
  INTEGER :: ik
  !! Counter on the k-point index 
  INTEGER :: ikk
  !! k-point index
  INTEGER :: ikq
  !! q-point index 
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: jbnd
  !! Counter on bands
  INTEGER :: imode
  !! Counter on mode
  INTEGER :: jmode
  !! Counter on mode
  INTEGER :: fermicount
  !! Number of states on the Fermi surface
  INTEGER :: ismear
  !! Upper bounds index after k or q paral
  !! Smearing for the Gaussian function 
  INTEGER :: n
  !! Counter on number of mode degeneracies
  ! 
  REAL(kind=DP) :: g2
  !! Electron-phonon matrix elements squared in Ry^2
  REAL(kind=DP) :: ekk
  !! Eigen energy on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level
  REAL(kind=DP) :: wq
  !! Phonon frequency on the fine grid
  REAL(kind=DP) :: wq_tmp
  !! Temporary Phonon frequency on the fine grid
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: wgkq
  !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
  REAL(kind=DP) :: weight
  !! Imaginary part of the phonhon self-energy factor 
  !!$$ \pi N_q \Im(\frac{f_{nk}(T) - f_{mk+q(T)}}{\varepsilon_{nk}-\varepsilon_{mk+q}-\omega_{q\nu}+i\delta}) $$
  !! In practice the imaginary is performed with a delta Dirac
  REAL(kind=DP) :: dosef
  !! Density of state N(Ef)
  REAL(kind=DP) :: w0g1
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: w0g2
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: inv_wq
  !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
  REAL(kind=DP) :: inv_eptemp0
  !! Inverse of temperature define for efficiency reasons
  REAL(kind=DP) :: g2_tmp
  !! If the phonon frequency is too small discart g
  REAL(kind=DP) :: gamma(nmodes)
  !! Gamma is the imaginary part of the phonon self-energy 
  REAL(kind=DP) :: gamma_v(nmodes)
  !! Gamma is the imaginary part of the phonon self-energy multiplied by (1-coskkq)
  REAL(kind=DP) :: coskkq(ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  !! $$(v_k \cdot v_{k+q}) / |v_k|^2$$
  REAL(kind=DP) :: DDOT
  !! Dot product function
  REAL(kind=DP) :: degaussw0
  !! degaussw0 = (ismear-1) * delta_smear + degaussw
  REAL(kind=DP) :: inv_degaussw0
  !! Inverse degaussw0 for efficiency reasons
  REAL(kind=DP) :: lambda_tot
  !! Integrated lambda function
  REAL(kind=DP) :: lambda_tr_tot
  !! Integrated transport lambda function
  REAL(kind=DP) :: wgkk
  !! Fermi-Dirac occupation factor $f_{nk}(T)$
  REAL(kind=DP) :: eptemp0
  !!eptemp0   = (ismear-1) * delta_smear + eptem
  REAL(kind=DP) :: vkk(3,ibndmax-ibndmin+1)
  !! Electronic velocity $v_{nk}$
  REAL(kind=DP) :: vkq(3,ibndmax-ibndmin+1)
  !! Electronic velocity $v_{nk+q}$
  REAL(kind=DP) :: tmp
  !! Temporary value of lambda for av.
  REAL(kind=DP) :: tmp2
  !! Temporary value of lambda_v for av.
  REAL(kind=DP) :: tmp3
  !! Temporary value of lambda_v for av.
  REAL(kind=DP) :: tmp4
  !! Temporary value of lambda_v for av.
  REAL(kind=DP) :: lambda_tmp(nmodes)
  !! Temporary value of lambda for av.  
  REAL(kind=DP) :: lambda_v_tmp(nmodes)
  !! Temporary value of lambda v for av.  
  REAL(kind=DP) :: gamma_tmp(nmodes)
  !! Temporary value of gamma for av.  
  REAL(kind=DP) :: gamma_v_tmp(nmodes)
  !! Temporary value of gamma v for av.  
  REAL(kind=DP), external :: dos_ef
  !! Function to compute the Density of States at the Fermi level
  REAL(kind=DP), external :: wgauss
  !! Fermi-Dirac distribution function (when -99)
  REAL(kind=DP), external :: w0gauss
  !! This function computes the derivative of the Fermi-Dirac function
  !! It is therefore an approximation for a delta function
  REAL(kind=DP), external :: efermig
  !! Return the fermi energy
  !  
  IF ( iq == 1 ) THEN 
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Phonon (Imaginary) Self-Energy in the Migdal Approximation")') 
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !
    IF ( fsthick.lt.1.d3 ) &
         WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
         'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
    !
    IF ( .not. ALLOCATED (lambda_all) )   ALLOCATE( lambda_all  (nmodes, totq, nsmear) )
    IF ( .not. ALLOCATED (lambda_v_all) ) ALLOCATE( lambda_v_all(nmodes, totq, nsmear) )
    lambda_all(:,:,:)   = zero
    lambda_v_all(:,:,:) = zero
    IF ( .not. ALLOCATED (gamma_all) )    ALLOCATE( gamma_all  (nmodes, totq, nsmear) )
    IF ( .not. ALLOCATED (gamma_v_all) )  ALLOCATE( gamma_v_all(nmodes, totq, nsmear) )
    gamma_all(:,:,:)   = zero
    gamma_v_all(:,:,:) = zero
    !
  ENDIF
  !
  DO ismear = 1, nsmear
    !
    degaussw0 = (ismear-1) * delta_smear + degaussw
    eptemp0   = (ismear-1) * delta_smear + eptemp
    ! 
    ! SP: Multiplication is faster than division ==> Important if called a lot
    !     in inner loops
    inv_degaussw0 = 1.0/degaussw0
    inv_eptemp0   = 1.0/eptemp0
    !
    ! Fermi level and corresponding DOS
    !
    IF ( efermi_read ) THEN
      !
      ef0 = fermi_energy
      !
    ELSE IF (nsmear > 1) THEN
      !
      ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw0, ngaussw, 0, isk_dummy)
      ! if some bands are skipped (nbndskip.neq.0), nelec has already been
      ! recalculated 
      ! in ephwann_shuffle
      !
    ELSE !SP: This is added for efficiency reason because the efermig routine is slow
      ef0 = efnew
    ENDIF
    !
    dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nkqf, nbndsub)
    !  N(Ef) in the equation for lambda is the DOS per spin
    dosef = dosef / two
    !
    IF ( iq == 1 ) THEN 
      WRITE (stdout, 100) degaussw0 * ryd2ev, ngaussw
      WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
      !WRITE (stdout, 101) dosef / ryd2ev, ef  * ryd2ev
    ENDIF
    !
    CALL start_clock('PH SELF-ENERGY')
    !
    fermicount = 0
    wgkk = 0.0_DP
    w0g1 = 0.0_DP
    gamma(:)   = zero
    gamma_v(:) = zero
    !
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      ! 
      coskkq = zero
      ! coskkq = (vk dot vkq) / |vk|^2  appears in Grimvall 8.20
      ! this is different from :   coskkq = (vk dot vkq) / |vk||vkq|
      ! In principle the only coskkq contributing to lambda_tr are both near the
      ! Fermi surface and the magnitudes will not differ greatly between vk and vkq
      ! we may implement the approximation to the angle between k and k+q
      ! vectors also listed in Grimvall
      !
      IF (vme ) THEN 
        DO ibnd = 1, ibndmax-ibndmin+1
          DO jbnd = 1, ibndmax-ibndmin+1
            !
            ! vmef is in units of Ryd * bohr
            !
            vkk(:, ibnd ) = REAL (vmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk ) )
            vkq(:, jbnd ) = REAL (vmef (:, ibndmin-1+jbnd, ibndmin-1+jbnd, ikq ) )
            IF ( abs ( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2) > eps4) &
                coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1)  / &
                DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
          ENDDO
        ENDDO
      ELSE
        DO ibnd = 1, ibndmax-ibndmin+1
          DO jbnd = 1, ibndmax-ibndmin+1
            !
            ! v_(k,i) = 1/m <ki|p|ki> = 2 * dmef (:, i,i,k)
            ! 1/m  = 2 in Rydberg atomic units
            !
            vkk(:, ibnd ) = 2.0 * REAL (dmef (:, ibndmin-1+ibnd, ibndmin-1+ibnd, ikk ) )
            vkq(:, jbnd ) = 2.0 * REAL (dmef (:, ibndmin-1+jbnd, ibndmin-1+jbnd, ikq ) )
            IF ( abs ( vkk(1,ibnd)**2 + vkk(2,ibnd)**2 + vkk(3,ibnd)**2) > eps4) &
                coskkq(ibnd, jbnd ) = DDOT(3, vkk(:,ibnd ), 1, vkq(:,jbnd),1)  / &
                DDOT(3, vkk(:,ibnd), 1, vkk(:,ibnd),1)
          ENDDO
        ENDDO
      ENDIF
      !
      !DBSP
      !if (ik==3) THEN
      !  print*,'vkk(:, 2)',vkk(:, 2)
      !  print*,'vkq(:, 2)',vkq(:, 2)
      !ENDIF         
      !
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
           ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
        !
        fermicount = fermicount + 1
        DO imode = 1, nmodes
          !
          ! the phonon frequency
          wq = wf (imode, iq)
          !
          ! SP : We should avoid branching statements (if statements) in
          !      innerloops. Therefore we do it here.
          inv_wq =  1.0/(two * wq)
          ! the coupling from Gamma acoustic phonons is negligible
          IF ( wq .gt. eps_acustic ) THEN
            g2_tmp = 1.0
          ELSE
            g2_tmp = 0.0
          ENDIF   
          !
          DO ibnd = 1, ibndmax-ibndmin+1
            !
            !  the fermi occupation for k
            ekk = etf (ibndmin-1+ibnd, ikk) - ef0
            IF (delta_approx) THEN
              w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
            ELSE
              wgkk = wgauss( -ekk*inv_eptemp0, -99)
            ENDIF
            !
            DO jbnd = 1, ibndmax-ibndmin+1
              !
              !  the fermi occupation for k+q
              ekq = etf (ibndmin-1+jbnd, ikq) - ef0
              !
              ! here we take into account the zero-point sqrt(hbar/2M\omega)
              ! with hbar = 1 and M already contained in the eigenmodes
              ! g2 is Ry^2, wkf must already account for the spin factor
              !
              IF ( shortrange .AND. ( abs(xqf (1, iq))> eps8 .OR. abs(xqf (2, iq))> eps8 &
                 .OR. abs(xqf (3, iq))> eps8 )) THEN              
                ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary 
                !     number, in which case its square will be a negative number. 
                g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp ) 
              ELSE
                g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
              ENDIF
              !
              IF (delta_approx) THEN 
                !
                w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                ! the expression below is positive-definite, but also an
                ! approximation which neglects some fine features
                weight = pi * wq * wkf (ikk) * w0g1 * w0g2
                !
              ELSE
                !
                wgkq = wgauss( -ekq*inv_eptemp0, -99)
                !
                ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q + id]
                ! This is the imaginary part of minus the phonon self-energy, sans
                ! the matrix elements
                !
                !weight = wkf (ikk) * (wgkk - wgkq) * &
                !     aimag ( cone / ( ekq - ekk - wq - ci * degaussw0 ) )
                !
                ! SP: The expression below (minus phonon self-energy) corresponds to
                !  = pi*k-point weight*[f(E_k) - f(E_k+q)]*delta[E_k+q - E_k - w_q]
                weight = pi * wkf (ikk) * (wgkk - wgkq)* &
                     w0gauss ( (ekq - ekk - wq) / degaussw0, 0) / degaussw0
                !
              ENDIF  
              !
              gamma   (imode) = gamma   (imode) + weight * g2 
              gamma_v (imode) = gamma_v (imode) + weight * g2 * (1-coskkq(ibnd, jbnd) ) 
              !
            ENDDO ! jbnd
            !
          ENDDO   ! ibnd
          !
        ENDDO ! loop on q-modes
        !
      ENDIF ! endif fsthick
      !
    ENDDO ! loop on k
    !
    CALL stop_clock('PH SELF-ENERGY')
    !
    ! collect contributions from all pools (sum over k-points)
    ! this finishes the integral over the BZ  (k)
    !
    CALL mp_sum(gamma,inter_pool_comm) 
    CALL mp_sum(gamma_v,inter_pool_comm) 
    CALL mp_sum(fermicount, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    ! 
    ! An average over degenerate phonon-mode is performed. 
    DO imode = 1, nmodes
      n = 0
      tmp = 0.0_DP
      tmp2 = 0.0_DP
      tmp3 = 0.0_DP
      tmp4 = 0.0_DP
      wq = wf (imode, iq)
      DO jmode = 1, nmodes
        wq_tmp = wf (jmode, iq)
        IF ( ABS(wq - wq_tmp) < eps6 ) THEN
          n = n + 1
          IF ( wq_tmp .gt. eps_acustic ) THEN 
            tmp  =  tmp  + gamma  ( jmode ) / pi / wq**two / dosef
            tmp2 =  tmp2 + gamma_v( jmode ) / pi / wq**two / dosef
          ENDIF
          tmp3 =  tmp3 + gamma(jmode)
          tmp4 =  tmp4 + gamma_v(jmode)
        ENDIF
      ENDDO ! jbnd
      lambda_tmp(imode)   = tmp / float(n)
      lambda_v_tmp(imode) = tmp2 / float(n)
      gamma_tmp(imode)    = tmp3 / float(n)
      gamma_v_tmp(imode)  = tmp4 / float(n)
    ENDDO
    lambda_all( :, iq, ismear )   = lambda_tmp(:)
    lambda_v_all( :, iq, ismear ) = lambda_v_tmp(:)
    gamma_all( :, iq, ismear )    = gamma_tmp(:)
    gamma_v_all( :, iq, ismear )  = gamma_v_tmp(:)
    lambda_tot = sum(lambda_all(:,iq,ismear))
    lambda_tr_tot = sum(lambda_v_all(:,iq,ismear))
    !
    WRITE(stdout,'(/5x,"ismear = ",i5," iq = ",i7," coord.: ", 3f9.5, " wt: ", f9.5)') ismear, iq, xqf(:,iq), wqf(iq)
    WRITE(stdout,'(5x,a)') repeat('-',67)
    !
    DO imode = 1, nmodes
      ! 
      wq = wf (imode, iq)
      WRITE(stdout, 102) imode, lambda_all(imode,iq,ismear),ryd2mev*gamma_all(imode,iq,ismear), ryd2mev*wq
      WRITE(stdout, 104) imode, lambda_v_all(imode,iq,ismear),ryd2mev*gamma_v_all(imode,iq,ismear), ryd2mev*wq
      !
    ENDDO
    !
    WRITE(stdout, 103) lambda_tot
    WRITE(stdout, 105) lambda_tr_tot
    ! 
    IF (.NOT. specfun_ph) THEN
      WRITE(stdout,'(5x,a/)') repeat('-',67)
      WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
          'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkqtotf/2
    ENDIF
    !
  ENDDO !smears
  !
100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'lambda___( ',i3,' )=',f15.6,'   gamma___=',f15.6,' meV','   omega=',f12.4,' meV')
103 FORMAT(5x,'lambda___( tot )=',f15.6)
104 FORMAT(5x,'lambda_tr( ',i3,' )=',f15.6,'   gamma_tr=',f15.6,' meV','   omega=',f12.4,' meV')
105 FORMAT(5x,'lambda_tr( tot )=',f15.6)
  !
  RETURN
  !
END SUBROUTINE selfen_phon_q
! 
!-----------------------------------------------------------------------
FUNCTION dos_ef_seq (ngauss, degauss, ef, et, wk, nks, nbnd)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE mp,        ONLY : mp_sum
  IMPLICIT NONE
  REAL(DP) :: dos_ef_seq
  INTEGER :: ngauss, nbnd, nks
  REAL(DP) :: et (nbnd, nks), wk (nks), ef, degauss
  !
  INTEGER :: ik, ibnd
  REAL(DP), EXTERNAL :: w0gauss
  !
  !     Compute DOS at E_F (states per Ry per unit cell)
  !
  dos_ef_seq = 0.0d0
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        dos_ef_seq = dos_ef_seq + wk (ik) * w0gauss ( (et (ibnd, ik) - ef) &
             / degauss, ngauss) / degauss
     ENDDO
  ENDDO
  !
  RETURN
END FUNCTION dos_ef_seq

