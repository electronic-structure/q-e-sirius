  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine nesting_fn_q(iqq, iq)
  !-----------------------------------------------------------------------
  !!
  !!  compute the imaginary part of the phonon self energy due to electron-
  !!  phonon interaction in the Migdal approximation. This corresponds to 
  !!  the phonon linewidth (half width). The phonon frequency is taken into
  !!  account in the energy selection rule.
  !!
  !!  Use matrix elements, electronic eigenvalues and phonon frequencies
  !!  from ep-wannier interpolation. 
  !!
  !-----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : nbndsub, fsthick, &
                        eptemp, ngaussw, degaussw,     &
                        nsmear, delta_smear, efermi_read, fermi_energy
  USE pwcom,     ONLY : nelec, ef
  USE klist_epw, ONLY : isk_dummy
  USE elph2,     ONLY : ibndmax, ibndmin, etf, &
                        wkf, xqf, wqf, nkqf, &
                        nkf, nkqtotf, xqf
  USE constants_epw, ONLY : ryd2ev, two
  USE mp,        ONLY : mp_barrier,mp_sum
  USE mp_global, ONLY : inter_pool_comm
  !
  implicit none
  !
  INTEGER, INTENT (in) :: iqq
  !! Current q-point index from selecq
  INTEGER, INTENT (in) :: iq
  !! Current q-point index
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
  INTEGER :: fermicount
  !! Number of states on the Fermi surface
  INTEGER :: ismear
  !! Upper bounds index after k or q paral
  !! Smearing for the Gaussian function 
  ! 
  REAL(kind=DP) :: ekk
  !! Eigen energy on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: weight
  !! Imaginary part of the phonhon self-energy factor 
  REAL(kind=DP) :: dosef
  !! Density of state N(Ef)
  REAL(kind=DP) :: w0g1
  !! Dirac delta for the imaginary part of $\Sigma$
  REAL(kind=DP) :: w0g2
  !! Dirac delta for the imaginary part of $\Sigma$
  real(kind=DP) :: w0gauss, dos_ef, gamma, degaussw0
  real(kind=DP), external :: efermig
  !
  !
  IF (iqq == 1) then 
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Nesting Function in the double delta approx")')
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !
    IF ( fsthick < 1.d3 ) &
      WRITE(stdout, '(/5x,a,f10.6,a)' ) &
      'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
      'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
  ENDIF
  !
  ! SP: The Gamma function needs to be put to 0 for each q
  gamma = 0.0
  ! 
  ! Here we loop on smearing values
  DO ismear = 1, nsmear
    !
    degaussw0 = (ismear-1)*delta_smear+degaussw
    !
    ! Fermi level and corresponding DOS
    !
    !   Note that the weights of k+q points must be set to zero here
    !   no spin-polarized calculation here
    IF ( efermi_read ) THEN
      ef0 = fermi_energy 
    ELSE
      ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw0, ngaussw, 0, isk_dummy)
    ENDIF
    !
    dosef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nkqf, nbndsub)
    !  N(Ef) in the equation for lambda is the DOS per spin
    dosef = dosef / two
    !
    IF (iqq == 1) then
      WRITE (stdout, 100) degaussw0 * ryd2ev, ngaussw
      WRITE (stdout, 101) dosef / ryd2ev, ef0 * ryd2ev
    ENDIF
    !
    !
    CALL start_clock('nesting')
    !
    fermicount = 0
    !
    DO ik = 1, nkf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      ! 
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      IF ( ( minval ( abs(etf (:, ikk) - ef) ) < fsthick ) .and. &
          ( minval ( abs(etf (:, ikq) - ef) ) < fsthick ) ) then
        !
        fermicount = fermicount + 1
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          !
          ekk = etf (ibndmin-1+ibnd, ikk) - ef0
          w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
          !
          DO jbnd = 1, ibndmax-ibndmin+1
            !
            ekq = etf (ibndmin-1+jbnd, ikq) - ef0
            w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
            !
            ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
            ! This is the imaginary part of the phonon self-energy, sans the matrix elements
            !
            ! weight = wkf (ikk) * (wgkk - wgkq) * &
            !      aimag ( cone / ( ekq - ekk  - ci * degaussw ) ) 
            !
            ! the below expression is positive-definite, but also an approximation
            ! which neglects some fine features
            !
            weight = wkf (ikk) * w0g1 * w0g2
            !
            gamma  = gamma  + weight  
            !
          ENDDO ! jbnd
        ENDDO   ! ibnd
        !
      ENDIF ! endif fsthick
      !
    ENDDO ! loop on k
    !
    ! collect contributions from all pools (sum over k-points)
    ! this finishes the integral over the BZ  (k)
    !
    CALL mp_sum(gamma,inter_pool_comm) 
    CALL mp_sum(fermicount, inter_pool_comm)
    !
    WRITE(stdout,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:,iq) , wqf(iq)
    WRITE(stdout,'(5x,a)') repeat('-',67)
       ! 
    WRITE(stdout, 102)  gamma
    WRITE(stdout,'(5x,a/)') repeat('-',67)
    !
    WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
      'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkqtotf/2
    !
    !
    CALL stop_clock('nesting')
  ENDDO !smears
  !
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' eV, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/eV/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,' Nesting function (q)=',e15.6,' [Adimensional]')
  !
  END SUBROUTINE nesting_fn_q
  !
