!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_q0correction()
!-----------------------------------------------------------------------
   USE crpacom,              ONLY : v_coul_bare, v_coul_scrd, w_freq
   USE crpa_pert,            ONLY : pert_basis, npert_tot, pert_iwlist, pert_jwlist, pert_Rlist, &
                                   nmels_tot, mels_iwlist, mels_jwlist, mels_Rlist
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE control_lr,           ONLY : lgamma
   !
   IMPLICIT NONE
   ! 
   LOGICAL :: exst
   INTEGER :: n_freq
   REAL(DP), ALLOCATABLE :: omega(:), epsr_file(:), epsi_file(:)
   REAL(DP) :: epsr_int, epsi_int, delta_bare
   COMPLEX(DP) :: delta_scrd
   !
   IF (lgamma) THEN
      ! Set the dielectric tensor
      ! If a file eps.dat exist, read and set up the dielectric function
      INQUIRE(file="eps.dat", exist=exst)
      IF (exst) THEN  
         !
         WRITE (stdout,'(/,5X, "INFO: Dielectric tensor read from file eps.dat")')  
         CALL read_eps_dyn(n_freq, omega, epsr_file, epsi_file)
         !
         CALL eps_interpolate(w_freq, n_freq, omega, epsr_file, epsi_file, epsr_int, epsi_int)
         !
         WRITE (stdout,'(  5X, "INFO: interpolated macroscopic eps", 1F12.6)')  epsr_int, epsi_int
         DEALLOCATE(omega, epsr_file, epsi_file)
         !
         CALL divergence(epsr_int, epsi_int, delta_bare, delta_scrd)
         !
      ELSE 
         WRITE (stdout,'(/,5X, "INFO: Dielectric function = NOT SPECIFIED")')
         WRITE (stdout,'(  5X, "      NO Correction for Screened Coulomb, bare only")')
         !
         epsr_int = 1.D0
         epsi_int = 0.D0
         !
         CALL divergence(epsr_int, epsi_int, delta_bare, delta_scrd)
         delta_scrd = CMPLX(0.D0,0.D0)
         !
      ENDIF
      !
      CALL correction(delta_bare, delta_scrd)
      !  
   ENDIF
!
CONTAINS
!
SUBROUTINE read_eps_dyn (n_freq, omega, eps_re, eps_im) 
    !
    USE kinds,               ONLY : DP
    USE io_global,           ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(OUT)                   :: n_freq
    REAL(DP), ALLOCATABLE, INTENT(OUT)     :: omega(:), eps_re(:), eps_im(:)
    INTEGER            :: ios, i
    CHARACTER(LEN=256) :: line
    !
    OPEN (765, file = 'eps.dat', status='old', iostat=ios)
    IF (ios /= 0) THEN 
       WRITE(stdout,'(/, 5X,A, I5)') "ERROR: Cannot open eps.dat", ios
       STOP
    ENDIF
    !
    n_freq = 0
    DO 
        READ (765, '(A)', IOSTAT=ios) line
        IF (ios < 0) EXIT 
        IF (ios > 0) THEN
            WRITE(stdout,'(/, 5X,A, I5)') "ERROR: Reading eps.dat during count", ios
            STOP
        ENDIF
        
        line = ADJUSTL(line)
        IF (line(1:1) == '#' .OR. TRIM(line) == '') CYCLE
        n_freq = n_freq + 1
    ENDDO
    !
    ALLOCATE(omega(n_freq), eps_re(n_freq), eps_im(n_freq))
    !
    REWIND(765)
    !
    i = 0
    DO 
        READ (765, '(A)', IOSTAT=ios) line
        IF (ios < 0) EXIT
        ! 
        line = ADJUSTL(line)
        IF (line(1:1) == '#' .OR. TRIM(line) == '') CYCLE
        !
        i = i + 1
        READ (line, *, IOSTAT=ios) omega(i), eps_re(i), eps_im(i)
        !
        IF (ios /= 0) THEN 
            WRITE(stdout,'(/, 5X,A, I5)') "ERROR: Reading data on line", i
            STOP
        ENDIF
    ENDDO

    CLOSE (765)
    !
    WRITE(stdout,'(/, 5X,A, I5)') "Successfully read eps.dat. Number of frequencies: ", n_freq
    !
    RETURN
    !
END SUBROUTINE read_eps_dyn
 !
SUBROUTINE eps_interpolate(w_freq, n_freq, omega, eps_re_data, eps_im_data, eps_re_out, eps_im_out)
    !
    USE kinds, ONLY : DP
    USE constants,     ONLY : RYTOEV
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: w_freq
    INTEGER,  INTENT(IN)  :: n_freq
    REAL(DP), INTENT(IN)  :: omega(n_freq), eps_re_data(n_freq), eps_im_data(n_freq)
    REAL(DP), INTENT(OUT) :: eps_re_out, eps_im_out
    INTEGER  :: i
    REAL(DP) :: t, dw, w, target_w
    !
    target_w = DBLE(w_freq)*RYTOEV
    IF (target_w <= omega(1)) THEN
        eps_re_out = eps_re_data(1)
        eps_im_out = eps_im_data(1)
        RETURN
    ENDIF
    !
    IF (target_w >= omega(n_freq)) THEN
        eps_re_out = eps_re_data(n_freq)
        eps_im_out = eps_im_data(n_freq)
        RETURN
    ENDIF
    !
    DO i = 1, n_freq - 1
        IF (target_w >= omega(i) .AND. target_w <= omega(i+1)) THEN
            dw = omega(i+1) - omega(i)
            IF (dw < 1.0d-12) THEN
                eps_re_out = eps_re_data(i)
                eps_im_out = eps_im_data(i)
                RETURN
            ENDIF
            t = (target_w - omega(i)) / dw
            eps_re_out = eps_re_data(i) + t * (eps_re_data(i+1) - eps_re_data(i))
            eps_im_out = eps_im_data(i) + t * (eps_im_data(i+1) - eps_im_data(i))
            RETURN
        ENDIF
    ENDDO
    !
END SUBROUTINE eps_interpolate
 !
SUBROUTINE divergence(epsr_int, epsi_int, delta_bare, delta_scrd)
   USE kinds,             ONLY : DP
   USE crpacom,           ONLY : q0div_treatment
   !
   IMPLICIT NONE
   !
   REAL(DP) :: epsr_int, epsi_int, delta_bare
   COMPLEX(DP) :: delta_scrd
   !
   IF (TRIM(q0div_treatment) == 'HL') THEN
      CALL divergence_HL(epsr_int, epsi_int, delta_bare, delta_scrd)
   ELSEIF (TRIM(q0div_treatment) == 'GB') THEN 
      CALL divergence_GB(epsr_int, epsi_int, delta_bare, delta_scrd)
   ENDIF
   !
   RETURN
   !
END SUBROUTINE divergence
 !
 !-----------------------------------------------------------------------
SUBROUTINE divergence_HL(epsr_int, epsi_int, delta_bare, delta_scrd)
   ! Hybertsen Louie correction scheme of the q->0 divergence  
   !-----------------------------------------------------------------------
   USE kinds,             ONLY : DP
   USE constants,         ONLY : pi, RYTOEV
   USE cell_base,         ONLY : omega
   USE qpoint,            ONLY : nqs
   USE io_global,         ONLY : stdout
   !
   IMPLICIT NONE
   !
   REAL(DP) :: epsr_int, epsi_int, qsz, eps2, delta_bare
   COMPLEX(DP) :: epsinv, delta_scrd
   !
   eps2 = epsr_int**2.d0+epsi_int**2.d0
   epsinv = CMPLX(epsr_int/eps2, -epsi_int/eps2)
   qsz=(6*pi**2/(nqs*omega))**(1.0/3.0)
   !
   delta_bare = 2.0d0*nqs*(2.0d0/pi)*qsz
   delta_scrd = 2.0d0*nqs*(2.0d0/pi)*qsz*epsinv
   !
   WRITE (stdout,'(/,5X, "INFO: Hybertsen-Louie correction scheme")')
   WRITE (stdout,'(  5X, "INFO: qsz [Ha]               ", 3X, 1F12.6 )') qsz
   WRITE (stdout,'(  5X, "INFO: (2/pi)*qsz [Ha]        ", 3X, 1F12.6 )') (2.0d0/pi)*qsz
   WRITE (stdout,'(  5X, "INFO: Delta_bare q=0 [eV]        ", 3X, 1F12.6 )') delta_bare*RYTOEV
   WRITE (stdout,'(  5X, "INFO: Delta_scrd q=0 [eV]        ", 3X, 1F12.6 )') delta_scrd*RYTOEV
   !
   RETURN
   !
END SUBROUTINE divergence_HL
 !
 !-----------------------------------------------------------------------
SUBROUTINE divergence_GB(epsr_int, epsi_int, delta_bare, delta_scrd)
    ! Gygi-Baldereschi correction scheme of the q->0 divergence
    !-----------------------------------------------------------------------
    USE kinds,              ONLY : DP
    USE constants,          ONLY : fpi, e2, RYTOEV
    USE cell_base,          ONLY : bg, at, alat, omega
    USE gvect,              ONLY : ngm, g
    USE gvecw,              ONLY : gcutw
    USE control_flags,      ONLY : gamma_only
    USE mp_global,          ONLY : intra_pool_comm
    USE mp,                 ONLY : mp_sum
    USE qpoint,             ONLY : nq1, nq2, nq3, nqs
    USE io_global,          ONLY : stdout
    !
    !
    IMPLICIT NONE
    !
    INTEGER  :: iq1,iq2,iq3, ig
    REAL(DP) :: div,dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)
    INTEGER  :: nqq, iq
    REAL(DP) :: aa, dq
    !
    COMPLEX(DP) :: q_eps_q
    COMPLEX(DP) :: aa_eps, div_eps, det_eps, alpha_eps, dq_eps
    !
    REAL(DP) :: epsr_int, epsi_int, eps2, delta_bare
    COMPLEX(DP) :: delta_scrd, epsinv
    COMPLEX(DP) :: eps_mat(3,3)
    !
    REAL(DP)      :: grid_factor = 1.d0
    REAL(DP)      :: yukawa = 0.d0
    !
    CALL start_clock( 'exx_div' )
    !
    tpiba2 = ( fpi / 2.d0 / alat ) ** 2
    alpha  = 10.d0 / gcutw
    !
    eps2 = epsr_int**2.d0+epsi_int**2.d0
    epsinv = CMPLX(epsr_int/eps2, -epsi_int/eps2)
    alpha_eps = alpha*epsinv
    !
    dq1= 1.d0/DBLE(nq1)
    dq2= 1.d0/DBLE(nq2) 
    dq3= 1.d0/DBLE(nq3)
    !
    div = 0.d0
    div_eps = CMPLX(0.d0, 0.d0)
    !
    ! TODO: Full tensor case
    eps_mat(:,:) = CMPLX(0.D0,0.D0)
    eps_mat(1,1) = CMPLX(epsr_int, epsi_int) 
    eps_mat(2,2) = CMPLX(epsr_int, epsi_int) 
    eps_mat(3,3) = CMPLX(epsr_int, epsi_int)
    !
    DO iq1 = 1, nq1
      DO iq2 = 1, nq2
        DO iq3 = 1, nq3
          !
          xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                  bg(:,2) * (iq2-1) * dq2 + &
                  bg(:,3) * (iq3-1) * dq3
          !
          DO ig = 1, ngm
            !
            q(1) = xq(1) + g(1,ig)
            q(2) = xq(2) + g(2,ig)
            q(3) = xq(3) + g(3,ig)
            qq = q(1)*q(1) + q(2)*q(2) + q(3)*q(3)
            q_eps_q = DOT_PRODUCT( q(:), MATMUL( eps_mat(:,:), q(:) ) )
            !
            IF ( qq > 1.d-8 ) THEN
                div = div + EXP( - alpha * qq ) / ( qq + yukawa / tpiba2 ) &
                                                    * grid_factor
                div_eps = div_eps + EXP( - alpha_eps * q_eps_q ) / ( q_eps_q + yukawa / tpiba2 ) &
                                                    * grid_factor
            ENDIF
            !
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum( div, intra_pool_comm )
    CALL mp_sum( div_eps, intra_pool_comm )
    !
    IF ( gamma_only ) div = 2.d0 * div
    IF ( gamma_only ) div_eps = 2.d0 * div_eps
    !
    IF ( yukawa < 1.d-8) THEN
      div = div - alpha
      div_eps = div_eps - alpha_eps
    ELSE
      div = div + tpiba2 / yukawa
      div_eps = div_eps + tpiba2 / yukawa
    ENDIF
    !
    div = div * e2 * fpi / tpiba2 / nqs
    div_eps = div_eps * e2 * fpi / tpiba2 / nqs
    !
    alpha = alpha / tpiba2
    alpha_eps = alpha_eps / tpiba2
    !
    nqq = 100000
    dq = 5.0d0 / SQRT( alpha ) / nqq
    dq_eps = 5.0d0 / SQRT( alpha_eps ) / nqq
    aa = 0.d0
    aa_eps = 0.d0
    !
    DO iq = 0, nqq
      !
      q_ = dq * ( iq + 0.5d0 )
      qq = q_ * q_
      aa = aa - EXP( - alpha * qq ) * yukawa / ( qq + yukawa ) * dq
      q_ = dq_eps * ( iq + 0.5d0 )
      qq = q_ * q_
      aa_eps = aa_eps - EXP( - alpha_eps * q_eps_q ) * yukawa / ( q_eps_q + yukawa ) * dq_eps
      !
    ENDDO
    !
    det_eps = eps_mat(1,1)*(eps_mat(2,2)*eps_mat(3,3)-eps_mat(2,3)*eps_mat(3,2)) & 
             -eps_mat(2,2)*(eps_mat(2,1)*eps_mat(3,3)-eps_mat(2,3)*eps_mat(3,1)) &
             +eps_mat(3,3)*(eps_mat(2,1)*eps_mat(3,2)-eps_mat(2,2)*eps_mat(3,1)) 
    !
    aa = aa * 8.d0 / fpi
    aa_eps = aa_eps * 8.d0 / fpi
    aa = aa + 1.d0 / SQRT( alpha * 0.25d0 * fpi )
    aa_eps = aa_eps + 1.d0 / SQRT( alpha_eps * 0.25d0 * fpi * det_eps)
    !  
    div = div - e2 * omega * aa
    div_eps = div_eps - e2 * omega * aa_eps
    !
    delta_bare = -(nqs*div)/omega
    delta_scrd = -(nqs*div_eps)/omega
    !
    WRITE (stdout,'(/,5X, "INFO: Gygi-Baldereschi correction scheme")')
    WRITE (stdout,'(  5X, "INFO: Delta_bare q=0 [eV]        ", 3X, 1F12.6 )') delta_bare*RYTOEV
    WRITE (stdout,'(  5X, "INFO: Delta_scrd q=0 [eV]        ", 3X, 1F12.6 )') delta_scrd*RYTOEV
    !
    CALL stop_clock( 'exx_div' )
    !
    RETURN
    !
END SUBROUTINE divergence_GB 
 !
SUBROUTINE correction(delta_bare, delta_scrd)
   !
   USE crpa_pert,         ONLY : pert_basis
   USE kinds,             ONLY : DP
   !
   IMPLICIT NONE
   !
   REAL(DP) :: delta_bare
   COMPLEX(DP) :: delta_scrd
   !
   IF (pert_basis == 'bands') THEN
      CALL correction_bands(delta_bare, delta_scrd)
   ELSEIF (pert_basis == 'wannier') THEN
      CALL correction_wannier(delta_bare, delta_scrd)
   ENDIF
   !
END SUBROUTINE correction
   !
SUBROUTINE correction_bands(delta_bare, delta_scrd)
   USE crpa_pert,         ONLY : npert_tot, nmels_tot, pert_nbnd
   USE crpacom,           ONLY : v_coul_bare, v_coul_scrd
   USE kinds,             ONLY : DP
   !
   IMPLICIT NONE
   !
   INTEGER :: ib1, ib2, jb1, jb2, ik1, ik2, ipert1, ipert2
   REAL(DP) :: delta_bare
   COMPLEX(DP) :: delta_scrd
   !
   DO ipert2 = 1, npert_tot
      !
      ib2 = MOD(ipert2 - 1, pert_nbnd) + 1
      jb2 = MOD((ipert2 - 1) / pert_nbnd, pert_nbnd) + 1
      ik2 = (ipert2 - 1) / (pert_nbnd * pert_nbnd) + 1
      !
      DO ipert1 = 1, npert_tot
         !
         ib1 = MOD(ipert1 - 1, pert_nbnd) + 1
         jb1 = MOD((ipert1 - 1) / pert_nbnd, pert_nbnd) + 1
         ik1 = (ipert1 - 1) / (pert_nbnd * pert_nbnd) + 1
         !
         IF ( (ib2 == jb2) .AND. (ib1 == jb1) ) THEN
            !
            v_coul_bare(ipert1, ipert2)  = v_coul_bare(ipert1, ipert2)  + delta_bare
            v_coul_scrd(ipert1, ipert2)  = v_coul_scrd(ipert1, ipert2)  + delta_scrd
            !
         ENDIF
      ENDDO
   ENDDO
   !
END SUBROUTINE correction_bands
   !
SUBROUTINE correction_wannier(delta_bare, delta_scrd)
   !
   USE crpa_pert,         ONLY : npert_tot, pert_iwlist, pert_jwlist, pert_Rlist, &
                                 nmels_tot, mels_iwlist, mels_jwlist, mels_Rlist
   USE crpacom,           ONLY : v_coul_bare, v_coul_scrd
   USE kinds,             ONLY : DP
   !
   IMPLICIT NONE
   !
   INTEGER :: iw1, iw2, jw1, jw2, R1(3), R2(3), ipert1, ipert2
   REAL(DP) :: delta_bare
   COMPLEX(DP) :: delta_scrd
   !
   DO ipert2 = 1, npert_tot
      !
      iw2 = pert_iwlist(ipert2)
      jw2 = pert_jwlist(ipert2)
      R2  = pert_Rlist(:, ipert2)
      !
      DO ipert1 = 1, nmels_tot
         !
         iw1 = mels_iwlist(ipert1)
         jw1 = mels_jwlist(ipert1)
         R1  = mels_Rlist(:, ipert1)
         !
         IF ( (iw2 == jw2) .AND. ALL(R2 == 0) .AND. &
              (iw1 == jw1) .AND. ALL(R1 == 0) ) THEN
            !
            v_coul_bare(ipert1, ipert2)  = v_coul_bare(ipert1, ipert2)  + delta_bare
            v_coul_scrd(ipert1, ipert2)  = v_coul_scrd(ipert1, ipert2)  + delta_scrd
            !
         ENDIF
      ENDDO
   ENDDO
  !
END SUBROUTINE correction_wannier
  !
END SUBROUTINE crpa_q0correction
