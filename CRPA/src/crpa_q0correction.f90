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
   !! q=0 correction for the screened Coulomb matrix elements. 
   !! If epsilon_w.dat is present the macroscopic dielectric tensor is 
   !! read from it, otherwise it is computed (if crpa_epsil=.true.) or
   !! the correction is set to zero (if crpa_epsil=.false.).
   !! The correction is then computed and applied either with the 
   !! Hybertsen-Louie (HL) or Gygi-Baldereschi (GB) scheme.
   !
   USE kinds,             ONLY : DP
   USE io_global,         ONLY : stdout
   USE constants,         ONLY : RYTOEV
   USE crpacom,           ONLY : crpa_llanczos, crpa_epsil, crpa_nfreq, crpa_freq, &
                                  epsilon_w, epsilon_w_diel, crpa_done_epsil,      &
                                  v_coul_bare, v_coul_scrd, q0div_treatment
   USE crpa_pert,         ONLY : pert_basis
   !
   IMPLICIT NONE
   !
   LOGICAL :: exst, save_done_epsil
   CHARACTER(LEN=32) :: filename
   INTEGER :: ifreq
   COMPLEX(DP), ALLOCATABLE :: eps_mat(:,:,:), delta_scrd(:)
   REAL(DP) :: delta_bare, bare_now
   LOGICAL :: lfreq
   !! epsilon is frequency resolved, in epsilon_w_diel (Lanczos, or the
   !! frequency grid of the CG path)
   !
   lfreq = crpa_llanczos .OR. crpa_nfreq > 1
   !
   ALLOCATE(eps_mat(3,3,crpa_nfreq), delta_scrd(crpa_nfreq))
   !
   filename = 'epsilon_w.dat'
   INQUIRE(file=TRIM(filename), exist=exst)
   IF (exst) THEN
      !
      WRITE(stdout,'(/,5X, "INFO: Dielectric tensor read from file ",A)') TRIM(filename)
      CALL read_eps_file(TRIM(filename), crpa_nfreq, crpa_freq, eps_mat)
      IF (.NOT. crpa_llanczos) crpa_epsil = .false.
      !
   ELSE IF (crpa_epsil) THEN
      !
      WRITE(stdout,'(/,5X, "INFO: Calculation of the dielectric tensor")')
      !
      IF (crpa_llanczos) THEN
         save_done_epsil = crpa_done_epsil
         crpa_done_epsil = .false.
      ENDIF
      !
      CALL crpa_phescf_w()
      !
      IF (crpa_llanczos) crpa_done_epsil = save_done_epsil
      IF (lfreq) THEN
         eps_mat = epsilon_w_diel
      ELSE
         eps_mat(:,:,1) = epsilon_w
      ENDIF
      !
   ELSE
      !
      WRITE(stdout,'(/,5X, "INFO: Dielectric function = NOT SPECIFIED")')
      WRITE(stdout,'(  5X, "      NO Correction for Screened Coulomb, bare only")')
      eps_mat = (0.0_DP, 0.0_DP)
      DO ifreq = 1, crpa_nfreq
         eps_mat(1,1,ifreq) = (1.0_DP, 0.0_DP)
         eps_mat(2,2,ifreq) = (1.0_DP, 0.0_DP)
         eps_mat(3,3,ifreq) = (1.0_DP, 0.0_DP)
      ENDDO
      !
   ENDIF
   !
   ! Evaluation of the correction
   !
   IF (TRIM(q0div_treatment) == 'HL') THEN
      CALL divergence_HL(crpa_nfreq, eps_mat, delta_bare, delta_scrd)
   ELSEIF (TRIM(q0div_treatment) == 'GB') THEN
      CALL divergence_GB(crpa_nfreq, eps_mat, delta_bare, delta_scrd)
   ENDIF
   IF (.NOT. exst .AND. .NOT. crpa_epsil) delta_scrd(:) = (0.0_DP, 0.0_DP)
   !
   WRITE (stdout,'(  5X, "INFO: Delta_bare q=0 [eV]        ", 3X, 1F12.6 )') delta_bare*RYTOEV
   WRITE (stdout,'(  5X, "INFO: Delta_scrd q=0 [eV] vs frequency (Ry)")')
   DO ifreq = 1, crpa_nfreq
      WRITE (stdout,'(  5X, "      w = (",F12.6,",",F12.6,") : Re = ",F12.6,3X,"Im = ",F12.6 )') &
         DBLE(crpa_freq(ifreq)), AIMAG(crpa_freq(ifreq)), &
         DBLE(delta_scrd(ifreq))*RYTOEV, AIMAG(delta_scrd(ifreq))*RYTOEV
   ENDDO
   !
   ! Correction is applied to the Coulomb matrix elements
   !
   DO ifreq = 1, crpa_nfreq
      !
      bare_now = delta_bare
      IF (ifreq > 1) bare_now = 0.0_DP
      !
      IF (pert_basis == 'bands') THEN
         CALL correction_bands(bare_now, delta_scrd(ifreq), v_coul_bare, v_coul_scrd(:,:,ifreq))
      ELSE
         CALL correction_wannier(bare_now, delta_scrd(ifreq), v_coul_bare, v_coul_scrd(:,:,ifreq))
      ENDIF
      !
   ENDDO
   !
   DEALLOCATE(eps_mat, delta_scrd)
   !
   RETURN
   !
CONTAINS
!
 !-----------------------------------------------------------------------
SUBROUTINE divergence_HL(nfreq, eps_mat, delta_bare, delta_scrd)
   ! Hybertsen Louie correction scheme of the q->0 divergence. 
   !-----------------------------------------------------------------------
   USE kinds,             ONLY : DP
   USE constants,         ONLY : pi, RYTOEV
   USE cell_base,         ONLY : omega
   USE qpoint,            ONLY : nq1, nq2, nq3
   USE io_global,         ONLY : stdout
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: nfreq
   COMPLEX(DP), INTENT(IN) :: eps_mat(3,3,nfreq)
   REAL(DP), INTENT(OUT) :: delta_bare
   COMPLEX(DP), INTENT(OUT) :: delta_scrd(nfreq)
   !
   INTEGER :: nqs, ifreq
   REAL(DP) :: epsr_avg, epsi_avg, qsz, eps2
   COMPLEX(DP) :: epsinv
   !
   nqs = nq1 * nq2 * nq3
   qsz = (6*pi**2/(nqs*omega))**(1.0/3.0)
   delta_bare = 2.0d0*nqs*(2.0d0/pi)*qsz
   !
   DO ifreq = 1, nfreq
      epsr_avg = ( DBLE(eps_mat(1,1,ifreq)) + DBLE(eps_mat(2,2,ifreq)) + DBLE(eps_mat(3,3,ifreq)) ) / 3.0d0
      epsi_avg = ( AIMAG(eps_mat(1,1,ifreq)) + AIMAG(eps_mat(2,2,ifreq)) + AIMAG(eps_mat(3,3,ifreq)) ) / 3.0d0
      eps2 = epsr_avg**2.d0 + epsi_avg**2.d0
      epsinv = CMPLX(epsr_avg/eps2, -epsi_avg/eps2)
      delta_scrd(ifreq) = 2.0d0*nqs*(2.0d0/pi)*qsz*epsinv
   ENDDO
   !
   WRITE (stdout,'(/,5X, "INFO: Hybertsen-Louie correction scheme")')
   WRITE (stdout,'(  5X, "INFO: qsz [Ha]               ", 3X, 1F12.6 )') qsz
   WRITE (stdout,'(  5X, "INFO: (2/pi)*qsz [Ha]        ", 3X, 1F12.6 )') (2.0d0/pi)*qsz
   !
   RETURN
   !
END SUBROUTINE divergence_HL
 !
 !-----------------------------------------------------------------------
SUBROUTINE divergence_GB(nfreq, eps_mat, delta_bare, delta_scrd)
    ! Gygi-Baldereschi correction scheme of the q->0 divergence. 
    !-----------------------------------------------------------------------
    USE kinds,              ONLY : DP
    USE constants,          ONLY : fpi, e2, RYTOEV
    USE cell_base,          ONLY : bg, alat, omega
    USE gvect,              ONLY : ngm, g
    USE gvecw,              ONLY : gcutw
    USE control_flags,      ONLY : gamma_only
    USE mp_global,          ONLY : intra_pool_comm
    USE mp,                 ONLY : mp_sum
    USE qpoint,             ONLY : nq1, nq2, nq3
    USE io_global,          ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nfreq
    COMPLEX(DP), INTENT(IN) :: eps_mat(3,3,nfreq)
    REAL(DP), INTENT(OUT) :: delta_bare
    COMPLEX(DP), INTENT(OUT) :: delta_scrd(nfreq)
    !
    INTEGER  :: iq1, iq2, iq3, ig, nqs, ifreq
    REAL(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, q(3)
    INTEGER  :: nqq, iq
    REAL(DP) :: aa, dq
    !
    COMPLEX(DP) :: q_eps_q
    COMPLEX(DP) :: aa_eps, div_eps, det_eps, alpha_eps, dq_eps
    COMPLEX(DP) :: eps_trace_inv
    !
    REAL(DP) :: grid_factor = 1.d0
    REAL(DP) :: yukawa = 0.d0
    !
    CALL start_clock( 'exx_div' )
    !
    nqs = nq1 * nq2 * nq3
    tpiba2 = ( fpi / 2.d0 / alat ) ** 2
    alpha  = 10.d0 / gcutw
    !
    dq1= 1.d0/DBLE(nq1)
    dq2= 1.d0/DBLE(nq2)
    dq3= 1.d0/DBLE(nq3)
    !
    ! Frequency-independent part: div, alpha, aa, delta_bare
    !
    div = 0.d0
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
            !
            IF ( qq > 1.d-8 ) THEN
                div = div + EXP( - alpha * qq ) / ( qq + yukawa / tpiba2 ) &
                                                    * grid_factor
            ENDIF
            !
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL mp_sum( div, intra_pool_comm )
    !
    IF ( gamma_only ) div = 2.d0 * div
    !
    IF ( yukawa < 1.d-8) THEN
      div = div - alpha
    ELSE
      div = div + tpiba2 / yukawa
    ENDIF
    !
    div = div * e2 * fpi / tpiba2 / nqs
    !
    alpha = alpha / tpiba2
    !
    nqq = 100000
    dq = 5.0d0 / SQRT( alpha ) / nqq
    aa = 0.d0
    !
    DO iq = 0, nqq
      q_ = dq * ( iq + 0.5d0 )
      qq = q_ * q_
      aa = aa - EXP( - alpha * qq ) * yukawa / ( qq + yukawa ) * dq
    ENDDO
    !
    aa = aa * 8.d0 / fpi
    aa = aa + 1.d0 / SQRT( alpha * 0.25d0 * fpi )
    !
    div = div - e2 * omega * aa
    delta_bare = -(nqs*div)/omega
    !
    ! Frequency-dependent part: div_eps, alpha_eps, aa_eps, det_eps, delta_scrd
    !
    DO ifreq = 1, nfreq
       !
       eps_trace_inv = 3.d0 / ( eps_mat(1,1,ifreq) + eps_mat(2,2,ifreq) + eps_mat(3,3,ifreq) )
       alpha_eps = (10.d0 / gcutw) * eps_trace_inv
       !
       div_eps = CMPLX(0.d0, 0.d0)
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
               q_eps_q = DOT_PRODUCT( q(:), MATMUL( eps_mat(:,:,ifreq), q(:) ) )
               !
               IF ( qq > 1.d-8 ) THEN
                   div_eps = div_eps + EXP( - alpha_eps * q_eps_q ) / ( q_eps_q + yukawa / tpiba2 ) &
                                                       * grid_factor
               ENDIF
               !
             ENDDO
           ENDDO
         ENDDO
       ENDDO
       !
       CALL mp_sum( div_eps, intra_pool_comm )
       !
       IF ( gamma_only ) div_eps = 2.d0 * div_eps
       !
       IF ( yukawa < 1.d-8) THEN
         div_eps = div_eps - alpha_eps
       ELSE
         div_eps = div_eps + tpiba2 / yukawa
       ENDIF
       !
       div_eps = div_eps * e2 * fpi / tpiba2 / nqs
       !
       alpha_eps = alpha_eps / tpiba2
       dq_eps = 5.0d0 / SQRT( alpha_eps ) / nqq
       aa_eps = 0.d0
       !
       DO iq = 0, nqq
         q_ = dq_eps * ( iq + 0.5d0 )
         qq = q_ * q_
         aa_eps = aa_eps - EXP( - alpha_eps * q_eps_q ) * yukawa / ( q_eps_q + yukawa ) * dq_eps
       ENDDO
       !
       det_eps = eps_mat(1,1,ifreq)*(eps_mat(2,2,ifreq)*eps_mat(3,3,ifreq)-eps_mat(2,3,ifreq)*eps_mat(3,2,ifreq)) &
                -eps_mat(2,2,ifreq)*(eps_mat(2,1,ifreq)*eps_mat(3,3,ifreq)-eps_mat(2,3,ifreq)*eps_mat(3,1,ifreq)) &
                +eps_mat(3,3,ifreq)*(eps_mat(2,1,ifreq)*eps_mat(3,2,ifreq)-eps_mat(2,2,ifreq)*eps_mat(3,1,ifreq))
       !
       aa_eps = aa_eps * 8.d0 / fpi
       aa_eps = aa_eps + 1.d0 / SQRT( alpha_eps * 0.25d0 * fpi * det_eps)
       !
       div_eps = div_eps - e2 * omega * aa_eps
       delta_scrd(ifreq) = -(nqs*div_eps)/omega
       !
    ENDDO ! ifreq
    !
    WRITE (stdout,'(/,5X, "INFO: Gygi-Baldereschi correction scheme")')
    !
    CALL stop_clock( 'exx_div' )
    !
    RETURN
    !
END SUBROUTINE divergence_GB
 !
SUBROUTINE correction_bands(delta_bare, delta_scrd, v_coul_bare, v_coul_scrd)
   USE crpa_pert,         ONLY : npert_tot, nmels_tot, pert_nbnd
   USE kinds,             ONLY : DP
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: delta_bare
   COMPLEX(DP), INTENT(IN) :: delta_scrd
   COMPLEX(DP), INTENT(INOUT) :: v_coul_bare(nmels_tot, npert_tot)
   COMPLEX(DP), INTENT(INOUT) :: v_coul_scrd(nmels_tot, npert_tot)
   !
   INTEGER :: ib1, ib2, jb1, jb2, ipert1, ipert2
   !
   DO ipert2 = 1, npert_tot
      !
      ib2 = MOD(ipert2 - 1, pert_nbnd) + 1
      jb2 = MOD((ipert2 - 1) / pert_nbnd, pert_nbnd) + 1
      !
      DO ipert1 = 1, nmels_tot
         !
         ib1 = MOD(ipert1 - 1, pert_nbnd) + 1
         jb1 = MOD((ipert1 - 1) / pert_nbnd, pert_nbnd) + 1
         !
         IF ( (ib2 == jb2) .AND. (ib1 == jb1) ) THEN
            !
            v_coul_bare(ipert1, ipert2) = v_coul_bare(ipert1, ipert2) + delta_bare
            v_coul_scrd(ipert1, ipert2) = v_coul_scrd(ipert1, ipert2) + delta_scrd
            !
         ENDIF
      ENDDO
   ENDDO
   !
END SUBROUTINE correction_bands
   !
SUBROUTINE correction_wannier(delta_bare, delta_scrd, v_coul_bare, v_coul_scrd)
   USE crpa_pert,         ONLY : npert_tot, pert_iwlist, pert_jwlist, pert_Rlist, &
                                 nmels_tot, mels_iwlist, mels_jwlist, mels_Rlist
   USE kinds,             ONLY : DP
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: delta_bare
   COMPLEX(DP), INTENT(IN) :: delta_scrd
   COMPLEX(DP), INTENT(INOUT) :: v_coul_bare(nmels_tot, npert_tot)
   COMPLEX(DP), INTENT(INOUT) :: v_coul_scrd(nmels_tot, npert_tot)
   !
   INTEGER :: iw1, iw2, jw1, jw2, R1(3), R2(3), ipert1, ipert2
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
            v_coul_bare(ipert1, ipert2) = v_coul_bare(ipert1, ipert2) + delta_bare
            v_coul_scrd(ipert1, ipert2) = v_coul_scrd(ipert1, ipert2) + delta_scrd
            !
         ENDIF
      ENDDO
   ENDDO
  !
END SUBROUTINE correction_wannier
  !
END SUBROUTINE crpa_q0correction
!
!-----------------------------------------------------------------------
SUBROUTINE write_eps_file(filename, nfreq, freq, eps_mat)
   !-----------------------------------------------------------------------
   USE kinds,      ONLY : DP
   USE io_global,  ONLY : stdout, ionode
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=*), INTENT(IN) :: filename
   INTEGER, INTENT(IN) :: nfreq
   COMPLEX(DP), INTENT(IN) :: freq(nfreq)
   COMPLEX(DP), INTENT(IN) :: eps_mat(3,3,nfreq)
   !
   INTEGER :: iunit, ios, ifreq
   !
   IF (.NOT. ionode) RETURN
   !
   OPEN (NEWUNIT=iunit, FILE=TRIM(filename), STATUS='replace', IOSTAT=ios)
   IF (ios /= 0) CALL errore('write_eps_file', 'cannot open '//TRIM(filename), 1)
   !
   WRITE(iunit, '("# w0 (Ry)    eta (Ry)     Re(e11)        Im(e11)        &
        &Re(e22)        Im(e22)        Re(e33)        Im(e33)        &
        &Re(e12)        Im(e12)        Re(e13)        Im(e13)        &
        &Re(e23)        Im(e23)")')
   !
   DO ifreq = 1, nfreq
      WRITE(iunit, '(2f12.6,12e15.7)') DBLE(freq(ifreq)), AIMAG(freq(ifreq)), &
           DBLE(eps_mat(1,1,ifreq)), AIMAG(eps_mat(1,1,ifreq)), &
           DBLE(eps_mat(2,2,ifreq)), AIMAG(eps_mat(2,2,ifreq)), &
           DBLE(eps_mat(3,3,ifreq)), AIMAG(eps_mat(3,3,ifreq)), &
           DBLE(eps_mat(1,2,ifreq)), AIMAG(eps_mat(1,2,ifreq)), &
           DBLE(eps_mat(1,3,ifreq)), AIMAG(eps_mat(1,3,ifreq)), &
           DBLE(eps_mat(2,3,ifreq)), AIMAG(eps_mat(2,3,ifreq))
   ENDDO
   !
   CLOSE(iunit)
   !
   WRITE(stdout,'(5X,A)') "INFO: Dielectric tensor written to "//TRIM(filename)
   !
   RETURN
   !
END SUBROUTINE write_eps_file
!
!-----------------------------------------------------------------------
SUBROUTINE read_eps_file(filename, nfreq, freq, eps_mat)
   !-----------------------------------------------------------------------
   !! Read the complex dielectric tensor from filename (same row format as
   !! write_eps_file above, one header line then nfreq rows: e21=e12,
   !! e31=e13, e32=e23 by symmetry). Errors out if the file's frequency
   !! grid does not match the caller-supplied freq(nfreq) (fewer/more rows,
   !! or a row whose (w0, eta) differ from freq(ifreq) beyond a small
   !! tolerance).
   !
   USE kinds,      ONLY : DP
   USE io_global,  ONLY : stdout
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=*), INTENT(IN) :: filename
   INTEGER, INTENT(IN) :: nfreq
   COMPLEX(DP), INTENT(IN) :: freq(nfreq)
   COMPLEX(DP), INTENT(OUT) :: eps_mat(3,3,nfreq)
   !
   INTEGER :: iunit, ios, ifreq
   REAL(DP) :: w0, eta
   REAL(DP) :: re11, im11, re22, im22, re33, im33, re12, im12, re13, im13, re23, im23
   CHARACTER(LEN=512) :: line
   REAL(DP), PARAMETER :: freq_tol = 1.0D-4
   !
   OPEN(NEWUNIT=iunit, FILE=TRIM(filename), STATUS='old', IOSTAT=ios)
   IF (ios /= 0) CALL errore('read_eps_file', 'cannot open '//TRIM(filename), 1)
   !
   READ(iunit, '(A)', IOSTAT=ios) line
   !
   DO ifreq = 1, nfreq
      !
      READ(iunit, *, IOSTAT=ios) w0, eta, re11, im11, re22, im22, re33, im33, &
                                  re12, im12, re13, im13, re23, im23
      IF (ios /= 0) CALL errore('read_eps_file', &
         'error reading '//TRIM(filename)//' (fewer frequencies than expected?)', ifreq)
      !
      IF (ABS(w0 - DBLE(freq(ifreq))) > freq_tol .OR. &
          ABS(eta - AIMAG(freq(ifreq))) > freq_tol) THEN
         WRITE(stdout,'(5X,"ERROR: frequency mismatch at entry ",I6,": file has (", &
              &F12.6,",",F12.6,"), expected (",F12.6,",",F12.6,")")') &
            ifreq, w0, eta, DBLE(freq(ifreq)), AIMAG(freq(ifreq))
         CALL errore('read_eps_file', &
            TRIM(filename)//' frequency grid does not match the expected one', 1)
      ENDIF
      !
      eps_mat(1,1,ifreq) = CMPLX(re11, im11, KIND=DP)
      eps_mat(2,2,ifreq) = CMPLX(re22, im22, KIND=DP)
      eps_mat(3,3,ifreq) = CMPLX(re33, im33, KIND=DP)
      eps_mat(1,2,ifreq) = CMPLX(re12, im12, KIND=DP)
      eps_mat(2,1,ifreq) = eps_mat(1,2,ifreq)
      eps_mat(1,3,ifreq) = CMPLX(re13, im13, KIND=DP)
      eps_mat(3,1,ifreq) = eps_mat(1,3,ifreq)
      eps_mat(2,3,ifreq) = CMPLX(re23, im23, KIND=DP)
      eps_mat(3,2,ifreq) = eps_mat(2,3,ifreq)
      !
   ENDDO
   !
   CLOSE(iunit)
   !
   RETURN
   !
END SUBROUTINE read_eps_file
