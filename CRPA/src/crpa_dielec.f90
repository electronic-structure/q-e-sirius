!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine crpa_dielec()
  !-----------------------------------------------------------------------
  !! Calculates the dielectric tensor.
  !
  USE kinds, only : DP
  USE io_global,  ONLY : stdout
  USE constants, ONLY: fpi
  USE cell_base, ONLY: omega
  USE klist, ONLY: wk
  USE wvfct, ONLY: npwx
  USE buffers, ONLY : get_buffer
  USE noncollin_module, ONLY : npol
  USE units_ph, ONLY : lrebar, iuebar
  USE units_lr, ONLY : lrdwf, iudwf
  USE eqv, ONLY : dpsi, dvpsi
  USE qpoint, ONLY : nksq, ikks
  USE control_lr, ONLY : nbnd_occ
  USE mp_pools,   ONLY : inter_pool_comm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE symme,      ONLY : symmatrix, crys_to_cart
  USE crpacom,    ONLY : epsilon_w, crpa_epsil, crpa_done_epsil, current_w, &
                          plasma_freq2, crpa_metal, crpa_diel_intraband, &
                          w_drude_min, crpa_nfreq

  implicit none

  integer :: ibnd, ipol, jpol, nrec, ik, ikk
  real(DP) :: w, weight, repsilon_w(3,3), iepsilon_w(3,3), wdr
  LOGICAL :: finite_freq

  complex(DP), external :: zdotc

  IF (.NOT.crpa_epsil.OR.crpa_done_epsil) RETURN

  finite_freq = (ABS(current_w) > 1.0d-10)

  call start_clock ('dielec')
  !
  ! Intraband (Drude) plasma-frequency tensor: zero, and crpa_metal=.false.,
  ! for insulators. See crpa_plasma_frequency.f90 for why this piece is
  ! missing from (and must not be added to) the interband sum below.
  !
  CALL crpa_plasma_frequency()
  IF (crpa_metal .AND. crpa_diel_intraband .AND. .NOT. finite_freq) CALL errore('crpa_dielec', &
     'a metal requires a finite complex frequency (with Im(w) > 0): '// &
     'the intraband/Drude term diverges at w_freq = 0', 1)
  !
  epsilon_w(:,:) = 0.d0
  do ik = 1, nksq
     ikk=ikks(ik)
     weight = wk (ikk)
     w = fpi * weight / omega
     do ipol = 1, 3
        nrec = (ipol - 1) * nksq + ik
        call get_buffer(dvpsi, lrebar, iuebar, nrec)
        do jpol = 1, 3
           nrec = (jpol - 1) * nksq + ik
           call get_buffer (dpsi, lrdwf, iudwf, nrec)
           do ibnd = 1, nbnd_occ (ikk)
              !
              !  this is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
              !
              epsilon_w(ipol,jpol)=epsilon_w(ipol,jpol)-2.d0*w* ( &
                   zdotc(npwx*npol, dvpsi(1,ibnd), 1, dpsi(1,ibnd), 1) )
           enddo
           IF (finite_freq) THEN
              call get_buffer (dpsi, lrdwf, iudwf, nrec + 3 * nksq)
              do ibnd = 1, nbnd_occ (ikk)
                 epsilon_w(ipol,jpol)=epsilon_w(ipol,jpol)-2.d0*w* ( &
                   zdotc(npwx*npol, dvpsi(1,ibnd), 1, dpsi(1,ibnd), 1) )
              enddo
           ENDIF
        enddo
     enddo
  enddo
  !
  IF (.NOT.finite_freq) epsilon_w(:,:) = epsilon_w(:,:) * 2.d0
  !
  call mp_sum ( epsilon_w, intra_bgrp_comm )
  call mp_sum ( epsilon_w, inter_pool_comm )
  !
  !      symmetrize
  !
  !WRITE( stdout,'(/,10x,"Unsymmetrized in cartesian axis ",/)')
  !do ipol = 1, 3
  !   WRITE( stdout,'(10x,3("(",f12.6,",",f12.6,")",2x))') &
  !        (DBLE(epsilon_w(ipol,jpol)), AIMAG(epsilon_w(ipol,jpol)), jpol=1,3)
  !enddo
  !
  repsilon_w = DBLE(epsilon_w)
  call crys_to_cart (repsilon_w)
  call symmatrix(repsilon_w)
  iepsilon_w = AIMAG(epsilon_w)
  call crys_to_cart (iepsilon_w)
  call symmatrix(iepsilon_w)
  epsilon_w = CMPLX(repsilon_w, iepsilon_w, KIND=DP)
  !
  !WRITE( stdout,'(/,10x,"Symmetrized in cartesian axis ",/)')
  !do ipol = 1, 3
  !   WRITE( stdout,'(10x,3("(",f12.6,",",f12.6,")",2x))') &
  !        (DBLE(epsilon_w(ipol,jpol)), AIMAG(epsilon_w(ipol,jpol)), jpol=1,3)
  !enddo
  !
  ! Intraband (Drude) contribution, metals only: - omega_p^2 / ( w * (w+i*eta) ),
  ! w = Re(current_w), eta = Im(current_w)
  !
  IF (crpa_metal .AND. crpa_diel_intraband) THEN
     wdr = DBLE(current_w)
     IF (ABS(wdr) < w_drude_min) THEN
        wdr = w_drude_min
        WRITE(stdout,'(5x,A,ES9.2,A)') 'Drude term: Re(w) = 0; denominator floored to Re(w) = ', w_drude_min,' Ry.'
     ENDIF
     do ipol = 1, 3
        do jpol = 1, 3
           epsilon_w(ipol,jpol) = epsilon_w(ipol,jpol) &
              - CMPLX(plasma_freq2(ipol,jpol), 0.d0, KIND=DP) &
                / ( CMPLX(wdr, 0.d0, KIND=DP) * CMPLX(wdr, AIMAG(current_w), KIND=DP) )
        enddo
     enddo
  ENDIF
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     epsilon_w(ipol,ipol) = epsilon_w(ipol,ipol) + 1.d0
  enddo
  !
  crpa_done_epsil=.TRUE.
  CALL crpa_summarize_epsilon()
  !
  ! frequency grid: crpa_phescf_w writes all the frequencies at the end
  !
  IF (crpa_nfreq == 1) THEN
     CALL crpa_write_epsilon_w()
     CALL write_eps_file('epsilon_w.dat', 1, (/current_w/), RESHAPE(epsilon_w, (/3,3,1/)))
  ENDIF

  call stop_clock ('dielec')

  return

end subroutine crpa_dielec
