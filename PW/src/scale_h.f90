!
! Copyright (C) 2001-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE scale_h
  !-----------------------------------------------------------------------
  !! When variable cell calculation are performed this routine scales the
  !! quantities needed in the calculation of the hamiltonian using the
  !! new and old cell parameters.
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout
  USE cell_base,      ONLY : bg, omega, set_h_ainv, tpiba
  USE cellmd,         ONLY : at_old, omega_old
  USE constants,      ONLY : eps8
  USE gvect,          ONLY : g, gg, ngm
  USE klist,          ONLY : xk, wk, nkstot, qnorm
  USE uspp_data,      ONLY : nqxq, dq, scale_uspp_data
  USE control_flags,  ONLY : iverbosity
  USE start_k,        ONLY : nks_start, xk_start, nk1,nk2,nk3
  USE exx_base,       ONLY : exx_grid_init, exx_mp_init
  USE exx,            ONLY : exx_gvec_reinit
  USE xc_lib,         ONLY : xclib_dft_is
  USE rism_module,    ONLY : lrism, rism_reinit3d
  USE mp,             ONLY : mp_max
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE beta_mod,       ONLY : scale_tab_beta
  USE vloc_mod,       ONLY : scale_tab_vloc
  USE rhoc_mod,       ONLY : scale_tab_rhc
  USE rhoat_mod,      ONLY : scale_tab_rhoat
  USE qrad_mod,       ONLY : scale_tab_qrad, init_tab_qrad
  USE mod_sirius
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, ik, ipol, ierr
  ! counters
  REAL(DP) :: gg_max, qmax
  !
  ! scale the k points
  !
  CALL cryst_to_cart( nkstot, xk, at_old, - 1 )
  CALL cryst_to_cart( nkstot, xk, bg, + 1 )
  IF (nks_start>0) THEN
    CALL cryst_to_cart( nks_start, xk_start, at_old, - 1 )
    CALL cryst_to_cart( nks_start, xk_start, bg, + 1 )
  ENDIF
  !
  ! Print new k-points only if given in input and if not Gamma
  !
  IF ( nk1==0 .AND. nk2==0 .AND. nk3 == 0 .AND. &
       ( nkstot > 1 .OR. ABS(xk(1,1)**2+xk(2,1)**2+xk(3,1)**2) > eps8 ) ) THEN
     IF ( iverbosity > 0 .OR. nkstot < 100 ) THEN
        WRITE( stdout, '(5x,a)' ) 'NEW k-points:'
        DO ik = 1, nkstot
           WRITE( stdout, '(3f12.7,f12.7)') (xk(ipol,ik) , ipol=1,3), wk(ik)
        ENDDO
     ELSE
        WRITE( stdout, '(5x,a)' ) "NEW k-points: use verbosity='high' to print them"
     ENDIF
  ENDIF
  !
  ! scale the g vectors (as well as gg and gl arrays)
  !
  CALL cryst_to_cart( ngm, g, at_old, - 1 )
  CALL cryst_to_cart( ngm, g, bg, + 1 )
  gg_max = 0.0_dp
  !
  DO ig = 1, ngm
     gg (ig) = g(1,ig) * g(1,ig) + g(2,ig) * g(2,ig) + g(3,ig) * g(3,ig)
     gg_max = MAX(gg(ig), gg_max)
  ENDDO
  !$acc update device(g,gg)
  !
  CALL mp_max( gg_max, intra_bgrp_comm )
  qmax = SQRT(gg_max)*tpiba
  ! qmax is the largest |G| actually needed in interpolation tables
  IF ( nqxq < INT(qmax/dq)+4 ) THEN
     CALL errore( 'scale_h', 'Not enough space allocated for radial FFT: '//&
                             'try restarting with a larger cell_factor.', 1 )
  ENDIF
  !
  ! scale the non-local pseudopotential tables
  !
#if defined(__SIRIUS)
  !IF ((use_sirius_scf.OR.use_sirius_nlcg).AND..NOT.use_veff_callback) THEN
  !  CONTINUE
  !ELSE
#endif
  !
  call scale_uspp_data( omega_old/omega )
  call scale_tab_beta( omega_old/omega )
  CALL scale_tab_rhc( omega_old/omega )
  CALL scale_tab_rhoat( omega_old/omega )
  CALL scale_tab_qrad( omega_old/omega )
  !
#if defined(__SIRIUS)
  !END IF
#endif
  !
  ! for hybrid functionals
  !
  IF ( xclib_dft_is('hybrid') ) THEN
     CALL exx_grid_init( reinit=.TRUE. )
     ! not sure next line is needed
     CALL exx_mp_init( )
     CALL exx_gvec_reinit( at_old )
     qmax = qmax + qnorm
     ! For USPP + hybrid, qmax is the largest |q+G| needed
     ! Note that qnorm is recomputed in exx_grid_init
  END IF
  !
  ! Check interpolation table, re-allocate if needed
  !
  CALL init_tab_qrad ( qmax, omega, intra_bgrp_comm, ierr)
  !
  ! recalculate the local part of the pseudopotential
  !
  CALL scale_tab_vloc( omega_old/omega )
  CALL init_vloc( )
  !
  ! for ts-vdw
  !
  CALL set_h_ainv()
  !
  ! for 3D-RISM
  !
  IF ( lrism ) CALL rism_reinit3d()
  !
  RETURN
  !
END SUBROUTINE scale_h

