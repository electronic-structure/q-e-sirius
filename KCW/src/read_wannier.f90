!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!-----------------------------------------------------------------------
subroutine read_wannier ()
  !-----------------------------------------------------------------------
  !
  !! This routine read the relevant infos from a previous Wannier90 run
  !! The unitary matrix U are stored in unimatrx, unimatrx, unimatrx_opt.
  !! U(i,j) = <KS_i|wann_j> 
  !! Two cases are addressed:
  !! 1) A standard wannierization (unique manifold)
  !! 2) A separate wannierization where occ and emp states are treated independently. 
  !! In both cases I define a unique U whose dimension is num_wann x num_wann 
  !! and eventually a U_opt whose dimension in num_bands x num_wann
  !
  USE control_kcw,          ONLY : l_unique_manifold
  !
  IMPLICIT NONE 
  !
  IF (l_unique_manifold) THEN 
    !
    ! ... standard wannierization (occ + emp all together) 
    CALL read_wannier_unique_manifold ( )
    !
  ELSE
   !
   ! ... Separate wannierization (No occ-emp mixing)
   CALL read_wannier_two_manifold ( )
   !
  ENDIF
  !
END subroutine read_wannier
  !
  !-----------------------------------------------------------------------
  subroutine read_wannier_unique_manifold ()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : unimatrx, seedname, has_disentangle, &
                                   unimatrx_opt, num_wann, kcw_iverbosity
  USE mp_global,            ONLY : intra_image_comm
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nkstot, xk
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : nbnd
  USE wannier90_readback,   ONLY : w90rb_read_umat
  !
  IMPLICIT NONE 
  !
  INTEGER :: i, nkstot_eff, nrows, ncols
  !
  !! Here and in W90 We treat always one spin component at the time.
  !! While in PWscf nkstot contain k point for spin up (the first nkstot/nspin 
  !! and spin down (the last nkstot.nspin). See also rotate_ks
  !
  IF (nspin == 4) THEN
      nkstot_eff = nkstot
  ELSE IF (nspin ==2)  THEN
      nkstot_eff = nkstot/nspin
  ELSE
      nkstot_eff = nkstot/nspin
  ENDIF
  !
  ! ... The unitary matrix
  !
  CALL w90rb_read_umat (TRIM(seedname)//'_u.mat', xk(:,1:nkstot_eff), intra_image_comm, &
                        unimatrx, nrows, num_wann)
  !
  ! ... The Optimal subspace Matrix (if any)
  !
  IF (has_disentangle) THEN
     !
     CALL w90rb_read_umat (TRIM(seedname)//'_u_dis.mat', xk(:,1:nkstot_eff), intra_image_comm, &
                           unimatrx_opt, nrows, ncols)
     !
     IF (ncols /= num_wann) &
           CALL errore ('read_wannier', 'Mismatch in num_wann from U vs U Optimal matrix', 1)
     IF (nrows /=  nbnd ) &
          CALL errore ('read_wannier', 'Mismatch between num KS state from PW and Wann90', 1)
     !
  ELSE
     !
     ! ... If no disentangle than the U-opt is a square matrix of dimension num_wan x num_wann
     ! ... set it to the identity (in apply_u_mat we anyway apply both Ui_opt and U)
     !
     ALLOCATE (unimatrx_opt(num_wann, num_wann, nkstot_eff))
     unimatrx_opt = CMPLX(0.D0, 0.D0, kind=DP)
     DO i = 1, num_wann; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind = DP); ENDDO
     !
  ENDIF
  !
  if ( kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Optimal Matrix READ")') 
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  !
  RETURN
  !
END subroutine read_wannier_unique_manifold 
  !
  !
  !-----------------------------------------------------------------------
  subroutine read_wannier_two_manifold ()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : unimatrx, have_empty, num_wann_occ, seedname, &
                                   has_disentangle, have_empty, num_wann_emp, & 
                                   unimatrx_opt, num_wann, kcw_iverbosity, spin_component
  USE mp_global,            ONLY : intra_image_comm
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nkstot, xk, nelec, nelup, neldw
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : nbnd
  USE wannier90_readback,   ONLY : w90rb_read_umat
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), ALLOCATABLE :: unimatrx_occ(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: unimatrx_emp(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: unimatrx_emp_opt(:,:,:)
  INTEGER :: i, j, nkstot_eff, nbnd_pw, num_ks_file, nrows, ncols
  INTEGER :: ieff, jeff
  !
  !! Here and in W90 We treat always one spin component at the time.
  !! While in PWscf nkstot contain k point for spin up (the first nkstot/nspin 
  !! and spin down (the last nkstot.nspin). See also rotate_ks
  !
  IF (nspin == 4) THEN
      nkstot_eff = nkstot
      nbnd_pw = nint(nelec)
  ELSE IF (nspin ==2)  THEN
      nkstot_eff = nkstot/nspin
      nbnd_pw = nint(nelup)
      IF (spin_component == 2) nbnd_pw = nint(neldw)
  ELSE
      nkstot_eff = nkstot/nspin
      nbnd_pw= nint(nelec) / 2
  ENDIF 
  !
  !#### In the case occupied and empty manifolds are treated separately
  !#### we assume an insulating systems and NO disentanglement for occupied manifold
  !  
  CALL w90rb_read_umat (TRIM(seedname)//'_u.mat', xk(:,1:nkstot_eff), intra_image_comm, &
                        unimatrx_occ, nrows, ncols)
  !
  IF (ncols /= num_wann_occ .AND. num_wann_occ /= 0) &
           CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U matrix', 1)
  num_wann_occ = ncols
  !! Store the number of occupied wannier in a global variable
  !
  IF (num_wann_occ /= nbnd_pw) &
       CALL errore ('read_wannier', 'Mismatch between  num occ bands and num wann', 1)
  !
  num_wann = num_wann_occ
  if (.NOT. have_empty .AND. kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  !
  IF (.NOT. have_empty) THEN
    !
    ! Store the unitary matrix in a global vabiable
    ALLOCATE (unimatrx (num_wann_occ, num_wann_occ, nkstot_eff)) 
    unimatrx = unimatrx_occ
    !
    ! The optimal subspace Matrix
    ALLOCATE (unimatrx_opt (num_wann_occ, num_wann_occ, nkstot_eff)) 
    unimatrx_opt=CMPLX(0.D0,0.D0, kind=DP)
    !
    ! The optimal matrix is simply the identity
    DO i = 1, num_wann_occ; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
    DEALLOCATE (unimatrx_occ)
    !
    RETURN   ! Nothing more to do
    !
  ENDIF
  !
  !####################
  !   EMPTY STATES 
  !###################
  !
  ! ... read unitary matrix empty states ...
  !
  CALL w90rb_read_umat (TRIM(seedname)//'_emp_u.mat', xk(:,1:nkstot_eff), intra_image_comm, &
                        unimatrx_emp, nrows, ncols)
  !
  IF (ncols /= num_wann_emp .AND. num_wann_emp /= 0) &
        CALL errore ('read_wannier', 'Mismatch in num_wann from input vs from U Empty matrix', 1)
  num_wann_emp = ncols
  !! Store the numebr of empty wannier in a global variable
  !
  ! ... The optimal subspace matrix for the empty states (if any)
  !
  IF (has_disentangle) THEN
     !
     CALL w90rb_read_umat (TRIM(seedname)//'_emp_u_dis.mat', xk(:,1:nkstot_eff), intra_image_comm, &
                           unimatrx_emp_opt, nrows, ncols)
     num_ks_file = nrows
     !
     IF (ncols /= num_wann_emp) &
           CALL errore ('read_wannier', 'Mismatch in num_wann from U Empty vs U Optimal matrix', 1)
     IF (num_ks_file /= ( nbnd - num_wann_occ ) ) &
          CALL errore ('read_wannier', 'Mismatch between num empty bands and num empty wann', 1)
     !
  ELSE
     !
     num_ks_file = num_wann_emp
     !
  ENDIF
  !
  if ( kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: Optimal Matrix READ")') 
  !
  num_wann = num_wann_occ + num_wann_emp
  if (kcw_iverbosity .gt. 1) WRITE(stdout,'(/,5X, "INFO: total number of Wannier functions", i5)') num_wann
  !
  ! Store the result in a unique matrix
  ALLOCATE (unimatrx (num_wann, num_wann, nkstot_eff))
  ALLOCATE (unimatrx_opt (nbnd, num_wann, nkstot_eff))
  !
  unimatrx=CMPLX(0.D0,0.D0, kind=DP)
  unimatrx_opt=CMPLX(0.D0,0.D0, kind=DP)
  !
  ! Build the unique disentangle matrix (identity for the occupied state and 
  ! unimatrx_emp_opt for empty states
  !
  ! 1)occ
  DO i = 1, num_wann_occ; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
  !
  ! 2)emp
  IF (.not. has_disentangle) THEN  
    !
    ! ...In this case the optimal matrix is the ideintity (nbnd=num_wann)
    DO i = 1, num_wann; unimatrx_opt(i,i,:) = CMPLX(1.D0, 0.D0, kind=DP); ENDDO
    !
  ELSE
    !
    DO i = 1,num_ks_file
      !
      ieff=i+num_wann_occ
      DO j = 1, num_wann_emp
        jeff=j+num_wann_occ
        unimatrx_opt(ieff,jeff,:) = unimatrx_emp_opt(i,j,:)
      ENDDO
      !
    ENDDO
    !
  ENDIF 
  !
  ! Build the unique unitary matrix
  !
  ! 1)occ
  unimatrx(1:num_wann_occ, 1:num_wann_occ, :) = unimatrx_occ(1:num_wann_occ, 1:num_wann_occ, :)
  !
  ! 2)emp
  DO i = 1, num_wann_emp
    !
    ieff=i+num_wann_occ
    DO j = 1, num_wann_emp
      ! 
      jeff=j+num_wann_occ
      unimatrx(ieff,jeff,:) = unimatrx_emp(i,j,:)
      !
    ENDDO
    !
  ENDDO
  !
  DEALLOCATE (unimatrx_occ)
  DEALLOCATE (unimatrx_emp)
  IF (has_disentangle) DEALLOCATE (unimatrx_emp_opt )
  !
  RETURN
  !
END subroutine read_wannier_two_manifold
