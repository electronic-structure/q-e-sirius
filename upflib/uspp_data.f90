!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE uspp_data
  !
  !! These parameters are needed with the US pseudopotentials.
  !
  USE upf_kinds,      ONLY : DP
  !
  SAVE
  PRIVATE
  !
  PUBLIC :: nqxq, nqx, dq
  PUBLIC :: qrad,   tab,   tab_at
  PUBLIC :: tab_d
  !
  PUBLIC :: allocate_uspp_data
  PUBLIC :: deallocate_uspp_data
  PUBLIC :: scale_uspp_data
#if defined(__SIRIUS)
  PUBLIC :: beta_ri_tab, aug_ri_tab, wfc_ri_tab
#endif
  !
  INTEGER :: nqxq
  !! size of interpolation table
  INTEGER :: nqx
  !! number of interpolation points
  REAL(DP), PARAMETER:: dq = 0.01D0
  !! space between points in the pseudopotential tab.
  REAL(DP), ALLOCATABLE :: qrad(:,:,:,:)
  !! interpolation table for radial FT of Q functions
  REAL(DP), ALLOCATABLE :: tab(:,:,:)
  !! interpolation table for PP projectorss
  REAL(DP), ALLOCATABLE :: tab_at(:,:,:)
  !! interpolation table for atomic wfc
  !
#if defined(__SIRIUS)
  REAL(DP), ALLOCATABLE :: beta_ri_tab(:,:,:)
  !! radial integrals of beta projectors without unit-cell volume (Omega) factor
  !
  REAL(DP), ALLOCATABLE :: aug_ri_tab(:,:,:,:)
  !! radial integrals of augmentation charge without unit-cell volume (Omega) factor
  !
  REAL(DP), ALLOCATABLE :: wfc_ri_tab(:,:,:)
  !! radial integrals of atomic wave-functions
  !
#endif
  !! GPUs variables - only those tables that is useful to have on GPUss
  !
  REAL(DP), ALLOCATABLE :: tab_d(:,:,:)
  !
#if defined(__CUDA)
  attributes (DEVICE) :: tab_d
#endif
  !
contains
  !
  subroutine allocate_uspp_data(use_gpu,nqxq_,nqx_,nbetam,nwfcm,lmaxq,nsp)
     implicit none
     logical, intent(in) :: use_gpu
     integer, intent(in) :: nqxq_,nqx_,nbetam,nwfcm,lmaxq,nsp
     !
     if (nqxq_/=nqxq) call upf_error("allocate_uspp_data","invalid nqxq_",1)
     if (nqx_/=nqx)   call upf_error("allocate_uspp_data","invalid nqx_",1)
     !
     if (lmaxq>0) allocate(qrad(nqxq_,nbetam*(nbetam+1)/2, lmaxq, nsp))
     allocate(tab(nqx_,nbetam,nsp))
     allocate(tab_at(nqx_,nwfcm,nsp))
     !$acc enter data create(qrad,tab_at)
     !
     IF (use_gpu) then
        ! allocations with zero size protected
        ! since problematic with CUDAfor
        if (nbetam>0)  allocate(tab_d(nqx_,nbetam,nsp))
     endif
     !
  end subroutine allocate_uspp_data
  !
  subroutine deallocate_uspp_data()
     implicit none
     if( allocated( tab ) )       deallocate( tab )
     !$acc exit data delete(qrad,tab_at)
     if( allocated( qrad ) )      deallocate( qrad )
     if( allocated( tab_at ) )    deallocate( tab_at )
     !
     if( allocated( tab_d ) )     deallocate( tab_d )
  end subroutine
  !
  subroutine scale_uspp_data( vol_ratio_m1 )
     ! vol_ratio_m1 = omega_old / omega
     implicit none
     real(DP), intent(in) :: vol_ratio_m1
     !
     tab(:,:,:)    = tab(:,:,:) * SQRT(vol_ratio_m1)
     qrad(:,:,:,:) = qrad(:,:,:,:) * vol_ratio_m1
     tab_at(:,:,:) = tab_at(:,:,:) * SQRT(vol_ratio_m1)
#if defined __CUDA
!$acc update device (qrad, tab_at)
     ! CUDA Fortran safeguard
     if(size(tab) > 0) tab_d = tab
#endif
  end subroutine scale_uspp_data
  !
END MODULE uspp_data

