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
  PUBLIC :: qrad_d, tab_d, tab_at_d
  PUBLIC :: tab_rho, tab_rhc
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
  REAL(DP), ALLOCATABLE :: tab_rho(:,:)
  !! interpolation table for atomic charge density
  REAL(DP), ALLOCATABLE :: tab_rhc(:,:)
  !! interpolation table for atomic pseudo-core charge density
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
  REAL(DP), ALLOCATABLE :: qrad_d(:,:,:,:)
  REAL(DP), ALLOCATABLE :: tab_d(:,:,:)
  REAL(DP), ALLOCATABLE :: tab_at_d(:,:,:)
  !
#if defined(__CUDA)
  attributes (DEVICE) :: qrad_d, tab_d, tab_at_d
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
     allocate(tab_rho(nqxq_,nsp))
     allocate(tab_rhc(nqxq_,nsp))
     !
     IF (use_gpu) then
        ! allocations with zero size protected
        ! since problematic with CUDAfor
        if (lmaxq>0.and.nbetam>0)  &
                       allocate(qrad_d(nqxq_,nbetam*(nbetam+1)/2, lmaxq, nsp))
        if (nbetam>0)  allocate(tab_d(nqx_,nbetam,nsp))
        if (nwfcm>0)   allocate(tab_at_d(nqx_,nwfcm,nsp))
     endif
     !
  end subroutine allocate_uspp_data
  !
  subroutine deallocate_uspp_data()
     implicit none
     if( allocated( qrad ) )      deallocate( qrad )
     if( allocated( tab ) )       deallocate( tab )
     if( allocated( tab_at ) )    deallocate( tab_at )
!$acc exit data delete(tab_rho, tab_rhc)
     if( allocated( tab_rho) )    deallocate( tab_rho)
     if( allocated( tab_rhc) )    deallocate( tab_rhc)
     !
     if( allocated( qrad_d ) )    deallocate( qrad_d )
     if( allocated( tab_d ) )     deallocate( tab_d )
     if( allocated( tab_at_d ) )  deallocate( tab_at_d )
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
     tab_rho(:,:)  = tab_rho(:,:) * vol_ratio_m1
     tab_rhc(:,:)  = tab_rhc(:,:) * vol_ratio_m1
#if defined __CUDA
!$acc enter data copyin (tab_rho, tab_rhc)
     ! CUDA Fortran safeguard
     if(size(tab) > 0) tab_d = tab
     if(size(qrad) > 0) qrad_d = qrad
     if(size(tab_at) > 0) tab_at_d = tab_at
#endif
  end subroutine scale_uspp_data
  !
END MODULE uspp_data

