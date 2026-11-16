! Copyright(C) 2005-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE exx_std
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only, tqr
  USE io_global,            ONLY : stdout
  USE exx_base,             ONLY : dfftt, exxbuff, locbuff, exxmat, locmat, evc0, &
                                   xkq_collect, index_xk, npwt, exxalfa, eps_occ, &
                                   ibnd_start, ibnd_end, ibnd_buff_start, ibnd_buff_end
  USE klist,                ONLY : nks, igk_k
  USE noncollin_module,     ONLY : noncolin, npol
  !
 CONTAINS
#define _CX(A)  CMPLX(A,0._dp,kind=DP)
#define _CY(A)  CMPLX(0._dp,-A,kind=DP)
  !
  !------------------------------------------------------------------------
  SUBROUTINE exxinit_std( DoLoc )
    !------------------------------------------------------------------------
    !! This subroutine is run before the first H_psi() of each iteration. 
    !! It saves the wavefunctions for the right density matrix, in real space.
    !
    USE wavefunctions,        ONLY : evc, psic
    USE io_files,             ONLY : nwordwfc, iunwfc
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx
    USE klist,                ONLY : ngk, nkstot
    USE mp_exx,               ONLY : me_egrp, inter_egrp_comm
    USE mp,                   ONLY : mp_bcast
    USE scatter_mod,          ONLY : gather_grid, scatter_grid
    USE fft_interfaces,       ONLY : invfft
    USE uspp,                 ONLY : okvan
    USE paw_variables,        ONLY : okpaw
    USE mp_orthopools,        ONLY : intra_orthopool_comm
    USE exx_base,             ONLY : nkqs, index_sym, rir, working_pool, &
                                     x_nbnd_occ, nbndproj, local_thr, d_spin
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: DoLoc
    !! TRUE:  Real Array locbuff(ir, nbnd, nkqs);  
    !! FALSE: Complex Array exxbuff(ir, nbnd/2, nkqs).
    !
    ! ... local variables
    !
    INTEGER :: ik, ibnd, i, isym, ikq, ig
    INTEGER :: ibnd_loop_start
    INTEGER :: ipol, jpol
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: temppsic
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) temppsic
#endif
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: psic_exx(:)
    INTEGER :: nxxs, nrxxs
#if defined(__MPI)
    COMPLEX(DP),ALLOCATABLE  :: temppsic_all(:), psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    INTEGER :: npw, current_ik
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    INTEGER :: h_ibnd
    INTEGER :: evc_offset
    !
    IF ( Doloc ) THEN
        WRITE(stdout,'(/,5X,"Using localization algorithm with threshold: ",&
                & D10.2)') local_thr
        ! IF (.NOT.gamma_only) CALL errore('exxinit','SCDM with K-points NYI',1)
        IF (okvan .OR. okpaw) CALL errore( 'exxinit','SCDM with USPP/PAW not &
                                           &implemented', 1 )
    ENDIF 
    !
    ! Note that nxxs is not the same as nrxxs in parallel case
    nxxs = dfftt%nr1x * dfftt%nr2x * dfftt%nr3x
    nrxxs = dfftt%nnr
    !
#if defined(__MPI)
    IF (noncolin) THEN
       ALLOCATE( psic_all_nc(nxxs,npol), temppsic_all_nc(nxxs,npol) )
    ELSEIF ( .NOT. gamma_only ) THEN
       ALLOCATE( psic_all(nxxs), temppsic_all(nxxs) )
    ENDIF
#endif
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs, npol), psic_nc(nrxxs, npol) )
    ELSEIF ( .NOT. gamma_only ) THEN
       ALLOCATE( temppsic(nrxxs) )
       !$acc enter data create(temppsic(1:nrxxs))
    ENDIF
    !
    CALL divide( inter_egrp_comm, x_nbnd_occ, ibnd_start, ibnd_end )
    !
    IF ( gamma_only ) THEN
        ibnd_buff_start = ibnd_start/2
        IF(mod(ibnd_start,2)==1) ibnd_buff_start = ibnd_buff_start +1
        !
        ibnd_buff_end = ibnd_end/2
        IF(mod(ibnd_end,2)==1) ibnd_buff_end = ibnd_buff_end +1
    ELSE
        ibnd_buff_start = ibnd_start
        ibnd_buff_end   = ibnd_end
    ENDIF
    !
    IF (DoLoc) THEN
      !
      IF (gamma_only) THEN
        IF (.NOT. ALLOCATED(locbuff)) ALLOCATE( locbuff(nrxxs*npol,nbnd,nks) )
        IF (.NOT. ALLOCATED(locmat))  ALLOCATE( locmat(nbnd,nbnd,nks) )
        locbuff = 0.0d0
        locmat = 0.0d0
      ELSE 
        IF (.NOT. ALLOCATED(exxbuff)) ALLOCATE( exxbuff(nrxxs*npol,nbnd,nkqs) )
        IF (.NOT. ALLOCATED(exxmat) ) ALLOCATE( exxmat(nbnd,nkqs,nbnd,nks) )
        exxbuff = (0.0d0, 0.0d0)
        exxmat = 0.0d0
      ENDIF
      !
      IF (.NOT. ALLOCATED(evc0)) then 
        ALLOCATE( evc0(npwx*npol,nbndproj,nks) )
        evc0 = (0.0d0,0.0d0)
      ENDIF
      !
    ELSE
      !
      IF (.NOT. ALLOCATED(exxbuff)) THEN
        IF (gamma_only) THEN
           ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nks))
        ELSE
           ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nkqs))
        END IF
        !$acc enter data create(exxbuff)
      END IF
      !
    ENDIF
    !
    !assign buffer
    IF(DoLoc) THEN
      locbuff = 0.0_DP
    ELSE
      exxbuff = (0.0_DP,0.0_DP)
    ENDIF
    !
    ! ... This is parallelized over pools. Each pool computes only its k-points
    !
    KPOINTS_LOOP : &
    DO ik = 1, nks
       !
       IF ( nks > 1 ) THEN
         CALL get_buffer( evc, nwordwfc, iunwfc, ik )
         !$acc update device(evc)
       END IF
       !
       ! ik         = index of k-point in this pool
       ! current_ik = index of k-point over all pools
       !
       current_ik = global_kpoint_index( nkstot, ik )
       !
       IF_GAMMA_ONLY : &
       IF (gamma_only) THEN
          !
          h_ibnd = ibnd_start/2
          !
          IF(mod(ibnd_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start-1
          ELSE
             ibnd_loop_start=ibnd_start
          ENDIF
          ! 
          DO ibnd = ibnd_loop_start, ibnd_end, 2
             h_ibnd = h_ibnd + 1
             !
             psic(:) = ( 0._DP, 0._DP )
             !
             IF ( ibnd < ibnd_end ) THEN
                DO ig = 1, npwt
                   psic(dfftt%nl(ig))  = evc(ig,ibnd) + ( 0._dp, 1._dp ) * evc(ig,ibnd+1)
                   psic(dfftt%nlm(ig)) = conjg( evc(ig,ibnd) ) + ( 0._dp, 1._dp ) * conjg( evc(ig,ibnd+1) )
                ENDDO
             ELSE
                DO ig = 1, npwt
                   psic(dfftt%nl(ig))  = evc(ig,ibnd)
                   psic(dfftt%nlm(ig)) = CONJG( evc(ig,ibnd) ) 
                ENDDO
             ENDIF
             !
             CALL invfft( 'Wave', psic, dfftt )
             !
             IF (DoLoc) THEN
               locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+1,ik) = DBLE(  psic_exx(1:nrxxs) )
               IF (ibnd-ibnd_loop_start+evc_offset+2 <= nbnd) &
                  locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+2,ik) = AIMAG( psic_exx(1:nrxxs) )
             ELSE
               exxbuff(1:nrxxs,h_ibnd,ik)=psic(1:nrxxs) 
             ENDIF
             !
          ENDDO
          !
       ELSE IF_GAMMA_ONLY
          !
          npw = ngk (ik)
          IBND_LOOP_K : &
          DO ibnd = ibnd_start, ibnd_end 
             !
             IF (noncolin) THEN
                temppsic_nc(:,:) = ( 0._DP, 0._DP )
                temppsic_nc(dfftt%nl(igk_k(1:npw,ik)),1) = evc(1:npw,ibnd)
                CALL invfft( 'Wave', temppsic_nc(:,1), dfftt )
                temppsic_nc(dfftt%nl(igk_k(1:npw,ik)),2) = evc(npwx+1:npwx+npw,ibnd)
                CALL invfft( 'Wave', temppsic_nc(:,2), dfftt )
             ELSE
                !$acc kernels
                temppsic(:) = ( 0._DP, 0._DP )
                !$acc end kernels
                !$acc parallel loop 
                do i = 1, npw
                  temppsic(dfftt%nl(igk_k(i,ik))) = evc(i,ibnd)
                end do 
                !$acc host_data use_device(temppsic)
                CALL invfft( 'Wave', temppsic, dfftt )
                !$acc end host_data
             ENDIF
             !
             DO ikq = 1, nkqs
                !
                IF (index_xk(ikq) /= current_ik) CYCLE
                isym = ABS(index_sym(ikq) )
                !
                IF (noncolin) THEN ! noncolinear
#if defined(__MPI)
                   DO ipol = 1, npol
                      CALL gather_grid( dfftt, temppsic_nc(:,ipol), temppsic_all_nc(:,ipol) )
                   ENDDO
                   !
                   IF ( me_egrp == 0 ) THEN
                      psic_all_nc(:,:) = (0.0_DP, 0.0_DP)
                      DO ipol = 1, npol
                         DO jpol = 1, npol
                            psic_all_nc(:,ipol) = psic_all_nc(:,ipol) + &
                                          CONJG(d_spin(jpol,ipol,isym)) * &
                                          temppsic_all_nc(rir(:,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDIF
                   !
                   DO ipol = 1, npol
                      CALL scatter_grid( dfftt, psic_all_nc(:,ipol), psic_nc(:,ipol) )
                   ENDDO
#else
                   DO ipol = 1, npol
                      psic_nc(:,ipol) = (0._DP,0._DP)
                      DO jpol = 1, npol
                         psic_nc(:,ipol) = psic_nc(:,ipol) + CONJG(d_spin(jpol,ipol,isym))* &
                                            temppsic_nc(rir(:,isym),jpol)
                      ENDDO
                   ENDDO
#endif
                   !
!civn: not sure about these two cases are really different in this new(old) implementation
                   IF (index_sym(ikq) > 0 ) THEN
                      ! sym. op. without time reversal: normal case
                      exxbuff(1:nrxxs,ibnd,ikq) = psic_nc(:,1)
                      exxbuff(nrxxs+1:nrxxs+nrxxs,ibnd,ikq) = psic_nc(:,2)
                   ELSE
                      ! sym. op. with time reversal: spin 1->2*, 2->-1*
                      exxbuff(1:nrxxs,ibnd,ikq) = CONJG(psic_nc(:,2))
                      exxbuff(nrxxs+1:nrxxs+nrxxs,ibnd,ikq) = -CONJG(psic_nc(:,1))
                   ENDIF
                ELSE ! noncolinear
                   !$acc update host(temppsic)
#if defined(__MPI)
                   CALL gather_grid( dfftt, temppsic, temppsic_all )
                   IF ( me_egrp == 0 ) psic_all(:) = temppsic_all(rir(:,isym))
                   CALL scatter_grid( dfftt, psic_all, psic )
#else
                   psic(:) = temppsic(rir(:,isym))
#endif
                   IF (index_sym(ikq) < 0 ) psic(1:nrxxs) = CONJG(psic(1:nrxxs))
                   exxbuff(1:nrxxs,ibnd,ikq) = psic(1:nrxxs)
                   !
                ENDIF ! noncolinear
                !
             ENDDO
             !
          ENDDO&
          IBND_LOOP_K
          !
       ENDIF&
       IF_GAMMA_ONLY
    ENDDO&
    KPOINTS_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc, psic_nc )
#if defined(__MPI)
       DEALLOCATE( temppsic_all_nc, psic_all_nc )
#endif
    ELSE IF ( .NOT. gamma_only ) THEN
       !$acc exit data delete(temppsic)
       DEALLOCATE( temppsic )
#if defined(__MPI)
       DEALLOCATE( temppsic_all, psic_all )
#endif
    ENDIF
    !
    ! Each wavefunction in exxbuff is computed by a single pool, collect among 
    ! pools in a smart way (i.e. without doing all-to-all sum and bcast)
    ! See also the initialization of working_pool in exx_mp_init
    ! Note that in Gamma-only LSDA can be parallelized over two pools, and there
    ! is no need to communicate anything: each pools deals with its own spin
    !
    IF ( .NOT. gamma_only ) THEN
       DO ikq = 1, nkqs
         CALL mp_bcast( exxbuff(:,:,ikq), working_pool(ikq), intra_orthopool_comm ) 
       ENDDO
    ENDIF
    !
    !$acc update device(exxbuff)
    !
  END SUBROUTINE exxinit_std
  !
#if defined(__CUDA)
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_std_gamma_outer(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... Gamma-specific version of vexx with batched outer loop over bands
    !     (this is a GPU-only version, howmany on CPU not implemented)
    !
    USE cell_base,      ONLY : omega
    USE wvfct,          ONLY : current_k
    USE klist,          ONLY : nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, nbgrp
    USE mp,             ONLY : mp_sum, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : npwt, gt, nqs, x_occupation, g2_convolution, batchsize
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: RESULT(:), temppsic(:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), deexx(:,:)
    REAL(DP),   ALLOCATABLE :: fac(:), facmask(:)
    INTEGER          :: ibnd
    INTEGER          :: ir, ig, idx, idxm
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP), PARAMETER :: xkg(3) = (/0.0_dp, 0.0_dp, 0.0_dp/)
    ! coordinates of Gamma point
    REAL(DP) :: x1, x2
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL  :: global_kpoint_index
    INTEGER :: ngmt
    !
    ! local variables to control FFT band batching
    !
    INTEGER :: nb, nb_last
    ! number of batched bands per fft 
    ! in principle nb=batchsize, but it may vary for the last batch (nb_last).
    ! using a fixed batchsize prevents allocating and deallocating batched variables inside the loop
    INTEGER :: nblocks, iblock
    ! loop over batches
    INTEGER :: i_global
    ! global index of bands ranging (1...m)
    INTEGER :: ib 
    ! local index of bands ranging (1...nb)
    INTEGER :: batchsize0, subbatchsize0
    ! backup values in case FFT internal batching is changed 
    LOGICAL, SAVE :: first_call = .true.
    !
    if ( mod(batchsize, 2) /= 0 ) Call errore('vexx_std_gamma_outer', 'odd batchsize values not allowed', 1)
    !
    ngmt = dfftt%ngm
    ALLOCATE( fac(ngmt) )
    nrxxs= dfftt%nnr
    !
    ! update internal FFT batching parameters 
    !
    batchsize0 = 0
    subbatchsize0 = 0
    if( batchsize > dfftt%batchsize ) then
      batchsize0 = dfftt%batchsize
      subbatchsize0 = dfftt%subbatchsize
      ! set new batchsize and subbatchsize
      if ( mod( batchsize, 4 ) == 0 ) then
        dfftt%subbatchsize = batchsize / 4
        dfftt%batchsize = batchsize
      else
        dfftt%subbatchsize = batchsize / 4 + 1
        dfftt%batchsize = dfftt%subbatchsize * 4
      endif
    endif
    IF( first_call ) THEN
      write(stdout, '(5X,3(A,I4))') &
        "EXX batchsize: ", batchsize, " FFT batchsize: ", dfftt%batchsize, " subbatchsize: ", dfftt%subbatchsize
      first_call = .false.
    END IF
    !
    nblocks = m / batchsize
    nb_last = mod(m, batchsize)             ! number of bands in the last block
    if( nb_last > 0 ) nblocks = nblocks + 1 ! add one more block for the last nb_last bands 
    !
    ALLOCATE( temppsic( batchsize/2 * nrxxs) ) ! Gamma trick here
    ALLOCATE( rhoc(batchsize*nrxxs), RESULT(batchsize*nrxxs), facmask(nrxxs) )
    !
    IF(okvan) ALLOCATE(deexx(nkb, batchsize))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=0.0_DP
       psi =0.0_DP
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    ENDIF
    !$acc update device(psi)
    !
    ! Here the loops start
    !
    !$acc data create(rhoc, fac, temppsic, result, facmask) copy(hpsi) present(exxbuff) 
    !
    ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
    CALL g2_convolution(ngmt, gt, xkg, xkg, fac)
    !$acc update device(fac)
    !
    ! facmask = fac inside the ngmt sphere, 0 otherwise 
    ! (can be used in loops over nrxxs)
    !
    !$acc kernels
    facmask = 0.0_dp
    !$acc end kernels 
    !
    !$acc parallel loop private(idx, idxm)
    do ig = 1, ngmt
      idx  = dfftt%nl(ig)
      idxm = dfftt%nlm(ig)
      facmask( idx ) = fac(ig) 
      facmask( idxm) = fac(ig) 
    end do
    !
    IF ( okvan .and..not.tqr ) CALL qvan_init (ngmt, xkg, xkg)
    !
    LOOP_ON_I_BLOCKS : &
    DO iblock = 1, nblocks
       !
       ! bring a block (temppsic) of psi functions to real space
       !
       ! set the actual number of bands in this block
       if(nb_last > 0 .and. iblock==nblocks) then
         nb = nb_last 
       else
         nb = batchsize
       endif
       ! set the global band starting index for this block
       i_global = (iblock-1) * batchsize 
       !
       !$acc kernels 
       temppsic = (0.0_dp, 0.0_dp)
       !$acc end kernels
       !
       !$acc parallel loop collapse(2) private(idx, idxm)
       DO ib = 1, nb-1, 2
         DO ig = 1, npwt
            idx  = ((ib-1)/2) * nrxxs + dfftt%nl(ig) 
            idxm = ((ib-1)/2) * nrxxs + dfftt%nlm(ig) 
            temppsic ( idx  ) =       psi(ig, i_global + ib) + (0._DP,1._DP) * psi(ig, i_global + ib+1)
            temppsic ( idxm ) = conjg(psi(ig, i_global + ib) - (0._DP,1._DP) * psi(ig, i_global + ib+1))
         ENDDO
       ENDDO 
       ! when m is odd, in the last block nb_last is also odd and the previous loop
       ! cannot increment the last iteration otherwise psi(ig,i_global+ib+1) goes out of bounds 
       ! so this is the last iteration: 
       if( mod(nb_last,2)==1 .and. iblock==nblocks ) then 
         !$acc parallel loop private(idx, idxm)
         DO ig = 1, npwt
            idx  = (nb/2) * nrxxs + dfftt%nl(ig) 
            idxm = (nb/2) * nrxxs + dfftt%nlm(ig) 
            temppsic ( idx  ) =       psi(ig, i_global + nb) 
            temppsic ( idxm ) = conjg(psi(ig, i_global + nb) )
         ENDDO
       end if 
       !
       !$acc host_data use_device(temppsic)
       CALL invfft ('Wave', temppsic, dfftt, howmany=batchsize/2)
       !$acc end host_data
       !
       ! now in temppsic there are batchsize wavefunctons, 
       ! half in the real part, half in the imaginary part.
       !
       IF(okvan) deexx = 0.0_DP
       !
       !$acc kernels
       RESULT = 0.0_DP
       !$acc end kernels
       !
       h_ibnd = ibnd_start/2
       IF(mod(ibnd_start,2)==0) THEN
          h_ibnd=h_ibnd-1
          ibnd_loop_start=ibnd_start-1
       ELSE
          ibnd_loop_start=ibnd_start
       ENDIF
       !
       IBND_LOOP_GAM : &
       DO ibnd=ibnd_loop_start,ibnd_end, 2 !for each band of psi
          !
          h_ibnd = h_ibnd + 1
          IF( ibnd < ibnd_start ) THEN
             x1 = 0.0_DP
          ELSE
             x1 = x_occupation(ibnd,  current_ik)
          ENDIF
          IF( ibnd == ibnd_end) THEN
             x2 = 0.0_DP
          ELSE
             x2 = x_occupation(ibnd+1,  current_ik)
          ENDIF
          IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
          !
          ! calculate rho in real space. Gamma tricks are used.
          ! temppsic is real; tempphic contains one band in the real part,
          ! another one in the imaginary part; the same applies to rhoc
          !
          !$acc parallel loop collapse(2) private(idx) 
          do ib = 1, batchsize, 2
             DO ir = 1, nrxxs
                idx = ib/2 *nrxxs + ir 
                rhoc((ib-1)*nrxxs + ir ) = exxbuff(ir,h_ibnd,current_k) *  dble(temppsic(idx)) / omega
                rhoc((ib)*nrxxs + ir )   = exxbuff(ir,h_ibnd,current_k) * aimag(temppsic(idx)) / omega
             ENDDO
          end do 
          !
          ! bring rho to G-space
          !
          !   >>>> add augmentation in REAL SPACE here
          IF(okvan .and. tqr) THEN
             !$acc update host(rhoc)
             do ib = 1, nb
               IF(ibnd>=ibnd_start) &
                 CALL addusxx_r(rhoc((ib-1)*nrxxs+1:ib*nrxxs),_CX(becxx(current_ik)%r(:,ibnd)), &
                                                                            _CX(becpsi%r(:,i_global+ib)))
               IF(ibnd<ibnd_end) &
                 CALL addusxx_r(rhoc((ib-1)*nrxxs+1:ib*nrxxs),_CY(becxx(current_ik)%r(:,ibnd+1)), &
                                                                            _CX(becpsi%r(:,i_global+ib)))
             end do
             !$acc update device(rhoc)
          ENDIF
          !
          !$acc host_data use_device(rhoc)
          CALL fwfft ('Rho', rhoc, dfftt, howmany=batchsize)
          !$acc end host_data
          !
          !   >>>> add augmentation in G SPACE here
          IF(okvan .and. .not. tqr) THEN
             !$acc update host(rhoc)
             ! contribution from one band added to real (in real space) part of rhoc
             do ib = 1, nb
               IF(ibnd>=ibnd_start) &
                  CALL addusxx_g(dfftt, rhoc((ib-1)*nrxxs+1), xkg,  xkg, 'r', &
                  becphi_r=becxx(current_ik)%r(:,ibnd), becpsi_r=becpsi%r(:,i_global+ib) )
               ! contribution from following band added to imaginary (in real space) part of rhoc
               IF(ibnd<ibnd_end) &
                  CALL addusxx_g(dfftt, rhoc((ib-1)*nrxxs+1), xkg,  xkg, 'i', &
                  becphi_r=becxx(current_ik)%r(:,ibnd+1), becpsi_r=becpsi%r(:,i_global+ib) )
             end do 
             !$acc update device(rhoc)
          ENDIF
          !   >>>> charge density done
          !
          !$acc parallel loop collapse(2) private(idx)
          do ib = 1, batchsize
            do ir = 1, nrxxs
              idx = (ib-1)*nrxxs + ir
              rhoc( idx ) = rhoc( idx ) * facmask(ir) 
            end do 
          end do  
          !
          !   >>>>  compute <psi|H_fock G SPACE here
          IF(okvan .and. .not. tqr) THEN
             !$acc update host(rhoc)
             do ib = 1, nb
               IF(ibnd>=ibnd_start) &
               CALL newdxx_g(dfftt, rhoc((ib-1)*nrxxs+1), xkg, xkg, 'r', deexx(:,ib), &
                             becphi_r=x1*becxx(current_ik)%r(:,ibnd))
               IF(ibnd<ibnd_end) &
               CALL newdxx_g(dfftt, rhoc((ib-1)*nrxxs+1), xkg, xkg, 'i', deexx(:,ib), &
                             becphi_r=x2*becxx(current_ik)%r(:,ibnd+1))
             end do
          ENDIF
          !
          !brings back v in real space
          !$acc host_data use_device(rhoc)
          CALL invfft ('Rho', rhoc, dfftt, howmany=batchsize)
          !$acc end host_data
          !
          !   >>>>  compute <psi|H_fock REAL SPACE here
          IF(okvan .and. tqr) THEN
             !$acc update host(rhoc)
             do ib = 1, nb 
               IF(ibnd>=ibnd_start) &
                 CALL newdxx_r(dfftt, rhoc((ib-1)*nrxxs+1:ib*nrxxs), _CX(x1*becxx(current_ik)%r(:,ibnd)), &
                                                                                                     deexx(:,ib))
               IF(ibnd<ibnd_end) &
                 CALL newdxx_r(dfftt, rhoc((ib-1)*nrxxs+1:ib*nrxxs), _CY(x2*becxx(current_ik)%r(:,ibnd+1)), &
                                                                                                     deexx(:,ib))
             end do 
          ENDIF
          !
          IF(okpaw) THEN
             do ib = 1, nb
               IF(ibnd>=ibnd_start) &
               CALL PAW_newdxx(x1/nqs, _CX(becxx(current_ik)%r(:,ibnd)),&
                                       _CX(becpsi%r(:,i_global+ib)), deexx(:,ib))
               IF(ibnd<ibnd_end) &
               CALL PAW_newdxx(x2/nqs, _CX(becxx(current_ik)%r(:,ibnd+1)), &
                                       _CX(becpsi%r(:,i_global+ib)), deexx(:,ib))
             end do 
          ENDIF
          !
          ! accumulates over bands and k points
          !
          !$acc parallel loop collapse(2) private(idx)
          do ib = 1, batchsize
            DO ir = 1, nrxxs
               idx  = (ib-1)*nrxxs + ir 
               RESULT(idx) = RESULT(idx)+x1* dble(rhoc(idx))* dble(exxbuff(ir,h_ibnd,current_k))&
                                        +x2*aimag(rhoc(idx))*aimag(exxbuff(ir,h_ibnd,current_k))
            ENDDO
          end do 
          !
       ENDDO &
       IBND_LOOP_GAM
       !
       IF(okvan) THEN
          CALL mp_sum(deexx,intra_bgrp_comm)
          CALL mp_sum(deexx,inter_bgrp_comm)
       ENDIF
       !
       !$acc host_data use_device(result)
       CALL mp_sum( RESULT, inter_bgrp_comm)
       !$acc end host_data
       !
       ! brings back result in G-space
       !
       !$acc host_data use_device(RESULT)
       CALL fwfft( 'Wave' , RESULT, dfftt, howmany=batchsize )
       !$acc end host_data
       !
       !$acc parallel loop collapse(2) private(idx)
       do ib = 1, nb
         DO ig = 1, n
            idx  = (ib-1)*nrxxs + dfftt%nl(ig)
            hpsi(ig,i_global+ib)=hpsi(ig,i_global+ib) - exxalfa*RESULT(idx)
         ENDDO
       end do 
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) THEN
         !$acc update host(hpsi(:,i_global+1:i_global+nb))
         do ib = 1, nb
           CALL add_nlxx_pot (lda, hpsi(:,i_global+ib), xkg, n, &
                        igk_k(1,current_k), deexx(:,ib), eps_occ, exxalfa)
         end do 
         !$acc update device(hpsi(:,i_global+1:i_global+nb))
       END IF
    ENDDO &
    LOOP_ON_I_BLOCKS
    !
    IF ( okvan .and..not.tqr ) CALL qvan_clean ()
    !
    !$acc end data
    !
    ! restore original FFT batching parameters
    if( batchsize0 > 0 ) then
      dfftt%batchsize = batchsize0
      dfftt%subbatchsize = subbatchsize0
    endif
    !
    DEALLOCATE( RESULT, temppsic)
    !
    DEALLOCATE(rhoc, fac, facmask )
    !
    IF(okvan) DEALLOCATE( deexx )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_std_gamma_outer
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_std_k_outer(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... generic, k-point version of vexx
    !     (this is a GPU-only version, howmany on CPU not implemented)
    !
    USE cell_base,      ONLY : omega
    USE wvfct,          ONLY : current_k
    USE klist,          ONLY : xk, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, nbgrp
    USE mp,             ONLY : mp_sum, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : gt, nqs, index_xkq, x_occupation, g2_convolution, batchsize
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
    COMPLEX(DP),ALLOCATABLE :: result(:)
    COMPLEX(DP),ALLOCATABLE :: rhoc(:)
    !
    COMPLEX(DP),ALLOCATABLE :: deexx(:,:)
    INTEGER, ALLOCATABLE :: nl(:), nl_k(:)
    REAL(DP),   ALLOCATABLE :: fac(:), facmask(:)
    INTEGER          :: ibnd, ik, ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: nrxxs
    REAL(DP) :: xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ngmt
    REAL(DP) :: x_occ
    !
    ! local variables to control FFT band batching
    !
    INTEGER :: nb, nb_last
    ! number of batched bands per fft 
    ! in principle nb=batchsize, but it may vary for the last batch (nb_last).
    ! using a fixed batchsize prevents allocating and deallocating batched variables inside the loop
    INTEGER :: nbatch, ibatch 
    ! loop over batches 
    INTEGER :: im_start
    ! global index of bands ranging (1...m)
    INTEGER :: ib
    ! local index of bands ranging (1...nb)
    INTEGER :: batchsize0, subbatchsize0
    !! backup values in case FFT internal batching is changed
    INTEGER :: idx, jdx
    !! index of batched quantities: 
    !!   temppsic, result -> npol*nrxxs*batchsize
    !!   rhoc -> nrxxs*batchsize
    !! Strides: 
    !!      ... (ipol-1) * batchsize * nrxxs
    !!      ... (ib-1)   * nrxxs
    LOGICAL, SAVE :: first_call = .true.
    !
    ngmt = dfftt%ngm 
    ALLOCATE( fac(ngmt) )
    nrxxs= dfftt%nnr
    !
    ! update internal FFT batching parameters 
    !
    batchsize0 = 0 
    subbatchsize0 = 0
    if( batchsize > dfftt%batchsize ) then
      batchsize0 = dfftt%batchsize
      subbatchsize0 = dfftt%subbatchsize
      ! set new batchsize and subbatchsize
      if ( mod( batchsize, 4 ) == 0 ) then
        dfftt%subbatchsize = batchsize / 4
        dfftt%batchsize = batchsize 
      else
        dfftt%subbatchsize = batchsize / 4 + 1 
        dfftt%batchsize = dfftt%subbatchsize * 4 
      endif
    endif
    IF( first_call ) THEN
      write(stdout, '(5X,3(A,I4))') &
        "EXX batchsize: ", batchsize, " FFT batchsize: ", dfftt%batchsize, " subbatchsize: ", dfftt%subbatchsize
      first_call = .false.
    END IF
    !
    nbatch = m / batchsize 
    !
    nb_last = mod( m, batchsize )  ! bands remaining in the last batch
    if(nb_last > 0 ) nbatch = nbatch + 1 
    !
    ALLOCATE( nl(ngmt), nl_k(n) )
    !
    ALLOCATE( temppsic(npol*nrxxs*batchsize) )
    ALLOCATE( result(npol*nrxxs*batchsize) )
    !
    ALLOCATE( facmask(nrxxs) )
    ALLOCATE( rhoc(nrxxs*batchsize) )
    !
    IF(okvan) ALLOCATE(deexx(nkb, batchsize))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=0.0_DP
       psi =0.0_DP
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    ENDIF
    !$acc update device(psi)
    !
    !$acc data create( rhoc, fac, facmask, nl, nl_k, temppsic, result) copy(hpsi) present(exxbuff) 
    !
    !$acc parallel loop
    do ig = 1, ngmt 
      nl(ig) = dfftt%nl(ig)
    end do
    !
    !$acc parallel loop
    do ig = 1, n  
      nl_k(ig) = dfftt%nl(igk_k(ig,current_k))
    end do
    !
    EXTERNAL_LOOP_ON_BATCHES: &
    do ibatch = 1, nbatch
       !
       CALL start_clock( 'vexx_k_setup' )
       !
       ! set the actual number of bands in this batch
       if(nb_last > 0 .and. ibatch==nbatch ) then
         nb = nb_last 
       else
         nb = batchsize
       endif
       ! set the global band starting index for this batch
       im_start = (ibatch-1) * batchsize 
       !
       ! inizializations on the whole allocation dimension (batchsize)
       !$acc parallel loop collapse(3) private(idx)
       do ipol = 1, npol
         do ib = 1, batchsize
           do ir = 1, nrxxs
             idx = (ipol-1) * batchsize * nrxxs + (ib-1) * nrxxs + ir
             temppsic( idx ) = 0._dp
             result( idx ) = 0._dp
           end do   
         end do 
       end do 
       !
       IF(okvan) deexx = 0._DP
       !
       !$acc parallel loop collapse(3) private(idx)
       do ipol = 1, npol
         do ib = 1, nb
           do ig = 1, n
              idx = (ipol-1) * batchsize * nrxxs + (ib-1) * nrxxs + nl_k(ig)
              temppsic( idx ) = psi( (ipol-1)*lda + ig, im_start + ib)
           enddo
         end do 
       end do 
       !
       !$acc host_data use_device(temppsic)
       CALL invfft ('Wave', temppsic, dfftt, howmany=nb)
       if(noncolin) CALL invfft ('Wave', temppsic(batchsize*nrxxs+1:), dfftt, howmany=nb)
       !$acc end host_data
       !
       Call stop_clock('vexx_k_setup')
       !
       Call start_clock('vexx_k_main')
          !
          INTERNAL_LOOP_ON_Q : &
          DO iq=1,nqs
             !
             ikq  = index_xkq(current_ik,iq)
             ik   = index_xk(ikq)
             xkq  = xkq_collect(:,ikq)
             !
             ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
             CALL g2_convolution(ngmt, gt, xkp, xkq, fac)
             !$acc update device(fac)
             !
             ! facmask = fac inside the ngmt sphere, 0 otherwise 
             ! (can be used in loops over nrxxs)
             !
             !$acc kernels
             facmask = 0.0_dp
             !$acc end kernels 
             !
             !$acc parallel loop private(idx)
             do ig = 1, ngmt
               idx  = nl(ig)
               facmask( idx ) = fac(ig) 
             end do
             !
             IF ( okvan .and..not.tqr ) CALL qvan_init (ngmt, xkq, xkp)
             !
             IBND_LOOP_K : &
             DO ibnd=ibnd_start,ibnd_end !for each band of psi
                !
                IF ( abs(x_occupation(ibnd,ik)) < eps_occ) CYCLE IBND_LOOP_K
                !
                x_occ = x_occupation(ibnd,ik) / nqs
                !
                !loads the phi from file
                !
                !   >>>> calculate rho in real space
                !$acc parallel loop collapse(2) private(idx)
                DO ib = 1,nb 
                  DO ir = 1, nrxxs
                     idx = (ib-1) * nrxxs + ir
                     rhoc( idx )=conjg(exxbuff(ir,ibnd,ikq))*temppsic( idx ) / omega
                     IF (noncolin) THEN
                       jdx = batchsize * nrxxs + (ib-1) * nrxxs + ir 
                       rhoc( idx )= rhoc( idx ) + conjg(exxbuff(nrxxs+ir,ibnd,ikq))*temppsic(jdx) / omega
                     END IF
                  ENDDO
                ENDDO
                !
                !   >>>> add augmentation in REAL space HERE
                IF(okvan .and. tqr) THEN ! augment the "charge" in real space
                   !$acc update host(rhoc)
                   do ib = 1,nb
                     idx = (ib-1) * nrxxs + 1 
                     CALL addusxx_r( rhoc(idx:idx+nrxxs-1), becxx(ikq)%k(:,ibnd), becpsi%k(:,im_start+ib))
                   end do 
                   !$acc update device(rhoc)
                ENDIF
                !
                !   >>>> brings it to G-space
                !$acc host_data use_device(rhoc)
                CALL fwfft('Rho', rhoc, dfftt, howmany=nb)
                !$acc end host_data
                !
                !   >>>> add augmentation in G space HERE
                IF(okvan .and. .not. tqr) THEN
                   !$acc update host(rhoc)
                   do ib = 1,nb
                     idx = (ib-1) * nrxxs + 1
                     CALL addusxx_g(dfftt, rhoc(idx), xkq, xkp, 'c', &
                        becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,im_start+ib))
                   end do 
                   !$acc update device(rhoc)
                ENDIF
                !   >>>> charge done
                !
                ! inizializations on the whole allocation dimension (batchsize)
                !
                !$acc parallel loop collapse(2) private(idx)
                do ib = 1, nb
                  do ir = 1, nrxxs
                     idx = (ib-1)*nrxxs + ir 
                     rhoc( idx ) = rhoc(idx) * facmask(ir) * x_occ 
                  enddo
                enddo
                !
                ! Add ultrasoft contribution (RECIPROCAL SPACE)
                ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
                IF(okvan .and. .not. tqr) THEN
                   !$acc update host(rhoc)
                   do ib = 1, nb
                     idx = (ib-1)*nrxxs + 1
                     CALL newdxx_g(dfftt, rhoc(idx), xkq, xkp, 'c', deexx(:,ib), &
                                 becphi_c=becxx(ikq)%k(:,ibnd))
                   end do 
                ENDIF
                !
                !brings back v in real space
                !$acc host_data use_device(rhoc)
                CALL invfft('Rho', rhoc, dfftt, howmany=nb)
                !$acc end host_data
                !
                ! Add ultrasoft contribution (REAL SPACE)
                IF(okvan .and. tqr) THEN
                  !$acc update host(rhoc)
                  do ib = 1, nb
                    idx = (ib-1)*nrxxs + 1
                    CALL newdxx_r(dfftt,rhoc(idx:idx+nrxxs-1), becxx(ikq)%k(:,ibnd), deexx(:,ib))
                  end do 
                END IF
                !
                ! Add PAW one-center contribution
                IF(okpaw) THEN
                   do ib = 1, nb
                     CALL PAW_newdxx(x_occ , becxx(ikq)%k(:,ibnd), becpsi%k(:,im_start+ib), deexx(:,ib))
                   end do 
                ENDIF
                !
                !accumulates over bands and k points
                !
                !$acc parallel loop collapse(2) private(idx)
                do ib = 1, nb
                  do ir = 1, nrxxs
                    idx = (ib-1)*nrxxs + ir 
                    result( idx ) = result( idx ) + rhoc( idx )*exxbuff(ir,ibnd,ikq)
                    IF (noncolin) THEN
                      jdx = batchsize * nrxxs + (ib-1)*nrxxs + ir 
                      result( jdx ) = result( jdx ) + rhoc( idx )*exxbuff(ir+nrxxs,ibnd,ikq)
                    END IF
                  enddo
                end do 
                !
             ENDDO &
             IBND_LOOP_K
             IF ( okvan .and..not.tqr ) CALL qvan_clean ()
             !
          ENDDO &
          INTERNAL_LOOP_ON_Q
          !
          IF(okvan) THEN
            CALL mp_sum(deexx,intra_bgrp_comm)
            CALL mp_sum(deexx,inter_bgrp_comm)
          ENDIF
          Call stop_clock('vexx_k_main')
          !
          CALL start_clock( 'vexx_k_fin' )
          !
          !$acc host_data use_device(result)
          CALL fwfft('Wave', result, dfftt, howmany=nb)
          if(noncolin) CALL fwfft('Wave', result(batchsize*nrxxs+1:), dfftt, howmany=nb)
          !$acc end host_data
          !
          !$acc host_data use_device(result)
          CALL mp_sum( result, inter_bgrp_comm)
          !$acc end host_data
          !
          !$acc parallel loop collapse(3) private(idx)
          do ipol = 1, npol
            do ib = 1, nb
              do ig = 1, n
                 idx = (ipol-1) * batchsize * nrxxs + (ib-1) * nrxxs + nl_k(ig) 
                 hpsi((ipol-1)*lda + ig,im_start+ib) = hpsi((ipol-1)*lda + ig,im_start+ib) - exxalfa*result( idx )
              enddo
            end do 
          end do
          !
          ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
          IF(okvan) THEN
            !$acc update host(hpsi(:,im_start+1:im_start+nb))
            do ib = 1, nb
              CALL add_nlxx_pot(lda, hpsi(:,im_start+ib), xkp, n, igk_k(:,current_k),&
                                            deexx(:,ib), eps_occ, exxalfa)
            end do 
            !$acc update device(hpsi(:,im_start+1:im_start+nb))
          END IF
          !
          CALL stop_clock( 'vexx_k_fin' )
          !
    END DO &
    EXTERNAL_LOOP_ON_BATCHES
    !
    !$acc end data
    !
    ! restore original FFT batching parameters
    if( batchsize0 > 0 ) then
      dfftt%batchsize = batchsize0
      dfftt%subbatchsize = subbatchsize0 
    endif
    !
    DEALLOCATE(temppsic, result )
    !
    DEALLOCATE(rhoc, fac, facmask )
    !
    IF(okvan) DEALLOCATE( deexx)
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_std_k_outer
  !-----------------------------------------------------------------------
  !
#else
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_std_gamma_one(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... Gamma-specific version of vexx
    !     (this is ported to GPU, but called only in CPU execution)
    !
    USE cell_base,      ONLY : omega
    USE wvfct,          ONLY : current_k
    USE klist,          ONLY : xk, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, nbgrp
    USE mp,             ONLY : mp_sum, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : npwt, gt, nqs, index_xkq, x_occupation, g2_convolution
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: RESULT(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL  :: global_kpoint_index
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    INTEGER :: ngmt
    !
    ngmt = dfftt%ngm
    ALLOCATE( fac(ngmt) )
    nrxxs= dfftt%nnr
    !
    ALLOCATE( RESULT(nrxxs), temppsic_dble(nrxxs), temppsic_aimag(nrxxs) )
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=0.0_DP
       psi =0.0_DP
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    ENDIF
    !$acc update device(psi)
    !
    ! Here the loops start
    !
    !$acc data create(rhoc, vc, fac, temppsic_dble, temppsic_aimag, result) copy(hpsi) present(exxbuff) 
    !
    INTERNAL_LOOP_ON_Q : &
    DO iq=1,nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution(ngmt, gt, xkp, xkq, fac)
       !$acc update device(fac)
       !
       IF ( okvan .and..not.tqr ) CALL qvan_init (ngmt, xkq, xkp)
       !
       LOOP_ON_PSI_BANDS : &
       DO im = 1,m !for each band of psi (the k cycle is outside band)
          !
          IF(okvan) deexx(:) = 0.0_DP
          !
          !$acc kernels
          RESULT = 0.0_DP
          !$acc end kernels
          !
          l_fft_doubleband = .false.
          l_fft_singleband = .false.
          !
          IF ( mod(im,2)==1 .and. (im+1)<=m ) l_fft_doubleband = .true.
          IF ( mod(im,2)==1 .and. im==m )     l_fft_singleband = .true.
          !
          IF( l_fft_doubleband ) THEN
#if defined(__CUDA)
             !$acc parallel loop 
#else
             !$omp parallel do  default(shared), private(ig)
#endif
             DO ig = 1, npwt
                RESULT( dfftt%nl(ig) )  =       psi(ig, im) + (0._DP,1._DP) * psi(ig, im+1)
                RESULT( dfftt%nlm(ig) ) = conjg(psi(ig, im) - (0._DP,1._DP) * psi(ig, im+1))
             ENDDO
#if !defined(__CUDA)
             !$omp end parallel do
#endif
          ENDIF
          !
          IF( l_fft_singleband ) THEN
#if defined(__CUDA)
             !$acc parallel loop
#else
             !$omp parallel do  default(shared), private(ig)
#endif
             DO ig = 1, npwt
                RESULT( dfftt%nl(ig) )  =       psi(ig,im)
                RESULT( dfftt%nlm(ig) ) = conjg(psi(ig,im))
             ENDDO
#if !defined(__CUDA)
             !$omp end parallel do
#endif
          ENDIF
          !
          IF( l_fft_doubleband.or.l_fft_singleband) THEN
             !$acc host_data use_device(RESULT)
             CALL invfft ('Wave', RESULT, dfftt)
             !$acc end host_data
#if defined(__CUDA) 
             !$acc parallel loop
#else
             !$omp parallel do default(shared), private(ir)
#endif
             DO ir = 1, nrxxs
                temppsic_dble(ir)  = dble ( RESULT(ir) )
                temppsic_aimag(ir) = aimag( RESULT(ir) )
             ENDDO
#if !defined(__CUDA) 
             !$omp end parallel do
#endif
          ENDIF
          !
          !$acc kernels
          RESULT = 0.0_DP
          !$acc end kernels
          !
          h_ibnd = ibnd_start/2
          IF(mod(ibnd_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start-1
          ELSE
             ibnd_loop_start=ibnd_start
          ENDIF
          !
          IBND_LOOP_GAM : &
          DO ibnd=ibnd_loop_start,ibnd_end, 2 !for each band of psi
             !
             h_ibnd = h_ibnd + 1
             IF( ibnd < ibnd_start ) THEN
                x1 = 0.0_DP
             ELSE
                x1 = x_occupation(ibnd,  ik)
             ENDIF
             IF( ibnd == ibnd_end) THEN
                x2 = 0.0_DP
             ELSE
                x2 = x_occupation(ibnd+1,  ik)
             ENDIF
             IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
             !
             ! calculate rho in real space. Gamma tricks are used.
             ! temppsic is real; tempphic contains one band in the real part,
             ! another one in the imaginary part; the same applies to rhoc
             !
             IF( mod(im,2) == 0 ) THEN
#if defined(__CUDA)
                !$acc kernels present(exxbuff, rhoc, temppsic_aimag)
#else
                !$omp parallel do default(shared), private(ir)
#endif
                DO ir = 1, nrxxs
                   rhoc(ir) = exxbuff(ir,h_ibnd,current_k) * temppsic_aimag(ir) / omega
                ENDDO
#if defined(__CUDA)
                !$acc end kernels
#else
                !$omp end parallel do
#endif
             ELSE
#if defined(__CUDA)
                !$acc kernels
#else
                !$omp parallel do default(shared), private(ir)
#endif
                DO ir = 1, nrxxs
                   rhoc(ir) = exxbuff(ir,h_ibnd,current_k) * temppsic_dble(ir) / omega
                ENDDO
#if defined(__CUDA)
                !$acc end kernels
#else
                !$omp end parallel do
#endif
             ENDIF
             !
             ! bring rho to G-space
             !
             !   >>>> add augmentation in REAL SPACE here
             IF(okvan .and. tqr) THEN
                !$acc update host(rhoc)
                IF(ibnd>=ibnd_start) &
                CALL addusxx_r(rhoc,_CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,im)))
                IF(ibnd<ibnd_end) &
                CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),_CX(becpsi%r(:,im)))
                !$acc update device(rhoc)
             ENDIF
             !
             !$acc host_data use_device(rhoc)
             CALL fwfft ('Rho', rhoc, dfftt)
             !$acc end host_data
             !   >>>> add augmentation in G SPACE here
             IF(okvan .and. .not. tqr) THEN
                !$acc update host(rhoc)
                ! contribution from one band added to real (in real space) part of rhoc
                IF(ibnd>=ibnd_start) &
                   CALL addusxx_g(dfftt, rhoc, xkq,  xkp, 'r', &
                   becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,im) )
                ! contribution from following band added to imaginary (in real space) part of rhoc
                IF(ibnd<ibnd_end) &
                   CALL addusxx_g(dfftt, rhoc, xkq,  xkp, 'i', &
                   becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,im) )
                !$acc update device(rhoc)
             ENDIF
             !   >>>> charge density done
             !
             !$acc kernels
             vc = 0._DP
             !$acc end kernels
             !
#if defined(__CUDA) 
             !$acc kernels present(vc, fac, rhoc)
#else 
             !$omp parallel do default(shared), private(ig)
#endif
             DO ig = 1, ngmt
                !
                vc(dfftt%nl(ig))  = fac(ig) * rhoc(dfftt%nl(ig))
                vc(dfftt%nlm(ig)) = fac(ig) * rhoc(dfftt%nlm(ig))
                !
             ENDDO
#if defined(__CUDA)
             !$acc end kernels
#else
             !$omp end parallel do
#endif
             !
             !   >>>>  compute <psi|H_fock G SPACE here
             IF(okvan .and. .not. tqr) THEN
                !$acc update host(vc)
                IF(ibnd>=ibnd_start) &
                CALL newdxx_g(dfftt, vc, xkq, xkp, 'r', deexx, &
                              becphi_r=x1*becxx(ikq)%r(:,ibnd))
                IF(ibnd<ibnd_end) &
                CALL newdxx_g(dfftt, vc, xkq, xkp, 'i', deexx, &
                              becphi_r=x2*becxx(ikq)%r(:,ibnd+1))
             ENDIF
             !
             !brings back v in real space
             !$acc host_data use_device(vc)
             CALL invfft ('Rho', vc, dfftt)
             !$acc end host_data
             !
             !   >>>>  compute <psi|H_fock REAL SPACE here
             IF(okvan .and. tqr) THEN
                !$acc update host(vc)
                IF(ibnd>=ibnd_start) &
                CALL newdxx_r(dfftt, vc, _CX(x1*becxx(ikq)%r(:,ibnd)), deexx)
                IF(ibnd<ibnd_end) &
                CALL newdxx_r(dfftt, vc, _CY(x2*becxx(ikq)%r(:,ibnd+1)), deexx)
             ENDIF
             !
             IF(okpaw) THEN
                IF(ibnd>=ibnd_start) &
                CALL PAW_newdxx(x1/nqs, _CX(becxx(ikq)%r(:,ibnd)),&
                                        _CX(becpsi%r(:,im)), deexx)
                IF(ibnd<ibnd_end) &
                CALL PAW_newdxx(x2/nqs, _CX(becxx(ikq)%r(:,ibnd+1)), &
                                        _CX(becpsi%r(:,im)), deexx)
             ENDIF
             !
             ! accumulates over bands and k points
             !
#if defined(__CUDA) 
             !$acc kernels present(RESULT, vc, exxbuff)
#else
             !$omp parallel do default(shared), private(ir)
#endif
             DO ir = 1, nrxxs
                RESULT(ir) = RESULT(ir)+x1* dble(vc(ir))* dble(exxbuff(ir,h_ibnd,current_k))&
                                       +x2*aimag(vc(ir))*aimag(exxbuff(ir,h_ibnd,current_k))
             ENDDO
#if defined(__CUDA) 
             !$acc end kernels 
#else
             !$omp end parallel do
#endif
             !
          ENDDO &
          IBND_LOOP_GAM
          !
          IF(okvan) THEN
             CALL mp_sum(deexx,intra_bgrp_comm)
             CALL mp_sum(deexx,inter_bgrp_comm)
          ENDIF
          !
          !$acc host_data use_device(result)
          CALL mp_sum( RESULT(1:nrxxs), inter_bgrp_comm)
          !$acc end host_data
          !
          ! brings back result in G-space
          !
          !$acc host_data use_device(RESULT)
          CALL fwfft( 'Wave' , RESULT, dfftt )
          !$acc end host_data
          !
#if defined(__CUDA) 
          !$acc kernels present(hpsi, RESULT, dfftt%nl)
#else  
          !$omp parallel do default(shared), private(ig)
#endif
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*RESULT(dfftt%nl(ig))
          ENDDO
#if defined(__CUDA)
          !$acc end kernels
#else
          !$omp end parallel do
#endif
          ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
          IF(okvan) THEN
            !$acc update host(hpsi(:,im))
            CALL add_nlxx_pot (lda, hpsi(:,im), xkp, n, &
                           igk_k(1,current_k), deexx, eps_occ, exxalfa)
            !$acc update device(hpsi(:,im))
          END IF
       ENDDO &
       LOOP_ON_PSI_BANDS
       IF ( okvan .and..not.tqr ) CALL qvan_clean ()
       !
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !
    !$acc end data
    !
    DEALLOCATE( RESULT, temppsic_dble, temppsic_aimag)
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_std_gamma_one
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_std_k_one(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... generic, k-point version of vexx
    !     (this is ported to GPU, but called only in CPU execution)
    !
    USE cell_base,      ONLY : omega
    USE wvfct,          ONLY : npwx, current_k
    USE klist,          ONLY : xk, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, nbgrp
    USE mp,             ONLY : mp_sum, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : gt, nqs, index_xkq, x_occupation, g2_convolution
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: temppsic(:), RESULT(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:),result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: result_g(:), result_nc_g(:,:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: nrxxs
    REAL(DP) :: xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ngmt
    REAL(DP) :: x_occ
    !
    ngmt = dfftt%ngm 
    ALLOCATE( fac(ngmt) )
    nrxxs= dfftt%nnr
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol), result_nc(nrxxs,npol) )
       ALLOCATE( result_nc_g(n,npol) )
    ELSE
       ALLOCATE( temppsic(nrxxs), RESULT(nrxxs) )
       ALLOCATE( result_g(n) )
    ENDIF
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=0.0_DP
       psi =0.0_DP
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    ENDIF
    !$acc update device(psi)
    !
    !$acc data create(rhoc, vc, fac, temppsic, result, result_g, temppsic_nc, result_nc, result_nc_g) copy(hpsi) present(exxbuff) 
    !
    LOOP_ON_PSI_BANDS : &
    DO im = 1,m !for each band of psi (the k cycle is outside band)
       IF(okvan) deexx = 0._DP
       !
       IF (noncolin) THEN
          !$acc kernels
          temppsic_nc = 0._DP
          !$acc end kernels
       ELSE
          !$acc kernels
          temppsic    = 0._DP
          !$acc end kernels
       ENDIF
       !
       IF (noncolin) THEN
          !
#if defined(__CUDA)
          !$acc parallel loop
#else
          !$omp parallel do  default(shared), private(ig)
#endif
          DO ig = 1, n
             temppsic_nc(dfftt%nl(igk_k(ig,current_k)),1) = psi(ig,im)
             temppsic_nc(dfftt%nl(igk_k(ig,current_k)),2) = psi(npwx+ig,im)
          ENDDO
#if !defined(__CUDA)
          !$omp end parallel do
#endif
          !
          !$acc host_data use_device(temppsic_nc)
          CALL invfft ('Wave', temppsic_nc(:,1), dfftt)
          CALL invfft ('Wave', temppsic_nc(:,2), dfftt)
          !$acc end host_data
          !
       ELSE
          !
#if defined(__CUDA)
          !$acc parallel loop
#else
          !$omp parallel do  default(shared), private(ig)
#endif
          DO ig = 1, n
             temppsic( dfftt%nl(igk_k(ig,current_k)) ) = psi(ig,im)
          ENDDO
#if !defined(__CUDA)
          !$omp end parallel do
#endif
          !
          !$acc host_data use_device(temppsic)
          CALL invfft ('Wave', temppsic, dfftt)
          !$acc end host_data
          !
       ENDIF
       !
       IF (noncolin) THEN
          !$acc kernels
          result_nc = 0.0_DP
          !$acc end kernels
       ELSE
          !$acc kernels
          RESULT    = 0.0_DP
          !$acc end kernels
       ENDIF
       !
       INTERNAL_LOOP_ON_Q : &
       DO iq=1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          xkq  = xkq_collect(:,ikq)
          !
          ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
          CALL g2_convolution(ngmt, gt, xkp, xkq, fac)
          !$acc update device(fac)
          IF ( okvan .and..not.tqr ) CALL qvan_init (ngmt, xkq, xkp)
          !
          IBND_LOOP_K : &
          DO ibnd=ibnd_start,ibnd_end !for each band of psi
             !
             IF ( abs(x_occupation(ibnd,ik)) < eps_occ) CYCLE IBND_LOOP_K
             !
             x_occ = x_occupation(ibnd,ik) / nqs
             !
             !loads the phi from file
             !
             !   >>>> calculate rho in real space
             IF (noncolin) THEN
#if defined(__CUDA)
                !$acc kernels 
#else
                !$omp parallel do default(shared), private(ir)
#endif
                DO ir = 1, nrxxs
                   rhoc(ir) = ( conjg(exxbuff(ir,ibnd,ikq))*temppsic_nc(ir,1) + &
                                conjg(exxbuff(nrxxs+ir,ibnd,ikq))*temppsic_nc(ir,2) )/omega
                ENDDO
#if defined(__CUDA)
                !$acc end kernels
#else
                !$omp end parallel do
#endif
             ELSE
#if defined(__CUDA)
                !$acc kernels 
#else
                !$omp parallel do default(shared), private(ir)
#endif
                DO ir = 1, nrxxs
                   rhoc(ir)=conjg(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
                ENDDO
#if defined(__CUDA)
                !$acc end kernels
#else
                !$omp end parallel do
#endif
             ENDIF
             !   >>>> add augmentation in REAL space HERE
             IF(okvan .and. tqr) THEN ! augment the "charge" in real space
                !$acc update host(rhoc)
                CALL addusxx_r( rhoc, becxx(ikq)%k(:,ibnd), becpsi%k(:,im))
                !$acc update device(rhoc)
             ENDIF
             !
             !   >>>> brings it to G-space
             !$acc host_data use_device(rhoc)
             CALL fwfft('Rho', rhoc, dfftt)
             !$acc end host_data
             !
             !   >>>> add augmentation in G space HERE
             IF(okvan .and. .not. tqr) THEN
                !$acc update host(rhoc)
                CALL addusxx_g(dfftt, rhoc, xkq, xkp, 'c', &
                   becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,im))
                !$acc update device(rhoc)
             ENDIF
             !   >>>> charge done
             !
             !$acc kernels
             vc = 0._DP
             !$acc end kernels
             !
#if defined(__CUDA)
             !$acc parallel loop present(dfftt, dfftt%nl)
#else
             !$omp parallel do default(shared), private(ig)
#endif
             DO ig = 1, ngmt
                vc(dfftt%nl(ig)) = fac(ig) * rhoc(dfftt%nl(ig)) * x_occ 
             ENDDO
#if !defined(__CUDA)
             !$omp end parallel do
#endif
             !
             ! Add ultrasoft contribution (RECIPROCAL SPACE)
             ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
             IF(okvan .and. .not. tqr) THEN
                !$acc update host(vc)
                CALL newdxx_g(dfftt, vc, xkq, xkp, 'c', deexx, &
                              becphi_c=becxx(ikq)%k(:,ibnd))
             ENDIF
             !
             !brings back v in real space
             !$acc host_data use_device(vc)
             CALL invfft ('Rho', vc, dfftt)
             !$acc end host_data
             !
             ! Add ultrasoft contribution (REAL SPACE)
             IF(okvan .and. tqr) THEN
               !$acc update host(vc)
               CALL newdxx_r(dfftt,vc, becxx(ikq)%k(:,ibnd), deexx)
             END IF
             !
             ! Add PAW one-center contribution
             IF(okpaw) THEN
                CALL PAW_newdxx(x_occ , becxx(ikq)%k(:,ibnd), becpsi%k(:,im), deexx)
             ENDIF
             !
             !accumulates over bands and k points
             !
             IF (noncolin) THEN
#if defined(__CUDA)
                !$acc kernels
#else
                !$omp parallel do default(shared), private(ir)
#endif
                DO ir = 1, nrxxs
                   result_nc(ir,1)= result_nc(ir,1) + vc(ir) * exxbuff(ir,ibnd,ikq)
                   result_nc(ir,2)= result_nc(ir,2) + vc(ir) * exxbuff(ir+nrxxs,ibnd,ikq)
                ENDDO
#if defined(__CUDA)
                !$acc end kernels
#else
                !$omp end parallel do
#endif
             ELSE
#if defined(__CUDA)
                !$acc kernels
#else
                !$omp parallel do default(shared), private(ir)
#endif
                DO ir = 1, nrxxs
                   RESULT(ir) = RESULT(ir) + vc(ir)*exxbuff(ir,ibnd,ikq)
                ENDDO
#if defined(__CUDA)
                !$acc end kernels
#else
                !$omp end parallel do
#endif
             ENDIF
             !
          ENDDO &
          IBND_LOOP_K
          IF ( okvan .and..not.tqr ) CALL qvan_clean ()
          !
       ENDDO &
       INTERNAL_LOOP_ON_Q
       !
       IF(okvan) THEN
         CALL mp_sum(deexx,intra_bgrp_comm)
         CALL mp_sum(deexx,inter_bgrp_comm)
       ENDIF
       !
       ! bring result back to G-space
       !
       IF (noncolin) THEN
          !
          !$acc host_data use_device(result_nc)
          CALL fwfft ('Wave', result_nc(:,1), dfftt)
          CALL fwfft ('Wave', result_nc(:,2), dfftt)
          !$acc end host_data
          !
          !communicate result
          !$acc parallel loop collapse(2)
          DO ipol = 1, npol
            DO ig = 1, n
               result_nc_g(ig,ipol) = result_nc(dfftt%nl(igk_k(ig,current_k)),ipol)
            ENDDO
          ENDDO
          !$acc host_data use_device(result_nc_g)
          CALL mp_sum( result_nc_g(1:n,1:npol), inter_bgrp_comm)
          !$acc end host_data
          !
          !adds it to hpsi
#if defined(__CUDA)
          !$acc parallel loop
#else
          !$omp parallel do default(shared), private(ig)
#endif
          DO ig = 1, n
             hpsi(ig,im)    = hpsi(ig,im)     - exxalfa*result_nc_g(ig,1)
             hpsi(lda+ig,im)= hpsi(lda+ig,im) - exxalfa*result_nc_g(ig,2)
          ENDDO
#if !defined(__CUDA)
          !$omp end parallel do
#endif
          !
       ELSE
          !
          !$acc host_data use_device(result)
          CALL fwfft ('Wave', RESULT, dfftt)
          !$acc end host_data
          !
          !communicate result
          !$acc parallel loop
          DO ig = 1, n
             result_g(ig) = RESULT(dfftt%nl(igk_k(ig,current_k)))
          ENDDO
          !$acc host_data use_device(result_g)
          CALL mp_sum( result_g(1:n), inter_bgrp_comm)
          !$acc end host_data
          !
          !adds it to hpsi
#if defined(__CUDA) 
          !$acc parallel loop 
#else
          !$omp parallel do default(shared), private(ig)
#endif
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*result_g(ig)
          ENDDO
#if !defined(__CUDA) 
          !$omp end parallel do
#endif
       ENDIF
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) THEN
         !$acc update host(hpsi(:,im))
         CALL add_nlxx_pot(lda, hpsi(:,im), xkp, n, igk_k(:,current_k),&
                                       deexx, eps_occ, exxalfa)
         !$acc update device(hpsi(:,im))
       END IF
       !
    ENDDO &
    LOOP_ON_PSI_BANDS
    !
    !$acc end data
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, result_nc, result_nc_g )
    ELSE
       DEALLOCATE(temppsic, RESULT, result_g )
    ENDIF
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx)
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_std_k_one
  !-----------------------------------------------------------------------
  !
#endif
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_stress_std()
    !-----------------------------------------------------------------------
    !
    ! This is Eq.(10) of PRB 73, 125120 (2006).
    !
    USE constants,            ONLY : fpi, e2, pi
    USE io_files,             ONLY : iunwfc, nwordwfc
    USE buffers,              ONLY : get_buffer
    USE cell_base,            ONLY : omega, bg, at, tpiba
    USE symm_base,            ONLY : s
    USE wvfct,                ONLY : nbnd, wg, current_k
    USE wavefunctions,        ONLY : evc
    USE klist,                ONLY : xk, ngk, igk_k
    USE lsda_mod,             ONLY : lsda, current_spin, isk
    USE gvect,                ONLY : g
    USE mp_pools,             ONLY : npool, inter_pool_comm
    USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE uspp,                 ONLY : okvan
    USE exx_base,             ONLY : exxdiv, erfc_scrlen, gau_scrlen, grid_factor, eps, &
                                     nq1, nq2, nq3, nqs, on_double_grid, use_coulomb_vcut_spheric, &
                                     use_coulomb_vcut_ws, vcut, x_gamma_extrapolation, yukawa, index_xkq, &
                                     index_sym, x_occupation, exx_divergence_ew, exx_divergence_stress, &
                                     exxbeta
    USE coulomb_vcut_module,  ONLY : vcut_get,  vcut_spheric_get
    USE xc_lib,               ONLY : xclib_get_id
    !
    ! ---- local variables -------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! local variables
    REAL(DP)   :: exx_stress_std(3,3), exx_stress_(3,3)
    !
    COMPLEX(DP),ALLOCATABLE :: tempphic(:), temppsic(:)
    COMPLEX(DP),ALLOCATABLE :: rhoc(:)
    REAL(DP),    ALLOCATABLE :: fac(:), fac_tens(:,:,:), fac_stress(:), fac_stress_div(:)
    INTEGER  :: npw, jbnd, ibnd, ik, ikk, ig, ir, ikq, iq, isym
    INTEGER  :: h_ibnd, nqi, iqi, beta, nrxxs, ngm
    INTEGER  :: ibnd_loop_start, igcx
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    REAL(DP) :: cam_beta, scr_exp, qscr
    ! temp array for vcut_spheric
    REAL(DP) :: delta(3,3)
    ! q=0 divergence: Ewald part of exxdiv and strain derivative of its G-sum part
    REAL(DP) :: exxdivergence_ew, exxdiv_stress(3,3)
    INTEGER :: ngmt
    !
    IF (npool>1) CALL errore('exx_stress1','stress not available with pools',1)
    IF (noncolin) CALL errore('exx_stress1','noncolinear stress not implemented',1)
    IF (okvan) CALL infomsg('exx_stress1','USPP stress not tested')
    igcx = xclib_get_id('GGA','EXCH')
    cam_beta = 0._DP
    IF (igcx == 52) THEN
        IF (exxalfa == 0._DP) CALL errore('exx_stress1','CAM stress requires nonzero exx_fraction',1)
        cam_beta = exxbeta / exxalfa
    ENDIF
    !
    nrxxs = dfftt%nnr
    ngmt = dfftt%ngm
    ngm   = ngmt
    delta = reshape( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    exxdivergence_ew = exx_divergence_ew()
    exxdiv_stress    = exx_divergence_stress()
    ALLOCATE( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    ALLOCATE( fac_tens(3,3,ngm), fac_stress(ngm), fac_stress_div(ngm) )
    !
    nqi=nqs
    !
    ! loop over k-points
    DO ikk = 1, nks
        current_k = ikk
        IF (lsda) current_spin = isk(ikk)
        npw = ngk(ikk)

        IF (nks > 1) &
            CALL get_buffer(evc, nwordwfc, iunwfc, ikk)

        ! loop over bands
        DO jbnd = 1, nbnd
            !
            temppsic(:) = ( 0._dp, 0._dp )
!$omp parallel do default(shared), private(ig)
            DO ig = 1, npw
                temppsic(dfftt%nl(igk_k(ig,ikk))) = evc(ig,jbnd)
            ENDDO
!$omp end parallel do
            !
            IF(gamma_only) THEN
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw
                    temppsic(dfftt%nlm(igk_k(ig,ikk))) = conjg(evc(ig,jbnd))
                ENDDO
!$omp end parallel do
            ENDIF

            CALL invfft ('Wave', temppsic, dfftt)

            DO iqi = 1, nqi
                !
                iq=iqi
                !
                ikq  = index_xkq(current_k,iq)
                ik   = index_xk(ikq)
                isym = abs(index_sym(ikq))

                ! FIXME: use cryst_to_cart and company as above..
                xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
                IF (index_sym(ikq) < 0) xk_cryst = -xk_cryst
                sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                         s(:,2,isym)*xk_cryst(2) + &
                         s(:,3,isym)*xk_cryst(3)
                xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

                !CALL start_clock ('exxen2_ngmloop')

!$omp parallel do default(shared), private(ig, beta, q, qq, on_double_grid, x, qscr, scr_exp)
                DO ig = 1, ngm
                  q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                  q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                  q(3)= xk(3,current_k) - xkq(3) + g(3,ig)

                  q = q * tpiba
                  qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )

                  DO beta = 1, 3
                      fac_tens(1:3,beta,ig) = q(1:3)*q(beta)
                  ENDDO

                  IF (x_gamma_extrapolation) THEN
                      on_double_grid = .true.
                      x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                  ELSE
                      on_double_grid = .false.
                  ENDIF

                  IF (use_coulomb_vcut_ws) THEN
                      fac(ig) = vcut_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      fac_stress_div(ig) = 0._dp
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSEIF ( use_coulomb_vcut_spheric ) THEN
                      fac(ig) = vcut_spheric_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      fac_stress_div(ig) = 0._dp
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSEIF (gau_scrlen > 0) THEN
                      fac(ig)=e2*((pi/gau_scrlen)**(1.5d0))* &
                            exp(-qq/4.d0/gau_scrlen) * grid_factor
                      fac_stress(ig) =  e2*2.d0/4.d0/gau_scrlen * &
                            exp(-qq/4.d0/gau_scrlen) *((pi/gau_scrlen)**(1.5d0))* &
                                                                     grid_factor
                      fac_stress_div(ig) = 0._dp
                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSEIF (qq > 1.d-8) THEN
                      fac_stress_div(ig) = 0._dp
                      IF ( erfc_scrlen > 0 ) THEN
                        qscr = qq / (4.d0 * erfc_scrlen**2)
                        scr_exp = exp(-qscr)
                        IF (igcx == 52) THEN
                          fac(ig) = e2*fpi/qq*(1._dp + cam_beta*scr_exp) * grid_factor
                          fac_stress(ig) = e2*fpi * 2.d0/qq**2 * &
                              (1._dp + cam_beta*scr_exp*(1._dp+qscr)) * grid_factor
                        ELSE
                          fac(ig)=e2*fpi/qq*(1._dp-scr_exp) * grid_factor
                          fac_stress(ig) = -e2*fpi * 2.d0/qq**2 * ( &
                              (1._dp+qscr)*scr_exp - 1._dp) * &
                              grid_factor
                        ENDIF
                      ELSE
                        fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
                        fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2 * grid_factor
                      ENDIF

                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSE
                      ! exxdiv = D_sum + D_ew, with D_ew proportional to omega:
                      ! -delta*fac gives delta*D_sum, fac_stress_div adds -dD_sum/deps
                      fac(ig) = -exxdiv + exxdivergence_ew
                      fac_stress(ig) = 0._dp
                      fac_stress_div(ig) = 1._dp
                      IF ( yukawa> 0._dp .and. .not. x_gamma_extrapolation) THEN
                        fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
                        fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2
                      ENDIF
                      IF (erfc_scrlen > 0._dp .and. .not. x_gamma_extrapolation) THEN
                        IF (igcx == 52) THEN
                          fac(ig) = fac(ig) - cam_beta * e2*fpi / (4.d0*erfc_scrlen**2)
                          fac_stress(ig) = -cam_beta * e2*fpi / (8.d0*erfc_scrlen**4)
                        ELSE
                          fac(ig) = fac(ig) + e2*fpi / (4.d0*erfc_scrlen**2)
                          fac_stress(ig) = e2*fpi / (8.d0*erfc_scrlen**4)
                        ENDIF
                      ENDIF
                  ENDIF
                ENDDO
!$omp end parallel do
                !CALL stop_clock ('exxen2_ngmloop')

                IF (gamma_only) THEN
                    !
                    h_ibnd = ibnd_start/2
                    !
                    IF(mod(ibnd_start,2)==0) THEN
                      h_ibnd=h_ibnd-1
                      ibnd_loop_start=ibnd_start-1
                    ELSE
                      ibnd_loop_start=ibnd_start
                    ENDIF
                    !
                    DO ibnd = ibnd_loop_start, ibnd_end, 2     !for each band of psi
                        !
                        h_ibnd = h_ibnd + 1
                        !
                        IF( ibnd < ibnd_start ) THEN
                            x1 = 0._dp
                        ELSE
                            x1 = x_occupation(ibnd,  ik)
                        ENDIF

                        IF( ibnd == ibnd_end) THEN
                            x2 = 0._dp
                        ELSE
                            x2 = x_occupation(ibnd+1,  ik)
                        ENDIF
                        IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
                        !
                        ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                        DO ir = 1, nrxxs
                            tempphic(ir) = exxbuff(ir,h_ibnd,ikq)
                            rhoc(ir)     = conjg(tempphic(ir))*temppsic(ir) / omega
                        ENDDO
!$omp end parallel do
                        ! bring it to G-space
                        CALL fwfft ('Rho', rhoc, dfftt)

                        vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                        DO ig = 1, ngm
                            !
                            vc(:,:) = vc(:,:) + x1 * 0.25_dp * &
                                      abs( rhoc(dfftt%nl(ig)) + &
                                      conjg(rhoc(dfftt%nlm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig) - delta(:,:)*fac(ig) + &
                                       exxdiv_stress(:,:)*fac_stress_div(ig))
                            vc(:,:) = vc(:,:) + x2 * 0.25_dp * &
                                      abs( rhoc(dfftt%nl(ig)) - &
                                      conjg(rhoc(dfftt%nlm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig) - delta(:,:)*fac(ig) + &
                                       exxdiv_stress(:,:)*fac_stress_div(ig))
                        ENDDO
!$omp end parallel do
                        vc = vc / nqs / 2.d0
                        exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                    ENDDO

                ELSE

                    DO ibnd = ibnd_start, ibnd_end    !for each band of psi
                      !
                      IF ( abs(x_occupation(ibnd,ik)) < 1.d-6) CYCLE
                      !
                      ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                          tempphic(ir) = exxbuff(ir,ibnd,ikq)
                          rhoc(ir)     = conjg(tempphic(ir))*temppsic(ir) / omega
                      ENDDO
!$omp end parallel do

                      ! bring it to G-space
                      CALL fwfft ('Rho', rhoc, dfftt)

                      vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm
                          vc(:,:) = vc(:,:) + rhoc(dfftt%nl(ig))  * &
                                        conjg(rhoc(dfftt%nl(ig)))* &
                                    (fac_tens(:,:,ig)*fac_stress(ig) - delta(:,:)*fac(ig) + &
                                       exxdiv_stress(:,:)*fac_stress_div(ig))
                      ENDDO
!$omp end parallel do
                      vc = vc * x_occupation(ibnd,ik) / nqs / 2.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)

                    ENDDO

                ENDIF ! gamma or k-points

            ENDDO ! iqi
        ENDDO ! jbnd
    ENDDO ! ikk

    DEALLOCATE(tempphic, temppsic, rhoc, fac, fac_tens, fac_stress, fac_stress_div )
    !
    CALL mp_sum( exx_stress_, intra_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    exx_stress_std = exx_stress_
    !
    !-----------------------------------------------------------------------
  END FUNCTION exx_stress_std
  !-----------------------------------------------------------------------
  !
END MODULE exx_std
