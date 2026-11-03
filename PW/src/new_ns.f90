!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE new_ns_hubbard ( noncolin, rho )
!-----------------------------------------------------------------------
  !! Wrapper for all DFT+U+V+J+... cases
  !
  USE ions_base, ONLY : ntyp => nsp
  USE ldaU,      ONLY : lda_plus_u_kind, is_hubbard_back
  USE scf,       ONLY : scf_type
  !
  IMPLICIT NONE
  TYPE(scf_type), INTENT(INOUT) :: rho
  LOGICAL, INTENT(in) :: noncolin
  INTEGER :: nt
  !
  CALL start_clock( 'new_ns' )
  !
  IF (lda_plus_u_kind==0) THEN
     !
     IF (noncolin) THEN
        CALL new_ns_nc(rho%ns_nc)
     ELSE
        CALL new_ns(rho%ns)
     ENDIF
     !
     DO nt = 1, ntyp
        IF (is_hubbard_back(nt)) CALL new_nsb( rho%nsb )
     ENDDO
     !
  ELSEIF (lda_plus_u_kind==1) THEN
     !
     IF (noncolin) THEN
        CALL new_ns_nc( rho%ns_nc )
     ELSE
        CALL new_ns( rho%ns )
     ENDIF
     !
  ELSEIF (lda_plus_u_kind==2) THEN
     !
     IF (noncolin) THEN
        CALL new_nsg_nc( rho%nsg )
     ELSE
        CALL new_nsg( rho%nsg )
     ENDIF
     !
  ENDIF
  !
  CALL stop_clock( 'new_ns' )
  !
END SUBROUTINE new_ns_hubbard

!-----------------------------------------------------------------------
SUBROUTINE new_ns( ns )
  !-----------------------------------------------------------------------
  !! This routine computes the new value for ns (the occupation numbers of
  !! orthogonalized atomic wfcs).  
  !! These quantities are defined as follows:
  !!
  !! ns_{I,s,m1,m2} = \sum_{k,v}
  !! f_{kv} <\fi^{at}_{I,m1}|\psi_{k,v,s}><\psi_{k,v,s}|\fi^{at}_{I,m2}>
  !!
  !! ns is symmetric (and real for collinear case, due to time reversal),
  !! hence the order of m1 and m2 does not matter.  
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : ldmx, Hubbard_l, q_ae, wfcU, &
                                   Hubbard_projectors, is_hubbard, nwfcU, offsetU
  USE symm_base,            ONLY : d1, d2, d3
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE symm_base,            ONLY : nsym, irt, t_rev
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE control_flags,        ONLY : gamma_only, offload_type
  USE wavefunctions,        ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
  USE buffers,              ONLY : get_buffer
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE becmod,               ONLY : bec_type, calbec, allocate_bec_type_acc, &
                                   deallocate_bec_type_acc
  USE noncollin_module,     ONLY : colin_mag
#if defined (__OSCDFT)
  USE plugin_flags,         ONLY : use_oscdft
  USE oscdft_base,          ONLY : oscdft_ctx
#endif  
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: ns(ldmx,ldmx,nspin,nat)
  !! the occupation numbers of orthogonalized atomic wfcs
  !
  ! ... local variables
  !
  TYPE(bec_type) :: proj
  ! proj(nwfcU,nbnd)
  INTEGER :: ik, ibnd, is, i, na, nb, nt, isym, m1, m2, m0, m00, npw, is2
  ! counters
  INTEGER :: off, ldim
  REAL(DP), ALLOCATABLE :: nr(:,:,:,:), r1(:,:), r2(:,:)
  COMPLEX(DP), ALLOCATABLE :: c1(:,:), c2(:,:)
  !$acc declare device_resident(r1,r2,c1,c2)
  REAL(DP) :: psum
  !
  CALL allocate_bec_type_acc( nwfcU, nbnd, proj ) 
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  ALLOCATE( nr(ldmx, ldmx, nspin, nat) )
  IF ( gamma_only ) THEN
     ALLOCATE ( r1(nbnd,ldmx), r2(nbnd,ldmx) )
  ELSE
     ALLOCATE ( c1(nbnd,ldmx), c2(nbnd,ldmx) )
  ENDIF
  !
  ! we start a loop on k points
  !$acc data create(nr) present(evc, wfcU, proj) copyin(wg)
  !$acc kernels
  nr(:,:,:,:) = 0.0_dp
  !$acc end kernels
  DO ik = 1, nks
     !
     IF (lsda) current_spin = isk(ik)
     !
     npw = ngk(ik)
     !
     IF (nks > 1) THEN
        CALL get_buffer (evc, nwordwfc, iunwfc, ik)
        !$acc update device(evc)
     END IF
     !
     ! make the projection
     !
     IF ( Hubbard_projectors == 'pseudo' ) THEN
        CALL compute_pproj( ik, q_ae, proj )
     ELSE
        IF (nks > 1)  THEN
           CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )
           !$acc update device(wfcU)
        END IF
        CALL calbec( offload_type, npw, wfcU, evc, proj )
     ENDIF
     !
     ! compute the occupation matrix (ns_{I,s,m1,m2}) of the
     ! atomic orbitals
     !
     DO na = 1, nat 
        nt = ityp(na)  
        IF ( is_hubbard(nt) ) THEN 
           off = offsetU(na)
           ldim = 2*Hubbard_l(nt) + 1
           ! Next lines compute nr using matrix-matrix multiplication
           ! summing over the band index. For complex proj%k, only the
           ! real part of the result is computed, using real algebra
           IF ( gamma_only ) THEN
              !$acc parallel loop collapse(2) present(proj%r,wg)
              DO m1 = 1, ldim
                 DO ibnd = 1, nbnd  
                    r1(ibnd,m1) = proj%r(off+m1,ibnd)
                    r2(ibnd,m1) = proj%r(off+m1,ibnd) * wg(ibnd,ik) 
                 ENDDO
              ENDDO
              !$acc host_data use_device(r1,r2,nr)
              CALL MYDGEMM( 'T','N', ldim, ldim, nbnd, 1.0_dp, r1, nbnd, &
                            r2, nbnd, 1.0_dp, nr(1,1,current_spin,na), ldmx )
              !$acc end host_data 
           ELSE
              !$acc parallel loop collapse(2) present(proj%k,wg)
              DO m1 = 1, ldim
                 DO ibnd = 1, nbnd
                    c1(ibnd,m1) = proj%k(off+m1,ibnd)
                    c2(ibnd,m1) = proj%k(off+m1,ibnd) * wg(ibnd,ik)
                 ENDDO
              ENDDO
              !$acc host_data use_device(c1,c2,nr)
              CALL MYDGEMM( 'C','N', ldim, ldim, 2*nbnd, 1.0_dp, &
                c1, 2*nbnd, c2, 2*nbnd, 1.0_dp, nr(1,1,current_spin,na), ldmx )
              !$acc end host_data 
           ENDIF
        ENDIF
     ENDDO
     !
  ENDDO
  !$acc update host(nr)
  IF ( gamma_only ) THEN
     DEALLOCATE ( r1, r2 )
  ELSE
     DEALLOCATE ( c1, c2 )
  ENDIF
  !
  CALL deallocate_bec_type_acc( proj ) 
  !$acc end data 
  !
  CALL mp_sum( nr, inter_pool_comm )
  !
  IF (nspin == 1) nr = 0.5d0 * nr
  !
  ! impose hermiticity of n_{m1,m2}
  !
  DO na = 1, nat  
     nt = ityp(na)
     DO is = 1, nspin  
        DO m1 = 1, 2*Hubbard_l(nt)+1
           DO m2 = m1+1, 2*Hubbard_l(nt)+1 
              nr(m2,m1,is,na) = nr(m1,m2,is,na)  
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2) .AND. &
                  .NOT.oscdft_ctx%inp%constraint_diag) THEN
     IF (.NOT.oscdft_ctx%conv) THEN
        ns(:,:,:,:) = nr(:,:,:,:)
        DEALLOCATE(nr)
        RETURN
     ENDIF
  ENDIF
#endif
  !
  ! symmetrize the quantities nr -> ns
  !
  ns (:,:,:,:) = 0.d0
  !
  DO na = 1, nat  
     nt = ityp(na)
     IF ( is_hubbard(nt) ) THEN
        DO is = 1, nspin
           !
           DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                 !
                 DO isym = 1, nsym
                    nb = irt (isym, na)  
                    ! flip spin for time-reversal in collinear systems
                    IF ( (colin_mag==2) .AND. (t_rev(isym) == 1) ) THEN
                       is2 = 3-is
                    ELSE
                       is2 = is
                    END IF     
                    !
                    DO m0 = 1, 2*Hubbard_l(nt)+1  
                       DO m00 = 1, 2*Hubbard_l(nt)+1  
                          !
                          IF (Hubbard_l(nt) == 0) THEN
                             ns(m1,m2,is,na) = ns(m1,m2,is,na) +  &
                                               nr(m0,m00,is2,nb) / nsym
                          ELSEIF (Hubbard_l(nt) == 1) THEN
                             ns(m1,m2,is,na) = ns(m1,m2,is,na) +  &
                                               d1(m0 ,m1,isym) * nr(m0,m00,is2,nb) * &
                                               d1(m00,m2,isym) / nsym
                          ELSEIF (Hubbard_l(nt) == 2) THEN
                             ns(m1,m2,is,na) = ns(m1,m2,is,na) +  &
                                               d2(m0 ,m1,isym) * nr(m0,m00,is2,nb) * &
                                               d2(m00,m2,isym) / nsym
                          ELSEIF (Hubbard_l(nt) == 3) THEN
                             ns(m1,m2,is,na) = ns(m1,m2,is,na) +  &
                                               d3(m0 ,m1,isym) * nr(m0,m00,is2,nb) * &
                                               d3(m00,m2,isym) / nsym
                          ELSE
                             CALL errore( 'new_ns', &
                                          'angular momentum not implemented', &
                                          ABS(Hubbard_l(nt)) )
                          ENDIF
                          !
                       ENDDO
                    ENDDO
                    !
                 ENDDO
                 !
              ENDDO
           ENDDO
           !
        ENDDO
     ENDIF
  ENDDO
  !
  DEALLOCATE( nr )
  !
  ! Now we make the matrix ns(m1,m2) strictly hermitean
  !
  DO na = 1, nat  
     nt = ityp (na)  
     IF ( is_hubbard(nt) ) THEN  
        DO is = 1, nspin  
           DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = m1, 2*Hubbard_l(nt)+1
                 psum = ABS( ns(m1,m2,is,na) - ns(m2,m1,is,na) )  
                 IF (psum > 1.d-10) THEN  
                    WRITE( stdout, * ) na, is, m1, m2  
                    WRITE( stdout, * ) ns(m1,m2,is,na)  
                    WRITE( stdout, * ) ns(m2,m1,is,na)  
                    CALL errore( 'new_ns', 'non hermitean matrix', 1 )
                 ELSE  
                    ns(m1,m2,is,na) = 0.5d0 * (ns(m1,m2,is,na) + &
                                               ns(m2,m1,is,na) )
                    ns(m2,m1,is,na) = ns(m1,m2,is,na)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF 
  ENDDO
  !
  RETURN
  !
END SUBROUTINE new_ns

!--------------------------------------------------------------------
SUBROUTINE compute_pproj( ik, q, p )
    !------------------------------------------------------------------
    !! Here we compute DFT+U projections using the <beta|psi> overlaps.
    !
    USE kinds,                ONLY : DP
    USE ions_base,            ONLY : nat, ityp, ntyp => nsp
    USE klist,                ONLY : xk, igk_k, ngk
    USE becmod,               ONLY : becp
    USE uspp,                 ONLY : nkb, vkb, ofsbeta
    USE uspp_param,           ONLY : nhm, nh
    USE wvfct,                ONLY : nbnd
    USE wavefunctions,        ONLY : evc
    USE control_flags,        ONLY : gamma_only, offload_type
    USE ldaU,                 ONLY : is_hubbard, nwfcU
    USE becmod,               ONLY : bec_type, calbec, allocate_bec_type_acc, & 
                                     deallocate_bec_type_acc
    USE uspp_init,            ONLY : init_us_2
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik
    !! k point
    REAL(DP), INTENT(IN) :: q(nwfcU,nhm,nat)
    !! coefficients for projecting onto beta functions
    TYPE(bec_type), INTENT(INOUT) :: p
    !! proj(nwfcU,nbnd)
    !
    ! ... local variables
    !
    INTEGER :: ib, iw, nt, na, ih, npw, off
    !
    IF ( nkb == 0 ) RETURN
    !
    ! Number of plane waves at a given k point
    !
    npw = ngk(ik)
    !
    ! Compute <beta|psi>
    !
    CALL allocate_bec_type_acc( nkb, nbnd, becp )
    CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .true. )
    CALL calbec( offload_type, npw, vkb, evc, becp )
    ! does not need mp_sum intra-pool, since it is already done in calbec 
    !$acc data present (p,becp) copyin(q)
    IF ( gamma_only ) THEN 
       !$acc kernels present(p%r)
       p%r(:,:) = 0.0_DP
       !$acc end kernels
    ELSE
       !$acc kernels present(p%k)
       p%k(:,:) = (0.0_DP,0.0_DP)
       !$acc end kernels 
    ENDIF
    !
    DO nt = 1, ntyp
       DO na = 1, nat
          IF ( ityp(na) == nt ) THEN
             IF ( is_hubbard(nt) ) THEN
                off = ofsbeta(na)
                !$acc parallel loop collapse(2)
                DO ib = 1, nbnd
                   DO iw = 1, nwfcU
                      !$acc loop seq
                      DO ih = 1, nh(nt)
                         IF ( gamma_only ) THEN
                            p%r(iw,ib) = p%r(iw,ib) + q(iw,ih,na)*becp%r(off+ih,ib)
                         ELSE
                            p%k(iw,ib) = p%k(iw,ib) + q(iw,ih,na)*becp%k(off+ih,ib)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !$acc end data
    CALL deallocate_bec_type_acc( becp )
    !
    RETURN
    !
END SUBROUTINE compute_pproj
!
!---------------------------------------------------------------------------
SUBROUTINE new_ns_nc( ns )
  !-----------------------------------------------------------------------
  !! Noncollinear version (A. Smogunov). 
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, wfcU, &
                                   d_spin_ldau, is_hubbard, nwfcU, offsetU
  USE symm_base,            ONLY : d1, d2, d3
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE noncollin_module,     ONLY : npol
  USE symm_base,            ONLY : nsym, irt, time_reversal, t_rev, invs
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE gvect,                ONLY : gstart
  USE io_files,             ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
  USE buffers,              ONLY : get_buffer
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  IMPLICIT NONE
  !
  COMPLEX(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! the occupation numbers of orthogonalized atomic wfcs
  !
  ! ... local variables
  !
  INTEGER :: ik, ibnd, is, js, i, j, sigmay2, na, nb, nt, isym,  &
             m1, m2, m3, m4, is1, is2, is3, is4, m0, m00, ldim, npw
  COMPLEX(DP) , ALLOCATABLE :: nr(:,:,:,:,:), nr1(:,:,:,:,:), proj(:,:)
  COMPLEX(DP) :: z  
  REAL(DP) :: psum
  !
  ldim = 2 * Hubbard_lmax + 1
  !
  ALLOCATE( nr(ldim,ldim,npol,npol,nat), nr1(ldim,ldim,npol,npol,nat) )  
  ALLOCATE( proj(nwfcU,nbnd) )
  !
  nr(:,:,:,:,:)  = 0.d0
  nr1(:,:,:,:,:) = 0.d0
  ns(:,:,:,:)    = 0.d0
  !
  !--
  !  loop on k points
  DO ik = 1, nks
     !
     npw = ngk (ik)
     IF (nks > 1) THEN
        CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )
     ENDIF
     !
     ! make the projection - FIXME: use ZGEMM or calbec instead
     !
     DO ibnd = 1, nbnd
        DO i = 1, nwfcU
           proj(i,ibnd) = dot_product(wfcU(1:npwx*npol,i), evc(1:npwx*npol,ibnd))
        ENDDO
     ENDDO
     !
     CALL mp_sum( proj, intra_bgrp_comm )
     !
     ! compute the occupation matrix
     !
     DO na = 1, nat  
        nt = ityp (na)  
        IF ( is_hubbard(nt) ) THEN
           ldim = 2 * Hubbard_l(nt) + 1
           !
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 DO is1 = 1, npol
                    DO is2 = 1, npol
                       !
                       DO ibnd = 1, nbnd
                        nr(m1,m2,is1,is2,na) = nr(m1,m2,is1,is2,na) + &
                            wg(ibnd,ik) * CONJG( proj(offsetU(na)+m1+ldim*(is1-1),ibnd) ) * &
                                          proj(offsetU(na)+m2+ldim*(is2-1),ibnd)
                       ENDDO
                       !
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           !
        ENDIF
     ENDDO
     !
  ENDDO
  !
  !--
  !
  CALL mp_sum( nr, inter_pool_comm )
  !
  !--  symmetrize: nr  -->  nr1
  !
  DO na = 1, nat  
    nt = ityp(na)
    IF ( is_hubbard(nt) ) THEN
      !
      DO m1 = 1, 2*Hubbard_l(nt)+1  
        DO m2 = 1, 2*Hubbard_l(nt)+1  
          DO is1 = 1, npol
            DO is2 = 1, npol
              !
loopisym:     DO isym = 1, nsym  
                ! nr1(na) must be built from the atom that isym maps INTO na, i.e. the image of na
                ! under the inverse operation.  With nb = irt(isym,na) (the image of na) the rotation
                ! below is applied in the wrong orientation; harmless for operations equal to their
                ! inverse (E, i, C2, mirrors) but wrong for C3, C4, C6, S4, S6 whenever atoms are permuted.
                nb = irt (invs(isym), na)  
                !
                DO m3 = 1, 2*Hubbard_l(nt)+1
                  DO m4 = 1, 2*Hubbard_l(nt)+1
                    DO is3 = 1, npol
                      DO is4 = 1, npol
                        !
                        IF (Hubbard_l(nt) == 0) THEN
                          IF (t_rev(isym) == 1) THEN
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*                &
                                     nr(m4,m3,is4,is3,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  
                          ELSE
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*                &
                                     nr(m3,m4,is3,is4,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  
                          ENDIF
                        ELSEIF (Hubbard_l(nt) == 1) THEN
                          IF (t_rev(isym) == 1) THEN
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d1(m1,m3,isym)* &
                                     nr(m4,m3,is4,is3,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d1(m2,m4,isym)
                          ELSE
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d1(m1,m3,isym)* &
                                     nr(m3,m4,is3,is4,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d1(m2,m4,isym)
                          ENDIF
                        ELSEIF (Hubbard_l(nt) == 2) THEN
                          IF (t_rev(isym) == 1) THEN
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d2(m1,m3,isym)* &
                                     nr(m4,m3,is4,is3,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d2(m2,m4,isym)
                          ELSE
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d2(m1,m3,isym)* &
                                     nr(m3,m4,is3,is4,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d2(m2,m4,isym)
                          ENDIF
                        ELSEIF (Hubbard_l(nt) == 3) THEN
                          !
                          IF (t_rev(isym) == 1) THEN
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d3(m1,m3,isym)* &
                                     nr(m4,m3,is4,is3,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d3(m2,m4,isym)
                          ELSE
                            nr1(m1,m2,is1,is2,na) = nr1(m1,m2,is1,is2,na) +      &
                              CONJG( d_spin_ldau(is1,is3,isym) )*d3(m1,m3,isym)* &
                                     nr(m3,m4,is3,is4,nb)/nsym  *                &
                                     d_spin_ldau(is2,is4,isym)  *d3(m2,m4,isym)
                          ENDIF
                          !
                        ELSE
                          !
                          CALL errore( 'new_ns', &
                                       'angular momentum not implemented', &
                                       ABS(Hubbard_l(nt)) )
                          !
                        ENDIF
                        !
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
                !
              ENDDO  loopisym 
              !
            ENDDO
          ENDDO
        ENDDO
      ENDDO 
      !
    ENDIF
  ENDDO
  !--
  !
  !-- Setup the output matrix ns with combined spin index 
  !
  DO na = 1, nat
     nt = ityp (na)
     IF ( is_hubbard(nt) ) THEN
        DO is1 = 1, npol
          DO is2 = 1, npol
            i = npol*(is1-1) + is2
            DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                ns(m1,m2,i,na) = nr1(m1,m2,is1,is2,na)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
     ENDIF
  ENDDO
  !--
  !
  !-- make the matrix ns strictly hermitean
  !
  DO na = 1, nat  
     nt = ityp (na)  
     IF ( is_hubbard(nt) ) THEN  
        DO is1 = 1, npol  
          DO is2 = 1, npol
            i = npol*(is1-1) + is2
            j = is1 + npol*(is2-1)
            DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                 psum = ABS( ns(m1,m2,i,na) - CONJG(ns(m2,m1,j,na)) )  
                 IF (psum.GT.1.d-10) THEN  
                    WRITE( stdout, * ) na, m1, m2, is1, is2  
                    WRITE( stdout, * ) ns(m1,m2,i,na)  
                    WRITE( stdout, * ) ns(m2,m1,j,na)  
                    CALL errore( 'new_ns', 'non hermitean matrix', 1 )
                 ELSE  
                   ns(m2,m1,j,na) = CONJG( ns(m1,m2,i,na) )
                 ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
     ENDIF 
  ENDDO
  !--
  DEALLOCATE( nr, nr1 )
  !
  RETURN
  !
END SUBROUTINE new_ns_nc
