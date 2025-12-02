!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusddens (drhop, dbecsum, mode0, npe)
  !----------------------------------------------------------------------
  !! This routine adds to the change of the charge and of the
  !! magnetization densities the part due to the US augmentation.
  !! It assumes that the array dbecsum has already accumulated the
  !! change of the becsum term. It calculates Eq. B31 of Ref [1].
  !! dbecsum and drhop contain the orthogonalization contribution to the
  !! change of the wavefunctions and the terms with alphasum and becsum are added.
  !! The contribution of the change of the Fermi energy is not calculated here
  !! but added later by ef_shift.
  !! [1] PRB 64, 235118 (2001).
  !
  !
  USE kinds, only : DP
  use fft_base,  only: dfftp
  use fft_interfaces, only: invfft
  USE gvect,  ONLY : ngm, g, eigts1, eigts2, eigts3, mill
  USE uspp,     ONLY : okvan, becsum
  USE cell_base, ONLY : tpiba
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE wavefunctions,  ONLY: psic
  USE uspp_param, ONLY: upf, lmaxq, nh, nhm
  USE paw_variables, ONLY : okpaw
  USE modes,     ONLY : u
  USE phus,    ONLY : becsumort, alphasum
  USE noncollin_module, ONLY : nspin_mag
  USE qpoint,     ONLY : xq, eigqts

  implicit none
  !
  !   the dummy variables
  !

  integer :: npe
  !! input: the number of perturbations
  complex(DP) :: drhop (dfftp%nnr, nspin_mag, npe)
  !! inp/out: change of the charge density
  complex(DP) :: dbecsum (nhm*(nhm+1)/2, nat, nspin_mag, npe)
  !! input: sum over kv of bec
  integer ::  mode0
  !! input:the mode of the representation
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, mu, mode, ipert, is, ijh
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on r vectors
  ! pointer on modes
  ! pointer on the mode
  ! counter on perturbations
  ! counter on spin
  ! counter on combined beta functions

  real(DP), allocatable  :: qmod (:), qpg (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(DP) :: fact, zsum, bb, alpha, u1, u2, u3
  ! auxiliary variables
  complex(DP), allocatable ::  sk (:), qgm(:), aux (:,:,:)
  ! the structure factor
  ! q_lm(G)
  ! auxiliary variable for drho(G)
  !
  IF (.not.okvan) return
  !
  CALL start_clock ('addusddens')
  !
  allocate (aux(  ngm , nspin_mag , npe))
  allocate (sk (  ngm))
  allocate (ylmk0(ngm , lmaxq * lmaxq))
  allocate (qgm(  ngm))
  allocate (qmod( ngm))
  allocate (qpg( 3  , ngm))
  !      WRITE( stdout,*) aux, ylmk0, qmod
  !
  !  And then we compute the additional charge in reciprocal space
  !
  call setqmod (ngm, xq, g, qmod, qpg)
  call ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (qmod (ig) ) * tpiba
  enddo
  !
  fact = cmplx (0.d0, - tpiba, kind=DP)
  aux(:,:,:) = (0.d0, 0.d0)
  do nt = 1, ntyp
     if (upf(nt)%tvanp  ) then
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              ijh = ijh + 1
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    mu = 3 * (na - 1)
                    !
                    ! calculate the structure factor
                    !
                    do ig = 1, ngm
                       sk (ig) = eigts1 (mill(1,ig), na) * &
                                 eigts2 (mill(2,ig), na) * &
                                 eigts3 (mill(3,ig), na) * &
                                 eigqts (na) * qgm (ig)
                    enddo
                    !
                    !  And qgmq and becp and dbecq
                    !
                    do ipert = 1, npe
                       do is = 1, nspin_mag
                          mode = mode0 + ipert
                          zsum = dbecsum (ijh, na, is, ipert)
                          !
                          u1 = u (mu + 1, mode)
                          u2 = u (mu + 2, mode)
                          u3 = u (mu + 3, mode)
                          if (abs(u1) + abs(u2) + abs(u3) > 1d-12) then
                             bb = becsum (ijh, na, is)
                             zsum = zsum + &
                                  ( alphasum (ijh, 1, na, is) * u1 &
                                  + alphasum (ijh, 2, na, is) * u2 &
                                  + alphasum (ijh, 3, na, is) * u3)
                             !
                             do ig = 1, ngm
                                alpha = qpg(1,ig)*u1 + qpg(2,ig)*u2 + qpg(3,ig)*u3
                                aux(ig,is,ipert) = aux(ig,is,ipert) + fact * alpha * bb * sk(ig)
                             enddo
                          endif
                          call zaxpy (ngm, zsum, sk, 1, aux(1,is,ipert), 1)
                          IF (okpaw) becsumort(ijh,na,is,mode) = zsum
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  !     convert aux to real space
  !
  do ipert = 1, npe
     do is = 1, nspin_mag
        psic(:) = (0.d0, 0.d0)
        do ig = 1, ngm
           psic (dfftp%nl (ig) ) = aux (ig, is, ipert)
        enddo
        CALL invfft('Rho', psic, dfftp)
        call zaxpy(dfftp%nnr, (1.d0, 0.d0), psic, 1, drhop(1,is,ipert), 1)
     enddo
  enddo
  deallocate (qpg)
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (sk)
  deallocate (aux)
  !
  call stop_clock ('addusddens')
  !
end subroutine addusddens
