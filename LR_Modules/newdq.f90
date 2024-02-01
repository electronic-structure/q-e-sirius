!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine newdq (dvscf, npe)
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system and adds it to dvpsi.
  !
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, mill, eigts1, eigts2, eigts3
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE paw_variables,        ONLY : okpaw
  USE mp_bands,             ONLY: intra_bgrp_comm
  USE mp,                   ONLY: mp_sum
  USE lrus,                 ONLY : int3, int3_paw
  USE qpoint,               ONLY : xq, eigqts
  USE control_lr,           ONLY : lgamma
  USE mod_sirius

  implicit none
  !
  !   The dummy variables
  !
  integer, intent(in) :: npe
  ! input: the number of perturbations

  complex(DP), intent(in) :: dvscf (dfftp%nnr, nspin_mag, npe)
  ! input: the change of the selfconsistent pot.
  !
  !   And the local variables
  !
  integer :: na, ig, nt, ir, ipert, is, ih, jh, ijh, nij, N_nt, na_
  ! countera

  complex(DP), allocatable :: aux2 (:,:), veff (:), tmp(:,:), res2(:,:)
  ! work space

  if (.not.okvan) return
  !
  call start_clock ('newdq')
  !
  int3 (:,:,:,:,:) = (0.d0, 0.0d0)
  allocate (aux2 (ngm , nspin_mag))
  allocate (veff (dfftp%nnr))
  !
  !     and for each perturbation of this irreducible representation
  !     integrate the change of the self consistent potential and
  !     the Q functions
  !
  do ipert = 1, npe

     do is = 1, nspin_mag
        do ir = 1, dfftp%nnr
           veff (ir) = dvscf (ir, is, ipert)
        enddo
        CALL fwfft ('Rho', veff, dfftp)
        do ig = 1, ngm
           aux2 (ig, is) = veff (dfftp%nl (ig) )
        enddo
     enddo

     do nt = 1, ntyp ! loop over atom types
        if (upf(nt)%tvanp ) then
           ! composite index for ih and jh (ksi and ksi')
           nij = nh(nt)*(nh(nt)+1)/2 ! max number of (ih,jh) pairs per atom type nt
           N_nt = 0 ! number of atoms of type nt
           DO na = 1, nat
              IF ( ityp(na) == nt ) N_nt = N_nt + 1
           ENDDO

           allocate (tmp(ngm, N_nt))
           allocate (res2(nij, N_nt))

           do is = 1, nspin_mag ! loop over spins
               na_ = 0 ! count atoms of type nt
               !=============== compute potential (aux2) * phase factors
               call start_clock ('aux2_x_phases')
               do na = 1, nat ! loop over all atoms
                  if (ityp(na) == nt) then
                     na_ = na_ + 1
                     !$omp parallel do default(shared) private(ig)
                     do ig = 1, ngm ! loop over G-vectors
                         tmp(ig, na_) = aux2(ig, is) * CONJG( eigts1(mill(1,ig),na) * &
                                                              eigts2(mill(2,ig),na) * &
                                                              eigts3(mill(3,ig),na) * &
                                                              eigqts(na) )
                     enddo
                     !$omp end parallel do
                  endif
               enddo
               call stop_clock ('aux2_x_phases')
               !=============== compute Q*V for all atoms of type nt
               call start_clock ('newdq_ZGEMM')
               ! qpw is a complex array of dimension (ngm, nij)
               ! tmp is a complex array of dimension (ngm, N_nt)
               ! --- (qpw)^H * tmp : (nij,  ngm) x (ngm, N_nt) = (nij, N_nt), dimensions of res2

               call ZGEMM('C', 'N', nij, N_nt, ngm, dcmplx(1.d0, 0.d0), atom_type(nt)%qpw, &
                          ngm, tmp, ngm, dcmplx(0.d0, 0.d0), res2, nij)
               call stop_clock ('newdq_ZGEMM')

               na_ = 0
               do na = 1, nat ! loop over all atoms
                  if (ityp(na) == nt) then
                     na_ = na_ + 1
                     ijh = 0
                     do ih = 1, nh(nt) ! loop over ksi
                        do jh = ih, nh(nt) ! loop over ksi'
                           ijh = ijh + 1
                           int3(ih,jh,na,is,ipert) = omega * res2(ijh, na_)
                           !                 lower triangle            upper triangle
                           IF (jh > ih) int3(jh,ih,na,is,ipert) = int3(ih,jh,na,is,ipert)
                        enddo
                     enddo
                  endif
               enddo
           enddo ! loop over spins
           !
           deallocate (tmp)
           deallocate (res2)
        endif ! if US-PP
     enddo ! nt
  enddo ! ipert
#if defined(__MPI)
  call mp_sum ( int3, intra_bgrp_comm )
#endif
  !
  IF (noncolin) CALL set_int3_nc(npe)
  !
  ! Sum of the USPP and PAW terms 
  ! (see last two terms in Eq.(12) in PRB 81, 075123 (2010))
  !
  IF (okpaw) int3 = int3 + int3_paw
  !
  deallocate (veff)
  deallocate (aux2)
  !
  call stop_clock ('newdq')
  !
  return
  !
end subroutine newdq
