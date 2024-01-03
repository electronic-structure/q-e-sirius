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
  integer :: na, ig, nt, ir, ipert, is, ih, jh, ijh
  ! countera

  complex(DP), allocatable :: aux1 (:), aux2 (:,:), veff (:)
  ! work space
  complex(DP) z1(nspin_mag), z2

  if (.not.okvan) return
  !
  call start_clock ('newdq')
  !
  int3 (:,:,:,:,:) = (0.d0, 0.0d0)
  allocate (aux1 (ngm))
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
          ijh = 0
           do ih = 1, nh (nt)
              do jh = ih, nh (nt)
                 ijh = ijh + 1
                 !call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
                 do na = 1, nat
                    if (ityp (na) == nt) then
                       z1 = (0.d0, 0.d0)
!$omp parallel do default(shared) private(z2, is) reduction(+:z1)
                       do ig = 1, ngm
                          z2 = atom_type(nt)%qpw(ig, ijh) * &
                                eigts1(mill(1,ig),na) * &
                                eigts2(mill(2,ig),na) * &
                                eigts3(mill(3,ig),na) * &
                                eigqts(na)
                          do is = 1, nspin_mag
                             z1(is) = z1(is) + conjg(z2) * aux2(ig, is)
                          enddo
                       enddo !ig
!$omp end parallel do
                       do is = 1, nspin_mag
                          int3(ih,jh,na,is,ipert) = omega * z1(is)
                       enddo
                    endif
                 enddo !na
              enddo !jh
           enddo !ih
           do na = 1, nat
              if (ityp(na) == nt) then
                 !
                 !    We use the symmetry properties of the ps factor
                 !
                 do ih = 1, nh (nt)
                    do jh = ih, nh (nt)
                       do is = 1, nspin_mag
                          !    lower triangle            upper triangle   
                          int3(jh,ih,na,is,ipert) = int3(ih,jh,na,is,ipert)
                       enddo
                    enddo
                 enddo
              endif
           enddo

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
  deallocate (aux1)
  !
  call stop_clock ('newdq')
  !
  return
  !
end subroutine newdq
