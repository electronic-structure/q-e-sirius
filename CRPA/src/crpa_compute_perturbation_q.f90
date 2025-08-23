!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
SUBROUTINE crpa_compute_perturbation_q(ipert, drho, dfpt_data)
   !----------------------------------------------------------------------
   !! Compute the bare charge density perturbation for the ipert-th perturbation
   !! and write it at drho.
   !! For a metal at q=0, compute bare charge perturbation and write to dfpt_data%dn0.
   !----------------------------------------------------------------------
   !
   USE kinds,                ONLY : DP
   USE mp,                   ONLY : mp_sum
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE io_global,            ONLY : stdout
   USE fft_base,             ONLY : dfftp
   USE cell_base,            ONLY : omega
   USE noncollin_module,     ONLY : nspin_mag
   USE klist,                ONLY : ltetra, lgauss
   USE control_lr,           ONLY : lgamma
   USE dfpt_type,            ONLY : dfpt_data_type
   USE crpa_pert,            ONLY : pert_basis
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: ipert
   !! Perturbation index
   COMPLEX(DP), INTENT(OUT) :: drho(dfftp%nnr, nspin_mag)
   !! Wavefunction pair density in real space
   TYPE(dfpt_data_type) :: dfpt_data
   !! Data that describes linear response quantities
   !
   COMPLEX(DP) :: drho_sum
   !! Sum of drho
   !
   CALL start_clock ('crpa_perturbation_q')
   !
   IF (dfpt_data%npert /= 1) CALL errore('crpa_compute_perturbation_q', 'npert is not 1', 1)
   !
   IF (pert_basis == 'bands') THEN
      CALL crpa_compute_perturbation_q_bands(ipert, drho)
      !
   ELSEIF (pert_basis == 'wannier') THEN
      CALL crpa_compute_perturbation_q_wannier(ipert, drho)
      !
   ELSE
      CALL errore('crpa_compute_perturbation_q', 'Wrong pert_basis', 1)
   ENDIF
   !
   ! Print the sum of the perturbation for sanity check
   drho_sum = SUM(drho)
   CALL mp_sum(drho_sum, intra_bgrp_comm)
   drho_sum = drho_sum * omega / DBLE(dfftp%nr1) / DBLE(dfftp%nr2) / DBLE(dfftp%nr3)
   WRITE(stdout, '(5x,A,2ES20.10)') 'Sum of drho (normalized to 1) = ', drho_sum
   !
   ! Set delta_n_ext (monopole charge of external perturbation).
   ! Needed for the calculation of the Fermi energy shift for q=0 with metals.
   !
   IF (ALLOCATED(dfpt_data%dn0)) THEN
      dfpt_data%dn0(1) = drho_sum
   ENDIF
   !
   CALL stop_clock ('crpa_perturbation_q')
   !
END SUBROUTINE crpa_compute_perturbation_q
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
SUBROUTINE crpa_compute_perturbation_q_bands(ipert, drho)
   !----------------------------------------------------------------------
   !
   ! This routine computes the density that corresponds to the pair of wavefunctions.
   !
   ! On output: drho(r) = psi_{ibnd, k+q}(r) * psi_{jbnd, k}(r)^*
   !
   USE kinds,                ONLY : DP
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE io_files,             ONLY : nwordwfcU
   USE io_global,            ONLY : stdout
   USE ions_base,            ONLY : nat, ityp
   USE cell_base,            ONLY : omega
   USE wavefunctions,        ONLY : evc
   USE klist,                ONLY : ngk, igk_k
   USE buffers,              ONLY : get_buffer
   USE wvfct,                ONLY : npwx, nbnd
   USE gvecs,                ONLY : doublegrid
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_wave,             ONLY : invfft_wave
   USE fft_interfaces,       ONLY : fft_interpolate
   USE mp_pools,             ONLY : intra_pool_comm
   USE mp,                   ONLY : mp_sum
   USE eqv,                  ONLY : dvpsi
   USE qpoint,               ONLY : nksq, ikks, ikqs, ikmks, ikmkmqs
   USE units_lr,             ONLY : iuatswfc, lrwfc, iuwfc
   USE control_lr,           ONLY : lgamma
   USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, offsetU, nwfcU
   USE ldaU_lr,              ONLY : swfcatomk, swfcatomkpq
   USE noncollin_module,     ONLY : noncolin, npol, domag, nspin_mag
   USE eqv,                  ONLY : evq
   USE crpa_pert,            ONLY : pert_nbnd, pert_ibnd_min
   USE mp_global, ONLY : npool
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: ipert
   !! Perturbation index
   COMPLEX(DP), INTENT(INOUT) :: drho(dfftp%nnr, nspin_mag)
   !! Wavefunction pair density in real space : pert_nbnd, pert_ibnd_min
   !
   INTEGER :: ik, ikk, ikq, na, nt, m, ibnd, jbnd, ig, ldim, ir, ipol, is
   !! Counters
   COMPLEX (DP), ALLOCATABLE :: evc_k_r(:, :)
   !! Wavefunction in real space
   COMPLEX (DP), ALLOCATABLE :: evc_kq_r(:, :)
   !! Wavefunction in real space
   COMPLEX (DP), ALLOCATABLE :: drhos(:, :)
   !! Wavefunction in real space in dffts grid
   !
   IF (npool > 1) CALL errore("crpa_compute_perturbation_q_wannier", "npool > 1 not implemented", 1)
   !
   ALLOCATE(evc_k_r(dffts%nnr, npol))
   ALLOCATE(evc_kq_r(dffts%nnr, npol))
   ALLOCATE(drhos(dffts%nnr, nspin_mag))
   !
   drho = (0.d0, 0.d0)
   drhos = (0.d0, 0.d0)
   !
   ! Compute the indices of the pair of wavefunctions to be perturbed
   !
   ! Perturbation by band pairs: ipert = (ibnd, jbnd, ik)
   ibnd = MOD(ipert - 1, pert_nbnd) + pert_ibnd_min
   jbnd = MOD((ipert - 1) / pert_nbnd, pert_nbnd) + pert_ibnd_min
   ik = (ipert - 1) / (pert_nbnd * pert_nbnd) + 1
   !
   ! FIXME: Pool parallelization. Would need some communication between pools.
   !
   ! Compute the wavefunctions in real space
   !
   ikk = ikks(ik) ! points to k+G indices
   ikq = ikqs(ik) ! points to k+q+G indices
   !
   CALL get_buffer(evc, lrwfc, iuwfc, ikk)
   CALL get_buffer(evq, lrwfc, iuwfc, ikq)
   !
   CALL invfft_wave(npwx, ngk(ikk), igk_k(1, ikk), evc(:, jbnd), evc_k_r)
   CALL invfft_wave(npwx, ngk(ikq), igk_k(1, ikq), evq(:, ibnd), evc_kq_r)
   !
   DO ipol = 1, npol
      DO ir = 1, dffts%nnr
         drhos(ir, 1) = drhos(ir, 1) + evc_kq_r(ir, ipol) * CONJG(evc_k_r(ir, ipol))
      ENDDO
   ENDDO
   !
   IF (doublegrid) THEN
      DO is = 1, nspin_mag
         CALL fft_interpolate(dffts, drhos(:, is), dfftp, drho(:, is))
      ENDDO
   ELSE
      CALL zcopy(dffts%nnr * nspin_mag, drhos(1, 1), 1, drho(1, 1), 1)
   ENDIF
   !
   DEALLOCATE(evc_k_r)
   DEALLOCATE(evc_kq_r)
   DEALLOCATE(drhos)
   !
   ! Normalization
   !
   drho = drho / omega
   !
END SUBROUTINE crpa_compute_perturbation_q_bands
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
SUBROUTINE crpa_compute_perturbation_q_wannier(ipert, drho)
   !----------------------------------------------------------------------
   !
   ! This routine computes the density that corresponds to the pair of Wannier functions.
   !
   ! On output: drho(r) = 1/N_k \sum_k w_{iw, k+q}(r) * w_{jw, k)(r)^* * exp(iqr) * exp(ikR)
   !
   ! w_{iw, k+q} = \sum_{ibnd} psi_{ibnd, k+q} U_{ibnd, iw}(k+q)
   ! w_{jw, k}   = \sum_{jbnd} psi_{jbnd, k}   U_{jbnd, jw}(k)
   !
   USE kinds,                ONLY : DP
   USE constants,            ONLY : tpi
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE io_files,             ONLY : nwordwfcU
   USE io_global,            ONLY : stdout
   USE ions_base,            ONLY : nat, ityp
   USE cell_base,            ONLY : omega, at
   USE wavefunctions,        ONLY : evc
   USE klist,                ONLY : ngk, igk_k, xk
   USE buffers,              ONLY : save_buffer, get_buffer
   USE wvfct,                ONLY : npwx, nbnd
   USE gvecs,                ONLY : doublegrid
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_wave,             ONLY : invfft_wave
   USE fft_interfaces,       ONLY : fft_interpolate
   USE mp_pools,             ONLY : intra_pool_comm, inter_pool_comm
   USE mp,                   ONLY : mp_sum
   USE eqv,                  ONLY : dvpsi
   USE qpoint,               ONLY : nksq, nksqtot, ikks, ikqs, ikmks, ikmkmqs
   USE control_lr,           ONLY : lgamma
   USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, offsetU, nwfcU
   USE ldaU_lr,              ONLY : swfcatomk, swfcatomkpq
   USE noncollin_module,     ONLY : noncolin, npol, domag, nspin_mag
   USE eqv,                  ONLY : evq
   USE crpacom,              ONLY : iuwan, lrwan
   USE w90_interface,        ONLY : num_wann, u_matrix
   USE crpa_pert,            ONLY : pert_iwlist, pert_jwlist, pert_Rlist
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: ipert
   !! Perturbation index
   COMPLEX(DP), INTENT(INOUT) :: drho(dfftp%nnr, nspin_mag)
   !! Wavefunction pair density in real space : pert_nbnd, pert_ibnd_min
   !
   INTEGER :: ik, ikk, ikq, ipol, iR, is, iw, jw
   !! Counters
   REAL(DP) :: arg
   !! k dot r
   REAL(DP) :: xkk(3)
   !! k point in Cartesian coordinates
   REAL(DP) :: R(3)
   !! lattice vector in Cartesian coordinates
   COMPLEX(DP) :: phase
   !! exp(ikR) phase
   COMPLEX (DP), ALLOCATABLE :: wan_k_G(:)
   !! Wannier function in G vector basis
   COMPLEX (DP), ALLOCATABLE :: wan_kq_G(:)
   !! Wannier function in G vector basis
   COMPLEX (DP), ALLOCATABLE :: wan_k_r(:, :)
   !! Wannier function in real space
   COMPLEX (DP), ALLOCATABLE :: wan_kq_r(:, :)
   !! Wannier function in real space
   COMPLEX (DP), ALLOCATABLE :: drhos(:, :)
   !! Wavefunction in real space in dffts grid
   !
   ALLOCATE(wan_k_G(npwx*npol))
   ALLOCATE(wan_kq_G(npwx*npol))
   ALLOCATE(wan_k_r(dffts%nnr, npol))
   ALLOCATE(wan_kq_r(dffts%nnr, npol))
   ALLOCATE(drhos(dffts%nnr, nspin_mag))
   !
   drho = (0.d0, 0.d0)
   drhos = (0.d0, 0.d0)
   !
   ! Compute the indices of the pair of wavefunctions to be perturbed
   !
   ! Perturbation by band pairs: ipert = (iw, jw, iR)
   iw = pert_iwlist(ipert)
   jw = pert_jwlist(ipert)
   R = REAL(pert_Rlist(:, ipert), DP)
   CALL cryst_to_cart(1, R, at, +1)
   !
   WRITE(stdout, '(5x, A, I6, A, I6, A, 3I6)') "iw = ", iw, " , jw = ", jw, &
                                               " , R (crystal) = ", pert_Rlist(:, ipert)
   !
   DO ik = 1, nksq
      !
      ikk = ikks(ik)
      ikq = ikqs(ik)
      xkk = xk(:, ikk)
      !
      arg = SUM(xkk * R) * tpi
      phase = CMPLX( COS(arg), SIN(arg), KIND = DP)
      !
      ! Read Wannier function from buffer
      !
      CALL get_buffer(evc, lrwan, iuwan, ikk)
      CALL get_buffer(evq, lrwan, iuwan, ikq)
      wan_k_G(:)  = evc(:, jw)
      wan_kq_G(:) = evq(:, iw)
      !
      ! Fourier transform the Wannier functions to real space
      !
      CALL invfft_wave(npwx, ngk(ikk), igk_k(1, ikk), wan_k_G, wan_k_r)
      CALL invfft_wave(npwx, ngk(ikq), igk_k(1, ikq), wan_kq_G, wan_kq_r)
      !
      ! Compute the density in real space
      !
      DO ipol = 1, npol
         DO ir = 1, dffts%nnr
            drhos(ir, 1) = drhos(ir, 1) + wan_kq_r(ir, ipol) * CONJG(wan_k_r(ir, ipol)) * phase
         ENDDO
      ENDDO
   ENDDO
   !
   CALL mp_sum(drhos, inter_pool_comm)
   !
   ! Normalization due to the sum over k points
   drhos = drhos / nksqtot
   !
   IF (doublegrid) THEN
      DO is = 1, nspin_mag
         CALL fft_interpolate(dffts, drhos(:, is), dfftp, drho(:, is))
      ENDDO
   ELSE
      CALL zcopy(dffts%nnr * nspin_mag, drhos(1, 1), 1, drho(1, 1), 1)
   ENDIF
   !
   DEALLOCATE(wan_k_G)
   DEALLOCATE(wan_kq_G)
   DEALLOCATE(wan_k_r)
   DEALLOCATE(wan_kq_r)
   DEALLOCATE(drhos)
   !
   ! Normalization
   !
   drho = drho / omega
   !
END SUBROUTINE crpa_compute_perturbation_q_wannier
