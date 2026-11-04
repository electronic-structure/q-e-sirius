!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_psym_diel(dvtosym, dfft)
  !-----------------------------------------------------------------------
  !! Symmetrizes the induced density/potential of the 3 Cartesian
  !! polarizations of the q=Gamma electric-field Lanczos chain
  !! (crpa_do_lanczos_epsilon.f90), reproducing exactly what
  !! LR_Modules/psymdvscf.f90 + symdvscf.f90 compute for this case
  !! (lr_npert=3, q=Gamma, no LSDA/noncolinear magnetism), but distributing
  !! the O(grid_points * nsymq) rotation-and-accumulation loop across MPI
  !! ranks - each rank symmetrizes only the points of its own local FFT
  !! slab - instead of redundantly repeating the full, non-distributed grid
  !! on every rank.
  !!
  !! This matters because symdvscf's per-call cost does NOT decrease with
  !! more MPI processes (confirmed empirically: ~33ms/call regardless of
  !! whether it runs on 1 or 10 processes for this test system), since every
  !! rank does the same full-grid computation independently. That is
  !! negligible for the ~10-30 calls of the self-consistent CG dielectric
  !! solve (LR_Modules/dfpt_kernels.f90), which is why it was never an
  !! issue there, but it dominates (>80% of the total runtime) when called
  !! thousands of times by the Lanczos chain driver.
  !!
  !! NOTE: an earlier version of this routine reused
  !! LR_Modules/lr_sym_mod.f90's rotate_mesh/rotate_mesh_1s/find_mesh_ijk to
  !! build the per-rank-local rotation table (mirroring how thermo_pw's own
  !! Lanczos implementation, thermo_pw/qe/lr_sym_mod_tpw.f90, avoids the same
  !! cost). That was WRONG: those 3 routines are hardcoded to
  !! "USE fft_base, ONLY : dfftp" internally instead of taking the FFT
  !! descriptor as an argument, so when called here for the smooth grid
  !! (dffts), they silently used dfftp's parallel distribution instead of
  !! dffts's - invisible with 1 MPI process (nothing to distribute
  !! differently), but silently wrong with more than one (confirmed: the
  !! Lanczos chain was measurably less stable, though not as catastrophically
  !! as with no symmetrization at all). This version instead builds the
  !! rotation table locally from scale_sym_ops/rotate_grid_point
  !! (PW/src/ruotaijk.f90), which - like symdvscf.f90 itself - take the grid
  !! dimensions as plain explicit arguments, so they are correct for either
  !! grid regardless of process count.
  !
  USE kinds,             ONLY : DP
  USE fft_base,          ONLY : fft_type_descriptor
  USE symm_base,         ONLY : s, ft, t_rev
  USE noncollin_module,  ONLY : nspin_lsda, nspin_mag
  USE control_lr,        ONLY : lgamma
  USE lr_symm_base,      ONLY : minus_q, nsymq, gi, lr_npert, upert
  USE lr_sym_mod,        ONLY : compute_phase
  USE scatter_mod,       ONLY : cgather_sym
  USE mp_pools,          ONLY : npool, my_pool_id, inter_pool_comm
  USE mp,                ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  !! FFT descriptor for dvtosym (dfftp for the hard grid, dffts for the smooth grid)
  COMPLEX(DP), INTENT(INOUT) :: dvtosym(dfft%nnr, nspin_mag, lr_npert)
  !! the potential/density to symmetrize, one slice per Cartesian polarization
  !
  INTEGER :: is, ir, i, j, k, ipert, jpert, isym, my_nrxx
  INTEGER :: nr1, nr2, nr3, nr1x, my_nr2p, my_nr3p, j0, k0, idx
  INTEGER :: ri, rj, rk
  INTEGER :: s_scaled(3,3,nsymq), ftau(3,nsymq)
  INTEGER, ALLOCATABLE :: rir(:,:), iir(:), jir(:), kir(:)
  COMPLEX(DP), ALLOCATABLE :: aux(:,:), dvsym(:,:), add_dvsym(:)
  COMPLEX(DP), ALLOCATABLE :: phase1(:,:), phase2(:,:), phase3(:,:)
  LOGICAL, ALLOCATABLE :: zb(:)
  COMPLEX(DP) :: phase
  !
  IF (nsymq == 1 .AND. (.NOT. minus_q)) RETURN
  !
  IF (minus_q .AND. .NOT. lgamma) CALL errore('crpa_psym_diel', &
     'minus_q with q /= Gamma is not supported: the electric-field &
     &Lanczos chain is only ever built at q = Gamma', 1)
  !
  CALL start_clock('crpa_psym_diel')
  !
  nr1 = dfft%nr1;  nr2 = dfft%nr2;  nr3 = dfft%nr3
  nr1x = dfft%nr1x
  my_nr2p = dfft%my_nr2p;  my_nr3p = dfft%my_nr3p
  j0 = dfft%my_i0r2p;  k0 = dfft%my_i0r3p
  my_nrxx = nr1x * my_nr2p * my_nr3p
  !
  CALL scale_sym_ops(nsymq, s, ft, nr1, nr2, nr3, s_scaled, ftau)
  !
  ALLOCATE(rir(my_nrxx, nsymq), iir(my_nrxx), jir(my_nrxx), kir(my_nrxx))
  ! only this pool's stride of isym is filled below; zero the rest so the
  ! untouched entries are never read as uninitialized memory
  rir(:,:) = 0
  ALLOCATE(aux(dfft%nr1x*dfft%nr2x*dfft%nr3x, lr_npert))
  ALLOCATE(dvsym(my_nrxx, lr_npert), add_dvsym(lr_npert))
  ALLOCATE(phase1(nr1, nsymq), phase2(nr2, nsymq), phase3(nr3, nsymq))
  ALLOCATE(zb(nsymq))
  !
  ! Build, once, the rotated-point table and (i,j,k) coordinates of this
  ! rank's own local slab (mirroring rotate_mesh_1s/find_mesh_ijk, but
  ! parameterized by the actual dfft passed in, not hardcoded to dfftp).
  !
  DO ir = 1, my_nrxx
     idx = ir - 1
     k   = idx / (nr1x*my_nr2p)
     idx = idx - (nr1x*my_nr2p)*k
     k   = k + k0
     j   = idx / nr1x
     idx = idx - nr1x*j
     j   = j + j0
     i   = idx
     i = i + 1;  j = j + 1;  k = k + 1
     iir(ir) = i;  jir(ir) = j;  kir(ir) = k
     IF (i > nr1 .OR. j > nr2 .OR. k > nr3) THEN
        rir(ir,:) = 0
        CYCLE
     ENDIF
     DO isym = my_pool_id + 1, nsymq, npool
        CALL rotate_grid_point(s_scaled(1,1,isym), ftau(1,isym), i, j, k, &
             nr1, nr2, nr3, ri, rj, rk)
        rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*dfft%nr2x
     ENDDO
  ENDDO
  !
  ! Phase e^{i*gi.r}, needed only if some symmetry in the small group of q
  ! maps q=Gamma to q+G with G /= 0 (not expected at q=Gamma, but computed
  ! the same way symdvscf.f90 does rather than assumed away). compute_phase
  ! also sets zb(isym), telling which symmetries need this phase at all.
  !
  CALL compute_phase(phase1, phase2, phase3, nr1, nr2, nr3, nsymq, gi, zb)
  !
  DO is = 1, nspin_lsda
     !
     ! Collect the quantity to symmetrize on every processor (as
     ! psymdvscf.f90 also does), then let each rank symmetrize only the
     ! points of its own local slab.
     !
     DO ipert = 1, lr_npert
        CALL cgather_sym(dfft, dvtosym(:,is,ipert), aux(:,ipert))
     ENDDO
     !
     ! q = Gamma special case of the minus_q correction: S(irotmq) = I,
     ! Gi = 0 (see symdvscf.f90), so it reduces to forcing the result real.
     !
     IF (minus_q) aux(:,:) = CMPLX(DBLE(aux(:,:)), 0.0_DP, KIND=DP)
     !
     dvsym(:,:) = (0.0_DP, 0.0_DP)
     !
     ! The symmetry sum is split over pools: pool p takes isym = p+1, p+1+npool,
     ! ... and the partial sums are combined by the mp_sum below. Without this
     ! every pool redundantly evaluated the FULL sum over its own copy of aux.
     ! This is valid because the FFT descriptor is distributed within a pool and
     ! is therefore identical in every pool.
     DO isym = my_pool_id + 1, nsymq, npool
        DO ir = 1, my_nrxx
           IF (rir(ir,isym) == 0) CYCLE
           IF (zb(isym)) THEN
              i = iir(ir); j = jir(ir); k = kir(ir)
              phase = phase1(i,isym) * phase2(j,isym) * phase3(k,isym)
           ELSE
              phase = (1.0_DP, 0.0_DP)
           ENDIF
           add_dvsym(:) = (0.0_DP, 0.0_DP)
           DO ipert = 1, lr_npert
              DO jpert = 1, lr_npert
                 add_dvsym(ipert) = add_dvsym(ipert) + upert(jpert,ipert,isym) * &
                      aux(rir(ir,isym), jpert) * phase
              ENDDO
           ENDDO
           IF (t_rev(isym) == 0) THEN
              dvsym(ir,:) = dvsym(ir,:) + add_dvsym(:)
           ELSE
              dvsym(ir,:) = dvsym(ir,:) + CONJG(add_dvsym(:))
           ENDIF
        ENDDO
     ENDDO
     !
     IF (npool > 1) CALL mp_sum(dvsym, inter_pool_comm)
     !
     DO ipert = 1, lr_npert
        dvtosym(1:my_nrxx, is, ipert) = dvsym(:,ipert) / DBLE(nsymq)
     ENDDO
     !
  ENDDO
  !
  DEALLOCATE(rir, iir, jir, kir, aux, dvsym, add_dvsym, phase1, phase2, phase3, zb)
  !
  CALL stop_clock('crpa_psym_diel')
  !
  RETURN
  !
END SUBROUTINE crpa_psym_diel
