!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------
SUBROUTINE crpa_skip_nscf_by_rotate()
!------------------------------------------------------------
   !
   ! Perform rotation and unfolding to setup eigenvalues and eigenstates at k and k+q
   ! from the eigenvalues at the irreducible BZ of the k-point grid.
   ! Used to skip any NSCF calculation in the DFPT calculation.
   ! Errors of the desired k or k+q point is not present in the unfolded k-point grid.
   !
   USE kinds,            ONLY : DP
   USE io_global,        ONLY : stdout
   USE mp,               ONLY : mp_sum
   USE mp_bands,         ONLY : intra_bgrp_comm
   USE klist,            ONLY : nks, xk, ngk, igk_k
   USE lsda_mod,         ONLY : lsda, isk, current_spin
   USE becmod,           ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, becp
   USE fft_base,         ONLY : dffts
   USE fft_interfaces,   ONLY : fwfft, invfft
   USE control_flags,    ONLY : gamma_only
   USE io_files,         ONLY : prefix, postfix
   USE gvect,            ONLY : g, gg
   USE cell_base,        ONLY : at, bg
   USE wvfct,            ONLY : nbnd, npwx, et, current_k
   USE wavefunctions,    ONLY : evc, psic
   USE noncollin_module, ONLY : npol
   USE units_lr,         ONLY : lrwfc, iuwfc
   USE control_lr,       ONLY : lgamma
   USE eqv,              ONLY : evq
   USE buffers,          ONLY : get_buffer, save_buffer
   USE uspp,             ONLY : nkb, vkb
   USE uspp_init,        ONLY : init_us_2
   USE crpacom,          ONLY : ik_to_ik_orig, tmp_dir_save, nks_orig, xk_orig
   !
   IMPLICIT NONE
   !
   LOGICAL :: found
   !! Flag to indicate if k' is found
   CHARACTER(LEN = 256) :: wfcdir
   !! Directory of the collecteed wavefunction files from the SCF calculation
   INTEGER :: ik, ik_orig
   !! k point index
   INTEGER :: ig, istart, iend, ig_k
   !! G vector index
   INTEGER :: ib
   !! Band index
   INTEGER :: npw
   !! Number of G vectors
   INTEGER :: ibnd
   !! Band index
   INTEGER :: ipol
   !! Spin polarization index
   INTEGER :: ierr
   !! Error number
   REAL(DP) :: gg_
   !! Square of the norm of G vector
   REAL(DP) :: xkk(3), xk_diff(3)
   !! k point vector
   REAL(DP) :: g_k(3)
   !! G vector such that xkk = k' + G. In Cartesian coordinates.
   COMPLEX(DP), ALLOCATABLE :: phase(:)
   !! Temporary storage for phase e^{-i * G_k * r}
   COMPLEX(DP), ALLOCATABLE :: h_evc(:, :)
   !! Temporary storage for H * evc
   !
   COMPLEX(DP) :: overlap(nbnd, nbnd) ! DEBUG
   CHARACTER(len=6), EXTERNAL :: int_to_char
   !
   CALL start_clock('crpa_nscf_rotate')
   !
   IF (lgamma) THEN
      !
      ! For q = 0, no rotation is needed. Set ik_to_ik_orig and return.
      !
      ALLOCATE(ik_to_ik_orig(nks))
      !
      DO ik = 1, nks
         !
         xkk = xk(:, ik)
         found = .FALSE.
         !
         DO ik_orig = 1, nks_orig
            xk_diff = xkk - xk_orig(:, ik_orig)
            CALL cryst_to_cart(1, xk_diff, at, -1)  ! Cartesian to crystal
            !
            IF ( ALL( ABS(REAL(NINT(xk_diff), DP) - xk_diff) < 1.d-6 ) ) THEN
               found = .TRUE.
               ik_to_ik_orig(ik) = ik_orig
            ENDIF
            !
         ENDDO ! ik_orig
         !
         IF (.NOT. found) THEN
            WRITE(stdout, *) "ik_orig not found for ik = ", ik, " xk = ", xkk
            CALL errore("crpa_skip_nscf_by_rotate", "ik_orig not found", 1)
         ENDIF
         !
      ENDDO
      !
      CALL stop_clock('crpa_nscf_rotate')
      !
      RETURN
      !
   ENDIF
   !
   ALLOCATE(phase(dffts%nnr))
   ALLOCATE(ik_to_ik_orig(nks))
   ALLOCATE(h_evc(npwx * npol, nbnd))
   CALL allocate_bec_type_acc(nkb, nbnd, becp)
   !
   wfcdir = TRIM(tmp_dir_save) // TRIM(prefix) // postfix
   !
   ! Loop over both k and k+q points (nks instead of nksq)
   ! Recompute the wavefunction from those of the SCF calculation
   !
   DO ik = 1, nks
      !
      current_k = ik
      IF (lsda) current_spin = isk(ik)
      xkk = xk(:, ik)
      npw = ngk(ik)
      !
      ! Find k' in the original BZ such that k = k' + G
      !
      found = .FALSE.
      !
      DO ik_orig = 1, nks_orig
         xk_diff = xkk - xk_orig(:, ik_orig)
         CALL cryst_to_cart(1, xk_diff, at, -1)  ! Cartesian to crystal
         !
         IF ( ALL( ABS(REAL(NINT(xk_diff), DP) - xk_diff) < 1.d-6 ) ) THEN
            found = .TRUE.
            g_k = xk_diff
            ik_to_ik_orig(ik) = ik_orig
            CALL cryst_to_cart(1, g_k, bg, 1)  ! Crystal to Cartesian
            EXIT
         ENDIF
         !
      ENDDO ! ik_orig
      !
      IF (.NOT. found) THEN
         WRITE(stdout, *) "ik_orig not found for ik = ", ik, " xk = ", xkk
         CALL errore("crpa_skip_nscf_by_rotate", "ik_orig not found", 1)
      ENDIF
      !
      ! Read SCF wavefunction at ik_orig from collected wavefunction file
      !
      CALL read_collected_wfc_new(wfcdir, ik, ik_orig, evc, "wfc", ierr)
      IF (ierr /= 0) CALL errore("crpa_skip_nscf_by_rotate", "error reading collected wavefunction", 1)
      !
      print*, "ik = ", ik, " ik_orig = ", ik_orig  ! DEBUG
      print*, "g_k = ", g_k  ! DEBUG
      !
      ! If G /= 0, multiply phase factor e^{-i * G * r} to evc.
      !
      IF ( ANY( ABS(g_k) > 1.d-6 ) ) THEN
         !
         ! Find index of -g_k in the list of G vectors. g(:, ig_k) = -g_k
         !
         ig_k = 0
         ig = 1
         gg_ = SUM(g_k * g_k)
         DO WHILE (gg(ig) <= gg_ + 1.d-6)
            IF ( ALL( ABS(g(:, ig) + g_k) < 1.d-6 ) ) THEN
               ig_k = ig
               EXIT
            ENDIF
            ig = ig + 1
         ENDDO
         !
         ! Compute the phase(ir) = e^{-i * gkq * r} if phase is not 1.
         ! Computed phase is used inside the loop over bands.
         !
         phase(:) = (0.0_DP, 0.0_DP)
         IF (ig_k > 0) phase( dffts%nl(ig_k) ) = (1.0_DP, 0.0_DP)
         CALL invfft('Wave', phase, dffts)
         !
         DO ibnd = 1, nbnd
            DO ipol = 1, npol
               !
               psic = (0.0_DP, 0.0_DP)
               !
               ! Copy evc to psic
               !
               istart = 1 + (ipol-1) * npwx
               iend = istart + npw - 1
               !
               psic(dffts%nl( igk_k(1:npw, ik) )) = evc(istart:iend, ibnd)
               !
               IF (gamma_only) psic(dffts%nlm( igk_k(1:npw, ik) )) = CONJG(evc(istart:iend, ibnd))
               !
               ! Multiply the phase factor e^{-i * g_kq * r}
               ! k + q = k' + G, u(k+q) = u(k') * e^{-i * G * r}
               !
               CALL invfft('Wave', psic, dffts)
               psic(1:dffts%nnr) = psic(1:dffts%nnr) * phase(1:dffts%nnr)
               CALL fwfft('Wave', psic, dffts)
               !
               ! Save psic * phase in evc
               !
               evc(istart:iend, ibnd) = psic(dffts%nl( igk_k(1:npw, ik) ))
               !
               IF (gamma_only) THEN
                  ! FIXME
                  CALL errore("crpa_skip_nscf_by_rotate", "gamma_only not implemented", 1)
                  ! IF (gstart == 2) psic(dffts%nlm(1)) = (0.0_DP, 0.0_DP)
                  ! evc(1:npw, n) = CONJG(psic(dffts%nlm( igk_k(1:npw, ik) )))
               ENDIF
               !
            ENDDO ! npol
         ENDDO ! nbnd
         !
      ENDIF ! zerophase
      !
      ! DEBUG: Check the overlap between the original (evq) and new wavefunctions (evc)
      !
      CALL get_buffer(evq, lrwfc, iuwfc, ik)  ! DEBUG
      CALL ZGEMM('C', 'N', nbnd, nbnd, npwx*npol, (1.d0, 0.d0), evc, npwx*npol, evq, npwx*npol, (0.d0, 0.d0), overlap, nbnd) ! DEBUG
      CALL mp_sum(overlap, intra_bgrp_comm) ! DEBUG
      print*, ABS(overlap(1, 1)) ! DEBUG
      ! print*, SUM(ABS(overlap(1:2, 1:2))**2) ! DEBUG
      !
      ! Check evc is an eigenstate of the Hamiltonian
      !
      CALL init_us_2(npw, igk_k(1, ik), xk(1, ik), vkb, .true.)
      CALL g2_kin(ik)
      CALL h_psi(npwx, npw, nbnd, evc, h_evc)
      CALL ZGEMM('C', 'N', nbnd, nbnd, npwx*npol, (1.d0, 0.d0), evc, npwx*npol, h_evc, &
                 npwx*npol, (0.d0, 0.d0), overlap, nbnd)
      CALL mp_sum(overlap, intra_bgrp_comm)
      DO ib = 1, nbnd
         overlap(ib, ib) = overlap(ib, ib) - et(ib, ik)
      ENDDO
      !
      IF (SUM(ABS(overlap)) > 1.d-10) THEN
         WRITE(stdout, '(5x,a)') "Rotated eigenvector is not an eigenstate of the Hamiltonian"
         WRITE(stdout, '(5x,a,I8)') "ik = ", ik
         WRITE(stdout, '(5x,a,ES20.10)') "<evc|H|evc> - E error = ", SUM(ABS(overlap))
       ENDIF
      !
      ! Write rotated evc to buffer
      !
      CALL save_buffer(evc, lrwfc, iuwfc, ik)
      !
   ENDDO ! ik
   !
   ! TODO: Setup et
   !
   ! TODO: Symmetry
   !
   ! call errore("aa", "aa", 1)
   !
   DEALLOCATE(phase)
   DEALLOCATE(h_evc)
   CALL deallocate_bec_type_acc (becp)
   !
   CALL stop_clock('crpa_nscf_rotate')
   !
CONTAINS
   !
   SUBROUTINE read_collected_wfc_new ( dirname, ik, ik_file, arr, label_, ierr_ )
      !------------------------------------------------------------------------
      !
      ! ... reads from directory "dirname" (new file format) for k-point "ik"
      ! ... wavefunctions from collected format into distributed array "arr"
      !
      ! Adapted from read_collected_wfc in pw_restart_new.f90
      ! The difference is the ik_file that gives the global index of the file to read.
      ! Allows to use G vectors at ik different from ik_file
      !
      USE kinds,                ONLY : dp
      USE io_files,             ONLY : iunpun
      USE control_flags,        ONLY : gamma_only
      USE lsda_mod,             ONLY : nspin, isk
      USE klist,                ONLY : nkstot, nks, ngk, igk_k
      USE wvfct,                ONLY : npwx, nbnd
      USE gvect,                ONLY : ig_l2g
      USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm
      USE mp_pools,             ONLY : me_pool, root_pool, intra_pool_comm
      USE mp,                   ONLY : mp_sum, mp_max
      USE io_base,              ONLY : read_wfc
      USE xc_lib,               ONLY : exx_is_active
      USE exx,                  ONLY : nbndproj
      USE io_global,            ONLY : stdout
      USE pw_restart_new,       ONLY : gk_l2gmap_kdip
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      INTEGER, INTENT(IN) :: ik
      INTEGER, INTENT(IN) :: ik_file
      COMPLEX(dp), INTENT(OUT) :: arr(:,:)
      CHARACTER(LEN=3), OPTIONAL, INTENT(IN) :: label_
      INTEGER, OPTIONAL, INTENT(OUT)  :: ierr_
      !
      CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
      CHARACTER(LEN=320)   :: filename, msg
      CHARACTER(LEN=3)     :: label
      LOGICAL              :: read_ace
      INTEGER              :: i, ik_g, ig
      INTEGER              :: npol_, nbnd_
      INTEGER              :: ike, iks, ngk_g, npw_g, ispin
      INTEGER, EXTERNAL    :: global_kpoint_index
      INTEGER, ALLOCATABLE :: mill_k(:,:)
      INTEGER, ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:)
      LOGICAL              :: ionode_k
      REAL(DP)             :: scalef, xk_(3), b1(3), b2(3), b3(3)
      CHARACTER(len=6), EXTERNAL :: int_to_char
      !
      ! ... decide whether to read wfc or ace
      !
      if(present(label_)) then
         label = label_
         if(label.eq."ace") then
            if(.not.exx_is_active()) CALL errore ('pw_restart-read_collected_wfc',&
               "ace but not exx_is_active", 1 )
            read_ace = .true.
         else if(label.eq."wfc") then
            read_ace = .false.
         else
            CALL errore ('pw_restart - read_collected_wfc', "wrong label", 1 )
         end if
      else
         label = "wfc"
         read_ace = .false.
      end if
      !
      ! ... the root processor of each pool reads
      !
      ionode_k = (me_pool == root_pool)
      !
      ik_g = ik_file
      !
      ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
      !
      ALLOCATE ( igk_l2g_kdip( npwx ) )
      !
      ! ... The igk_l2g array yields the correspondence between the
      ! ... local k+G index and the global G index - requires arrays
      ! ... igk_k (k+G indices) and ig_l2g (local to global G index map)
      !
      ALLOCATE ( igk_l2g( npwx ) )
      igk_l2g = 0
      DO ig = 1, ngk(ik)
         igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
      END DO
      !
      ! ... npw_g: the maximum G vector index among all processors
      ! ... ngk_g: global number of k+G vectors for all k points
      !
      npw_g = MAXVAL( igk_l2g(1:ngk(ik)) )
      CALL mp_max( npw_g, intra_pool_comm )
      ngk_g = ngk(ik)
      CALL mp_sum( ngk_g, intra_bgrp_comm)
      !
      ! ... now compute the igk_l2g_kdip local-to-global map
      !
      igk_l2g_kdip = 0
      CALL gk_l2gmap_kdip( npw_g, ngk_g, ngk(ik), igk_l2g, &
         igk_l2g_kdip )
      DEALLOCATE ( igk_l2g )
      !
      IF ( nspin == 2 ) THEN
         !
         ! ... LSDA: spin mapped to k-points, isk(ik) tracks up and down spin
         !
         ik_g = MOD ( ik_g-1, nkstot/2 ) + 1
         ispin = isk(ik)
         filename = TRIM(dirname) // label // updw(ispin) // &
         & TRIM(int_to_char(ik_g))
         !
      ELSE
         !
         filename = TRIM(dirname) // label // TRIM(int_to_char(ik_g))
         !
      ENDIF
      !
      ! ... Miller indices are read from file (but not used)
      !
      ALLOCATE( mill_k ( 3,npwx ) )
      !
      arr = (0.0_DP, 0.0_DP)
      !
      CALL read_wfc( iunpun, filename, root_bgrp, intra_bgrp_comm, &
         ik_g, xk_, ispin, npol_, arr, npw_g, gamma_only, nbnd_, &
         igk_l2g_kdip(:), ngk(ik), b1, b2, b3, mill_k, scalef, ierr_ )
      !
      DEALLOCATE ( mill_k )
      DEALLOCATE ( igk_l2g_kdip )
      !
      IF ( PRESENT (ierr_) ) THEN
         IF ( ierr_ /= 0 ) RETURN
      END IF
      !
      ! ... here one should check for consistency between what is read
      ! ... and what is expected
      !
      IF(read_ace) THEN
         !
         WRITE(stdout, '(5X,A,I8,A)') 'ACE potential read for ', nbnd_, ' bands'
         nbndproj = nbnd_
         !
      ELSE IF ( nbnd_ < nbnd .and..not. read_ace) THEN
         !
         WRITE (msg,'("The number of bands for this run is",I6,", but only",&
         & I6," bands were read from file")')  nbnd, nbnd_
         CALL errore ('read_collected_wfc_new', msg, 1 )
         !
      END IF
      !
      RETURN
      !
   END SUBROUTINE read_collected_wfc_new
   !
END SUBROUTINE crpa_skip_nscf_by_rotate
