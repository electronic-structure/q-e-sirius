!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
SUBROUTINE crpa_prepare_perturbations_q(iq)
  !----------------------------------------------------------------------------
  !
  ! This routine prepares the perturbations for the iq-th q point.
  ! The common prepration that applies to all q points are done in
  ! crpa_prepare_perturbations.
  !
  !----------------------------------------------------------------------------
  !
  USE io_global,         ONLY : stdout
  USE qpoint,            ONLY : nksqtot
  USE crpacom,           ONLY : v_coul_bare, v_coul_scrd, dist_thr, active_bands_min, &
                                active_bands_max, dist_thr_large
  USE crpa_pert,         ONLY : pert_basis, npert_tot, nmels_tot, pert_nbnd, pert_ibnd_max, &
                                pert_ibnd_min, pert_iwlist, pert_jwlist, pert_Rlist, &
                                mels_iwlist, mels_jwlist, mels_Rlist
  USE w90_interface,     ONLY : num_wann
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !! Index of the q point to prepare the perturbations for.
  !
  IF (pert_basis == 'bands') THEN
    !
    pert_ibnd_max = active_bands_max
    pert_ibnd_min = active_bands_min
    pert_nbnd = pert_ibnd_max - pert_ibnd_min + 1
    !
    npert_tot = pert_nbnd * pert_nbnd * nksqtot
    nmels_tot = npert_tot
    !
    ALLOCATE(v_coul_bare(nmels_tot, npert_tot))
    ALLOCATE(v_coul_scrd(nmels_tot, npert_tot))
    !
  ELSEIF (pert_basis == 'wannier') THEN
    !
    WRITE(stdout, '(/5x, A)') 'Wannier function pairs for perturbation'
    CALL crpa_find_close_WF_pairs(npert_tot, pert_iwlist, pert_jwlist, pert_Rlist, dist_thr)
    !
    WRITE(stdout, '(/5x, A)') 'Wannier function pairs for matrix element calculation'
    CALL crpa_find_close_WF_pairs(nmels_tot, mels_iwlist, mels_jwlist, mels_Rlist, dist_thr_large)
    !
    pert_nbnd = num_wann
    !
    ALLOCATE(v_coul_bare(nmels_tot, npert_tot))
    ALLOCATE(v_coul_scrd(nmels_tot, npert_tot))
    !
  ENDIF
  !
  v_coul_bare = (0.d0, 0.d0)
  v_coul_scrd = (0.d0, 0.d0)
  !
  RETURN
  !
  CONTAINS
  !
  !------------------------------------------------------------------------------
  SUBROUTINE crpa_find_close_WF_pairs(npairs, iwlist, jwlist, Rlist, dist_thr)
  !------------------------------------------------------------------------------
  !! Find the list of WF pairs (i0, jR) whose distance is below dist_thr.
  !! We consider all supercell periodic images of the WFs using w90_find_wigner_seitz.
  !------------------------------------------------------------------------------
    !
    USE kinds,             ONLY : DP
    USE io_global,         ONLY : stdout
    USE cell_base,         ONLY : at, alat
    USE w90_interface,     ONLY : num_wann, mp_grid, wann_centers
    USE w90_wigner,        ONLY : irvec, w90_find_wigner_seitz, w90_wigner_allocate, &
                                  w90_wigner_deallocate
    USE crpacom,           ONLY : crpa_mode
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(OUT) :: npairs
    !! Number of selected WF pairs
    INTEGER, ALLOCATABLE, INTENT(OUT) :: iwlist(:)
    !! List of selected R vectors in crystal coordinates.
    INTEGER, ALLOCATABLE, INTENT(OUT) :: jwlist(:)
    !! List of selected R vectors in crystal coordinates.
    INTEGER, ALLOCATABLE, INTENT(OUT) :: Rlist(:, :)
    !! List of selected R vectors in crystal coordinates.
    REAL(DP), INTENT(IN) :: dist_thr
    !! Cutoff distance for the WFs in bohr units.
    !
    INTEGER :: iR1, iR2, iR3
    !! Indices for the supercell periodic images.
    INTEGER :: iw, jw
    !! Indices for the WFs.
    INTEGER :: i
    !! Indices for the selected WF pairs.
    INTEGER :: nrr
    !! Wigner-Seitz degeneracy. (Needed to call w90_find_wigner_seitz, not used here.)
    INTEGER :: npairs_max
    !! Maximum possible number of selected WF pairs.
    REAL(DP) :: mindist
    !! Minimum distance between WF pairs.
    REAL(DP) :: R(3)
    !! R vector in Cartesian coordinates.
    REAL(DP) :: r_ijR(3)
    !! Relative distance in Cartesian coordinates.
    INTEGER(DP), ALLOCATABLE :: iwlist_tmp(:)
    !! Temporary list of i indices
    INTEGER(DP), ALLOCATABLE :: jwlist_tmp(:)
    !! Temporary list of j indices
    INTEGER(DP), ALLOCATABLE :: Rlist_tmp(:, :)
    !! Temporary list of R vectors.
    REAL(DP), ALLOCATABLE :: mindist_list(:)
    !! List of minimum distances for selected R vectors.
    !
    npairs_max = num_wann**2 * mp_grid(1) * mp_grid(2) * mp_grid(3)
    ALLOCATE(Rlist_tmp(3, npairs_max))
    ALLOCATE(iwlist_tmp(npairs_max))
    ALLOCATE(jwlist_tmp(npairs_max))
    ALLOCATE(mindist_list(npairs_max))
    !
    CALL w90_wigner_allocate(mp_grid(1), mp_grid(2), mp_grid(3))
    !
    npairs = 0
    !
    DO iR3 = 0, mp_grid(3) - 1
      DO iR2 = 0, mp_grid(2) - 1
        DO iR1 = 0, mp_grid(1) - 1
          R(1) = REAL(iR1, DP)
          R(2) = REAL(iR2, DP)
          R(3) = REAL(iR3, DP)
          CALL cryst_to_cart(1, R, at, +1)
          !
          ! Compute minimum distance |w_i - w_j - R - T| for all WF pairs (iw, jw) and
          ! Born-von Karman vectors T
          !
          DO jw = 1, num_wann
            DO iw = 1, num_wann
              !
              IF (TRIM(crpa_mode)=='dHP' .AND. (iw /= jw)) CYCLE
              !
              IF (TRIM(crpa_mode)=='debug' .AND. (jw /= 1 .OR. iw /= 1)) CYCLE
              !
              r_ijR = wann_centers(:, iw) - wann_centers(:, jw) - R
              CALL w90_find_wigner_seitz(r_ijR, nrr, mindist)
              !
              ! mindist is in alat units, dist_thr is in bohr units.
              !
              IF (mindist <= dist_thr / alat) THEN
                npairs = npairs + 1
                iwlist_tmp(npairs) = iw
                jwlist_tmp(npairs) = jw
                Rlist_tmp(:, npairs) = (/ iR1, iR2, iR3 /) + irvec(:, 1)  ! R + T vector
                mindist_list(npairs) = mindist
              ENDIF
              !
            ENDDO ! iw
          ENDDO ! jw
          !
        ENDDO ! iR3
      ENDDO ! iR2
    ENDDO ! iR1
    !
    CALL w90_wigner_deallocate()
    !
    ALLOCATE(iwlist(npairs))
    ALLOCATE(jwlist(npairs))
    ALLOCATE(Rlist(3, npairs))
    iwlist = iwlist_tmp(1:npairs)
    jwlist = jwlist_tmp(1:npairs)
    Rlist = Rlist_tmp(:, 1:npairs)
    !
    DEALLOCATE(iwlist_tmp)
    DEALLOCATE(jwlist_tmp)
    DEALLOCATE(Rlist_tmp)
    !
    WRITE(stdout, '(5x, A, I6)') 'Number of selected WF pairs : ', npairs
    WRITE(stdout, '(5x, A, F12.6)') 'Distance cutoff (bohr) : ', dist_thr
    WRITE(stdout, '(5x, A)') '# i j R (crystal) |r_i0 - r_jR| (bohr)'
    DO i = 1, npairs
      WRITE(stdout, '(5x, 5I8, F12.6)') iwlist(i), jwlist(i), Rlist(:, i), mindist_list(i) * alat
    ENDDO
    !
  !------------------------------------------------------------------------------
  END SUBROUTINE crpa_find_close_WF_pairs
  !------------------------------------------------------------------------------
  !
!------------------------------------------------------------------------------
END SUBROUTINE crpa_prepare_perturbations_q
!------------------------------------------------------------------------------
