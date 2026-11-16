!
! Copyright (C) 2001-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_write_coul_lanczos_chain(iq)
  !-----------------------------------------------------------------------
  !! Writes the RAW (pre-extrapolation) non-Hermitian Lanczos chains built
  !! by crpa_do_lanczos_coulomb() for q point iq to a formatted file
  !! (filU.coul_lanczos_chain#iq), so that the screened Coulomb matrix
  !! elements can later be reconstructed at a different frequency grid or
  !! extrapolation scheme (crpa_spectrum_only = .true., see
  !! crpa_run_spectrum_only.f90) without redoing the chains themselves.
  !! Self-contained: also stores v_coul_bare and the Wannier-pair index
  !! lists (pert_iwlist/jwlist/Rlist, mels_iwlist/jwlist/Rlist), since
  !! crpa_spectrum_only mode never runs crpa_setup_wannier/
  !! crpa_prepare_perturbations_q to rebuild them.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode, stdout
  USE lsda_mod,   ONLY : nspin
  USE crpacom,    ONLY : crpa_coul_lanczos_steps, filU, v_coul_bare
  USE crpa_pert,  ONLY : npert_tot, nmels_tot, pert_iwlist, pert_jwlist, pert_Rlist, &
                          mels_iwlist, mels_jwlist, mels_Rlist
  USE crpa_coul_lanczos, ONLY : coul_beta_store, coul_gamma_store, coul_zeta_store
  USE qpoint,     ONLY : xq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !! index of the q point
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  CHARACTER(LEN=256) :: filename
  INTEGER :: iun, ios, ipert, jmels
  REAL(DP) :: degspin
  !
  IF (.NOT. ionode) RETURN
  !
  degspin = 2.0_DP
  IF (nspin /= 1) degspin = 1.0_DP
  !
  filename = TRIM(filU) // ".coul_lanczos_chain" // TRIM(int_to_char(iq))
  OPEN(NEWUNIT=iun, FILE=TRIM(filename), STATUS='replace', IOSTAT=ios)
  IF (ios /= 0) CALL errore('crpa_write_coul_lanczos_chain', &
     'cannot open '//TRIM(filename), 1)
  !
  WRITE(iun, '(a)') "# CRPA screened-Coulomb Lanczos chains (pert_basis=wannier, sequential)"
  WRITE(iun, *) npert_tot, nmels_tot, crpa_coul_lanczos_steps
  WRITE(iun, *) xq(1), xq(2), xq(3)
  WRITE(iun, *) degspin
  !
  DO ipert = 1, npert_tot
     WRITE(iun, *) pert_iwlist(ipert), pert_jwlist(ipert), pert_Rlist(:, ipert)
  ENDDO
  DO jmels = 1, nmels_tot
     WRITE(iun, *) mels_iwlist(jmels), mels_jwlist(jmels), mels_Rlist(:, jmels)
  ENDDO
  !
  DO ipert = 1, npert_tot
     WRITE(iun, *) v_coul_bare(:, ipert)
  ENDDO
  !
  DO ipert = 1, npert_tot
     WRITE(iun, *) coul_beta_store(ipert, :)
     WRITE(iun, *) coul_gamma_store(ipert, :)
     DO jmels = 1, nmels_tot
        WRITE(iun, *) coul_zeta_store(ipert, jmels, :)
     ENDDO
  ENDDO
  !
  CLOSE(iun)
  !
  WRITE(stdout, '(5X,A)') "INFO: Coulomb Lanczos chain written to "//TRIM(filename)
  !
  RETURN
  !
END SUBROUTINE crpa_write_coul_lanczos_chain
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_read_coul_lanczos_chain(iq, exst, nmax, nfile)
  !-----------------------------------------------------------------------
  !! Reads back the chain(s) written by crpa_write_coul_lanczos_chain() for
  !! q point iq. Sets exst = .false. (and does nothing else) if the file for
  !! this q is not found, so the caller can loop over iq until exhausted.
  !!
  !! Two modes:
  !!  - nmax = 0 (crpa_spectrum_only): allocates crpa_pert's index lists,
  !!    crpacom's v_coul_bare, and crpa_coul_lanczos's chain-coefficient
  !!    arrays (and the _ext arrays used by crpa_extrapolate_coul/
  !!    crpa_calc_coulomb_lanczos) from the file's own npert_tot/nmels_tot/
  !!    chain length, overriding crpa_coul_lanczos_steps (and, if needed,
  !!    crpa_coul_lanczos_steps_ext) from the current namelist: only the
  !!    frequency grid and extrapolation scheme are taken fresh from the
  !!    current input.
  !!  - nmax > 0 (q-level restart, crpa_coul_q_from_chain_file): the
  !!    index lists, xq and v_coul_bare are those of the current calculation
  !!    (crpa_prepare_q/crpa_prepare_perturbations_q), whose input must be the
  !!    same. Only the first nmax (= crpa_coul_lanczos_steps) coefficients are
  !!    loaded, and no global is changed. If the file holds fewer than nmax
  !!    steps, nothing is loaded and exst = .false.; nfile returns the length.
  !!
  !! spectrum_only mode also allocates v_coul_scrd (crpa_allocate_coul_lanczos,
  !! which allocates the wavefunction-sized chain arrays, is not called here).
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE crpacom,    ONLY : crpa_coul_lanczos_steps, crpa_coul_lanczos_steps_ext, &
                          crpa_coul_extrapolation, filU, v_coul_bare, v_coul_scrd, &
                          crpa_nfreq
  USE crpa_pert,  ONLY : npert_tot, nmels_tot, pert_iwlist, pert_jwlist, pert_Rlist, &
                          mels_iwlist, mels_jwlist, mels_Rlist, pert_basis
  USE crpa_coul_lanczos, ONLY : coul_beta_store, coul_gamma_store, coul_zeta_store, &
                                 coul_beta_store_ext, coul_gamma_store_ext
  USE qpoint,     ONLY : xq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  !! index of the q point
  LOGICAL, INTENT(OUT) :: exst
  INTEGER, INTENT(IN) :: nmax
  !! restart mode if > 0: number of coefficients to load
  INTEGER, INTENT(OUT) :: nfile
  !! chain length stored on file (0 if there is no file)
  !
  CHARACTER(LEN=*), PARAMETER :: sub = 'crpa_read_coul_lanczos_chain'
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  CHARACTER(LEN=256) :: filename, line
  INTEGER :: iun, ios, ipert, jmels, nsteps, nload, npert_f, nmels_f
  REAL(DP) :: degspin, xq_f(3)
  REAL(DP), ALLOCATABLE :: rtmp(:)
  COMPLEX(DP), ALLOCATABLE :: ztmp(:)
  LOGICAL :: lrestart
  !
  lrestart = (nmax > 0)
  nfile = 0
  !
  filename = TRIM(filU) // ".coul_lanczos_chain" // TRIM(int_to_char(iq))
  INQUIRE(FILE=TRIM(filename), EXIST=exst)
  IF (.NOT. exst) RETURN
  !
  OPEN(NEWUNIT=iun, FILE=TRIM(filename), STATUS='old', IOSTAT=ios)
  IF (ios /= 0) CALL errore(sub, 'cannot open '//TRIM(filename), 1)
  !
  READ(iun, '(A)') line
  READ(iun, *) npert_f, nmels_f, nsteps
  READ(iun, *) xq_f(1), xq_f(2), xq_f(3)
  READ(iun, *) degspin
  nfile = nsteps
  !
  IF (lrestart) THEN
     !
     IF (nsteps < nmax) THEN
        CLOSE(iun)
        exst = .FALSE.
        RETURN
     ENDIF
     IF (npert_f /= npert_tot .OR. nmels_f /= nmels_tot) CALL errore(sub, TRIM(filename)// &
        ' does not match this calculation (number of perturbations or matrix elements)', 1)
     nload = nmax
     DO ipert = 1, npert_tot + nmels_tot
        READ(iun, *)          ! pair lists: those of the current calculation
     ENDDO
     !
  ELSE
     !
     pert_basis = 'wannier'
     npert_tot = npert_f
     nmels_tot = nmels_f
     xq(:) = xq_f(:)
     nload = nsteps
     !
     crpa_coul_lanczos_steps = nsteps
     IF (TRIM(crpa_coul_extrapolation) == 'no') THEN
        crpa_coul_lanczos_steps_ext = crpa_coul_lanczos_steps
     ELSE IF (crpa_coul_lanczos_steps_ext < crpa_coul_lanczos_steps) THEN
        CALL errore(sub, 'crpa_coul_lanczos_steps_ext must be >= the chain length stored on file', 1)
     ENDIF
     !
     IF (ALLOCATED(pert_iwlist)) DEALLOCATE(pert_iwlist, pert_jwlist, pert_Rlist)
     IF (ALLOCATED(mels_iwlist)) DEALLOCATE(mels_iwlist, mels_jwlist, mels_Rlist)
     IF (ALLOCATED(v_coul_bare)) DEALLOCATE(v_coul_bare)
     IF (ALLOCATED(v_coul_scrd)) DEALLOCATE(v_coul_scrd)
     !
     ALLOCATE(pert_iwlist(npert_tot), pert_jwlist(npert_tot), pert_Rlist(3, npert_tot))
     ALLOCATE(mels_iwlist(nmels_tot), mels_jwlist(nmels_tot), mels_Rlist(3, nmels_tot))
     ALLOCATE(v_coul_bare(nmels_tot, npert_tot))
     ALLOCATE(v_coul_scrd(nmels_tot, npert_tot, crpa_nfreq))
     !
     DO ipert = 1, npert_tot
        READ(iun, *) pert_iwlist(ipert), pert_jwlist(ipert), pert_Rlist(:, ipert)
     ENDDO
     DO jmels = 1, nmels_tot
        READ(iun, *) mels_iwlist(jmels), mels_jwlist(jmels), mels_Rlist(:, jmels)
     ENDDO
     !
  ENDIF
  !
  IF (ALLOCATED(coul_beta_store))  DEALLOCATE(coul_beta_store)
  IF (ALLOCATED(coul_gamma_store)) DEALLOCATE(coul_gamma_store)
  IF (ALLOCATED(coul_zeta_store))  DEALLOCATE(coul_zeta_store)
  IF (ALLOCATED(coul_beta_store_ext))  DEALLOCATE(coul_beta_store_ext)
  IF (ALLOCATED(coul_gamma_store_ext)) DEALLOCATE(coul_gamma_store_ext)
  ALLOCATE(coul_beta_store(npert_tot, nload))
  ALLOCATE(coul_gamma_store(npert_tot, nload))
  ALLOCATE(coul_zeta_store(npert_tot, nmels_tot, nload))
  ALLOCATE(coul_beta_store_ext(npert_tot, crpa_coul_lanczos_steps_ext))
  ALLOCATE(coul_gamma_store_ext(npert_tot, crpa_coul_lanczos_steps_ext))
  !
  DO ipert = 1, npert_tot
     READ(iun, *) v_coul_bare(:, ipert)
  ENDDO
  !
  ! Rows hold nsteps entries; only the first nload are kept (truncation).
  !
  ALLOCATE(rtmp(nsteps), ztmp(nsteps))
  DO ipert = 1, npert_tot
     READ(iun, *) rtmp
     coul_beta_store(ipert, :) = rtmp(1:nload)
     READ(iun, *) rtmp
     coul_gamma_store(ipert, :) = rtmp(1:nload)
     DO jmels = 1, nmels_tot
        READ(iun, *) ztmp
        coul_zeta_store(ipert, jmels, :) = ztmp(1:nload)
     ENDDO
  ENDDO
  DEALLOCATE(rtmp, ztmp)
  !
  CLOSE(iun)
  !
  WRITE(stdout, '(5X,A,I8,A,I8,A)') "INFO: Coulomb Lanczos chain read from "//TRIM(filename)// &
     " (", nsteps, " steps,", nload, " used)"
  !
  RETURN
  !
END SUBROUTINE crpa_read_coul_lanczos_chain
