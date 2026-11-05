!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_phescf_w()
  !-----------------------------------------------------------------------
  !! This is the main driver for the calculation of the
  !! response to an electric field and related quantities.
  !
  USE io_global,       ONLY : stdout
  USE paw_variables,   ONLY : okpaw
  USE uspp,            ONLY : okvan, nkb
  USE uspp_param,      ONLY : nhm
  USE ions_base,       ONLY : nat, ntyp => nsp
  USE noncollin_module,ONLY : noncolin, nspin_mag, lspinorb
  USE lsda_mod,        ONLY : nspin
  USE lsda_mod,        ONLY : nspin
  USE wvfct,           ONLY : nbnd
  USE becmod,          ONLY : becp, allocate_bec_type_acc, deallocate_bec_type_acc
  USE crpacom,         ONLY : crpa_epsil, crpa_done_epsil, crpa_diel_llanczos, &
                              crpa_nfreq, crpa_freq, current_w, epsilon_w, epsilon_w_diel
  USE output,          ONLY : fildrho
  USE lrus,            ONLY : int3, int3_nc, int3_paw, dpqq, dpqq_so
  USE control_lr,      ONLY : convt, rec_code, rec_code_read, where_rec
  USE ldaU,            ONLY : lda_plus_u, Hubbard_lmax
  USE ldaU_lr,         ONLY : dnsscf
  USE ldaU_ph,         ONLY : dnsscf_all_modes
  USE control_flags,   ONLY : iverbosity
  USE write_hub
  USE ph_restart,      ONLY : ph_writefile
  !
  IMPLICIT NONE
  !
  INTEGER :: ierr, ifreq
  !
  !
  ! The single-frequency restart mechanism (crpa_write/read_epsilon_w) only
  ! covers the CG path's epsilon_w(3,3). It is not wired up for the
  ! frequency-array Lanczos path (epsilon_w_diel), so that path must always
  ! recompute rather than risk silently short-circuiting on a stale
  ! "after_diel" status record from a previous Lanczos run. The same holds
  ! for the frequency grid of the CG path (crpa_nfreq > 1).
  !
  IF ( rec_code_read >  1 .AND. .NOT. crpa_diel_llanczos .AND. crpa_nfreq == 1 ) THEN
     CALL crpa_read_epsilon_w(ierr)
     IF (ierr == 0) THEN
        WRITE( stdout, '(/,5X,"Dielectric function read from previous run")' )
        crpa_done_epsil = .TRUE.
        CALL crpa_summarize_epsilon()
     ENDIF
     RETURN
  ENDIF
  !
  IF (okvan) THEN
     ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, 3))
     IF (okpaw) ALLOCATE (int3_paw ( nhm, nhm, nat, nspin_mag, 3))
     IF (noncolin) ALLOCATE(int3_nc( nhm, nhm, nat, nspin, 3))
     !
     ! USPP/PAW electric-field term: dvpsi_e -> adddvepsi_us needs dpqq, the
     ! dipole moment of the augmentation functions Q_nm (LR_Modules/lrcom.f90).
     ! PHonon builds it in phq_init.f90:
     !    IF ( (epsil .OR. zue .OR. l_head) .AND. okvan ) THEN
     !       CALL compute_qdipol(dpqq); ... ; CALL qdipol_cryst()
     ! with the array itself allocated in allocate_phq.f90:127. CRPA has its own
     ! init (crpa_init_q/crpa_allocate_q) and did NEITHER, so adddvepsi_us
     ! dereferenced an unallocated dpqq and SEGFAULTED for USPP (verified on
     ! Si + Si.pbe-nl-rrkjus_psl.1.0.0.UPF, 2026-09-02). qdipol_cryst() rotates
     ! dpqq to crystal axes, which is the convention dvpsi_e uses (at(:,ipol)).
     !
     ALLOCATE (dpqq( nhm, nhm, 3, ntyp))
     CALL compute_qdipol(dpqq)
     IF (lspinorb) THEN
        ALLOCATE (dpqq_so( nhm, nhm, nspin, 3, ntyp))
        CALL compute_qdipol_so(dpqq, dpqq_so)
     ENDIF
     CALL qdipol_cryst()
  ENDIF
  !
  ! Set symmetry representation in lr_symm_base
  !
  CALL ph_set_upert_e()
  !
  ! DFPT+U: dnsscf in the electric field calculation
  ! is the scf change of atomic occupations ns induced by the electric field.
  ! dnsscf_all_modes = dnsscf because nirr=1, number of perturbations = 3.
  !
  IF (lda_plus_u) THEN
     ALLOCATE (dnsscf(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,3))
     ALLOCATE (dnsscf_all_modes(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,3))
     dnsscf = (0.d0, 0.d0)
     dnsscf_all_modes = (0.d0, 0.d0)
  ENDIF
  !
  IF ((crpa_epsil.AND..NOT.crpa_done_epsil)) THEN
     !
     WRITE( stdout, '(/,5X,"Electric Fields Calculation")' )
     !
     IF (crpa_diel_llanczos) THEN
        !
        ! Lanczos-chain path (see crpa_diel_lanczos.f90): no self-consistent
        ! Sternheimer solve, so there is no convergence criterion to check.
        !
        CALL allocate_bec_type_acc (nkb, nbnd, becp)
        CALL crpa_do_lanczos_epsilon()
        CALL deallocate_bec_type_acc (becp)
        convt = .TRUE.
        !
     ELSE IF (crpa_nfreq == 1) THEN
        !
        CALL allocate_bec_type_acc (nkb, nbnd, becp)
        CALL crpa_solve_e_w()
        CALL deallocate_bec_type_acc (becp)
        !
     ELSE
        !
        ! CG frequency grid: one self-consistent solve and one dielectric
        ! tensor per frequency, collected in epsilon_w_diel as on the
        ! Lanczos path and written to epsilon_w.dat at the end
        !
        IF (ALLOCATED(epsilon_w_diel)) DEALLOCATE(epsilon_w_diel)
        ALLOCATE(epsilon_w_diel(3, 3, crpa_nfreq))
        !
        CALL allocate_bec_type_acc (nkb, nbnd, becp)
        DO ifreq = 1, crpa_nfreq
           current_w = crpa_freq(ifreq)
           WRITE(stdout, '(/,5X,"Frequency #",i4," of",i4,": w = (",f10.6,",",f10.6,") Ry")') &
                ifreq, crpa_nfreq, current_w
           crpa_done_epsil = .FALSE.
           CALL crpa_solve_e_w()
           IF (.NOT. convt) CALL errore('crpa_phescf_w', &
              'Dielectric tensor evaluation not converged.', ifreq)
           CALL crpa_dielec()
           epsilon_w_diel(:,:,ifreq) = epsilon_w(:,:)
        ENDDO
        CALL deallocate_bec_type_acc (becp)
        !
        CALL write_eps_file('epsilon_w.dat', crpa_nfreq, crpa_freq, epsilon_w_diel)
        !
     ENDIF
     !
     WRITE( stdout, '(/,5X,"End of electric fields calculation")' )
     !
     IF ( convt ) THEN
        !
        ! ... calculate the dielectric tensor epsilon
        !
        IF (.NOT. crpa_done_epsil) THEN
           IF (crpa_diel_llanczos) THEN
              CALL crpa_calc_chi_epsilon()
           ELSE
              CALL crpa_dielec()
           ENDIF
        ELSE
           CALL crpa_summarize_epsilon()
        ENDIF
        !
        IF ( fildrho /= ' ' .AND. .NOT. crpa_diel_llanczos ) CALL punch_plot_e()
        !
     ELSE
        !
        CALL errore('crpa_phescf_w', 'Dielectric tensor evaluation not converged.', 1)
        !
     END IF
     !
     where_rec='after_diel'
     rec_code=2
     CALL ph_writefile('status_ph', 0, 0, ierr)
  ENDIF
  !
  IF (okvan) THEN
     DEALLOCATE (int3)
     IF (okpaw) DEALLOCATE (int3_paw)
     IF (noncolin) DEALLOCATE(int3_nc)
     IF (ALLOCATED(dpqq))    DEALLOCATE (dpqq)
     IF (ALLOCATED(dpqq_so)) DEALLOCATE (dpqq_so)
  ENDIF
  !
  CALL ph_deallocate_upert()
  !
  ! DFPT+U
  !
  IF (lda_plus_u) THEN
     !
     ! Write dnsscf_all_modes in the cartesian coordinates
     ! to the standard output
     !
     IF (iverbosity==1) CALL write_dnsscf_e()
     !
     DEALLOCATE (dnsscf)
     DEALLOCATE (dnsscf_all_modes)
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE crpa_phescf_w
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_summarize_epsilon()
  !-----------------------------------------------------------------------
  !! Print the complex macroscopic dielectric tensor to stdout.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE crpacom,    ONLY : epsilon_w, crpa_done_epsil
  USE control_lr, ONLY : lrpa, lnoloc
  !
  IMPLICIT NONE
  INTEGER :: ipol, jpol
  !
  IF (.NOT. crpa_done_epsil) RETURN
  !
  IF (lnoloc) THEN
     WRITE(stdout, '(/,10x,"Dielectric constant in cartesian axis (DV_Hxc=0)",/)')
  ELSE IF (lrpa) THEN
     WRITE(stdout, '(/,10x,"RPA dielectric constant in cartesian axis (DV_xc=0)",/)')
  ELSE
     WRITE(stdout, '(/,10x,"Dielectric constant in cartesian axis ",/)')
  ENDIF
  !
  DO ipol = 1, 3
     WRITE(stdout, '(10x,3("(",f15.6,",",f15.6,")",2x))') &
          (DBLE(epsilon_w(ipol,jpol)), AIMAG(epsilon_w(ipol,jpol)), jpol=1,3)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE crpa_summarize_epsilon
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_write_epsilon_w()
  !-----------------------------------------------------------------------
  !! Write the complex dielectric tensor epsilon_w to a binary file
  !! in tmp_dir_epsw so it can be recovered in subsequent runs.
  !
  USE io_global,  ONLY : ionode, ionode_id
  USE mp,         ONLY : mp_bcast
  USE mp_images,  ONLY : intra_image_comm
  USE crpacom,    ONLY : epsilon_w
  USE control_ph, ONLY : tmp_dir_ph
  !
  IMPLICIT NONE
  INTEGER :: iunit = 98, ios
  CHARACTER(LEN=320) :: filename
  !
  filename = TRIM(tmp_dir_ph) // 'epsilon_w.dat'
  IF (ionode) THEN
     OPEN(UNIT=iunit, FILE=TRIM(filename), FORM='unformatted', &
          STATUS='unknown', IOSTAT=ios)
     IF (ios == 0) THEN
        WRITE(iunit) epsilon_w
        CLOSE(iunit, STATUS='keep')
     ENDIF
  ENDIF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  IF (ios /= 0) CALL errore('crpa_write_epsilon_w', &
       'cannot write '//TRIM(filename), 1)
  !
END SUBROUTINE crpa_write_epsilon_w
!
!-----------------------------------------------------------------------
SUBROUTINE crpa_read_epsilon_w(ierr)
  !-----------------------------------------------------------------------
  !! Read the complex dielectric tensor epsilon_w from the binary file
  !! written by crpa_write_epsilon_w. Returns ierr /= 0 if not found.
  !
  USE io_global,  ONLY : ionode, ionode_id
  USE mp,         ONLY : mp_bcast
  USE mp_images,  ONLY : intra_image_comm
  USE crpacom,    ONLY : epsilon_w
  USE control_ph, ONLY : tmp_dir_ph
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ierr
  INTEGER :: iunit = 97
  LOGICAL :: exst
  CHARACTER(LEN=320) :: filename
  !
  ierr = 0
  filename = TRIM(tmp_dir_ph) // 'epsilon_w.dat'
  IF (ionode) THEN
     INQUIRE(FILE=TRIM(filename), EXIST=exst)
     IF (exst) THEN
        OPEN(UNIT=iunit, FILE=TRIM(filename), FORM='unformatted', &
             STATUS='old', IOSTAT=ierr)
        IF (ierr == 0) THEN
           READ(iunit, IOSTAT=ierr) epsilon_w
           CLOSE(iunit, STATUS='keep')
        ENDIF
     ELSE
        ierr = 1
     ENDIF
  ENDIF
  CALL mp_bcast(ierr, ionode_id, intra_image_comm)
  CALL mp_bcast(epsilon_w, ionode_id, intra_image_comm)
  !
END SUBROUTINE crpa_read_epsilon_w

