  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Roxana Margine
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE io_eliashberg
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the IO part of the superconductivity part of EPW
  !!  
  IMPLICIT NONE
  ! 
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_read_aniso_iaxis( itemp )
    !-----------------------------------------------------------------------
    !!  
    !! This routine reads from file the anisotropic Delta and Znorm on the imaginary-axis
    !! 
    !! input
    !!
    !! itemp  - temperature point
    !!
    !---------------------------------------------------------------------- 
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nstemp, fsthick
    USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, Agap, wsi, NZnormi, Znormi, Deltai, & 
                              AZnormi, NAZnormi, ADeltai, nkfs, nbndfs, ef0, ekfs, &
                              dosef, wkfs, w0g
    USE constants_epw, ONLY : kelvin2eV, eps6, zero
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE superconductivity, ONLY : mem_size_eliashberg, free_energy
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on band
    INTEGER :: imelt
    !! Required allocation of memory
    INTEGER :: ios
    !! Status variables when reading a file
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: eband
    !! Temporary variable for eigenvalue
    REAL(KIND = DP) :: omega
    !! Temporary variable for frequency
    REAL(KIND = DP) :: weight
    CHARACTER(LEN = 256) :: name1, word
    !
    ! get the size of required allocated memory 
    imelt = ( 1 + nbndfs * nkfs ) * nstemp + ( 3 + 3 * nbndfs * nkfs ) * nsiw(itemp)
    CALL mem_size_eliashberg( imelt )
    !
    IF (.NOT. ALLOCATED(gap) )      ALLOCATE(gap(nstemp) )
    IF (.NOT. ALLOCATED(Agap) )     ALLOCATE(Agap(nbndfs,nkfs,nstemp) )
    IF (.NOT. ALLOCATED(Deltai) )   ALLOCATE(Deltai(nsiw(itemp)) )
    IF (.NOT. ALLOCATED(Znormi) )   ALLOCATE(Znormi(nsiw(itemp)) )
    IF (.NOT. ALLOCATED(NZnormi) )  ALLOCATE(NZnormi(nsiw(itemp)) )
    IF (.NOT. ALLOCATED(ADeltai) )  ALLOCATE(ADeltai(nbndfs,nkfs,nsiw(itemp)) )
    IF (.NOT. ALLOCATED(AZnormi) )  ALLOCATE(AZnormi(nbndfs,nkfs,nsiw(itemp)) )
    IF (.NOT. ALLOCATED(NAZnormi) ) ALLOCATE(NAZnormi(nbndfs,nkfs,nsiw(itemp)) )
    gap(:) = zero
    Agap(:, :, :) = zero
    Deltai(:) = zero
    Znormi(:) = zero
    NZnormi(:) = zero
    ADeltai(:, :, :) = zero
    AZnormi(:, :, :) = zero
    NAZnormi(:, :, :) = zero
    !
    IF (mpime == ionode_id) THEN     
      !   
      temp = estemp(itemp) / kelvin2eV
      ! anisotropic case
      IF (temp < 10.d0) THEN
         WRITE(name1,'(a,a14,f4.2)') TRIM(prefix),'.imag_aniso_00', temp
      ELSEIF (temp >= 10.d0) THEN
         WRITE(name1,'(a,a13,f5.2)') TRIM(prefix),'.imag_aniso_0', temp
      ELSEIF (temp >= 100.d0) THEN
         WRITE(name1,'(a,a12,f6.2)') TRIM(prefix),'.imag_aniso_', temp
      ENDIF 
      ! 
      OPEN(iufilgap, FILE = name1, FORM = 'formatted', err=100, IOSTAT = ios)
100 CALL errore('eliashberg_read_aniso_iaxis','opening file '//name1,ABS(ios))
      READ(iufilgap,'(a)') word
      DO iw = 1, nsiw(itemp) ! loop over omega
         DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
               IF (ABS(ekfs(ibnd,ik) - ef0 ) < fsthick) THEN
                  READ(iufilgap,'(5ES20.10)') omega, eband, AZnormi(ibnd,ik,iw), ADeltai(ibnd,ik,iw), NAZnormi(ibnd,ik,iw)
                  IF (iw == 1 ) & 
                     Agap(ibnd,ik,itemp) = ADeltai(ibnd,ik,1)
               ENDIF
            ENDDO ! ibnd
         ENDDO ! ik             
         IF (ABS(wsi(iw)-omega) > eps6 ) &
            CALL errore('eliashberg_read_aniso_iaxis','temperature not the same with the input',1)
      ENDDO ! iw
      CLOSE(iufilgap)
      !
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = 1, nkfs
           DO ibnd = 1, nbndfs
              IF (ABS(ekfs(ibnd,ik) - ef0 ) < fsthick) THEN
                 weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
                 Znormi(iw) = Znormi(iw) + weight * AZnormi(ibnd,ik,iw)
                 Deltai(iw) = Deltai(iw) + weight * ADeltai(ibnd,ik,iw)
                 NZnormi(iw) = NZnormi(iw) + weight * NAZnormi(ibnd,ik,iw)
              ENDIF
           ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! iw
      gap(itemp) = Deltai(1)
      gap0 = gap(itemp)
      !
      CALL gap_FS( itemp )
      !
      IF (iverbosity == 2 ) &
         CALL free_energy( itemp )
      !
    ENDIF
    CALL mp_bcast( Deltai, ionode_id, inter_pool_comm )
    CALL mp_bcast( Znormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( NZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( ADeltai, ionode_id, inter_pool_comm )
    CALL mp_bcast( AZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( NAZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap0, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap, ionode_id, inter_pool_comm )
    CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    RETURN
    !
    END SUBROUTINE eliashberg_read_aniso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_write_iaxis( itemp )
    !-----------------------------------------------------------------------
    !!
    !! This routine writes to files results from the solutions of the Eliashberg equations
    !! on the imaginary-axis
    !!
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : fsthick, laniso, liso
    USE eliashbergcom, ONLY : nsiw, estemp, Agap, wsi, & 
                              NAZnormi, AZnormi, ADeltai, NZnormi, Znormi, & 
                              Deltai, nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV 
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter for temperature
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    CHARACTER(LEN = 256) :: name1, cname
    !
    temp = estemp(itemp) / kelvin2eV
    !
    cname = 'imag'
    !
    IF (laniso) THEN 
       !
       IF (temp < 10.d0) THEN
          WRITE(name1,'(a,a1,a4,a9,f4.2)') TRIM(prefix), '.', cname, '_aniso_00', temp
       ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
          WRITE(name1,'(a,a1,a4,a8,f5.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
       ELSEIF (temp >= 100.d0) THEN
          WRITE(name1,'(a,a1,a4,a7,f6.2)') TRIM(prefix), '.', cname, '_aniso_',temp
       ENDIF     
       OPEN(iufilgap, FILE = name1, FORM = 'formatted')
       WRITE(iufilgap,'(5a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Znorm(w)', 'Delta(w) [eV]', 'NZnorm(w)'
       DO iw = 1, nsiw(itemp) ! loop over omega
          DO ik = 1, nkfs
             DO ibnd = 1, nbndfs
                IF (ABS(ekfs(ibnd,ik) - ef0 ) < fsthick) THEN
                   WRITE(iufilgap,'(5ES20.10)') wsi(iw), ekfs(ibnd,ik)-ef0,&
                         AZnormi(ibnd,ik,iw), ADeltai(ibnd,ik,iw), NAZnormi(ibnd,ik,iw)
                   IF (iw == 1 ) Agap(ibnd,ik,itemp) = ADeltai(ibnd,ik,iw)
                ENDIF
             ENDDO ! ibnd                   
          ENDDO ! ik
       ENDDO ! iw
       CLOSE(iufilgap)
       !
       CALL gap_distribution_FS ( itemp, cname )
       !
       CALL gap_FS ( itemp )
       !
    ENDIF
    !
    ! isotropic case
    ! SP: Only write isotropic for laniso if user really wants that
    IF (( laniso .AND. iverbosity == 2 ) .OR. liso) THEN
       IF (temp < 10.d0) THEN
          WRITE(name1,'(a,a1,a4,a7,f4.2)') TRIM(prefix), '.', cname, '_iso_00', temp
       ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
          WRITE(name1,'(a,a1,a4,a6,f5.2)') TRIM(prefix), '.', cname, '_iso_0', temp
       ELSEIF (temp >= 100.d0) THEN
          WRITE(name1,'(a,a1,a4,a5,f6.2)') TRIM(prefix), '.', cname, '_iso_', temp
       ENDIF
       OPEN(iufilgap, FILE = name1, FORM = 'formatted')
       WRITE(iufilgap,'(4a20)') 'w [eV]', 'Znorm(w)', 'Delta(w) [eV]', 'NZnorm(w)'
       DO iw = 1, nsiw(itemp) ! loop over omega
          WRITE(iufilgap,'(4ES20.10)') wsi(iw), Znormi(iw), Deltai(iw), NZnormi(iw)
       ENDDO
       CLOSE(iufilgap)
    ENDIF 
    !
    RETURN
    !
    END SUBROUTINE eliashberg_write_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_write_raxis( itemp, cname )
    !-----------------------------------------------------------------------
    !
    !
    ! This routine writes to files results from the solutions of the Eliashberg
    ! equations on the real-axis 
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nqstep, fsthick, laniso, liso
    USE eliashbergcom, ONLY : nsw, estemp, ws, gap, Agap, &
                              Delta, Znorm, ADelta, AZnorm, &
                              nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter for temperature
    CHARACTER(len=256), INTENT(in) :: cname
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    CHARACTER(LEN = 256) :: name1
    LOGICAL :: lgap
    !! True if gap found
    !
    temp = estemp(itemp) / kelvin2eV
    !
    IF (laniso) THEN 
       IF (iverbosity == 2) THEN
          IF (temp < 10.d0) THEN
             WRITE(name1,'(a,a1,a4,a9,f4.2)') TRIM(prefix), '.', cname, '_aniso_00', temp
          ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a8,f5.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
          ELSEIF (temp >= 100.d0) THEN
             WRITE(name1,'(a,a1,a4,a7,f6.2)') TRIM(prefix), '.', cname, '_aniso_', temp
          ENDIF
          OPEN(iufilgap, FILE = name1, FORM = 'formatted')
          WRITE(iufilgap,'(6a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Re[Znorm(w)]', 'Im[Znorm(w)]',&
                                                            'Re[Delta(w)] [eV]', 'Im[Delta(w)] [eV]'
       ENDIF
       !
       DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
             IF (ABS(ekfs(ibnd,ik) - ef0 ) < fsthick) THEN
                lgap = .TRUE.
                ! DO iw = 1, nsw
                DO iw = 1, nsw-1   ! FG: this change is to prevent segfault in ws(iw+1) and ADelta(*,*,iw+1)
                   IF (lgap .AND. iw < nqstep .AND. REAL(ADelta(ibnd,ik,iw)) > 0.d0 &
                        .AND. REAL(ADelta(ibnd,ik,iw+1)) > 0.d0 &
                        .AND. ( ws(iw) - REAL(ADelta(ibnd,ik,iw)) )*( ws(iw+1) - REAL(ADelta(ibnd,ik,iw+1)) ) < 0.d0) THEN
                      Agap(ibnd,ik,itemp) = (   ( REAL(ADelta(ibnd,ik,iw))   - ws(iw)   ) * ws(iw+1) &
                                              - ( REAL(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) * ws(iw) ) &
                                          / ( ( REAL(ADelta(ibnd,ik,iw)) - ws(iw) ) - ( REAL(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) )
                      lgap = .FALSE.
                   ENDIF
                   IF (iverbosity == 2) THEN
                      WRITE(iufilgap,'(6ES20.10)') ws(iw), ekfs(ibnd,ik)-ef0, &
                                     REAL(AZnorm(ibnd,ik,iw)), aimag(AZnorm(ibnd,ik,iw)), &
                                     REAL(ADelta(ibnd,ik,iw)), aimag(ADelta(ibnd,ik,iw))
                   ENDIF
                ENDDO ! iw
                IF (lgap ) & 
                   Agap(ibnd,ik,itemp) = REAL(ADelta(ibnd,ik,1))
             ENDIF
          ENDDO ! ibnd
       ENDDO ! ik
       IF (iverbosity == 2 ) CLOSE(iufilgap)
       !
       CALL gap_distribution_FS ( itemp, cname )
       !
    ENDIF
    !
    ! isotropic case
    ! SP: Only write isotropic for laniso if user really wants that
    IF (( laniso .AND. iverbosity == 2 ) .OR. liso) THEN
       IF (temp < 10.d0) THEN
          WRITE(name1,'(a,a1,a4,a7,f4.2)') TRIM(prefix), '.', cname, '_iso_00', temp
       ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
          WRITE(name1,'(a,a1,a4,a6,f5.2)') TRIM(prefix), '.', cname, '_iso_0', temp
       ELSEIF (temp >= 100.d0) THEN
          WRITE(name1,'(a,a1,a4,a5,f6.2)') TRIM(prefix), '.', cname, '_iso_', temp
       ENDIF
       OPEN(iufilgap, FILE = name1, FORM = 'formatted')
       WRITE(iufilgap,'(5a20)') 'w [eV]', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)] [eV]', 'Im[Delta(w)] [eV]'
       lgap = .TRUE.
       ! DO iw = 1, nsw
       DO iw = 1, nsw-1   ! this change is to prevent segfault in Delta(iw+1) and ws(iw+1)
          IF (lgap .AND. iw < nqstep .AND. REAL(Delta(iw)) > 0.d0 .AND. REAL(Delta(iw+1)) > 0.d0 .AND. &
               ( ws(iw) - REAL(Delta(iw)) )*( ws(iw+1) - REAL(Delta(iw+1)) ) < 0.d0) THEN
             gap(itemp) = ( ( REAL(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( REAL(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                        / ( ( REAL(Delta(iw)) - ws(iw) ) - ( REAL(Delta(iw+1)) - ws(iw+1) ) )
             lgap = .FALSE.
          ENDIF
          WRITE(iufilgap,'(5ES20.10)') ws(iw), REAL(Znorm(iw)), aimag(Znorm(iw)), &
                                       REAL(Delta(iw)), aimag(Delta(iw))
       ENDDO ! iw
       CLOSE(iufilgap)
       IF (lgap ) & 
          gap(itemp) = REAL(Delta(1))
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE eliashberg_write_raxis
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_write_cont_raxis( itemp, cname )
    !-----------------------------------------------------------------------
    !
    !
    ! This routine writes to files results from the solutions of the Eliashberg
    ! equations on the real-axis 
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nqstep, fsthick, laniso, liso
    USE eliashbergcom, ONLY : nsw, estemp, ws, gap, Agap, &
                              Delta, Znorm, ADelta, AZnorm, &
                              nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV
    !
    IMPLICIT NONE
    !
    INTEGER :: iw, itemp, ik, ibnd
    REAL(KIND = DP) :: temp
    LOGICAL :: lgap
    CHARACTER(len=256) :: name1, cname
    !
    temp = estemp(itemp) / kelvin2eV
    !
    IF (laniso) THEN
       IF (iverbosity == 2) THEN
          IF (temp < 10.d0) THEN
             WRITE(name1,'(a,a1,a4,a9,f4.2)') TRIM(prefix), '.', cname, '_aniso_00', temp
          ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a8,f5.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
          ELSEIF (temp >= 100.d0) THEN
             WRITE(name1,'(a,a1,a4,a7,f6.2)') TRIM(prefix), '.', cname, '_aniso_', temp
          ENDIF
          OPEN(iufilgap, FILE = name1, FORM = 'formatted')
          WRITE(iufilgap,'(6a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Re[Znorm(w)]', 'Im[Znorm(w)]',&
                                                            'Re[Delta(w)] [eV]', 'Im[Delta(w)] [eV]'
       ENDIF
       !
       DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
             IF (ABS(ekfs(ibnd,ik) - ef0 ) < fsthick) THEN
                lgap = .TRUE.
                ! DO iw = 1, nsw
                DO iw = 1, nsw-1   ! FG: this change is to prevent segfault in ws(iw+1) and ADelta(*,*,iw+1)
                   IF (lgap .AND. iw < nqstep .AND. REAL(ADelta(ibnd,ik,iw)) > 0.d0 &
                        .AND. REAL(ADelta(ibnd,ik,iw+1)) > 0.d0 &
                        .AND. ( ws(iw) - REAL(ADelta(ibnd,ik,iw)) )*( ws(iw+1) - REAL(ADelta(ibnd,ik,iw+1)) ) < 0.d0) THEN
                      Agap(ibnd,ik,itemp) = (   ( REAL(ADelta(ibnd,ik,iw))   - ws(iw)   ) * ws(iw+1) &
                                              - ( REAL(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) * ws(iw) ) &
                                          / ( ( REAL(ADelta(ibnd,ik,iw)) - ws(iw) ) - ( REAL(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) )
                      lgap = .FALSE.
                   ENDIF
                   IF (iverbosity == 2) THEN
                      WRITE(iufilgap,'(6ES20.10)') ws(iw), ekfs(ibnd,ik)-ef0, &
                                     REAL(AZnorm(ibnd,ik,iw)), aimag(AZnorm(ibnd,ik,iw)), &
                                     REAL(ADelta(ibnd,ik,iw)), aimag(ADelta(ibnd,ik,iw))
                   ENDIF
                ENDDO ! iw
                IF (lgap ) &
                   Agap(ibnd,ik,itemp) = REAL(ADelta(ibnd,ik,1))
             ENDIF
          ENDDO ! ibnd
       ENDDO ! ik
       IF (iverbosity == 2 ) CLOSE(iufilgap)
       !
       CALL gap_distribution_FS ( itemp, cname )
       !
    ENDIF
    !
    ! isotropic case
    ! SP: Only write isotropic for laniso if user really wants that
    IF (( laniso .AND. iverbosity == 2 ) .OR. liso) THEN
       IF (temp < 10.d0) THEN
          WRITE(name1,'(a,a1,a4,a7,f4.2)') TRIM(prefix), '.', cname, '_iso_00', temp
       ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
          WRITE(name1,'(a,a1,a4,a6,f5.2)') TRIM(prefix), '.', cname, '_iso_0', temp
       ELSEIF (temp >= 100.d0) THEN
          WRITE(name1,'(a,a1,a4,a5,f6.2)') TRIM(prefix), '.', cname, '_iso_', temp
       ENDIF
       OPEN(iufilgap, FILE = name1, FORM = 'formatted')
       WRITE(iufilgap,'(5a20)') 'w [eV]', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)] [eV]', 'Im[Delta(w)] [eV]'
       lgap = .TRUE.
       ! DO iw = 1, nsw
       DO iw = 1, nsw-1   ! this change is to prevent segfault in Delta(iw+1) and ws(iw+1)
          IF (lgap .AND. iw < nqstep .AND. REAL(Delta(iw)) > 0.d0 .AND. REAL(Delta(iw+1)) > 0.d0 .AND. &
               ( ws(iw) - REAL(Delta(iw)) )*( ws(iw+1) - REAL(Delta(iw+1)) ) < 0.d0) THEN
             gap(itemp) = ( ( REAL(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( REAL(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                        / ( ( REAL(Delta(iw)) - ws(iw) ) - ( REAL(Delta(iw+1)) - ws(iw+1) ) )
             lgap = .FALSE.
          ENDIF
          WRITE(iufilgap,'(5ES20.10)') ws(iw), REAL(Znorm(iw)), aimag(Znorm(iw)), &
                                       REAL(Delta(iw)), aimag(Delta(iw))
       ENDDO ! iw
       CLOSE(iufilgap)
       IF (lgap ) &
          gap(itemp) = REAL(Delta(1))
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE eliashberg_write_cont_raxis
    !-----------------------------------------------------------------------
    SUBROUTINE read_a2f
    !-----------------------------------------------------------------------
    !!
    !! Read the eliashberg spectral function from fila2f
    !!
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nqstep, fila2f
    USE eliashbergcom, ONLY : wsphmax, wsph, a2f_iso, memlt_pool
    USE constants_epw, ONLY : zero
    USE mp_global,     ONLY : npool
    USE io_var,        ONLY : iua2ffil 
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iwph
    !! Counter for the number of freq
    INTEGER :: ios
    !! Status when opening a2F file
    !
    IF (.NOT. ALLOCATED(a2f_iso) ) ALLOCATE(a2f_iso(nqstep))
    IF (.NOT. ALLOCATED(wsph) ) ALLOCATE(wsph(nqstep)) 
    a2f_iso(:) = zero
    wsph(:) = zero
    !
    IF (mpime == ionode_id) THEN
      OPEN(iua2ffil, FILE = fila2f, STATUS = 'unknown', err=100, IOSTAT = ios)
100   CALL errore('read_a2f','opening file'//fila2f,ABS(ios))
    !
      DO iwph = 1, nqstep
         READ(iua2ffil,*) wsph(iwph), a2f_iso(iwph) ! freq from meV to eV
         wsph(iwph) = wsph(iwph) / 1000.d0
      ENDDO
      wsphmax = wsph(nqstep) 
      CLOSE(iua2ffil)
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( a2f_iso, ionode_id, inter_pool_comm )
    CALL mp_bcast( wsph, ionode_id, inter_pool_comm )
    CALL mp_bcast( wsphmax, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading a2f file '
    !
    IF (.NOT. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(npool))
    memlt_pool(:) = 0.d0
    !
    RETURN
    !
    END SUBROUTINE read_a2f
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_frequencies
    !-----------------------------------------------------------------------
    !
    ! read the frequencies obtained from a previous epw run
    !
    USE io_global, ONLY : stdout, ionode_id
    USE io_var,    ONLY : iufilfreq
    USE io_files,  ONLY : prefix, tmp_dir
    USE phcom,     ONLY : nmodes
    USE elph2,   ONLY : nqtotf, wf, wqf, xqf
    USE eliashbergcom, ONLY : wsphmax
    USE constants_epw, ONLY : ryd2ev, zero
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: imode
    !! Counter on modes
    CHARACTER(LEN = 256) :: filfreq
    !
    ! read frequencies from file
    IF (mpime == ionode_id) THEN
      filfreq = TRIM(tmp_dir) // TRIM(prefix) // '.freq'
      !OPEN(iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'formatted', err=100, IOSTAT = ios)
      OPEN(iufilfreq, FILE = filfreq, STATUS = 'unknown', FORM = 'unformatted', err=100, IOSTAT = ios)
100   CALL errore('read_frequencies','opening file '//filfreq,ABS(ios))
      !READ(iufilfreq,'(2i7)') nqtotf, nmodes
      READ(iufilfreq) nqtotf, nmodes
    ENDIF
    CALL mp_bcast( nqtotf, ionode_id, inter_pool_comm )
    CALL mp_bcast( nmodes, ionode_id, inter_pool_comm )
    !
    IF (.NOT. ALLOCATED(wf) )  ALLOCATE(wf(nmodes,nqtotf))
    IF (.NOT. ALLOCATED(wqf) ) ALLOCATE(wqf(nqtotf))
    IF (.NOT. ALLOCATED(xqf) ) ALLOCATE(xqf(3,nqtotf))
    wf(:, :) = zero
    wqf(:) = 1.d0 / DBLE(nqtotf)
    xqf(:, :) = zero
    !
    IF (mpime == ionode_id) THEN
      DO iq = 1, nqtotf ! loop over q-points
         !READ(iufilfreq,'(3f15.9)') xqf(1,iq), xqf(2,iq), xqf(3,iq)
         READ(iufilfreq) xqf(1,iq), xqf(2,iq), xqf(3,iq)
         DO imode = 1, nmodes
            !READ(iufilfreq,'(ES20.10)') wf(imode,iq)
            READ(iufilfreq) wf(imode,iq)
         ENDDO
      ENDDO 
      CLOSE(iufilfreq)
      ! go from Ryd to eV
      wf(:, :) = wf(:, :) * ryd2ev ! in eV
      wsphmax = 1.1d0 * MAXVAL( wf(:, :) ) ! increase by 10%
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( wf, ionode_id, inter_pool_comm )
    CALL mp_bcast( xqf, ionode_id, inter_pool_comm )
    CALL mp_bcast( wsphmax, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading .freq file '
    !
    RETURN
    !
    END SUBROUTINE read_frequencies
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_eigenvalues
    !-----------------------------------------------------------------------
    !!
    !! read the eigenvalues obtained from a previous epw run
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_files,      ONLY : prefix, tmp_dir
    USE pwcom,         ONLY : ef
    USE epwcom,        ONLY : nkf1, nkf2, nkf3, degaussw, fsthick, mp_mesh_k
    USE eliashbergcom, ONLY : nkfs, nbndfs, dosef, ef0, ekfs, wkfs, xkfs, w0g
    USE constants_epw, ONLY : ryd2ev, zero
    USE io_var,        ONLY : iufilegnv
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER ::  nkftot
    !! Number of k-points
    INTEGER ::  n, nbnd_
    !! Band indexes
    !
    REAL(KIND = DP), ALLOCATABLE :: ekf_(:, :)
    !! Temporary eigenvalues on the k point grid
    !
    CHARACTER(LEN = 256) :: filegnv
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !
    IF (mpime == ionode_id) THEN
      !
      ! SP: Needs to be initialized
      nbnd_ = 0 
      nkfs = 0
      !
      ! read eigenvalues on the irreducible fine k-mesh
      !  
      filegnv = TRIM(tmp_dir) // TRIM(prefix) // '.egnv'
      !OPEN(iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'formatted', err=100, IOSTAT = ios)
      OPEN(iufilegnv, FILE = filegnv, STATUS = 'unknown', FORM = 'unformatted', err=100, IOSTAT = ios)
100   CALL errore('read_eigenvalues','opening file '//filegnv,ABS(ios))
      !
      !READ(iufilegnv,'(5i7)') nkftot, nkf1, nkf2, nkf3, nkfs 
      !READ(iufilegnv,'(i7,5ES20.10)') nbnd_, ef, ef0, dosef, degaussw, fsthick
      READ(iufilegnv) nkftot, nkf1, nkf2, nkf3, nkfs
      READ(iufilegnv) nbnd_, ef, ef0, dosef, degaussw, fsthick
      degaussw = degaussw * ryd2ev
      ef0 = ef0 * ryd2ev
      ef = ef * ryd2ev
      fsthick = fsthick * ryd2ev
      dosef = dosef / ryd2ev
      WRITE(stdout,'(5x,a32,ES20.10)') 'Fermi level (eV) = ', ef0
      WRITE(stdout,'(5x,a32,ES20.10)') 'DOS(states/spin/eV/Unit Cell) = ', dosef
      WRITE(stdout,'(5x,a32,ES20.10)') 'Electron smearing (eV) = ', degaussw
      WRITE(stdout,'(5x,a32,ES20.10)') 'Fermi window (eV) = ', fsthick
      IF (mp_mesh_k) THEN 
         WRITE(stdout,'(5x,a,i9,a,i9)') 'Nr irreducible k-points within the Fermi shell = ', nkfs, ' out of ', nkftot
      ELSE
         WRITE(stdout,'(5x,a,i9,a,i9)') 'Nr k-points within the Fermi shell = ', nkfs, ' out of ', nkftot
      ENDIF
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( nkf1, ionode_id, inter_pool_comm )
    CALL mp_bcast( nkf2, ionode_id, inter_pool_comm )
    CALL mp_bcast( nkf3, ionode_id, inter_pool_comm )
    CALL mp_bcast( nkfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( degaussw, ionode_id, inter_pool_comm )
    CALL mp_bcast( ef0, ionode_id, inter_pool_comm )
    CALL mp_bcast( dosef, ionode_id, inter_pool_comm )
    CALL mp_bcast( fsthick, ionode_id, inter_pool_comm )
    CALL mp_bcast( ef, ionode_id, inter_pool_comm )
    !
    IF (.NOT. ALLOCATED(wkfs) ) ALLOCATE(wkfs(nkfs))
    IF (.NOT. ALLOCATED(xkfs) ) ALLOCATE(xkfs(3,nkfs))
    wkfs(:) = zero
    xkfs(:, :) = zero
    !
    IF (mpime == ionode_id) THEN
      !
      ! at each k-point keep only the bands within the Fermi shell
      !
      ALLOCATE(ekf_(nbnd_,nkfs))
      ekf_(:, :) = zero
      !
      ! nbndfs - nr of bands within the Fermi shell
      !
      nbndfs = 0
      DO ik = 1, nkfs ! loop over irreducible k-points
         !READ(iufilegnv,'(4f15.9)') wkfs(ik), xkfs(1,ik), xkfs(2,ik), xkfs(3,ik)
         READ(iufilegnv) wkfs(ik), xkfs(1,ik), xkfs(2,ik), xkfs(3,ik)
         DO ibnd = 1, nbnd_
            !READ(iufilegnv,'(ES20.10)') ekf_(ibnd,ik)
            READ(iufilegnv) ekf_(ibnd,ik)
         ENDDO
         n = 0
         DO ibnd = 1, nbnd_
            ! go from Ryd to eV
            ekf_(ibnd,ik) = ekf_(ibnd,ik) * ryd2ev
            IF (ABS(ekf_(ibnd,ik) - ef0 ) < fsthick) THEN
               n = n + 1
               IF (nbndfs < n ) nbndfs = n
            ENDIF
         ENDDO
      ENDDO
      WRITE(stdout,'(5x,i7,a/)') nbndfs, ' bands within the Fermi window'
      CLOSE(iufilegnv)
      ! 
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( nbndfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( wkfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( xkfs, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF (.NOT. ALLOCATED(ekfs) ) ALLOCATE(ekfs(nbndfs,nkfs))
    IF (.NOT. ALLOCATED(w0g) )  ALLOCATE(w0g(nbndfs,nkfs))
    ! sanity choice
    ekfs(:, :) = ef0 - 10.d0 * fsthick
    w0g(:, :) = zero
    IF (mpime == ionode_id) THEN
      DO ik = 1, nkfs ! loop over k-points
         n = 0
         DO ibnd = 1, nbnd_
            IF (ABS(ekf_(ibnd,ik) - ef0 ) < fsthick) THEN
               n = n + 1
               ekfs(n,ik) = ekf_(ibnd,ik)
               w0g(n,ik) = w0gauss( ( ekfs(n,ik) - ef0 ) / degaussw, 0 ) / degaussw
            ENDIF
         ENDDO
      ENDDO
      DEALLOCATE(ekf_)
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( ekfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( w0g, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading .egnv file '
    !
    RETURN
    !
    END SUBROUTINE read_eigenvalues
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_kqmap
    !-----------------------------------------------------------------------
    !
    ! read the map index of k+(sign)q on the k-mesh
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, ionode_id
    USE io_var,    ONLY : iufilikmap
    USE io_files,  ONLY : prefix, tmp_dir
    USE symm_base, ONLY : t_rev, time_reversal, s, set_sym_bl
    USE phcom,     ONLY : nmodes
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, mp_mesh_k
    USE elph2,     ONLY : nqtotf, xqf
    USE eliashbergcom, ONLY : ixkff, xkff, ixkf, xkfs, nkfs, ixkqf, ixqfs, nbndfs, nqfs, memlt_pool
    USE superconductivity, ONLY : mem_size_eliashberg, mem_integer_size_eliashberg
    USE constants_epw, ONLY : eps5, zero
    USE symm_base, ONLY : nrot
    USE mp_global, ONLY : inter_pool_comm, npool
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,  ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    INTEGER :: i, j, k, ik, nk, n
    !! Counter on k points
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: nkftot
    !! Total number of k points
    INTEGER :: nkf_mesh
    !! Nr. of k points read from .ikmap file
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: nks
    !! Number of non-equivalent k points
    INTEGER :: ns
    !! Counter on rotation operations
    INTEGER :: ios
    !! Integer variable for I/O control
    INTEGER :: imelt
    !! Memory allocated
    INTEGER, ALLOCATABLE :: equiv_(:)
    !! Index of equivalence of k points
    INTEGER, ALLOCATABLE :: index_(:, :)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    !
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3)
    !! coordinates of q points
    REAL(KIND = DP) :: xkr(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xx, yy, zz
    !! Temporary variables
    !
    LOGICAL :: in_the_list
    !! Check if k point is in the list
    CHARACTER(LEN = 256) :: filikmap
    !! Name of the file
    !
    IF (.NOT. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(npool))
    memlt_pool(:) = zero
    !
    ! get the size of arrays for frequency and eigenvalue variables allocated in 
    ! read_frequencies and read_eigenvalues
    imelt = ( nmodes + 4 ) * nqtotf + ( 4 + 2 * nbndfs ) * nkfs
    CALL mem_size_eliashberg( imelt )
    !
    nkftot = nkf1 * nkf2 * nkf3
    !
    ! get the size of required memory for ixkff  
    imelt = nkftot
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF (.NOT. ALLOCATED(ixkff) ) ALLOCATE(ixkff(nkftot))
    ixkff(:) = 0
    !
    IF (mpime == ionode_id) THEN
      !
      filikmap = TRIM(tmp_dir) // TRIM(prefix) // '.ikmap'
      !OPEN(iufilikmap, FILE = filikmap, STATUS = 'old', FORM = 'formatted', err=100, IOSTAT = ios)
      OPEN(iufilikmap, FILE = filikmap, STATUS = 'old', FORM = 'unformatted', err=100, IOSTAT = ios)
100   CALL errore('read_kqmap','opening file '//filikmap,ABS(ios))
      !
      ! nkf_mesh - Total number of k points
      !          - These are irreducible k-points if mp_mesh_k = .TRUE.
      READ(iufilikmap) nkf_mesh
      !
      IF (.NOT. ALLOCATED(ixkf) ) ALLOCATE(ixkf(nkf_mesh))
      ixkf(:) = 0
      !
      DO ik = 1, nkf_mesh
         !READ(iufilikmap,'(i9)') ixkf(ik)
         READ(iufilikmap) ixkf(ik)
      ENDDO
      CLOSE(iufilikmap)
      !
      IF (.NOT. ALLOCATED(xkff) )  ALLOCATE(xkff(3,nkftot))
      xkff(:, :) = zero
      !
      DO i = 1, nkf1
         DO j = 1, nkf2
            DO k = 1, nkf3
               ik = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + k
               xkff(1,ik) = DBLE(i-1) / DBLE(nkf1)
               xkff(2,ik) = DBLE(j-1) / DBLE(nkf2)
               xkff(3,ik) = DBLE(k-1) / DBLE(nkf3)
            ENDDO
         ENDDO
      ENDDO
      !
      IF (.NOT. ALLOCATED(equiv_) )  ALLOCATE(equiv_(nkftot))
      !  equiv_(nk) =nk : k-point nk is not equivalent to any previous k-point
      !  equiv_(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
      !
      DO nk = 1, nkftot
         equiv_(nk) = nk
      ENDDO
      !
      IF (mp_mesh_k) THEN 
         CALL set_sym_bl( ) 
         DO nk = 1, nkftot
            !  check if this k-point has already been found equivalent to another
            IF (equiv_(nk) == nk) THEN
               !  check if there are equivalent k-point to this in the list
               !  (excepted those previously found to be equivalent to another)
               !  check both k and -k
               DO ns = 1, nrot
                  DO i = 1, 3
                     xkr(i) = s(i,1,ns) * xkff(1,nk) &
                            + s(i,2,ns) * xkff(2,nk) &
                            + s(i,3,ns) * xkff(3,nk)
                     xkr(i) = xkr(i) - NINT( xkr(i) )
                  ENDDO
                  IF (t_rev(ns) == 1 ) xkr = -xkr
                  xx = xkr(1)*nkf1
                  yy = xkr(2)*nkf2
                  zz = xkr(3)*nkf3
                  in_the_list = ABS(xx-NINT(xx) ) <= eps5 .AND. &
                                ABS(yy-NINT(yy) ) <= eps5 .AND. &
                                ABS(zz-NINT(zz) ) <= eps5
                  IF (in_the_list) THEN
                     i = MOD( NINT( xkr(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                     j = MOD( NINT( xkr(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                     k = MOD( NINT( xkr(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                     n = (k-1) + (j-1)*nkf3 + (i-1)*nkf2*nkf3 + 1
                     IF (n > nk .AND. equiv_(n) == n) THEN
                        equiv_(n) = nk
                     ELSE
                        IF (equiv_(n) /= nk .OR. n < nk ) CALL errore('kmesh_fine', &
                           'something wrong in the checking algorithm',1)
                     ENDIF
                  ENDIF
                  IF (time_reversal) THEN
                     xx = -xkr(1)*nkf1
                     yy = -xkr(2)*nkf2
                     zz = -xkr(3)*nkf3
                     in_the_list = ABS(xx-NINT(xx) ) <= eps5 .AND. &
                                   ABS(yy-NINT(yy) ) <= eps5 .AND. &
                                   ABS(zz-NINT(zz) ) <= eps5
                     IF (in_the_list) THEN
                        i = MOD( NINT( -xkr(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                        j = MOD( NINT( -xkr(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                        k = MOD( NINT( -xkr(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                        n = (k-1) + (j-1)*nkf3 + (i-1)*nkf2*nkf3 + 1
                        IF (n > nk .AND. equiv_(n) == n) THEN
                           equiv_(n) = nk
                        ELSE
                           IF (equiv_(n) /= nk .OR. n < nk ) CALL errore('kmesh_fine', &
                              'something wrong in the checking algorithm',2)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      !
      !  define index of k on the full mesh (ixkff) using index of k-point within the
      !  Fermi shell (ixkf)
      !
      nks = 0
      DO nk = 1, nkftot
         IF (equiv_(nk) == nk) THEN
            nks = nks + 1
            ixkff(nk) = ixkf(nks)
         ELSE
            ixkff(nk) = ixkff(equiv_(nk))
         ENDIF
      ENDDO
      IF (nks /= nkf_mesh) CALL errore('read_kmap_mp', 'something wrong with the mesh',1)
      !
      IF (ALLOCATED(equiv_) ) DEALLOCATE(equiv_)
      IF (ALLOCATED(ixkf) )   DEALLOCATE(ixkf)
      !
    ENDIF
    CALL mp_bcast( ixkff, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get the size of required memory for ixkqf, nqfs, index_
    imelt = ( nqtotf + 1 ) * nkfs + ( upper_bnd - lower_bnd + 1 ) * nqtotf
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF (.NOT. ALLOCATED(ixkqf) ) ALLOCATE(ixkqf(nkfs,nqtotf))
    IF (.NOT. ALLOCATED(nqfs) )  ALLOCATE(nqfs(nkfs))
    IF (.NOT. ALLOCATED(index_) ) ALLOCATE(index_(lower_bnd:upper_bnd,nqtotf))
    ixkqf(:, :) = 0
    nqfs(:) = 0
    index_(:, :) = 0
    !
    !
    ! find the index of k+sign*q on the fine k-mesh
    ! nkfs - total nr. of k-points within the Fermi shell (fine mesh)
    !      - these are irreducible k-points if mp_mesh_k=.TRUE.
    ! nqtotf - total nr of q-points on the fine mesh
    !
    DO ik = lower_bnd, upper_bnd
       DO iq = 1, nqtotf
          xk(:) = xkfs(:,ik)
          xq(:) = xqf(:,iq)
          !
          !  nkq - index of k+sign*q on the full fine k-mesh.
          !
          CALL kpmq_map( xk, xq, +1, nkq )
          !
          !  ixkqf(ik,iq) - index of k+sign*q on the fine k-mesh within the Fermi shell
          !
          ixkqf(ik,iq) = ixkff(nkq)
          !
          ! nqfs(ik) - nr of q-points at each k-point for which k+sign*q is within the Fermi shell 
          ! index_   - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell
          !
          IF (ixkqf(ik,iq) > 0) THEN
             nqfs(ik) = nqfs(ik) + 1
             index_(ik,nqfs(ik)) = iq
          ENDIF
       ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( ixkqf, inter_pool_comm )
    CALL mp_sum( nqfs,  inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    ! get the size of required memory for ixqfs
    imelt = nkfs * MAXVAL(nqfs(:))
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF (.NOT. ALLOCATED(ixqfs) ) ALLOCATE(ixqfs(nkfs,MAXVAL(nqfs(:))))
    ixqfs(:, :) = 0
    !
    DO ik = lower_bnd, upper_bnd
       DO iq = 1, nqfs(ik)
          !
          ! ixqfs - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell 
          !
          ixqfs(ik,iq) = index_(ik,iq)
       ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( ixqfs, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF (ALLOCATED(index_) ) DEALLOCATE(index_)
    IF (ALLOCATED(xqf) )    DEALLOCATE(xqf)
    !
    ! remove memory allocated for index_
    imelt = nqtotf * ( upper_bnd - lower_bnd + 1 ) 
    CALL mem_integer_size_eliashberg( -imelt )
    !
    ! remove memory allocated for xqf
    imelt = 3 * nqtotf
    CALL mem_size_eliashberg( -imelt )
    !
    WRITE(stdout,'(/5x,a,i9/)') 'Max nr of q-points = ', MAXVAL(nqfs(:))  
    WRITE(stdout,'(/5x,a/)') 'Finish reading .ikmap files'
    !
    RETURN
    !
    END SUBROUTINE read_kqmap
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_ephmat
    !-----------------------------------------------------------------------
    !!
    !! Read the electron-phonon matrix elements 
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_var,        ONLY : iufileph
    USE io_files,      ONLY : prefix, tmp_dir
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : nqtotf, wf
    USE epwcom,        ONLY : eps_acustic, fsthick
    USE eliashbergcom, ONLY : nkfs, nbndfs, ef0, ekfs, g2, ixkqf, nqfs
    USE superconductivity, ONLY : mem_size_eliashberg
    USE constants_epw, ONLY : ryd2ev, zero
    USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, npool
    USE division,      ONLY : fkbounds
    USE low_lvl,       ONLY : set_ndnmbr
    !  
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! INTEGER variable for I/O control
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: nnk
    !! Number of k-points within the Fermi shell
    INTEGER :: nnq(nkfs)
    !! Number of k+q points within the Fermi shell for a given k-point
    INTEGER :: ipool
    !! Counter on pools
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: tmp_pool_id
    !! Pool index read from file
    INTEGER :: nkpool(npool)
    !! nkpool(ipool) - sum of nr. of k points from pool 1 to pool ipool
    INTEGER :: nmin
    !! Lower bound index for .ephmat file read in current pool
    INTEGER :: nmax
    !! Lower bound index for .ephmat file read in current pool
    INTEGER :: nks
    !! Counter on k points within the Fermi shell
    INTEGER :: imelt
    !! Memory allocated
    !
    REAL(KIND = DP) :: gmat
    !! Electron-phonon matrix element square
    !
    CHARACTER(LEN = 256) :: filephmat
    CHARACTER (len=3) :: filelab
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get the size of the e-ph matrices that need to be stored in each pool
    imelt = ( upper_bnd - lower_bnd + 1 ) * MAXVAL(nqfs(:)) * nbndfs**2 * nmodes
    CALL mem_size_eliashberg( imelt ) 
    !
    IF (.NOT. ALLOCATED(g2) ) ALLOCATE(g2(lower_bnd:upper_bnd,MAXVAL(nqfs(:)),nbndfs,nbndfs,nmodes))
    g2(:, :, :, :, :) = zero
    !
    ! go from Ryd to eV
    ! eps_acustic is given in units of cm-1 in the input file and converted to Ryd in epw_readin
    eps_acustic = eps_acustic * ryd2ev
    !
    WRITE(stdout,'(/5x,a/)') 'Start reading .ephmat files'
    !
    DO ipool = 1, npool ! nr of pools 
       CALL set_ndnmbr(0,ipool,1,npool,filelab)
#if defined(__MPI)
       filephmat = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat' // filelab
#else
       filephmat = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
#endif
       !OPEN(iufileph, FILE = filephmat, STATUS = 'old', FORM = 'formatted', err=100, IOSTAT = ios)
       OPEN(iufileph, FILE = filephmat, STATUS = 'old', FORM = 'unformatted', err=100, IOSTAT = ios)
100 CALL errore('read_ephmat','opening file '//filephmat,ABS(ios))
       !READ(iufileph,'(2i7)') tmp_pool_id, nkpool(ipool)
       READ(iufileph) tmp_pool_id, nkpool(ipool)
       IF (ipool /= tmp_pool_id )  CALL errore('read_ephmat', &
           'npool should be equal to the number of .ephmat files',1)
       IF (ipool > 1 ) & 
          nkpool(ipool) = nkpool(ipool) + nkpool(ipool-1)
       !WRITE(stdout,'(2i7)') tmp_pool_id, nkpool(ipool)
       CLOSE(iufileph)
    ENDDO
    CALL mp_barrier(inter_pool_comm)
    !
    ! since the nkfs k-points within the Fermi shell are not evenly distrubed
    ! among the .ephmat files, we re-distribute them here among the npool-pools
    nmin = npool
    nmax = npool
    DO ipool = npool, 1, -1
       IF (lower_bnd <= nkpool(ipool)) THEN
          nmin = ipool
       ENDIF
       IF (upper_bnd <= nkpool(ipool)) THEN
          nmax = ipool
       ENDIF
    ENDDO
    !
    nnk = 0
    nnq(:) = 0
    DO ipool = 1, npool ! nr of pools 
       CALL set_ndnmbr(0,ipool,1,npool,filelab)
#if defined(__MPI)
       filephmat = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat' // filelab
#else
       filephmat = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
#endif     
       OPEN(iufileph, FILE = filephmat, STATUS = 'old', FORM = 'unformatted')
       READ(iufileph) tmp_pool_id, nks
       IF (ipool >= nmin .AND. ipool <= nmax) THEN
          DO iq = 1, nqtotf ! loop over q-points 
             DO ik = 1, nks ! loop over k-points in the pool
                IF (ixkqf(ik+nnk,iq) > 0) THEN 
                   nnq(ik+nnk) = nnq(ik+nnk) + 1
                   DO imode = 1, nmodes ! loop over phonon modes
                      DO ibnd = 1, nbndfs ! loop over iband's 
                         IF (ABS(ekfs(ibnd,ik+nnk) - ef0 ) < fsthick) THEN
                            DO jbnd = 1, nbndfs ! loop over jband's 
                               IF (ABS(ekfs(jbnd,ixkqf(ik+nnk,iq)) - ef0 ) < fsthick) THEN
                                  !READ(iufileph,'(ES20.10)') gmat
                                  READ(iufileph) gmat
                                  IF (ik+nnk >= lower_bnd .AND. ik+nnk <= upper_bnd) THEN
                                     ! go from Ryd to eV
                                     IF (wf(imode,iq) > eps_acustic) THEN
                                        g2(ik+nnk,nnq(ik+nnk),ibnd,jbnd,imode) = gmat * ryd2ev * ryd2ev
                                     ELSE
                                        g2(ik+nnk,nnq(ik+nnk),ibnd,jbnd,imode) = zero
                                     ENDIF
                                  ENDIF
                               ENDIF ! ekq
                            ENDDO ! jbnd
                         ENDIF ! ekk
                      ENDDO ! ibnd
                   ENDDO ! imode
                ENDIF ! ekk and ekq
             ENDDO ! ik
          ENDDO ! iq
          CLOSE(iufileph)
       ENDIF ! ipool
       nnk = nnk + nks
       IF (ipool == npool .AND. nnk /= nkfs )  CALL errore('read_ephmat', &
           'nnk should be equal to nkfs',1)
    ENDDO ! ipool
    !
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading .ephmat files '
    !
    RETURN
    !
    END SUBROUTINE read_ephmat
    !
    !-----------------------------------------------------------------------
    SUBROUTINE write_ephmat(iq)
    !-----------------------------------------------------------------------
    !!
    !!  This SUBROUTINE writes the elph matrix elements in a format required 
    !!  by Eliashberg equations
    !! 
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,      ONLY : DP
    USE io_global,  ONLY : stdout
    USE io_var,     ONLY : iufilfreq, iufilegnv, iufileph
    USE io_files,   ONLY : prefix, tmp_dir
    USE phcom,      ONLY : nmodes
    USE epwcom,     ONLY : nbndsub, fsthick, ngaussw, degaussw, shortrange, & 
                           nkf1, nkf2, nkf3, efermi_read, fermi_energy
    USE pwcom,      ONLY : ef 
    USE elph2,      ONLY : etf, nkqf, epf17, wkf, nkf, &
                           nqtotf, wf, xqf, efnew, nbndfst, nktotf 
    USE eliashbergcom, ONLY : nkfs, ekfs, wkfs, xkfs, dosef, ixkf, ixkqf, nbndfs
    USE superconductivity, ONLY : mem_size_eliashberg, mem_integer_size_eliashberg
    USE constants_epw, ONLY : ryd2ev, ryd2mev, two, eps8
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE division,      ONLY : fkbounds
    USE low_lvl,       ONLY : set_ndnmbr
    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT(in) :: iq
    !! Current q-points
    !
    ! Local variables
    !
    CHARACTER(LEN = 256) :: filfreq
    !! 
    CHARACTER(LEN = 256) :: filegnv
    !! 
    CHARACTER(LEN = 256) :: filephmat
    !! 
    CHARACTER(LEN = 3) :: filelab
    !! 
    INTEGER :: ik
    !! Counter on the k-point index 
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index 
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: nkftot
    !! Total number of k+q points 
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: nks
    !! Number of k-point on the current pool
    INTEGER :: imelt
    !! Memory allocated
    !
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: wq
    !! phonon freq on the fine grid
    REAL(KIND = DP) :: inv_wq
    !! 
    REAL(KIND = DP):: g2
    !! Electron-phonon matrix element square
    REAL(KIND = DP), EXTERNAL :: dos_ef
    !! Function to compute the density of states at the Fermi level
    REAL(KIND = DP), EXTERNAL :: efermig
    !! Return the fermi energy
    !
    ! write phonon frequencies to file
    IF (my_pool_id == 0) THEN
      filfreq = TRIM(tmp_dir) // TRIM(prefix) // '.freq'
      IF (iq == 1) THEN
        OPEN(iufilfreq, FILE = filfreq, FORM = 'unformatted')
        WRITE(iufilfreq) nqtotf, nmodes
      ELSE
        OPEN(iufilfreq, FILE = filfreq, POSITION = 'append', FORM = 'unformatted')
      ENDIF
      WRITE(iufilfreq) xqf(1, iq), xqf(2, iq), xqf(3, iq)
      DO imode = 1, nmodes
        WRITE(iufilfreq) wf(imode, iq)
      ENDDO
      CLOSE(iufilfreq)
    ENDIF
    ! 
    ! Fermi level and corresponding DOS
    !  
    ! since wkf(:,ikq) = 0 these bands do not bring any contribution to ef0 or dosef
    ! 
    IF (efermi_read) THEN
      ef0 = fermi_energy 
    ELSE
      ef0 = efnew 
      !ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
      ! if some bands are skipped (nbndskip /= 0), nelec has already been recalculated 
      ! in ephwann_shuffle
    ENDIF
    !     
    dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
    ! N(Ef) in the equation for lambda is the DOS per spin
    dosef = dosef / two
    !
    ! find the bounds of k-dependent arrays in the parallel case
    nkftot = nktotf 
    CALL fkbounds(nkftot, lower_bnd, upper_bnd)
    !
    IF (iq == 1) THEN
      !
      ! find fermicount - nr of k-points within the Fermi shell per pool
      ! for mp_mesh_k=true. femicount is the nr of irreducible k-points within the Fermi shell per pool
      ! 
      fermicount = 0
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        IF (MINVAL(ABS(etf(:,ikk) - ef)) < fsthick) THEN
          fermicount = fermicount + 1
        ENDIF
        !
      ENDDO
      !
      ! nks = irr nr of k-points within the Fermi shell (fine mesh)
      nks = fermicount
      !
      ! collect contributions from all pools (sum over k-points)
      CALL mp_sum(nks, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      ! write eigenvalues to file
      IF (my_pool_id == 0) THEN
        filegnv = TRIM(tmp_dir) // TRIM(prefix) // '.egnv'
        OPEN(iufilegnv, FILE = filegnv, FORM = 'unformatted')
        IF (nks /= nkfs ) CALL errore('write_ephmat', &
          'nks should be equal to nr. of irreducible k-points within the Fermi shell on the fine mesh',1)
        WRITE(iufilegnv) nkftot, nkf1, nkf2, nkf3, nks
        WRITE(iufilegnv) nbndfst, ef, ef0, dosef, degaussw, fsthick
        DO ik = 1, nks
          WRITE(iufilegnv) wkfs(ik), xkfs(1, ik), xkfs(2, ik), xkfs(3, ik) 
          DO ibnd = 1, nbndfst
            WRITE(iufilegnv) ekfs(ibnd, ik)
          ENDDO
        ENDDO
        CLOSE(iufilegnv)
      ENDIF
      !
    ENDIF ! iq
    !
    ! write the e-ph matrix elements in the Bloch representation on the fine mesh
    ! in .ephmat files (one for each pool)
    !
#if defined(__MPI)
    CALL set_ndnmbr(0, my_pool_id + 1, 1, npool, filelab)
    filephmat = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat' // filelab
#else
    filephmat = TRIM(tmp_dir) // TRIM(prefix) // '.ephmat'
#endif
    IF (iq == 1) THEN 
       !OPEN(iufileph, file = filephmat, form = 'formatted')
       OPEN(iufileph, file = filephmat, form = 'unformatted')
    ELSE
       !OPEN(iufileph, file = filephmat, position='append', form = 'formatted')
       OPEN(iufileph, file = filephmat, position='append', form = 'unformatted')
    ENDIF
    !
    !IF (iq == 1 ) WRITE(iufileph,'(2i7)') my_pool_id+1, fermicount
    IF (iq == 1 ) WRITE(iufileph) my_pool_id+1, fermicount
    !
    ! nkf - nr of k-points in the pool (fine mesh)
    ! for mp_mesh_k = true nkf is nr of irreducible k-points in the pool 
    !
    DO ik = 1, nkf
      !  
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      ! go only over irreducible k-points
      !
      !
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      !
      !   IF (ixkf(lower_bnd+ik-1) > 0 .AND. ixkqf(ixkf(lower_bnd+ik-1),iq) > 0) THEN
      ! FG: here it can happen that ixkf is 0 and this leads to ixqf(0,iq) after .AND.
      !     modified to prevent crash
      IF (ixkf(lower_bnd+ik-1) > 0) THEN
        IF (ixkqf(ixkf(lower_bnd+ik-1),iq) > 0) THEN
          !
          ! 
          DO imode = 1, nmodes ! phonon modes
            wq = wf(imode, iq)
            inv_wq =  1.0/(two * wq) 
            !
            DO ibnd = 1, nbndfst
              IF (ABS(ekfs(ibnd,ixkf(lower_bnd+ik-1)) - ef0 ) < fsthick) THEN
                DO jbnd = 1, nbndfst
                  IF (ABS(ekfs(jbnd,ixkqf(ixkf(lower_bnd+ik-1),iq)) - ef0 ) < fsthick) THEN
                    !
                    ! here we take into account the zero-point DSQRT(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    IF (shortrange .AND. ( ABS(xqf (1, iq))> eps8 .OR. ABS(xqf (2, iq))> eps8 &
                         .OR. ABS(xqf (3, iq))> eps8 )) THEN
                      ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                      !     number, in which case its square will be a negative number.
                      g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two) *inv_wq )
                    ELSE
                      g2 = ABS(epf17(jbnd, ibnd, imode, ik) )**two * inv_wq 
                    ENDIF
                    !WRITE(iufileph,'(ES20.10)') g2
                    WRITE(iufileph) g2
                  ENDIF
                ENDDO ! jbnd
              ENDIF
            ENDDO ! ibnd
          ENDDO ! imode
          !
        ENDIF
      ENDIF ! fsthick
      !
    ENDDO ! ik's
    CLOSE(iufileph)
    !
    IF (iq == nqtotf) THEN 
       IF (ALLOCATED(ekfs) )   DEALLOCATE(ekfs)
       IF (ALLOCATED(wkfs) )   DEALLOCATE(wkfs)
       IF (ALLOCATED(xkfs) )   DEALLOCATE(xkfs)
       IF (ALLOCATED(ixkqf) )  DEALLOCATE(ixkqf)
       IF (ALLOCATED(ixkf) )   DEALLOCATE(ixkf)
       !
       ! remove memory allocated for ekfs, wkfs, xkfs 
       imelt = ( nbndfs + 4 ) * nkfs
       CALL mem_size_eliashberg( -imelt )
       !
       ! remove memory allocated for ixkqf 
       imelt = nqtotf * nkfs
       CALL mem_integer_size_eliashberg( -imelt )
       !
       ! remove memory allocated for ixkf
       imelt = nkftot
       CALL mem_integer_size_eliashberg( -imelt )
       !
       WRITE(stdout,'(5x,a32,d24.15)') 'Fermi level (eV) = ', ef0 * ryd2ev
       WRITE(stdout,'(5x,a32,d24.15)') 'DOS(states/spin/eV/Unit Cell) = ', dosef / ryd2ev
       WRITE(stdout,'(5x,a32,d24.15)') 'Electron smearing (eV) = ', degaussw * ryd2ev
       WRITE(stdout,'(5x,a32,d24.15)') 'Fermi window (eV) = ', fsthick * ryd2ev
       WRITE(stdout,'(5x,a)')          ' '
       WRITE(stdout,'(5x,a)')          'Finished writing .ephmat files'
       !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE write_ephmat
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE count_kpoints(iq)
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE epwcom,    ONLY : nbndsub, fsthick, ngaussw, degaussw, & 
                          efermi_read, fermi_energy, mp_mesh_k
    USE pwcom,     ONLY : nelec, ef
    USE klist_epw, ONLY : isk_dummy
    USE elph2,     ONLY : etf, nkqf, wkf, nkf, nktotf
    USE constants_epw, ONLY : two
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iq
    !! Current q-points
    !
    ! Local variables
    !
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: nks
    !! Number of k-point on the current pool
    !
    REAL(KIND = DP) :: ef0
    !! Fermi energy level
    REAL(KIND = DP) :: dosef
    !! density of states at the Fermi level
    !
    REAL(KIND = DP), EXTERNAL :: efermig, dos_ef
    ! 
    IF (iq == 1) THEN
       ! 
       ! Fermi level and corresponding DOS
       !  
       ! since wkf(:,ikq) = 0 these bands do not bring any contribution to ef0 or dosef
       !
       IF (efermi_read) THEN
          ef0 = fermi_energy 
       ELSE
          ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk_dummy)
       ENDIF  
       !     
       dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
       ! N(Ef) in the equation for lambda is the DOS per spin
       dosef = dosef / two
       !
       ! fermicount = nr of k-points within the Fermi shell per pool
       !
       fermicount = 0
       DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF (minval( ABS(etf(:,ikk) - ef  ) ) < fsthick ) &
             fermicount = fermicount + 1 
          !
       ENDDO
       !
       ! nks =  nr of k-points within the Fermi shell (fine mesh)
       nks = fermicount
       !
       ! collect contributions from all pools (sum over k-points)
       CALL mp_sum( nks, inter_pool_comm )
       CALL mp_barrier(inter_pool_comm)
       !
       IF (mp_mesh_k) THEN
          WRITE(stdout,'(5x,a,i9,a,i9)') 'Nr irreducible k-points within the Fermi shell = ', nks, ' out of ', nktotf
       ELSE
          WRITE(stdout,'(5x,a,i9,a,i9)') 'Nr k-points within the Fermi shell = ', nks, ' out of ', nktotf
       ENDIF
    ENDIF ! iq
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE count_kpoints
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE kmesh_fine
    !-----------------------------------------------------------------------
    !!
    !!   This routine defines the nr. of k-points on the fine k-mesh 
    !!   within the Fermi shell
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : ionode_id, stdout
    USE io_files,  ONLY : prefix, tmp_dir
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, fsthick, mp_mesh_k
    USE pwcom,     ONLY : ef
    USE io_var,    ONLY : iufilikmap
    USE elph2,     ONLY : xkf, wkf, etf, nkf, ibndmin, ibndmax, nktotf, nbndfst
    USE eliashbergcom, ONLY : nkfs, ixkf, xkfs, wkfs, ekfs, nbndfs, memlt_pool
    USE superconductivity, ONLY : mem_size_eliashberg, mem_integer_size_eliashberg
    USE constants_epw, ONLY : zero
    USE mp_global, ONLY : inter_pool_comm, npool
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE division,  ONLY : fkbounds
    USE kfold,     ONLY : backtoBZ 
    !
    IMPLICIT NONE
    !
    INTEGER :: nk
    !! Counter on k points
    INTEGER :: nks
    !! Counter on k points within the Fermi shell
    INTEGER :: ikk
    !! k-point index
    INTEGER :: nkf_mesh
    !! Total number of k points
    !! These are irreducible k-points if mp_mesh_k = .TRUE.
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: imelt
    !! Memory allocated
    REAL(KIND = DP) :: xx, yy, zz
    !!
    REAL(KIND = DP), ALLOCATABLE :: xkf_(:, :)
    !! Temporary k point grid
    REAL(KIND = DP), ALLOCATABLE :: wkf_(:)
    !! Temporary weights on the k point grid
    REAL(KIND = DP), ALLOCATABLE :: ekf_(:, :)
    !! Temporary eigenvalues on the k point grid
    CHARACTER(LEN = 256) :: filikmap
    !! Name of the file
    !
    nkf_mesh = nktotf 
    nbndfs = nbndfst
    !
#if defined(__MPI)
    IF (.NOT. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(npool))
    memlt_pool(:) = zero
#else
    IF (.NOT. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(1))
    memlt_pool(1) = zero
#endif
    !
    ! get the size of required memory for ekf_, wkf_, xkf_
    imelt = ( nbndfs + 4 ) * nkf_mesh 
    CALL mem_size_eliashberg( imelt )
    !
    ! get the size of required memory for ixkf 
    imelt = nkf_mesh
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF (.NOT. ALLOCATED(ekf_) )   ALLOCATE(ekf_(nbndfs,nkf_mesh))
    IF (.NOT. ALLOCATED(wkf_) )   ALLOCATE(wkf_(nkf_mesh))
    IF (.NOT. ALLOCATED(xkf_) )   ALLOCATE(xkf_(3,nkf_mesh))
    IF (.NOT. ALLOCATED(ixkf) )   ALLOCATE(ixkf(nkf_mesh))
    xkf_(:, :) = zero
    ekf_(:, :) = zero
    wkf_(:) = zero
    ixkf(:) = 0
    !
    CALL fkbounds( nkf_mesh, lower_bnd, upper_bnd )
    !
    ! nkf - nr of k-blocks in the pool (fine grid)
    !
    DO nk = 1, nkf
       ikk = 2 * nk - 1
       xkf_(:,lower_bnd+nk-1) = xkf(:,ikk)
       wkf_(lower_bnd+nk-1)   = wkf(ikk)
       ekf_(:,lower_bnd+nk-1) = etf(ibndmin:ibndmax,ikk)
    ENDDO
       !
       ! collect contributions from all pools (sum over k-points)
       CALL mp_sum( ekf_, inter_pool_comm )
       CALL mp_sum( xkf_, inter_pool_comm )
       CALL mp_sum( wkf_, inter_pool_comm )
       CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime == ionode_id) THEN
      !
      IF (mp_mesh_k) THEN
         WRITE(stdout,'(/5x,a,i9/)') 'Nr. of irreducible k-points on the uniform grid: ', nkf_mesh
      ELSE
         WRITE(stdout,'(/5x,a,i9/)') 'Nr. of k-points on the uniform grid: ', nkf_mesh
      ENDIF
      !
      filikmap = TRIM(tmp_dir) // TRIM(prefix) // '.ikmap'
      !OPEN(iufilikmap, file = filikmap, form = 'formatted')
      !WRITE(iufilikmap,'(i9)') nkf_mesh
      OPEN(iufilikmap, file = filikmap, form = 'unformatted')
      WRITE(iufilikmap) nkf_mesh
      !
      ! nkfs - find nr of k-points within the Fermi shell (fine grid)
      ! only a fraction of nkf_mesh are contained in the Fermi shell
      !
      ! ixkf - find the index of k-point within the Fermi shell (fine grid)
      ! if the k-point lies outside the Fermi shell the index is 0
      !
      nkfs = 0  
      DO nk = 1, nkf_mesh
         IF (minval( ABS(ekf_(:,nk) - ef  ) ) < fsthick) THEN
            nkfs = nkfs + 1
            ixkf(nk) = nkfs
         ELSE
            ixkf(nk) = 0
         ENDIF
         !  bring back into to the first BZ
         xx = xkf_(1,nk) * nkf1
         yy = xkf_(2,nk) * nkf2
         zz = xkf_(3,nk) * nkf3
         CALL backtoBZ(xx, yy, zz, nkf1, nkf2, nkf3)
         xkf_(1,nk) = xx / DBLE(nkf1)
         xkf_(2,nk) = yy / DBLE(nkf2)
         xkf_(3,nk) = zz / DBLE(nkf3)
         !WRITE(iufilikmap,'(i9)') ixkf(nk)
         WRITE(iufilikmap) ixkf(nk)
      ENDDO
      CLOSE(iufilikmap)
      !
    ENDIF
    CALL mp_bcast( nkfs, ionode_id, inter_pool_comm )
    !
    ! get the size of required memory for ekfs, wkfs, xkfs 
    imelt = ( nbndfs + 4 ) * nkfs
    CALL mem_size_eliashberg( imelt )
    ! 
    IF (.NOT. ALLOCATED(ekfs) ) ALLOCATE(ekfs(nbndfs,nkfs))
    IF (.NOT. ALLOCATED(wkfs) ) ALLOCATE(wkfs(nkfs))
    IF (.NOT. ALLOCATED(xkfs) ) ALLOCATE(xkfs(3,nkfs))
    xkfs(:, :) = zero
    wkfs(:) = zero
    ekfs(:, :) = zero
    !
    IF (mpime == ionode_id) THEN
      nks = 0
      DO nk = 1, nkf_mesh
         IF (minval( ABS(ekf_(:,nk) - ef  ) ) < fsthick) THEN
            nks = nks + 1
            IF (nks > nkf_mesh ) CALL errore('kmesh_fine','too many k-points',1)
            wkfs(nks)   = wkf_(nk)
            xkfs(:,nks) = xkf_(:,nk)
            ekfs(:,nks) = ekf_(:,nk)
         ENDIF
      ENDDO
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( ixkf, ionode_id, inter_pool_comm )
    CALL mp_bcast( xkfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( wkfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( ekfs, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF (ALLOCATED(ekf_) ) DEALLOCATE(ekf_)
    IF (ALLOCATED(xkf_) ) DEALLOCATE(xkf_)
    IF (ALLOCATED(wkf_) ) DEALLOCATE(wkf_)
    !
    ! remove memory allocated for ekf_, xkf_, wkf_
    imelt = ( nbndfs + 4 ) * nkf_mesh
    CALL mem_size_eliashberg( -imelt )
    !
    WRITE(stdout,'(/5x,a/)') 'Finished writing .ikmap file '
    !
    RETURN
    !
    END SUBROUTINE kmesh_fine
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kqmap_fine
    !-----------------------------------------------------------------------
    !!
    !! this routine finds the index of k+sign*q on the fine k-mesh
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, mp_mesh_k
    USE elph2,     ONLY : nqtotf, xqf
    USE eliashbergcom, ONLY : ixkff, xkff, ixkf, xkfs, nkfs, ixkqf, ixqfs, nqfs
    USE superconductivity, ONLY : mem_size_eliashberg, mem_integer_size_eliashberg
    USE constants_epw, ONLY : eps5, zero
    USE symm_base, ONLY : nrot
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE division,  ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    INTEGER :: i, j, k, ik, nk, n
    !! Counter on k points
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: nkftot
    !! Total number of k points
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: nks
    !! Number of non-equivalent k points
    INTEGER :: ns
    !! Counter on rotation operations
    INTEGER :: imelt
    !! Memory allocated
    INTEGER, ALLOCATABLE :: equiv_(:)
    !! Index of equivalence of k points
    INTEGER, ALLOCATABLE :: index_(:, :)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    REAL(KIND = DP) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xq(3)
    !! coordinates of q points
    REAL(KIND = DP) :: xkr(3)
    !! coordinates of k points
    REAL(KIND = DP) :: xx, yy, zz
    !! Temporary variables
    LOGICAL :: in_the_list
    !! Check if k point is in the list
    !
    nkftot = nkf1 * nkf2 * nkf3
    !
    ! get the size of required memory for xkff
    imelt = 3 * nkftot 
    CALL mem_size_eliashberg( imelt )
    !
    ! get the size of required memory for ixkff and equiv_
    imelt = 2 * nkftot
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF (.NOT. ALLOCATED(xkff) )  ALLOCATE(xkff(3,nkftot))
    IF (.NOT. ALLOCATED(ixkff) ) ALLOCATE(ixkff(nkftot))
    xkff(:, :) = zero
    ixkff(:) = 0
    !
    ! to map k+q onto k we need to define the index of k on the full mesh (ixkff) 
    ! using index of the k-point within the Fermi shell (ixkf)
    !
    IF (mpime == ionode_id) THEN
      !
      DO i = 1, nkf1
         DO j = 1, nkf2
            DO k = 1, nkf3
               ik = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + k
               xkff(1,ik) = DBLE(i-1) / DBLE(nkf1)
               xkff(2,ik) = DBLE(j-1) / DBLE(nkf2)
               xkff(3,ik) = DBLE(k-1) / DBLE(nkf3)
            ENDDO
         ENDDO
      ENDDO
      !
      IF (.NOT. ALLOCATED(equiv_) ) ALLOCATE(equiv_(nkftot))
      !  equiv_(nk) =nk : k-point nk is not equivalent to any previous k-point
      !  equiv_(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
      !
      DO nk = 1, nkftot
         equiv_(nk) = nk
      ENDDO
      !
      IF (mp_mesh_k) THEN
        CALL set_sym_bl( )
        DO nk = 1, nkftot
          !  check if this k-point has already been found equivalent to another
          IF (equiv_(nk) == nk) THEN
            !  check if there are equivalent k-point to this in the list
            !  (excepted those previously found to be equivalent to another)
            !  check both k and -k
            DO ns = 1, nrot
              DO i = 1, 3
                xkr(i) = s(i,1,ns) * xkff(1,nk) &
                       + s(i,2,ns) * xkff(2,nk) &
                       + s(i,3,ns) * xkff(3,nk)
                xkr(i) = xkr(i) - NINT( xkr(i) )
              ENDDO
              IF (t_rev(ns) == 1 ) xkr = -xkr
              xx = xkr(1)*nkf1
              yy = xkr(2)*nkf2
              zz = xkr(3)*nkf3
              in_the_list = ABS(xx-NINT(xx) ) <= eps5 .AND. &
                            ABS(yy-NINT(yy) ) <= eps5 .AND. &
                            ABS(zz-NINT(zz) ) <= eps5
              IF (in_the_list) THEN
                i = MOD( NINT( xkr(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                j = MOD( NINT( xkr(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                k = MOD( NINT( xkr(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                n = (k-1) + (j-1)*nkf3 + (i-1)*nkf2*nkf3 + 1
                IF (n > nk .AND. equiv_(n) == n) THEN
                   equiv_(n) = nk
                ELSE
                   IF (equiv_(n) /= nk .OR. n < nk ) CALL errore('kmesh_fine', &
                      'something wrong in the checking algorithm',1)
                ENDIF
              ENDIF
              IF (time_reversal) THEN
                xx = -xkr(1)*nkf1
                yy = -xkr(2)*nkf2
                zz = -xkr(3)*nkf3
                in_the_list = ABS(xx-NINT(xx) ) <= eps5 .AND. &
                              ABS(yy-NINT(yy) ) <= eps5 .AND. &
                              ABS(zz-NINT(zz) ) <= eps5
                IF (in_the_list) THEN
                  i = MOD( NINT( -xkr(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                  j = MOD( NINT( -xkr(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                  k = MOD( NINT( -xkr(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                  n = (k-1) + (j-1)*nkf3 + (i-1)*nkf2*nkf3 + 1
                  IF (n > nk .AND. equiv_(n) == n) THEN
                    equiv_(n) = nk
                  ELSE
                    IF (equiv_(n) /= nk .OR. n < nk ) CALL errore('kmesh_fine', &
                       'something wrong in the checking algorithm',2)
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      !
      ! find index of k on the full mesh (ixkff) using index of k within the Fermi shell (ixkf)
      ! 
      nks = 0
      DO nk = 1, nkftot
        IF (equiv_(nk) == nk) THEN
          nks = nks + 1
          ixkff(nk) = ixkf(nks)
        ELSE
          ixkff(nk) = ixkff(equiv_(nk))
        ENDIF
      ENDDO
      !
      IF (ALLOCATED(equiv_) ) DEALLOCATE(equiv_)
      !
    ENDIF
    CALL mp_bcast( xkff, ionode_id, inter_pool_comm )
    CALL mp_bcast( ixkff, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF (ALLOCATED(xkff) ) DEALLOCATE(xkff)
    !
    ! remove memory allocated for xkff
    imelt = 3 * nkftot
    CALL mem_size_eliashberg( -imelt )
    !
    ! remove memory allocated for equiv_
    imelt = nkftot
    CALL mem_integer_size_eliashberg( -imelt )
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get the size of required memory for ixkqf, nqfs, index_
    imelt = ( nqtotf + 1 ) * nkfs + ( upper_bnd - lower_bnd + 1 ) * nqtotf
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF (.NOT. ALLOCATED(ixkqf) )  ALLOCATE(ixkqf(nkfs,nqtotf))
    IF (.NOT. ALLOCATED(nqfs) )   ALLOCATE(nqfs(nkfs))
    IF (.NOT. ALLOCATED(index_) ) ALLOCATE(index_(lower_bnd:upper_bnd,nqtotf))
    ixkqf(:, :) = 0
    nqfs(:) = 0
    index_(:, :) = 0
    !
    ! find the index of k+sign*q on the fine k-mesh
    ! nkfs - total nr. of k-points within the Fermi shell
    !      - these are irreducible k-points if mp_mesh_k = true
    ! nqtotf - total nr of q-points on the fine mesh
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqtotf
        xk(:) = xkfs(:,ik)
        xq(:) = xqf(:,iq)
        !
        ! find nkq - index of k+sign*q on the full fine k-mesh.
        !
        CALL kpmq_map( xk, xq, +1, nkq )
        !
        ! find ixkqf(ik,iq) - index of k+sign*q on the fine k-mesh within the Fermi shell
        !
        ixkqf(ik,iq) = ixkff(nkq) 
        !
        ! nqfs(ik) - nr of q-points at each k-point for which k+sign*q is within the Fermi shell 
        ! index_   - index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        IF (ixkqf(ik,iq) > 0) THEN
          nqfs(ik) = nqfs(ik) + 1
          index_(ik,nqfs(ik)) = iq
        ENDIF
      ENDDO ! loop over full set of q-points (fine mesh)
    ENDDO ! loop over k-points within the Fermi shell in each pool (fine mesh) 
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( ixkqf, inter_pool_comm )
    CALL mp_sum( nqfs,  inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    ! get the size of required memory for ixqfs
    imelt = nkfs * MAXVAL(nqfs(:))
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF (.NOT. ALLOCATED(ixqfs) ) ALLOCATE(ixqfs(nkfs,MAXVAL(nqfs(:))))
    ixqfs(:, :) = 0
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqfs(ik)
        !
        ! ixqfs - index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell 
        !
        ixqfs(ik,iq) = index_(ik,iq)   
      ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( ixqfs, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    ! remove memory allocated for ixkff, ixqfs, index_, nqfs
    imelt = nkftot + nkfs * MAXVAL(nqfs(:)) + nqtotf * ( upper_bnd - lower_bnd + 1 ) + nkfs
    CALL mem_integer_size_eliashberg( -imelt )
    !
    IF (ALLOCATED(ixkff) )  DEALLOCATE(ixkff)
    IF (ALLOCATED(ixqfs) )  DEALLOCATE(ixqfs)
    IF (ALLOCATED(index_) ) DEALLOCATE(index_)
    IF (ALLOCATED(nqfs) )   DEALLOCATE(nqfs)
    !
    IF (mp_mesh_k) THEN 
      WRITE(stdout,'(/5x,a/)') 'Finished mapping k+sign*q onto the fine irreducibe k-mesh'
    ELSE
      WRITE(stdout,'(/5x,a/)') 'Finished mapping k+sign*q onto the fine k-mesh'
    ENDIF
    ! 
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kqmap_fine
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kpmq_map(xk, xq, sign1, nkq)
    !-----------------------------------------------------------------------
    !!
    !! this routine finds the index of k+q or k-q point on the fine k-mesh
    !!
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : nkf1, nkf2, nkf3
    USE constants_epw, ONLY : eps5
    USE mp,        ONLY : mp_bcast, mp_barrier
    USE kfold,     ONLY : backtoBZ
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: sign1
    !! +1 for searching k+q, -1 for k-q
    INTEGER, INTENT(out) :: nkq
    !! the index of k+sign*q
    ! 
    REAL(KIND = DP), INTENT(in) :: xk(3)
    !! coordinates of k points
    REAL(KIND = DP), INTENT(in) :: xq(3)
    !! coordinates of q points
    ! 
    ! Local variables
    REAL(KIND = DP) :: xx, yy, zz, xxk(3)
    LOGICAL :: in_the_list
    !
    !
    xxk(:) = xk(:) + DBLE(sign1) * xq(:)
    xx = xxk(1) * nkf1
    yy = xxk(2) * nkf2
    zz = xxk(3) * nkf3
    in_the_list = ABS(xx-NINT(xx)) <= eps5 .AND. &
                  ABS(yy-NINT(yy)) <= eps5 .AND. &
                  ABS(zz-NINT(zz)) <= eps5
    IF (.NOT. in_the_list ) CALL errore('kpmq_map','k+q does not fall on k-grid',1)
    !
    !  find the index of this k+q or k-q in the k-grid
    !  make sure xx, yy, zz are in the 1st BZ
    !
    CALL backtoBZ(xx, yy, zz, nkf1, nkf2, nkf3)
    !
    ! since k- and q- meshes are commensurate, nkq can be easily found
    !
    nkq = NINT(xx) * nkf2 * nkf3 + NINT(yy) * nkf3 + NINT(zz) + 1
    !
    !  Now nkq represents the index of k+sign*q on the fine k-grid.
    !
    RETURN
    ! 
    END SUBROUTINE kpmq_map
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gap_distribution_FS ( itemp, cname )
    !-----------------------------------------------------------------------
    !
    ! This routine writes to files the distribution of the superconducting
    ! gap on the Fermi surface
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : fsthick
    USE eliashbergcom, ONLY : estemp, Agap, nkfs, nbndfs, ef0, ekfs, w0g
    USE constants_epw, ONLY : kelvin2eV, zero, eps5
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    CHARACTER(LEN = 256), INTENT(in) :: cname
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on band
    INTEGER :: ibin
    !! Counter on bins
    INTEGER :: nbin
    !! Number of bins
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: dbin
    !! Step size in nbin
    REAL(KIND = DP) :: delta_max
    !! Max value of superconducting gap
    REAL(KIND = DP) :: weight
    !! Variable for weight
    REAL(KIND = DP), ALLOCATABLE :: delta_k_bin(:)
    !! Histogram superconducting gap
    REAL(KIND = DP), EXTERNAL :: w0gauss
    CHARACTER(LEN = 256) :: name1
    !
    temp = estemp(itemp) / kelvin2eV
    !
    delta_max = 1.1d0 * MAXVAL(Agap(:,:,itemp))
    nbin = NINT(delta_max / eps5) + 1
    dbin = delta_max / DBLE(nbin)
    IF (.NOT. ALLOCATED(delta_k_bin) ) ALLOCATE(delta_k_bin(nbin) )
    delta_k_bin(:) = zero
    !
    DO ik = 1, nkfs
       DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd,ik) - ef0 ) < fsthick) THEN
            ibin = NINT( Agap(ibnd,ik,itemp) / dbin ) + 1
            weight = w0g(ibnd,ik)
            delta_k_bin(ibin) = delta_k_bin(ibin) + weight
          ENDIF
       ENDDO
    ENDDO
    !
    IF (temp < 10.d0) THEN
       WRITE(name1,'(a,a1,a4,a14,f4.2)') TRIM(prefix), '.', cname, '_aniso_gap0_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0 ) THEN
       WRITE(name1,'(a,a1,a4,a13,f5.2)') TRIM(prefix), '.', cname, '_aniso_gap0_0', temp
    ELSEIF (temp >= 100.d0) THEN
       WRITE(name1,'(a,a1,a4,a12,f6.2)') TRIM(prefix), '.', cname, '_aniso_gap0_', temp
    ENDIF
    !
    OPEN(iufilgap, FILE = name1, FORM = 'formatted')
    DO ibin = 1, nbin
       WRITE(iufilgap,'(2ES20.10)') temp + delta_k_bin(ibin)/MAXVAL(delta_k_bin(:)), dbin*DBLE(ibin)
    ENDDO
    CLOSE(iufilgap)
    !
    IF (ALLOCATED(delta_k_bin) ) DEALLOCATE(delta_k_bin)
    !
    RETURN
    !
    END SUBROUTINE gap_distribution_FS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gap_FS( itemp )
    !-----------------------------------------------------------------------
    !
    ! This routine writes to files the superconducting gap on the Fermi surface
    !
    USE kinds,         ONLY : DP
    USE io_var,        ONLY : iufilgapFS
    USE io_files,      ONLY : prefix
    USE cell_base,     ONLY : bg
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : fsthick, nkf1, nkf2, nkf3
    USE eliashbergcom, ONLY : estemp, Agap, nkfs, nbndfs, ef0, ekfs, ixkff
    USE constants_epw, ONLY : kelvin2eV, zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: i
    !! Counter on grid points nkf1
    INTEGER :: j
    !! Counter on grid points nkf2
    INTEGER :: k
    !! Counter on grid points nkf3
    !
    REAL(KIND = DP) :: temp
    !! Temperature in K
    REAL(KIND = DP) :: x1
    !! Cartesian coordinates of grid points nkf1
    REAL(KIND = DP) :: x2
    !! Cartesian coordinates of grid points nkf2
    REAL(KIND = DP) :: x3
    !! Cartesian coordinates of grid points nkf3
    REAL(KIND = DP), ALLOCATABLE :: Agap_tmp(:, :)
    !! Temporary array for superconducting gap at ik, ibnd
    CHARACTER(LEN = 256) :: name1, cname
    !
    temp = estemp(itemp) / kelvin2eV
    !
    cname = 'imag'
    !
    ! RM - If the k-point is outside the Fermi shell,
    ! ixkff(ik)=0 and Agap_tmp(:,0) = 0.0
    !
    IF (.NOT. ALLOCATED(Agap_tmp) ) ALLOCATE(Agap_tmp(nbndfs,0:nkfs))
    Agap_tmp(:,1:nkfs) = Agap(:,1:nkfs,itemp)
    Agap_tmp(:,0) = zero
    !
    ! SP & RM: .cube file for VESTA plotting (only if iverbosity = 2)
    !
    IF (iverbosity == 2) THEN
      !
      DO ibnd = 1, nbndfs
        !
        IF (ibnd < 10) THEN
          ! We make the assumption that there are no superconductor with Tc
          ! higher than 999 K.
          IF (temp < 10.d0) THEN
             WRITE(name1,'(a,a1,a4,a14,f4.2,a1,i1,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0) THEN
             WRITE(name1,'(a,a1,a4,a13,f5.2,a1,i1,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0) THEN
             WRITE(name1,'(a,a1,a4,a12,f6.2,a1,i1,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSEIF (ibnd < 100) THEN
          IF (temp < 10.d0) THEN
             WRITE(name1,'(a,a1,a4,a14,f4.2,a1,i2,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0 .AND. temp > 9.9999d0) THEN
             WRITE(name1,'(a,a1,a4,a13,f5.2,a1,i2,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0 .AND. temp > 99.9999d0) THEN
             WRITE(name1,'(a,a1,a4,a12,f6.2,a1,i2,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSEIF (ibnd < 1000) THEN
          IF (temp < 10.d0) THEN
             WRITE(name1,'(a,a1,a4,a14,f4.2,a1,i3,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 100.d0 .AND. temp > 9.9999d0 ) THEN
             WRITE(name1,'(a,a1,a4,a13,f5.2,a1,i3,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF (temp < 1000.d0 .AND. temp > 99.9999d0) THEN
             WRITE(name1,'(a,a1,a4,a12,f6.2,a1,i3,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSE
          CALL errore( 'eliashberg_write', ' Too many bands ',1)
        ENDIF
        !
        OPEN(iufilgapFS, FILE = name1, FORM = 'formatted')
        WRITE(iufilgapFS,*) 'Cubfile created from EPW calculation'
        WRITE(iufilgapFS,*) 'gap'
        WRITE(iufilgapFS,'(i5,3f12.6)') 1, 0.0d0, 0.0d0, 0.0d0
        WRITE(iufilgapFS,'(i5,3f12.6)') nkf1, (bg(i,1)/DBLE(nkf1),i = 1,3)
        WRITE(iufilgapFS,'(i5,3f12.6)') nkf2, (bg(i,2)/DBLE(nkf2),i = 1,3)
        WRITE(iufilgapFS,'(i5,3f12.6)') nkf3, (bg(i,3)/DBLE(nkf3),i = 1,3)
        WRITE(iufilgapFS,'(i5,4f12.6)') 1, 1.0d0, 0.0d0, 0.0d0, 0.0d0
        WRITE(iufilgapFS,'(6f12.6)') ( Agap_tmp(ibnd,ixkff(ik)),ik = 1,nkf1*nkf2*nkf3 )
        CLOSE(iufilgapFS)
      ENDDO
      !
    ENDIF
    !
    ! SP & RM : Write on file the superconducting gap close to the Fermi surface
    ! along with
    !     Cartesian coordinate, band index, energy distance from Fermi level and
    !     gap value.
    !
    IF (temp < 10.d0) THEN
       WRITE(name1,'(a,a1,a4,a16,f4.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_00', temp
    ELSEIF (temp >= 10.d0 .AND. temp < 100.d0) THEN
       WRITE(name1,'(a,a1,a4,a15,f5.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_0', temp
    ELSEIF (temp >= 100.d0) THEN
       WRITE(name1,'(a,a1,a4,a14,f6.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_', temp
    ENDIF
    OPEN(iufilgapFS, FILE = name1, FORM = 'formatted')
    WRITE(iufilgapFS,'(a78)') '#               k-point                  Band Enk-Ef [eV]        Delta(0) [eV]'
    DO i = 1, nkf1
      DO j = 1, nkf2
        DO k = 1, nkf3
          ik = k + (j-1)*nkf3 + (i-1)*nkf2*nkf3
          !IF (ixkff(ik) > 0) THEN
            DO ibnd = 1, nbndfs
              ! RM: Everything is in eV here.
              ! SP: Here take a 0.2 eV interval around the FS.
              IF (ABS(ekfs(ibnd,ixkff(ik)) - ef0 ) < fsthick) THEN
              !IF (ABS(ekfs(ibnd,ixkff(ik)) - ef0 ) < 0.2) THEN
                 x1 = bg(1,1)*(i-1)/nkf1+bg(1,2)*(j-1)/nkf2+bg(1,3)*(k-1)/nkf3
                 x2 = bg(2,1)*(i-1)/nkf1+bg(2,2)*(j-1)/nkf2+bg(2,3)*(k-1)/nkf3
                 x3 = bg(3,1)*(i-1)/nkf1+bg(3,2)*(j-1)/nkf2+bg(3,3)*(k-1)/nkf3
                 WRITE(iufilgapFS,'(3f12.6,i8,f12.6,f24.15)') x1, x2, x3, ibnd, &
                       ekfs(ibnd,ixkff(ik))-ef0, Agap_tmp(ibnd,ixkff(ik))
              ENDIF
            ENDDO ! ibnd
          !ENDIF
        ENDDO  ! k
      ENDDO ! j
    ENDDO ! i
    CLOSE(iufilgapFS)
    !
    IF (ALLOCATED(Agap_tmp) ) DEALLOCATE(Agap_tmp)
    !
    RETURN
    !
    END SUBROUTINE gap_FS
    !                                                                            
    !----------------------------------------------------------------------
    ! 
  END MODULE io_eliashberg


