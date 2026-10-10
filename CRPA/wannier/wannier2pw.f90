!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM wannier2pw
  !-----------------------------------------------------------------------
  !
  ! Sample code, showing how to read QE data and re-use QE variables:
  ! 1. reads the data directory of QE, then
  ! 2. fills the hamiltonian matrix and diagonalizes it for each k-point
  !    (conventional, not iterative diagonalization)
  !
  ! Input: namelist &inputpp [outdir=...] [prefix=...] / as in QE input
  ! (default values as in QE).
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : ionode
  USE mp_global,        ONLY : mp_startup, mp_global_end
  USE environment,      ONLY : environment_start, environment_end
  USE w90_interface,    ONLY : seedname, w90_read
  USE w90_interpolate,  ONLY : w90_interpolate_wfc
  USE w90_real_space,   ONLY : w90_write_real_space
  !
  IMPLICIT NONE
  !
  LOGICAL :: needwf
  CHARACTER(LEN=256) :: outdir
  CHARACTER(LEN=256) :: filkf
  !! File containing k points
  CHARACTER(LEN=256) :: folder_wan_Rr
  !! Directory for the real-space collected Wannier functions
  LOGICAL :: write_wan_Rr
  !! If .true. write the collected real-space Wannier functions to file.
  !! If .false. use the existing files in folder_wan_Rr.
  !
  ! initialise environment
  !
  CALL mp_startup ( )
  CALL environment_start ( 'Wannier2PW' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  CALL w90_readin(outdir)
  !
  !   Read xml file, allocate and initialize general variables
  !
  needwf = .true.
  CALL read_file_new ( needwf )
  !
  ! Read Wannier90 information from checkpoint file
  !
  CALL w90_read(seedname)
  !
  CALL w90_run_interpolate()
  !
  CALL print_clock_w90()
  !
  CALL environment_end ( 'Wannier2PW' )
  !
  CALL mp_global_end()
  !
  CONTAINS
  !
  !-------------------------------------------------------------------------
  SUBROUTINE print_clock_w90()
  !-------------------------------------------------------------------------
    USE io_global, ONLY : stdout
    !
    IMPLICIT NONE
    !
    WRITE(stdout, '()')
    CALL print_clock("w90_write_rs")
    CALL print_clock("w90_buf_save")
    CALL print_clock("w90_interpolate")
    !
    WRITE(stdout, '()')
    ! CALL print_clock("w90_wigner")
    CALL print_clock("w90_readwrite")
    CALL print_clock("w90_buf_get")
    CALL print_clock("w90_mult_iqr")
    CALL print_clock("w90_shift_wfc")
    !
    WRITE(stdout, '()')
    CALL print_clock("davcio")
    CALL print_clock("fft")
    CALL print_clock("ffts")
    CALL print_clock("fftw")
    !
  END SUBROUTINE print_clock_w90
  !-------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------
  SUBROUTINE w90_run_interpolate()
  !---------------------------------------------------------------------------
  !! Read input namelist
  !---------------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout
    USE io_files,         ONLY : tmp_dir
    USE constants,        ONLY : rytoev
    USE mp_bands,         ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE fft_base,         ONLY : dfftp
    USE wvfct,            ONLY : npwx, current_k
    USE noncollin_module, ONLY : npol
    USE becmod,           ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type
    USE klist,            ONLY : xk, nks, igk_k, ngk
    USE lsda_mod,         ONLY : nspin, current_spin
    USE scf,              ONLY : vrs, vltot, v, kedtau
    USE uspp,             ONLY : nkb, vkb
    USE uspp_init,        ONLY : init_us_2
    USE gvecs,            ONLY : doublegrid
    USE w90_interface,    ONLY : num_wann
    USE w90_wan_Rr_buffer,ONLY : w90_wan_Rr_save_buffer, w90_wan_Rr_close_buffer
    !
    IMPLICIT NONE
    !
    INTEGER :: ik
    !! k point index
    INTEGER :: npw
    !! Number of plane waves at xkp
    REAL(DP) :: xkp(3)
    !! k point to interpolate the wavefunctions to (in Cartesian coordinates).
    COMPLEX(DP), ALLOCATABLE :: evc_kp(:, :)
    !! Interpolated wavefunction at kp in reciprocal space. (npwx*npol, num_wann)
    COMPLEX(DP), ALLOCATABLE :: aux(:,:)
    COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    REAL(DP),    ALLOCATABLE :: en(:)
    !
    INCLUDE 'laxlib.h'
    !
    ! 1. Read u_nk(G) wavefunctions from file
    ! 2. Multiply gauge matrix to get u_ik(G)
    ! 3. Fourier transform to u_ik(r)
    ! 4. Compute factors for k to R to k' transformation
    ! 5. Compute the wavefunction at k'
    !
    ! Step 1-3
    !
    ! Compute and write real-space Wannier functions to file in collected form
    !
    IF (write_wan_Rr) CALL w90_write_real_space(tmp_dir, folder_wan_Rr)
    !
    ! Read the collected real-space Wannier functions from file to buffer
    !
    CALL w90_wan_Rr_save_buffer(folder_wan_Rr, tmp_dir)
    !
    IF (TRIM(filkf) /= '') THEN
      CALL w90_read_kpoints(filkf)
    ENDIF
    !
    ALLOCATE(evc_kp(npwx*npol, num_wann))
    !
    ALLOCATE( aux(npwx, num_wann ) )
    ALLOCATE( hc( num_wann, num_wann) )
    ALLOCATE( sc( num_wann, num_wann) )
    ALLOCATE( vc( num_wann, num_wann) )
    ALLOCATE( en( num_wann ) )
    CALL allocate_bec_type (nkb, num_wann, becp)
    CALL set_vrs(vrs,vltot,v%of_r,kedtau,v%kin_r,dfftp%nnr,nspin,doublegrid)
    !
    ! Step 4-5
    !
    DO ik = 1, nks
      !
      current_k = ik
      current_spin = 1
      xkp = xk(:, ik)
      npw = ngk(ik)
      !
      CALL w90_interpolate_wfc(ik, evc_kp)
      !
      ! Compute Hamiltonian matrix in the evc_kp basis
      !
      CALL init_us_2 (npw, igk_k(:, ik), xkp, vkb)
      CALL calbec ( npw, vkb, evc_kp, becp)
      CALL g2_kin (ik)
      !
      CALL h_psi( npwx, npw, num_wann, evc_kp, aux )
      CALL calbec ( npw, evc_kp, aux, hc )
      CALL s_psi( npwx, npw, num_wann, evc_kp, aux )
      CALL calbec ( npw, evc_kp, aux, sc )
      !
      CALL diaghg( num_wann, num_wann, hc, sc, num_wann, en, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
      !
      WRITE(stdout, '(/,12x,"k =",3f7.4," (",i6," PWs)  bands (eV):",/)') xkp, npw
      WRITE(stdout, '(8f9.4)') en(:) * rytoev
      !
    ENDDO
    !
    CALL deallocate_bec_type (becp)
    CALL w90_wan_Rr_close_buffer()
    !
  !---------------------------------------------------------------------------
  END SUBROUTINE w90_run_interpolate
  !---------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------
  SUBROUTINE w90_read_kpoints(filename)
  !---------------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE io_global,        ONLY : stdout
    USE constants,        ONLY : rytoev
    USE cell_base,        ONLY : bg
    USE wvfct,            ONLY : g2kin
    USE becmod,           ONLY : calbec, allocate_bec_type, deallocate_bec_type
    USE klist,            ONLY : xk, nks, igk_k, ngk, init_igk
    USE uspp,             ONLY : vkb
    USE uspp_init,        ONLY : init_us_2
    USE wavefunctions,    ONLY : evc
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(in) :: filename
    !
    CHARACTER(LEN = 10) :: coordinate_type
    !! Coordinate type of k points. Cartesian or crystal
    INTEGER :: ik
    !! k point index
    INTEGER :: iun
    !! File unit
    INTEGER :: ios
    !! I/O status
    REAL(DP) :: dummy
    !! Dummy variable for reading k points
    LOGICAL, EXTERNAL :: imatches
    !
    ! Reset k point and G vectors
    !
    DEALLOCATE(evc)
    DEALLOCATE(vkb)
    DEALLOCATE(g2kin)
    DEALLOCATE(igk_k)
    DEALLOCATE(ngk)
    !
    ! Read k points from file
    !
    WRITE(stdout, '(5x,a,a)') 'Reading k-mesh file: ', TRIM(filename)
    !
    OPEN(NEWUNIT = iun, FILE = filename, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
    IF (ios /= 0) CALL errore('loadkmesh_para', 'opening file ' // filename, ABS(ios))
    READ(iun, *) nks, coordinate_type
    !
    ! coordinate_type: crystal (default) or cartesian
    IF (TRIM(coordinate_type) .EQ. '') coordinate_type = 'crystal'
    IF (.NOT. imatches("crystal", coordinate_type) .AND. .NOT. imatches("cartesian", coordinate_type)) THEN
      CALL errore('loadkmesh_para', 'ERROR: Specify either crystal or cartesian coordinates in the file', 1)
    ENDIF
    !
    DO ik = 1, nks
      READ(iun, *) xk(:, ik), dummy
    ENDDO
    !
    CLOSE(iun)
    !
    IF (imatches("crystal", coordinate_type)) THEN
      CALL cryst_to_cart(nks, xk, bg, +1)
    ENDIF
    !
    CALL allocate_wfc_k()
    !
  !---------------------------------------------------------------------------
  END SUBROUTINE w90_read_kpoints
  !---------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------
  SUBROUTINE w90_readin(outdir)
  !---------------------------------------------------------------------------
  !! Read input namelist
  !---------------------------------------------------------------------------
    !
    USE io_global,     ONLY : meta_ionode, meta_ionode_id
    USE io_files,      ONLY : tmp_dir, prefix
    USE mp_world,      ONLY : world_comm
    USE mp,            ONLY : mp_bcast
    USE control_flags, ONLY : io_level
    USE w90_interface, ONLY : seedname, w90_read
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(out) :: outdir
    !
    INTEGER :: ios
    LOGICAL :: reduce_io
    CHARACTER(LEN=256), EXTERNAL :: trimcheck
    !
    NAMELIST / wannierpw / outdir, prefix, seedname, filkf, folder_wan_Rr, reduce_io, write_wan_Rr
    !
    ! filkf : File containing k points.
    !
    ! folder_wan_Rr : Directory for the real-space collected Wannier functions.
    !                 Default: outdir
    !
    ! reduce_io : If .true. then the code will store the distributed real-space Wannier
    !             functions in memory. This reduces I/O but increases memory usage.
    !             If .false. then the code will write the Wannier functions buffer to file.
    !             Note in either case the collected Wannier functions are written to files
    !             in folder_wan_Rr.
    !             Default: .true.
    !
    ! write_wan_Rr : If .true. write the collected real-space Wannier functions to file.
    !                If .false. use the existing files in folder_wan_Rr. In this case,
    !                prefix.wan_Rr.xxx and prefix.wan_R.dat files must exist inside folder_wan_Rr.
    !
    !   set default values for variables in namelist
    !
    prefix = 'pwscf'
    CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
    IF ( trim( outdir ) == ' ' ) outdir = './'
    !
    seedname = 'wannier'
    filkf = ''
    write_wan_Rr = .true.
    !
    IF ( meta_ionode )  THEN
      !
      !     reading the namelist wannier
      !
      READ (5, wannierpw, iostat = ios)
      !
      tmp_dir = trimcheck ( outdir )
      !
    ENDIF
    !
    CALL mp_bcast (ios, meta_ionode_id, world_comm)
    IF ( ios /= 0) CALL errore ('wannier2pw', 'reading wannierpw namelist', abs(ios))
    !
    IF (TRIM(folder_wan_Rr) == '') folder_wan_Rr = TRIM(tmp_dir)
    !
    ! ... Broadcast variables
    !
    CALL mp_bcast(tmp_dir, meta_ionode_id, world_comm)
    CALL mp_bcast(prefix, meta_ionode_id, world_comm)
    CALL mp_bcast(seedname, meta_ionode_id, world_comm)
    CALL mp_bcast(filkf, meta_ionode_id, world_comm)
    CALL mp_bcast(folder_wan_Rr, meta_ionode_id, world_comm)
    CALL mp_bcast(reduce_io, meta_ionode_id, world_comm)
    CALL mp_bcast(write_wan_Rr, meta_ionode_id, world_comm)
    !
    IF (reduce_io) THEN
      io_level = 0  ! Store buffer in memory, no I/O
    ELSE
      io_level = 1  ! Write buffer to file
    ENDIF
    !
  !---------------------------------------------------------------------------
  END SUBROUTINE w90_readin
  !---------------------------------------------------------------------------
  !
!---------------------------------------------------------------------------
END PROGRAM wannier2pw
!---------------------------------------------------------------------------
