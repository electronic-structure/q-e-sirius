!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM elph_main
  !-----------------------------------------------------------------------
  !
  ! This is the main driver of the elph code.
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout
  USE check_stop,        ONLY : check_stop_init
  USE mp_global,         ONLY : mp_startup, mp_global_end
  USE environment,       ONLY : environment_start, environment_end
  USE ions_base,         ONLY : nat
  USE control_flags,     ONLY : use_para_diag, use_gpu
  USE qpoint,            ONLY : nqs
  USE ph_restart,        ONLY : allocate_grid_variables
  USE control_ph,        ONLY : tmp_dir_ph, trans
  USE modes,             ONLY : nirr
  USE save_ph,           ONLY : save_ph_input_variables
  USE el_phon,           ONLY : elph_mat
  USE elphcom,           ONLY : wannier_seedname, folder_wan_Rr, write_wan_Rr, tmp_dir_save
  USE w90_interface,     ONLY : w90_read
  USE w90_real_space,    ONLY : w90_write_real_space
  USE w90_wan_Rr_buffer, ONLY : w90_wan_Rr_save_buffer
  USE mp_pools, ONLY : npool
  USE control_flags, ONLY : iverbosity
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=4)   :: code = 'ELPH'
  LOGICAL :: do_band
  !! If .true., perform NSCF band structure calculation
  INTEGER :: iq
  LOGICAL :: do_iq, setup_pw
  LOGICAL,EXTERNAL :: check_gpu_support
  !
  use_gpu = check_gpu_support()
  !
  ! Initialize MPI, clocks, print initial messages
  !
  CALL mp_startup(start_images = .TRUE.)
  !
  CALL environment_start(code)
  iverbosity = 1
  !
  ! Print the preamble
  !
  CALL print_preamble()
  !
  ! Read the input parameters and the data produced by PWscf
  !
  CALL elph_readin()
  !
  IF (npool > 1) CALL errore('elph_main', 'npool not yet implemented', 1)
  !
  ! Initialization
  !
  CALL check_stop_init()
  !
  CALL elph_q_points()
  !
  ! Initialize PHonon variables
  !
  trans = .false.
  elph_mat = .false.
  nirr = 3 * nat
  call allocate_grid_variables()
  CALL init_representations()
  CALL allocate_part ( nat )
  CALL save_ph_input_variables()
  !
  WRITE (stdout, '(5x, a)') 'Reading Wannier90 checkpoint file'
  CALL w90_read(wannier_seedname)
  !
  IF (write_wan_Rr) THEN
    WRITE (stdout, '(5x, a)') 'Writing real-space Wannier functions in collected form'
    CALL w90_write_real_space(tmp_dir_save, folder_wan_Rr)
  ELSE
    WRITE (stdout, '(5x, a)') 'Using collected real-space Wannier functions in ' // TRIM(folder_wan_Rr)
  ENDIF
  !
  WRITE (stdout, '(5x, a)') 'Loading real-space Wannier functions to buffer in distributed form'
  CALL w90_wan_Rr_save_buffer(folder_wan_Rr, tmp_dir_ph)
  !
  DO iq = 1, nqs
    !
    CALL elph_prepare_q(iq, do_iq, setup_pw)
    !
    do_band = .FALSE.
    ! do_band = .TRUE.
    !
    CALL run_nscf(do_band, iq)
    !
    !  Initialize the quantities which do not depend on
    !  the linear response of the system
    !
    CALL elph_load_q()
    !
    CALL elph_compute_mel()
    !
    CALL clean_pw_ph(iq)
    !
  ENDDO
! !   !
! !   ! Deallocate some arrays
! !   !
! !   CALL hp_dealloc_2()
!   !
!   ! Print clocks
!   !
!   WRITE( stdout, * )
!   WRITE( stdout, * )  '    PRINTING TIMING FROM PWSCF ROUTINES: '
!   CALL print_clock_pw()
!   CALL crpa_print_clock()
  !
  CALL environment_end(code)
  !
  IF ( use_para_diag ) CALL laxlib_end()
  CALL mp_global_end()
  !
3336 FORMAT('     ',69('='))
  !
  STOP
  !
CONTAINS
  !
SUBROUTINE print_preamble()
  !
  IMPLICIT NONE
  !
  WRITE( stdout, '(/,5X,"=---------------------------------------------------------------------------=")')
  WRITE( stdout, '(/,5X,"   Calculation of the electron-phonon coupling using the ELPH code           ")')
  WRITE( stdout, '(/,5X,"      Please cite the following papers when using this program:              ")')
  WRITE( stdout, '(/,5X,"         - ELPH code : To be added                                           ")')
  WRITE( stdout, '(/,5X,"=---------------------------------------------------------------------------=")')
  !
  RETURN
  !
END SUBROUTINE print_preamble
  !
END PROGRAM elph_main
