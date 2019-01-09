! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qexsd_module
  !----------------------------------------------------------------------------
  !
  ! This module contains some common subroutines used to read and write
  ! in XML format the data produced by Quantum ESPRESSO package.
  !
  ! Written by Giovanni Borghi, A. Ferretti, ... (2015).
  !
  ! Based on the qexml.f90 routine:
  ! Written by Andrea Ferretti (2006).
  ! Initial work by Carlo Sbraccia (xml_io_base.f90)
  ! Modified by Simone Ziraldo (2013).
  !
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : input_xml_schema_file, title
  USE mp_world,         ONLY : nproc
  USE mp_images,         ONLY : nimage,nproc_image
  USE mp_pools,         ONLY : npool
  USE mp_bands,         ONLY : ntask_groups, nproc_bgrp, nbgrp
  USE global_version,   ONLY:  version_number, svn_revision
  !
  USE constants,        ONLY : e2
  USE qes_types_module
  USE qes_libs_module
  !
  USE FoX_wxml
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  ! definitions for the fmt
  !
  CHARACTER(5), PARAMETER :: fmt_name = "QEXSD"
  CHARACTER(5), PARAMETER :: fmt_version = "0.1.0"
  !
  ! some default for kinds
  !
  !INTEGER,   PARAMETER :: DP = selected_real_kind( 14, 200 )
  !
  ! internal data to be set
  !
  CHARACTER(256)   :: datadir_in, datadir_out
  INTEGER          :: iunit, ounit
  TYPE(xmlf_t)     :: qexsd_xf
  !
  ! vars to manage back compatibility
  !
  CHARACTER(10)    :: qexsd_current_version = " "
  CHARACTER(10)    :: qexsd_default_version = trim( fmt_version  )
  LOGICAL          :: qexsd_current_version_init = .FALSE.
  !
  LOGICAL          :: qexsd_use_large_indent = .FALSE.
  !
  ! 
  TYPE (input_type)                            :: qexsd_input_obj
  TYPE (general_info_type)                     :: general_info
  TYPE (parallel_info_type)                    :: parallel_info
  TYPE (berryPhaseOutput_type),TARGET          :: qexsd_bp_obj
  TYPE (k_points_IBZ_type)                     :: qexsd_start_k_obj 
  TYPE (occupations_type)                      :: qexsd_occ_obj
  TYPE (smearing_type)                         :: qexsd_smear_obj
  TYPE ( step_type),ALLOCATABLE                :: steps(:)
  INTEGER                                      :: exit_status
  TYPE ( closed_type )                         :: qexsd_closed_element
  INTEGER                                      :: step_counter
  !
  ! end of declarations
  !
  PUBLIC :: qexsd_current_version, qexsd_default_version
  PUBLIC :: qexsd_current_version_init
  PUBLIC :: qexsd_xf  
  !
  PUBLIC :: qexsd_input_obj, qexsd_start_k_obj, qexsd_occ_obj, qexsd_smear_obj
  ! 
  PUBLIC :: qexsd_init_schema,  qexsd_openschema, qexsd_closeschema
  !
  PUBLIC :: qexsd_init_convergence_info, qexsd_init_algorithmic_info, &
            qexsd_init_atomic_species, qexsd_init_atomic_structure, &
            qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
            qexsd_init_magnetization, qexsd_init_band_structure, & 
            qexsd_init_total_energy, qexsd_init_forces, qexsd_init_stress, &
            qexsd_init_dipole_info, qexsd_init_outputElectricField,   &
            qexsd_init_outputPBC, qexsd_init_gate_info  
  !
  PUBLIC :: qexsd_step_addstep, qexsd_set_status, qexsd_reset_steps    
  ! 
  PUBLIC :: qexsd_init_berryPhaseOutput, qexsd_bp_obj
CONTAINS
!
!-------------------------------------------
! ... basic (public) subroutines
!-------------------------------------------
!
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_schema( unit_in, unit_out )
      !------------------------------------------------------------------------
      !
      ! just init module data
      !
      IMPLICIT NONE
      INTEGER,                INTENT(in) :: unit_in
      INTEGER,      OPTIONAL, INTENT(in) :: unit_out
      !
      iunit       = unit_in
      ounit       = unit_in
      IF ( present( unit_out ) ) ounit  = unit_out
      !
      !
    END SUBROUTINE qexsd_init_schema
    !
    !
    !------------------------------------------------------------------------
    FUNCTION check_file_exst( filename )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL          :: check_file_exst
      CHARACTER(len=*) :: filename
      !
      LOGICAL :: lexists
      !
      INQUIRE( FILE = trim( filename ), EXIST = lexists )
      !
      check_file_exst = lexists
      RETURN
      !
    END FUNCTION check_file_exst
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_openschema( filename , prog)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      CHARACTER(len=*), INTENT(IN) :: filename, prog
      CHARACTER(len=16) :: subname = 'qexsd_openschema'
      INTEGER :: ierr, len_steps, i_step
      !
      ! we need a qes-version number here
      CALL xml_OpenFile(FILENAME = TRIM(filename), XF = qexsd_xf, UNIT = ounit, PRETTY_PRINT = .TRUE., &
                        REPLACE  = .TRUE., NAMESPACE = .TRUE., IOSTAT = ierr ) 
      !
      CALL xml_DeclareNamespace (XF=qexsd_xf, PREFIX = "xsi", nsURI ="http://www.w3.org/2001/XMLSchema-instance")
      CALL xml_DeclareNamespace (XF=qexsd_xf, PREFIX = "qes", nsURI ="http://www.quantum-espresso.org/ns/qes/qes-1.0")
      CALL xml_NewElement (XF=qexsd_xf, NAME = "qes:espresso")
      CALL xml_addAttribute(XF=qexsd_xf, NAME = "xsi:schemaLocation", &
                            VALUE = "http://www.quantum-espresso.org/ns/qes/qes-1.0 "//&
                                    "http://www.quantum-espresso.org/ns/qes/qes-1.0.xsd" )
      CALL xml_addAttribute(XF=qexsd_xf, NAME="Units", VALUE="Hartree atomic units")
      CALL xml_addComment(XF = qexsd_xf, &
              COMMENT = "If not explicitely indicated, all quantities are expressed in Hartree atomic units" ) 
      !
      IF (ierr /= 0) call errore(subname, 'opening xml output file', ierr)
      ! the input file is mandatory to have a validating schema 
      ! here an error should be issued, instead
      !
      CALL qexsd_init_general_info(general_info, prog(1:2) )
      CALL qes_write_general_info(qexsd_xf,general_info)
      CALL qes_reset_general_info(general_info)
      !
      CALL qexsd_init_parallel_info(parallel_info)
      CALL qes_write_parallel_info(qexsd_xf,parallel_info)
      CALL qes_reset_parallel_info(parallel_info) 
      IF ( check_file_exst(input_xml_schema_file) )  THEN
         CALL xml_addComment( XF = qexsd_xf, &
                              COMMENT= "")
         CALL qexsd_cp_line_by_line(ounit ,input_xml_schema_file, spec_tag="input")
      ELSE IF ( TRIM(qexsd_input_obj%tagname) == "input") THEN 
         CALL qes_write_input(qexsd_xf, qexsd_input_obj)
      END IF
      ! 
      IF (ALLOCATED(steps) ) THEN 
         len_steps= step_counter 
         IF (TRIM (steps(1)%tagname ) .EQ. 'step') THEN
            DO i_step = 1, len_steps
               CALL qes_write_step(qexsd_xf, steps(i_step) )
            END DO 
         END IF
      END IF
      ! 
    END SUBROUTINE qexsd_openschema
    !
    !
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_general_info(obj, prog )
    !---------------------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE( general_info_type )         ::  obj
      CHARACTER(LEN=*),INTENT(IN)       ::  prog
      CHARACTER(LEN=*),PARAMETER        ::  TAGNAME="general_info"
      TYPE( creator_type )              ::  creator_obj
      TYPE( created_type )              ::  created_obj
      TYPE( xml_format_type)            ::  xml_fmt_obj
      CHARACTER(LEN=256)                ::  version
      CHARACTER(9)                      ::  cdate, ctime
      CHARACTER(60)                     ::  timestamp
      !
      version=TRIM(version_number)
      IF (svn_revision .NE. "unknown") version=TRIM(version) // &
                                       " (svn rev. " // TRIM (svn_revision) // ")"
      SELECT CASE( prog(1:2))
         CASE ('pw','PW') 
            CALL qes_init_creator(creator_obj, "creator", "PWSCF", version, "XML file generated by PWSCF")
         CASE ('cp', 'CP') 
            CALL qes_init_creator(creator_obj, "creator", "CP", version, "XML file generated by CP") 
      END SELECT
      !
      CALL date_and_tim(cdate, ctime) 
      timestamp = 'This run was terminated on:  ' // ctime // ' ' // cdate(1:2) // & 
                  ' '//cdate(3:5) // ' '// cdate (6:9)
     
      CALL qes_init_created (created_obj, "created", cdate, ctime, timestamp ) 
      !
      CALL qes_init_xml_format(xml_fmt_obj, "xml_format", fmt_name, fmt_version, fmt_name//"_"//fmt_version)
      !
      CALL qes_init_general_info( obj, TAGNAME, xml_fmt_obj, creator = creator_obj, created = created_obj,&
                                  job=title)
      !
      CALL qes_reset_creator(creator_obj)
      CALL qes_reset_created(created_obj)
      CALL qes_reset_xml_format(xml_fmt_obj) 
    END SUBROUTINE qexsd_init_general_info    
    !
    !---------------------------------------------------------------------------------------------
    SUBROUTINE   qexsd_init_parallel_info(obj)
    !---------------------------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE ( parallel_info_type )           :: obj
      !
      INTEGER                               :: nthreads=1
#if defined(__OMP) 
      INTEGER,EXTERNAL                      :: omp_get_max
      !     
      nthreads = omp_get_max()
#endif      
      CALL qes_init_parallel_info(obj, "parallel_info", nproc, nthreads, ntask_groups, &
                                  nbgrp, npool, nproc_bgrp)
    END SUBROUTINE qexsd_init_parallel_info
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_closeschema()
      !------------------------------------------------------------------------
      USE mytime,    ONLY: nclock, clock_label
      IMPLICIT NONE
      REAL(DP),EXTERNAL    :: get_clock
      !
      CHARACTER(len=17) :: subname = 'qexsd_closeschema'
      INTEGER :: ierr
      !
      IF (exit_status .ge. 0 ) THEN 
         CALL xml_NewElement(qexsd_xf, "status")
         CALL xml_AddCharacters(qexsd_xf, exit_status)
         CALL xml_EndElement(qexsd_xf, "status")          
         CALL qexsd_set_closed()
         CALL xml_NewElement (qexsd_xf, "cputime")
         CALL xml_addCharacters(qexsd_xf, MAX(nint(get_clock('PWSCF')),nint(get_clock('CP'))) )
         CALL xml_EndElement ( qexsd_xf, "cputime")
         CALL qes_write_closed(qexsd_xf, qexsd_closed_element)
      END IF
         CALL xml_Close(qexsd_xf) 
      !
    END SUBROUTINE qexsd_closeschema
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_cp_line_by_line(iun_out,filename,spec_tag)
      !------------------------------------------------------------------------
      implicit none
      !
      integer,      intent(in) :: iun_out
      character(*), intent(in) :: filename
      character(*), optional, intent(in) :: spec_tag
      !
      integer :: iun, ierr
      character(256) :: str
      logical :: icopy, exists
      integer, external  :: find_free_unit

      iun =  find_free_unit()
      !
      INQUIRE(FILE=trim(filename), EXIST=exists)
      !
      IF(.not.exists) THEN
         CALL errore('qexsd_cp_line_by_line', 'input xml file "' // & 
        &             TRIM(filename) // '" not found', 1)
      ENDIF
      !
      open(iun,FILE=trim(filename),status="old", IOSTAT=ierr)
      !
      icopy=.false.
      copy_loop: do
         !
         read(iun,"(a256)",iostat=ierr) str
         if (ierr<0) exit copy_loop
         if (present(spec_tag)) then
            !
            if (index(str,"<"//trim(adjustl(spec_tag))//">")/=0) then
               !
               icopy=.true.
               !
            endif
            !
         else
            !
            icopy=.true.
            !
         endif
         ! 
         ! filtering
         ! 
         if ( index(str,"<Root>")/=0 .or. index(str,"<Root>")/=0 .or. &
              index(str,"<?")/=0     .or. .not.icopy) cycle copy_loop
         !
         write(iun_out,"(a)") trim(str)
         !
         if (present(spec_tag)) then
            if (index(str,"</input>")/=0) icopy=.false.
         endif
         ! 
      enddo copy_loop
      !
      close(iun)
      ! 
    END SUBROUTINE qexsd_cp_line_by_line
    !
!
!-------------------------------------------
! ... write subroutines
!-------------------------------------------
!
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_convergence_info(obj, n_scf_steps, scf_has_converged, scf_error, &
                                           optimization_has_converged, n_opt_steps, grad_norm )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(convergence_info_type)   :: obj
      INTEGER,           INTENT(IN) :: n_scf_steps
      LOGICAL,           INTENT(IN) :: scf_has_converged   
      REAL(DP),          INTENT(IN) :: scf_error
      LOGICAL, OPTIONAL, INTENT(IN) :: optimization_has_converged
      INTEGER, OPTIONAL, INTENT(in) :: n_opt_steps
      REAL(DP),OPTIONAL, INTENT(IN) :: grad_norm
      !
      CHARACTER(27)       :: subname="qexsd_init_convergence_info"
      TYPE(scf_conv_type) :: scf_conv
      TYPE(opt_conv_type),POINTER :: opt_conv => NULL() 
      !
      call qes_init_scf_conv(scf_conv, "scf_conv", scf_has_converged, n_scf_steps, scf_error)
      !
      IF ( PRESENT(optimization_has_converged ))  THEN
          !
          IF ( .NOT. PRESENT(n_opt_steps) ) CALL errore(subname,"n_opt_steps not present",10)
          IF ( .NOT. PRESENT(grad_norm) )   CALL errore(subname,"grad_norm not present",10)
          ALLOCATE ( opt_conv) 
          !
          call qes_init_opt_conv(opt_conv, "opt_conv", optimization_has_converged, n_opt_steps, grad_norm)
      ENDIF
      !
      call qes_init_convergence_info(obj, "convergence_info", scf_conv, opt_conv)
      !
      call qes_reset_scf_conv(scf_conv)
      IF (ASSOCIATED(opt_conv)) THEN
         CALL qes_reset_opt_conv(opt_conv)
         NULLIFY ( opt_conv) 
      END IF 
      !
    END SUBROUTINE qexsd_init_convergence_info
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_algorithmic_info(obj, real_space_beta, real_space_q, uspp, paw )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(algorithmic_info_type)   :: obj
      LOGICAL,           INTENT(IN) :: real_space_beta, real_space_q, uspp, paw
      !
      CALL qes_init_algorithmic_info(obj, "algorithmic_info", REAL_SPACE_Q = real_space_q, &
                                     REAL_SPACE_BETA = real_space_beta, USPP = uspp, PAW = paw)
      !
    END SUBROUTINE qexsd_init_algorithmic_info
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_atomic_species(obj, nsp, atm, psfile, amass, starting_magnetization,&
                                         angle1,angle2) 
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(atomic_species_type)    :: obj
      INTEGER,          INTENT(IN) :: nsp
      CHARACTER(len=*), INTENT(IN) :: atm(:)
      CHARACTER(len=*), INTENT(IN) :: psfile(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: amass(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: starting_magnetization(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: angle1(:),angle2(:)
      !
      TYPE(species_type), ALLOCATABLE :: species(:)
      REAL(DP)  :: amass_ = 0.0d0
      REAL(DP)  :: start_mag_ = 0.0d0
      REAL(DP)  :: spin_teta = 0.0d0
      REAL(DP)  :: spin_phi  = 0.0d0
      INTEGER   :: i
      
      ALLOCATE(species(nsp))
      !
      DO i = 1, nsp
          !
          IF ( PRESENT(amass) ) amass_=amass(i)
          IF ( PRESENT(starting_magnetization) ) start_mag_=starting_magnetization(i)
          IF ( PRESENT( angle1 ) )  spin_teta =angle1(i)
          IF ( PRESENT( angle2 ) )  spin_phi = angle2(i)
          !
          CALL qes_init_species( species(i), "species", TRIM(atm(i)),PRESENT(amass),amass_, &
                               TRIM(psfile(i)), PRESENT(starting_magnetization), start_mag_,&
                               PRESENT(angle1),spin_teta,PRESENT(angle2),spin_phi)
      ENDDO
      !
      CALL qes_init_atomic_species(obj, "atomic_species", nsp, SIZE(species), species)
      !
      DO i = 1, nsp
          CALL qes_reset_species(species(i))
      ENDDO
      DEALLOCATE(species)
      !
    END SUBROUTINE qexsd_init_atomic_species
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_atomic_structure(obj, nsp, atm, ityp, nat, tau, &
                                           alat, a1, a2, a3, ibrav)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(atomic_structure_type)  :: obj
      INTEGER,          INTENT(IN) :: nsp, nat
      INTEGER,          INTENT(in) :: ityp(:)
      CHARACTER(LEN=*), INTENT(in) :: atm(:)
      REAL(DP),         INTENT(IN) :: tau(3,*)! cartesian atomic positions, a.u.
      REAL(DP),         INTENT(IN) :: alat
      REAL(DP),         INTENT(IN) :: a1(:), a2(:), a3(:)
      INTEGER,          INTENT(IN) :: ibrav
      !
      INTEGER         :: ia
      TYPE(atom_type), ALLOCATABLE :: atom(:)
      TYPE(cell_type) :: cell
      TYPE(atomic_positions_type)  :: atomic_pos
      TYPE(wyckoff_positions_type) :: wyckoff_pos
      REAL(DP)                     :: new_alat
      LOGICAL                      :: ibrav_ispresent
      !
      ! atomic positions
      !
      IF ( ibrav .gt. 0 ) THEN 
         ibrav_ispresent = .TRUE.
      ELSE
         ibrav_ispresent = .FALSE.
      END IF
      !
      ALLOCATE(atom(nat))
      DO ia = 1, nat
          CALL qes_init_atom( atom(ia), "atom", name=trim(atm(ityp(ia))), &
                             position="", position_ispresent=.FALSE., &
                             atom=tau(1:3,ia), index_ispresent = .TRUE.,&
                             index = ia )
      ENDDO
      !
      CALL qes_init_atomic_positions(atomic_pos, "atomic_positions", SIZE(atom), atom)
      !
      DO ia = 1, nat
          CALL qes_reset_atom( atom(ia) )
      ENDDO
      DEALLOCATE(atom)
      !
      ! cell
      !
      CALL qes_init_cell(cell, "cell", a1, a2, a3)
      !
      ! global init
      !
      CALL qes_init_atomic_structure(obj, "atomic_structure", nat=nat, &
                     alat=alat, alat_ispresent=.TRUE., atomic_positions_ispresent=.TRUE., &
                     atomic_positions=atomic_pos, wyckoff_positions_ispresent=.FALSE., &
                     wyckoff_positions=wyckoff_pos, cell=cell ,& 
                     bravais_index_ispresent = ibrav_ispresent, bravais_index=ibrav)
      ! 
      ! cleanup 
      ! 
      CALL qes_reset_atomic_positions(atomic_pos)
      CALL qes_reset_cell(cell)
      !
    END SUBROUTINE qexsd_init_atomic_structure
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_symmetries(obj, nsym, nrot, space_group, s, ft, sname, t_rev, nat, irt, &
                                     class_names, verbosity, noncolin)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(symmetries_type)    :: obj
      INTEGER,          INTENT(IN) :: nsym, nrot, nat
      INTEGER,          INTENT(IN) :: space_group
      INTEGER,          INTENT(IN) :: s(:,:,:), irt(:,:)
      REAL(DP),         INTENT(IN) :: ft(:,:)
      INTEGER,          INTENT(IN) :: t_rev(:)
      CHARACTER(LEN=*), INTENT(IN) :: sname(:), verbosity
      CHARACTER(LEN=15),INTENT(IN) :: class_names(:)
      LOGICAL,INTENT(IN)           :: noncolin
      !
      TYPE(symmetry_type), ALLOCATABLE  :: symm(:)
      TYPE(equivalent_atoms_type)  :: equiv_atm
      TYPE(info_type)              :: info
      TYPE(matrix_type)            :: matrix
      CHARACTER(LEN=15)            :: classname
      CHARACTER(LEN=256)            :: la_info
      LOGICAL                      :: class_ispresent = .FALSE., time_reversal_ispresent = .FALSE.
      INTEGER                      :: i
      REAL(DP)                     :: mat_(3,3)
      
      ALLOCATE(symm(nrot))
      !
      IF ( TRIM(verbosity) .EQ. 'high' .OR. TRIM(verbosity) .EQ. 'medium')  class_ispresent= .TRUE.
      IF ( noncolin  ) time_reversal_ispresent = .TRUE.
      DO i = 1, nrot
          !
          classname = class_names(i)
          IF ( i .LE. nsym ) THEN 
             la_info = "crystal_symmetry"
          ELSE 
             la_info = "lattice_symmetry"
          END IF
          CALL qes_init_info(info, "info", name=sname(i), name_ispresent=.TRUE., &
                             class=classname, class_ispresent = class_ispresent,   &
                             time_reversal=(t_rev(i)==1), time_reversal_ispresent = time_reversal_ispresent, &
                             INFO= TRIM(la_info) )
          !
          mat_ = real(s(:,:,i),DP)
          CALL qes_init_matrix(matrix, "rotation", DIMS=[3,3], mat=mat_ )
          !
          IF ( i .LE. nsym ) THEN 
             CALL qes_init_equivalent_atoms(equiv_atm, "equivalent_atoms", nat=nat, index_list=irt(i,1:nat)  )
          !
             CALL qes_init_symmetry(symm(i),"symmetry", info=info, rotation=matrix, &
                                 fractional_translation_ispresent=.TRUE., fractional_translation=ft(:,i), &
                                 equivalent_atoms_ispresent=.TRUE., equivalent_atoms=equiv_atm)
          ELSE 
             CALL qes_init_symmetry ( symm(i), "symmetry", INFO = info, ROTATION = matrix, &
                                      FRACTIONAL_TRANSLATION_ISPRESENT = .FALSE., FRACTIONAL_TRANSLATION=ft(:,i), &
                                      EQUIVALENT_ATOMS_ISPRESENT = .FALSE.,  EQUIVALENT_ATOMS=equiv_atm) 
          END IF
          !
          CALL qes_reset_info(info)
          CALL qes_reset_matrix(matrix)
          IF ( i .LT. nsym ) THEN 
             CALL qes_reset_equivalent_atoms( equiv_atm )
          ELSE IF ( i .EQ. nrot ) THEN  
            CALL qes_reset_equivalent_atoms( equiv_atm )
          END IF
          !
      ENDDO
      !
      CALL qes_init_symmetries(obj,"symmetries",NSYM = nsym, NROT=nrot, SPACE_GROUP = space_group, &
                               NDIM_SYMMETRY=SIZE(symm), SYMMETRY=symm )
      !
      DO i = 1, nsym
         CALL qes_reset_symmetry(symm(i))
      ENDDO
      DEALLOCATE(symm)
      !
    END SUBROUTINE qexsd_init_symmetries
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_basis_set(obj, gamma_only, ecutwfc, ecutrho, &
                                    nr1, nr2, nr3, nr1s, nr2s, nr3s, &
                                    fft_box_ispresent, nr1b, nr2b, nr3b, &
                                    ngm, ngms, npwx, b1, b2, b3 )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(basis_set_type)    :: obj
      LOGICAL,          INTENT(IN) :: gamma_only
      INTEGER,          INTENT(IN) :: nr1, nr2, nr3
      INTEGER,          INTENT(IN) :: nr1s, nr2s, nr3s
      LOGICAL,          INTENT(IN) :: fft_box_ispresent
      INTEGER,          INTENT(IN) :: nr1b, nr2b, nr3b
      INTEGER,          INTENT(IN) :: ngm, ngms, npwx
      REAL(DP),         INTENT(IN) :: ecutwfc, ecutrho
      REAL(DP),         INTENT(IN) :: b1(3), b2(3), b3(3)
      !
      TYPE(basisSetItem_type) :: fft_grid
      TYPE(basisSetItem_type) :: fft_smooth
      TYPE(basisSetItem_type) :: fft_box
      TYPE(reciprocal_lattice_type) :: recipr_latt

      CALL qes_init_basisSetItem(fft_grid, "fft_grid", nr1, nr2, nr3, "")
      CALL qes_init_basisSetItem(fft_smooth, "fft_smooth", nr1s, nr2s, nr3s, "")
      CALL qes_init_basisSetItem(fft_box, "fft_box", nr1b, nr2b, nr3b, "" )
      CALL qes_init_reciprocal_lattice(recipr_latt, "reciprocal_lattice", b1, b2, b3)

      CALL qes_init_basis_set(obj, "basis_set", GAMMA_ONLY_ISPRESENT=.TRUE., GAMMA_ONLY=gamma_only, &
                              ECUTWFC=ecutwfc, ECUTRHO_ISPRESENT=.TRUE., ECUTRHO=ecutrho, FFT_GRID=fft_grid, &
                              FFT_SMOOTH_ISPRESENT=.TRUE., FFT_SMOOTH=fft_smooth, &
                              FFT_BOX_ISPRESENT=fft_box_ispresent, FFT_BOX=fft_box, NGM=ngm, &
                              NGMS_ISPRESENT=.TRUE., NGMS=ngms, NPWX=npwx, RECIPROCAL_LATTICE=recipr_latt )
      !
      CALL qes_reset_basisSetItem(fft_grid)
      CALL qes_reset_basisSetItem(fft_smooth)
      CALL qes_reset_basisSetItem(fft_box)
      CALL qes_reset_reciprocal_lattice(recipr_latt)
      !
    END SUBROUTINE qexsd_init_basis_set
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_dft (obj, functional, root_is_output, dft_is_hybrid,&
         nqx1, nqx2, nqx3, ecutfock, exx_fraction, screening_parameter, &
         exxdiv_treatment, x_gamma_extrapolation, ecutvcut,        &
         dft_is_vdW, vdw_corr, nonlocal_term, london_s6, london_c6, &
         london_rcut, xdm_a1, xdm_a2 ,ts_vdw_econv_thr, ts_vdw_isolated, &
         dft_is_lda_plus_U, lda_plus_U_kind, llmax, noncolin, nspin, nsp, &
         nat, species, ityp, Hubbard_U, Hubbard_J0, Hubbard_alpha,  &
         Hubbard_beta, Hubbard_J, starting_ns, U_projection_type, is_hubbard, &
         psd,  Hubbard_ns, Hubbard_ns_nc )
      !------------------------------------------------------------------------
      USE  constants,            ONLY:  eps16
      USE  parameters,           ONLY:  lqmax
      USE  input_parameters,     ONLY:  nspinx
      IMPLICIT NONE
      !
      TYPE(dft_type)    :: obj
      CHARACTER(len=*), INTENT(IN) :: functional, nonlocal_term 
      LOGICAL,          INTENT(IN) :: dft_is_hybrid
      LOGICAL,          INTENT(IN) :: root_is_output
      INTEGER,          INTENT(IN) :: nqx1, nqx2, nqx3
      REAL(DP),         INTENT(IN) :: ecutfock
      REAL(DP),         INTENT(IN) :: exx_fraction
      REAL(DP),         INTENT(IN) :: screening_parameter
      CHARACTER(len=*), INTENT(IN) :: exxdiv_treatment
      LOGICAL,          INTENT(IN) :: x_gamma_extrapolation
      REAL(DP),         INTENT(IN) :: ecutvcut
      !
      LOGICAL,          INTENT(IN) :: dft_is_lda_plus_U, noncolin 
      INTEGER,          INTENT(IN) :: lda_plus_U_kind
      INTEGER,          INTENT(IN) :: llmax, nspin, nsp, nat
      CHARACTER(len=*), INTENT(IN) :: species(nsp)
      INTEGER,          INTENT(IN) :: ityp(nat)
      REAL(DP),         INTENT(IN) :: Hubbard_U(nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_J0(nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_alpha(nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_beta(nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_J(3,nsp)
      REAL(DP),         INTENT(IN) :: starting_ns(lqmax,nspinx,nsp)
      REAL(DP),   OPTIONAL, ALLOCATABLE, INTENT(IN) :: Hubbard_ns(:,:,:,:)
      COMPLEX(DP),OPTIONAL, ALLOCATABLE, INTENT(IN) :: Hubbard_ns_nc(:,:,:,:)
      CHARACTER(len=*), INTENT(IN) :: U_projection_type
      LOGICAL,INTENT(IN)           :: is_hubbard(nsp)
      CHARACTER(LEN=2),INTENT(IN)  :: psd(nsp)
      !
      LOGICAL,          INTENT(IN) :: dft_is_vdW, ts_vdw_isolated
      CHARACTER(len=*), INTENT(IN) :: vdw_corr
      REAL(DP),         INTENT(IN) :: london_s6
      REAL(DP),         INTENT(IN) :: london_rcut
      REAL(DP),         INTENT(IN) :: xdm_a1
      REAL(DP),         INTENT(IN) :: xdm_a2
      REAL(DP),         INTENT(IN) :: london_c6(nsp), ts_vdw_econv_thr
      !
      INTEGER  :: i, is, isp, ind,hubb_l,hubb_n, ldim
      TYPE(hybrid_type) :: hybrid
      TYPE(qpoint_grid_type) :: qpoint_grid
      TYPE(dftU_type) :: dftU
      TYPE(vdW_type) :: vdW
      TYPE(HubbardCommon_type), ALLOCATABLE :: Hubbard_U_(:)
      TYPE(HubbardCommon_type), ALLOCATABLE :: Hubbard_J0_(:)
      TYPE(HubbardCommon_type), ALLOCATABLE :: Hubbard_alpha_(:)
      TYPE(HubbardCommon_type), ALLOCATABLE :: Hubbard_beta_(:)
      TYPE(HubbardJ_type),      ALLOCATABLE :: Hubbard_J_(:)
      TYPE(starting_ns_type),   ALLOCATABLE :: starting_ns_(:)
      TYPE(Hubbard_ns_type),    ALLOCATABLE :: Hubbard_ns_(:)
      TYPE(HubbardCommon_type), ALLOCATABLE :: london_c6_obj(:)
      REAL(DP),                 ALLOCATABLE :: Hubb_occ_aux(:,:) 
      INTEGER                               :: m1, m2
      LOGICAL  :: Hubbard_U_ispresent
      LOGICAL  :: Hubbard_J0_ispresent
      LOGICAL  :: Hubbard_alpha_ispresent
      LOGICAL  :: Hubbard_beta_ispresent
      LOGICAL  :: Hubbard_J_ispresent
      LOGICAL  :: starting_ns_ispresent
      LOGICAL  :: Hubbard_ns_ispresent
      LOGICAL  :: london_c6_ispresent, london_s6_ispresent, london_rvdw_ispresent, ts_vdw_econv_thr_ispresent, & 
                  london_rcut_ispresent, ts_vdw_isolated_ispresent, xdm_a1_ispresent, xdm_a2_ispresent, &
                  empirical_vdw = .FALSE. 
      INTEGER  :: ndim_london_c6, ndim_starting_ns                   
      CHARACTER(10), ALLOCATABLE :: label(:)
      CHARACTER                  :: hubbard_shell 
      INTEGER,EXTERNAL           :: set_hubbard_l,set_hubbard_n
      !
      !
      IF ( dft_is_hybrid ) THEN
          !
          CALL qes_init_qpoint_grid(qpoint_grid, "qpoint_grid", nqx1, nqx2, nqx3, "")
          !
          CALL qes_init_hybrid(hybrid, "hybrid", qpoint_grid, ecutfock, exx_fraction, &
                               screening_parameter, exxdiv_treatment, x_gamma_extrapolation, ecutvcut)
          !
          CALL qes_reset_qpoint_grid(qpoint_grid)
          !
      ENDIF
      !
      IF ( dft_is_lda_plus_U ) THEN
          !
          ALLOCATE(label(nsp))
          DO i = 1, nsp
             IF (is_hubbard(i)) THEN
                hubb_l=set_hubbard_l(psd(i))
                hubb_n=set_hubbard_n(psd(i))
                SELECT CASE ( hubb_l ) 
                   CASE ( 0)  
                      hubbard_shell='s'
                   CASE ( 1 ) 
                      hubbard_shell='p'
                   CASE( 2 ) 
                      hubbard_shell='d'
                   CASE( 3 ) 
                      hubbard_shell='f'
                END SELECT
                WRITE (label(i),'(I0,A)') hubb_n,hubbard_shell
             ELSE
                label(i)="no Hubbard"
             END IF
          END DO
          !
          Hubbard_U_ispresent = (SIZE(Hubbard_U)>0)
          Hubbard_J0_ispresent = (SIZE(Hubbard_J0)>0)
          Hubbard_alpha_ispresent = (SIZE(Hubbard_alpha)>0)
          Hubbard_beta_ispresent = (SIZE(Hubbard_beta)>0)
          Hubbard_J_ispresent = (SIZE(Hubbard_J)>0)
          starting_ns_ispresent = (SIZE(starting_ns)>0)
          !
          ALLOCATE( Hubbard_U_(nsp) )
          ALLOCATE( Hubbard_J0_(nsp) )
          ALLOCATE( Hubbard_alpha_(nsp) )
          ALLOCATE( Hubbard_beta_(nsp) )
          ALLOCATE( Hubbard_J_(nsp) )
          !
          IF (noncolin ) THEN 
             ALLOCATE (starting_ns_(nsp))
             ALLOCATE (Hubbard_ns_(nat))
          ELSE 
            ALLOCATE( starting_ns_(min(nspin,nspinx)*nsp) )
            ALLOCATE( Hubbard_ns_(nspin*nat) )
          END IF
          !
          DO i = 1, nsp
              CALL qes_init_HubbardCommon(Hubbard_U_(i),"Hubbard_U",TRIM(species(i)),TRIM(label(i)),Hubbard_U(i))
              CALL qes_init_HubbardCommon(Hubbard_J0_(i),"Hubbard_J0",TRIM(species(i)),TRIM(label(i)),Hubbard_J0(i))
              CALL qes_init_HubbardCommon(Hubbard_alpha_(i),"Hubbard_alpha",TRIM(species(i)),TRIM(label(i)),&
                                          Hubbard_alpha(i))
              CALL qes_init_HubbardCommon(Hubbard_beta_(i),"Hubbard_beta",TRIM(species(i)),TRIM(label(i)),&
                                          Hubbard_beta(i))
              CALL qes_init_HubbardJ(Hubbard_J_(i),"Hubbard_J",TRIM(species(i)),TRIM(label(i)),Hubbard_J(1:3,i))
          ENDDO
          !
          ind = 0
          IF (starting_ns_ispresent) THEN 
             IF (noncolin) THEN 
                DO i = 1, nsp 
                   IF (.NOT. ANY(starting_ns(1:2*llmax,1,i) > 0)) CYCLE
                   ind = ind + 1 
                   CALL qes_init_starting_ns(starting_ns_(ind), "starting_ns", TRIM (species(i)),TRIM (label(i)),&
                                             1, starting_ns(1:2*llmax, 1, i))
                END DO
             ELSE 
                DO is = 1, MIN(nspin,nspinx) 
                   DO i  = 1, nsp
                      IF (.NOT. ANY (starting_ns(1:llmax,is,i) > 0)) CYCLE 
                      ind = ind+1 
                      CALL qes_init_starting_ns(starting_ns_(ind),"starting_ns",TRIM(species(i)),TRIM(label(i)), &
                                                is, max(starting_ns(1:llmax,is,i),0._DP)  )
                  ENDDO 
                ENDDO 
             END IF
             ndim_starting_ns = ind
          END IF
          IF ( ndim_starting_ns == 0) starting_ns_ispresent = .FALSE.
          !
          ind = 0
          IF (noncolin .AND. PRESENT(Hubbard_ns_nc) ) THEN
             !
             IF (.NOT. ALLOCATED(Hubbard_ns_nc)) &
                CALL errore('qexsd_init_dft', 'Hubbard_ns_nc not alloc', 10)
             !
             Hubbard_ns_ispresent = .TRUE.
             ldim = SIZE(Hubbard_ns_nc,1)
             ALLOCATE (Hubb_occ_aux(2*ldim,2*ldim))
             DO i = 1, nat 
                Hubb_occ_aux = 0.d0
                DO m1 =1, ldim 
                   DO m2 = 1, ldim 
                      Hubb_occ_aux(     m1,     m2) = SQRT(DCONJG(Hubbard_ns_nc(m1,m2,1,i))*Hubbard_ns_nc(m1,m2,1,i))
                      Hubb_occ_aux(     m1,ldim+m2) = SQRT(DCONJG(Hubbard_ns_nc(m1,m2,2,i))*Hubbard_ns_nc(m1,m2,2,i)) 
                      Hubb_occ_aux(ldim+m1,     m2) = SQRT(DCONJG(Hubbard_ns_nc(m1,m2,3,i))*Hubbard_ns_nc(m1,m2,3,i))
                      Hubb_occ_aux(ldim+m1,ldim+m2) = SQRT(DCONJG(Hubbard_ns_nc(m1,m2,4,i))*Hubbard_ns_nc(m1,m2,4,i))
                   END DO
                END DO
                CALL qes_init_Hubbard_ns(Hubbard_ns_(i),"Hubbard_ns_mod", TRIM(species(ityp(i))),TRIM(label(ityp(i))), &
                                         1, i, Hubb_occ_aux(:,:))
             END DO 
             DEALLOCATE ( Hubb_occ_aux) 
          ELSE IF ( PRESENT(Hubbard_ns) ) THEN
             !
             IF (.NOT. ALLOCATED(Hubbard_ns)) &
                CALL errore('qexsd_init_dft', 'Hubbard_ns not alloc', 10)
             !
             Hubbard_ns_ispresent = .TRUE.
             ldim = SIZE(Hubbard_ns,1)
             DO i = 1, nat
                DO is = 1, nspin
                   ind = ind+1
                   CALL qes_init_Hubbard_ns(Hubbard_ns_(ind),"Hubbard_ns", TRIM(species(ityp(i))),TRIM(label(ityp(i))), &
                                       is, i, Hubbard_ns(:,:,is,i) )
                ENDDO
             ENDDO
          ELSE
             Hubbard_ns_ispresent = .FALSE.
          END IF
          !
          ! main init
          CALL qes_init_dftU(dftU, "dftU", .TRUE., lda_plus_u_kind, &
                              Hubbard_U_ispresent, SIZE(Hubbard_U_), Hubbard_U_, &
                              Hubbard_J0_ispresent, SIZE(Hubbard_J0_), Hubbard_J0_, &
                              Hubbard_alpha_ispresent, SIZE(Hubbard_alpha_), Hubbard_alpha_, &
                              Hubbard_beta_ispresent, SIZE(Hubbard_beta_), Hubbard_beta_, &
                              Hubbard_J_ispresent, SIZE(Hubbard_J_), Hubbard_J_, &
                              starting_ns_ispresent, ndim_starting_ns, starting_ns_, &
                              Hubbard_ns_ispresent, SIZE(Hubbard_ns_), Hubbard_ns_, &
                              .TRUE., U_projection_type)
          !
          DO i = 1, nsp
              CALL qes_reset_HubbardCommon(Hubbard_U_(i))
              CALL qes_reset_HubbardCommon(Hubbard_J0_(i))
              CALL qes_reset_HubbardCommon(Hubbard_alpha_(i))
              CALL qes_reset_HubbardCommon(Hubbard_beta_(i))
              CALL qes_reset_HubbardJ(Hubbard_J_(i))
          ENDDO
          !
          DEALLOCATE(Hubbard_U_)
          DEALLOCATE(Hubbard_J0_)
          DEALLOCATE(Hubbard_alpha_)
          DEALLOCATE(Hubbard_beta_)
          DEALLOCATE(Hubbard_J_)
          !
          DO i = 1, SIZE(starting_ns_)
              CALL qes_reset_starting_ns(starting_ns_(i))
          ENDDO
          DEALLOCATE(starting_ns_)
          !
          DO i = 1, SIZE(Hubbard_ns_)
              CALL qes_reset_Hubbard_ns(Hubbard_ns_(i))
          ENDDO
          DEALLOCATE(Hubbard_ns_)
          !
          DEALLOCATE(label)
          !
      ENDIF
      !
      SELECT CASE ( TRIM (vdw_corr )) 
      CASE ( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d') 
           empirical_vdw = .TRUE.
           london_s6_ispresent = .TRUE. 
           london_rcut_ispresent = .TRUE. 
           xdm_a1_ispresent = .TRUE. 
           xdm_a2_ispresent = .TRUE.
           IF ( ANY(london_c6 .GT.  -eps16 )) THEN ! -eps16 to allow london_c6(i) = 0.0 
              london_c6_ispresent = .TRUE.
              ndim_london_c6 = 0 
              DO isp = 1, nsp 
                 IF ( london_c6(isp) .GT. -eps16 ) THEN 
                    ndim_london_c6 = ndim_london_c6 + 1
                 END IF 
              END DO
              ALLOCATE (london_c6_obj(ndim_london_c6))
              ndim_london_c6 = 0 
              DO isp = 1, nsp
                 IF ( london_c6(isp) .GT. -eps16 ) THEN
                    ndim_london_c6 = ndim_london_c6 + 1  
                    CALL qes_init_hubbardcommon(london_c6_obj(ndim_london_c6), "london_c6", TRIM(species(isp)),"",&
                                                london_c6(isp))
                 END IF 
              END DO                        
           ELSE 
              london_c6_ispresent = .FALSE. 
              ALLOCATE ( london_c6_obj(1))
           END IF
           ts_vdw_econv_thr_ispresent = .FALSE. 
           ts_vdw_isolated_ispresent = .FALSE. 
      CASE ( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler')
           empirical_vdw = .TRUE.
           london_s6_ispresent =   .FALSE.
           london_c6_ispresent = .FALSE.
           ALLOCATE ( london_c6_obj(1)) 
           london_rcut_ispresent = .FALSE.
           xdm_a1_ispresent =      .FALSE.
           xdm_a2_ispresent =      .FALSE.
           ts_vdw_econv_thr_ispresent = .TRUE. 
           ts_vdw_isolated_ispresent  = .TRUE. 
      CASE default 
           empirical_vdw = .FALSE.
           ts_vdw_econv_thr_ispresent = .FALSE.
           ts_vdw_isolated_ispresent = .FALSE.
           london_s6_ispresent =   .FALSE.
           london_c6_ispresent = .FALSE.
           ALLOCATE (london_c6_obj(1))
           london_rcut_ispresent = .FALSE.
           xdm_a1_ispresent =      .FALSE.
           xdm_a2_ispresent =      .FALSE.
           london_c6_ispresent   = .FALSE.
      END SELECT 

      IF ( dft_is_vdW .OR. empirical_vdw ) THEN
          !
          CALL qes_init_vdW(vdW, "vdW", TRIM(vdw_corr), root_is_output,  TRIM(nonlocal_term), london_s6_ispresent, london_s6, &
                            ts_vdw_econv_thr_ispresent, ts_vdw_econv_thr, ts_vdw_isolated_ispresent, ts_vdw_isolated,& 
                            london_rcut_ispresent, london_rcut, xdm_a1_ispresent, xdm_a1, xdm_a2_ispresent, xdm_a2, &
                            london_c6_ispresent, ndim_london_c6, london_c6_obj )
          !
          IF (london_c6_ispresent )   THEN
             DO isp=1, ndim_london_c6
                CALL qes_reset_hubbardcommon(london_c6_obj(isp))
             END DO 
          END IF
          DEALLOCATE ( london_c6_obj) 
      ENDIF
        
      CALL qes_init_dft(obj, "dft", functional, dft_is_hybrid, hybrid, &
                             dft_is_lda_plus_U, dftU, (dft_is_vdW .OR. empirical_vdw) , vdW)
      !
      IF (dft_is_hybrid)      CALL qes_reset_hybrid(hybrid)
      IF (dft_is_lda_plus_U)  CALL qes_reset_dftU(dftU)
      IF (dft_is_vdW .OR. empirical_vdw )  CALL qes_reset_vdW(vdW)
      !
    END SUBROUTINE qexsd_init_dft
    !
    !--------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_outputPBC(obj,assume_isolated)
    !--------------------------------------------------------------------------------------------
    ! 
    IMPLICIT NONE
    ! 
    TYPE (outputPBC_type)                       :: obj
    CHARACTER(LEN=*),INTENT(IN)                  :: assume_isolated
    CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="boundary_conditions"
    !
    CALL qes_init_outputPBC(obj,TAGNAME,ASSUME_ISOLATED =assume_isolated)
    END SUBROUTINE qexsd_init_outputPBC
    !
    !
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_magnetization(obj, lsda, noncolin, spinorbit, total_mag, total_mag_nc, &
                                        absolute_mag, do_magnetization)
      !------------------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(magnetization_type)    :: obj
      LOGICAL,         INTENT(IN) :: lsda, noncolin, spinorbit
      REAL(DP),        INTENT(IN) :: total_mag, absolute_mag
      REAL(DP),        INTENT(IN) :: total_mag_nc(3)
      LOGICAL,         INTENT(IN) :: do_magnetization
      !
      CALL qes_init_magnetization(obj, "magnetization", lsda, noncolin, spinorbit, total_mag, absolute_mag, &
                                 do_magnetization)
      !
    END SUBROUTINE qexsd_init_magnetization 
    !
    ! 
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_band_structure(obj, lsda, noncolin, lspinorb, nbnd_up, nbnd_dw, nelec, n_wfc_at, occupations_are_fixed, & 
                                         fermi_energy, two_fermi_energies, ef_updw, et, wg, nks, xk, ngk, wk, & 
                                         starting_kpoints, occupation_kind, smearing, wf_collected)
    !----------------------------------------------------------------------------------------
    IMPLICIT NONE
    !
    TYPE(band_structure_type)               :: obj
    CHARACTER(LEN=*), PARAMETER             :: TAGNAME="band_structure"
    LOGICAL,INTENT(IN)                      :: lsda, noncolin, lspinorb, occupations_are_fixed
    INTEGER,INTENT(IN)                      :: nbnd_up, nbnd_dw, nks, n_wfc_at
    REAL(DP),INTENT(IN)                     :: nelec, fermi_energy
    REAL(DP),DIMENSION(:,:),INTENT(IN)      :: et, wg, xk
    REAL(DP),DIMENSION(:),INTENT(IN)        :: wk
    INTEGER,DIMENSION(:),INTENT(IN)         :: ngk      
    REAL(DP),DIMENSION(2),INTENT(IN)        :: ef_updw 
    LOGICAL,INTENT(IN)                      :: two_fermi_energies
    TYPE(k_points_IBZ_type),INTENT(IN)      :: starting_kpoints
    TYPE(occupations_type), INTENT(IN)      :: occupation_kind
    TYPE(smearing_type),OPTIONAL,INTENT(IN) :: smearing
    LOGICAL,INTENT(IN)                      :: wf_collected                    
    ! 
    LOGICAL                                 :: nbnd_up_ispresent, nbnd_dw_ispresent, &
                                               fermi_energy_ispresent, HOL_ispresent, & 
                                               n_wfc_at_ispresent = .TRUE.  
    INTEGER                                 :: ndim_ks_energies, nbnd, ik
    TYPE(k_point_type)                      :: kp_obj
    TYPE(ks_energies_type),ALLOCATABLE      :: ks_objs(:)
    TYPE (k_points_IBZ_type)                :: starting_k_points_
    TYPE ( occupations_type)                :: occupations_kind_ 
    REAL(DP),DIMENSION(:),ALLOCATABLE       :: eigenvalues, occupations
    TYPE (smearing_type)                    :: smearing_ 
    !
    !
    ndim_ks_energies=nks   
    !
    IF ( lsda ) THEN 
       ndim_ks_energies=ndim_ks_energies/2
       nbnd=nbnd_up+nbnd_dw
       nbnd_up_ispresent=.true.
       nbnd_dw_ispresent=.true.
    ELSE 
       nbnd=nbnd_up
       nbnd_up_ispresent=.false.
       nbnd_dw_ispresent=.false. 
    END IF 
    IF (fermi_energy.GT.-1.D6 .AND. ( .NOT. two_fermi_energies ) ) THEN
      IF ( occupations_are_fixed ) THEN 
         fermi_energy_ispresent = .FALSE. 
         HOL_ispresent = .TRUE. 
      ELSE 
         fermi_energy_ispresent = .TRUE.
         HOL_ispresent = .FALSE. 
      END IF 
    ELSE 
      fermi_energy_ispresent=.FALSE.
      HOL_ispresent = .FALSE.
    END IF  
    !
    !   
    ALLOCATE(eigenvalues(nbnd),occupations(nbnd))
    ALLOCATE(ks_objs(ndim_ks_energies))
    !  
    ks_objs%tagname="ks_energies"
    DO ik=1,ndim_ks_energies
       CALL qes_init_k_point(kp_obj,"k_point",wk(ik),.true.,"",.FALSE., xk(:,ik))
       IF ( lsda ) THEN 
          eigenvalues(1:nbnd_up)=et(1:nbnd_up,ik)/e2
          eigenvalues(nbnd_up+1:nbnd)=et(1:nbnd_dw,ndim_ks_energies+ik)/e2
       ELSE 
          eigenvalues(1:nbnd)= et(1:nbnd,ik)/e2
       END IF
       !
       !
       IF (lsda) THEN 
          IF ( ABS(wk(ik)).GT.1.d-10) THEN 
             occupations(1:nbnd_up)=wg(1:nbnd_up,ik)/wk(ik)
             occupations(nbnd_up+1:nbnd)=wg(1:nbnd_dw,ndim_ks_energies+ik)/wk(ndim_ks_energies+ik)
          ELSE 
             occupations(1:nbnd_up)=wg(1:nbnd_up,ik)
             occupations(nbnd_up+1:nbnd)=wg(1:nbnd_dw,ik) 
          END IF            
       ELSE 
          IF (ABS(wk(ik)).GT.1.d-10) THEN
              occupations(1:nbnd)=wg(1:nbnd,ik)/wk(ik)
          ELSE
              occupations(1:nbnd)=wg(1:nbnd,ik)
          END IF
       END IF
       !
       !
       ks_objs(ik)%k_point = kp_obj
       ks_objs(ik)%npw = ngk(ik)
       CALL  qes_init_vector(ks_objs(ik)%eigenvalues, "eigenvalues",eigenvalues)
       CALL  qes_init_vector(ks_objs(ik)%occupations, "occupations",occupations)
       !
       eigenvalues=0.d0
       occupations=0.d0
       CALL qes_reset_k_point(kp_obj)  
    END DO 
    ks_objs%lwrite = .TRUE.
    ks_objs%lread  = .TRUE.
    !
    IF ( PRESENT(smearing) ) smearing_ = smearing
!
    starting_k_points_ = starting_kpoints
    starting_k_points_%tagname = "starting_k_points"
!
    occupations_kind_  = occupation_kind
    occupations_kind_%tagname = "occupations_kind"
! 
    CALL qes_init_band_structure( obj,TAGNAME,lsda,noncolin,lspinorb, nbnd, nbnd_up_ispresent,&
                  nbnd_up,nbnd_dw_ispresent,nbnd_dw,nelec, n_wfc_at_ispresent, n_wfc_at, wf_collected, & 
                  fermi_energy_ispresent, fermi_energy/e2, HOL_ispresent, fermi_energy/e2,     &
                  two_fermi_energies, ef_updw/e2, starting_k_points_, ndim_ks_energies,      &
                  occupations_kind_, PRESENT(smearing), smearing_, ndim_ks_energies, ks_objs )
    DO ik=1,ndim_ks_energies
       CALL qes_reset_ks_energies(ks_objs(ik))
    END DO
    CALL qes_reset_k_points_IBZ ( starting_k_points_ ) 
    DEALLOCATE (ks_objs,eigenvalues,occupations)
    END SUBROUTINE qexsd_init_band_structure 
    !
    ! 
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_total_energy(obj, etot, eband, ehart, vtxc, etxc, &
         ewald, degauss, demet, electric_field_corr, potentiostat_contr, gate_contribution)
    !----------------------------------------------------------------------------------------
    !
    ! 
    IMPLICIT NONE
    ! 
    TYPE (total_energy_type)        :: obj
    REAL(DP),INTENT(IN)             :: etot, ehart,vtxc,etxc
    REAL(DP),OPTIONAL,INTENT(IN)    :: ewald,demet, eband, degauss 
    REAL(DP),OPTIONAL               :: electric_field_corr
    REAL(DP),OPTIONAL               :: potentiostat_contr
    REAL(DP),OPTIONAL               :: gate_contribution
    !
    LOGICAL                         :: demet_ispresent
    CHARACTER(LEN=*),PARAMETER      :: TAGNAME="total_energy"
    ! 
    CALL  qes_init_total_energy(obj,TAGNAME,etot, EBAND = eband , EHART = ehart, VTXC = vtxc,& 
                               ETXC = etxc , EWALD = ewald, DEMET = demet, &
                               EFIELDCORR=electric_field_corr, POTENTIOSTAT_CONTR = potentiostat_contr,  & 
                               GATE_CONTRIBUTION = gate_contribution )

    END SUBROUTINE qexsd_init_total_energy
    ! 
    !
    !-------------------------------------------------------------------------------------------------------- 
    SUBROUTINE qexsd_init_forces(obj,nat,forces,tprnfor)
    !-------------------------------------------------------------------------------------------------------- 
    !
    IMPLICIT NONE
    !
    TYPE(matrix_type)                            :: obj
    INTEGER,INTENT(IN)                           :: nat 
    REAL(DP),DIMENSION(:,:),INTENT(IN)           :: forces
    LOGICAL,INTENT(IN)                           :: tprnfor
    !
    CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="forces"
    REAL(DP),DIMENSION(:,:),ALLOCATABLE          :: forces_aux
    ! 
    IF (.NOT. tprnfor) THEN
       obj%lwrite=.FALSE.
       obj%lread =.FALSE.
       RETURN
    END IF 
    !
    ALLOCATE (forces_aux(3,nat))
    forces_aux(1:3,1:nat)=forces(1:3,1:nat)/e2
    !
    CALL qes_init_matrix(obj,TAGNAME,[3,nat],forces_aux )
    !
    DEALLOCATE (forces_aux)
    !
    END SUBROUTINE qexsd_init_forces
    ! 
    ! 
    !---------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_stress(obj,stress,tstress) 
    !---------------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE( matrix_type)                           :: obj
    REAL(DP),DIMENSION(3,3),INTENT(IN)           :: stress
    LOGICAL,INTENT(IN)                           :: tstress
    ! 
    CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="stress"
    REAL(DP),DIMENSION(3,3)                      :: stress_aux
    
    IF ( .NOT. tstress ) THEN 
       obj%lwrite = .FALSE.
       obj%lread  = .FALSE.
       stress_aux = 0.d0
       RETURN
    END IF
    ! 
    stress_aux=stress/e2
    CALL qes_init_matrix(obj,TAGNAME,[3,3],stress_aux )
    ! 
    END SUBROUTINE qexsd_init_stress
    !
    !
    !------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_dipole_info (dipole_info, el_dipole, ion_dipole, edir, eamp, emaxpos, eopreg) 
       !------------------------------------------------------------------------------------------------
       ! 
       USE kinds,           ONLY : DP
       USE constants,       ONLY : e2, fpi 
       USE qes_types_module,ONLY : dipoleOutput_type, scalarQuantity_type
       USE qes_libs_module, ONLY : qes_init_scalarQuantity, qes_reset_scalarQuantity
       USE cell_base,       ONLY : alat, at, omega
       ! 
       IMPLICIT NONE  
       ! 
       TYPE ( dipoleOutput_type ), INTENT(OUT)  :: dipole_info
       REAL(DP),INTENT(IN)                      :: el_dipole, ion_dipole, eamp, emaxpos, eopreg
       INTEGER , INTENT(IN)                     :: edir
       ! 
       REAL(DP)                                 :: tot_dipole, length, vamp, fac
       TYPE ( scalarQuantity_type)              :: temp_qobj
       ! 
       tot_dipole = -el_dipole+ion_dipole
       ! 
       dipole_info%idir = edir  
       fac=omega/fpi
       dipole_info%tagname = "dipoleInfo"
       dipole_info%lwrite  = .TRUE.
       dipole_info%lread   = .TRUE.
       CALL qes_init_scalarQuantity(dipole_info%ion_dipole,"ion_dipole" , units="Atomic Units", &
                                    scalarQuantity= ion_dipole*fac)
       CALL qes_init_scalarQuantity(dipole_info%elec_dipole,"elec_dipole" , units="Atomic Units",&
                                     scalarQuantity= el_dipole*fac)
       CALL qes_init_scalarQuantity(dipole_info%dipole,"dipole" , units="Atomic Units", &
                                    scalarQuantity= tot_dipole*fac)
       CALL qes_init_scalarQuantity(dipole_info%dipoleField,"dipoleField" , units="Atomic Units", &
                                    scalarQuantity= tot_dipole)
       ! 
       length=(1._DP-eopreg)*(alat*SQRT(at(1,edir)**2+at(2,edir)**2+at(3,edir)**2))
       vamp=e2*(eamp-tot_dipole)*length
       !
       CALL qes_init_scalarQuantity(dipole_info%potentialAmp,"potentialAmp" , units="Atomic Units",&
                                     scalarQuantity= vamp)
       CALL qes_init_scalarQuantity(dipole_info%totalLength, "totalLength", units = "Bohr",&
                                     scalarQuantity = length )
  
    END SUBROUTINE qexsd_init_dipole_info
    !---------------------------------------------------------------------------------------------
    SUBROUTINE  qexsd_init_outputElectricField(obj, lelfield, tefield, ldipole, lberry, bp_obj, el_pol, &
                                               ion_pol, dipole_obj , gateInfo)
    !---------------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    ! 
    TYPE(outputElectricField_type)                    :: obj 
    ! 
    LOGICAL,INTENT(IN)                                :: lberry, lelfield, tefield, ldipole
    REAL(DP),OPTIONAL,INTENT(IN)                      :: el_pol(:), ion_pol(:)
    TYPE(berryPhaseOutput_type),OPTIONAL,INTENT(IN)   :: bp_obj
    TYPE ( dipoleOutput_type ),OPTIONAL, INTENT(IN)   :: dipole_obj 
    TYPE ( gateInfo_type),OPTIONAL,INTENT(IN)         :: gateInfo
    ! 
    CHARACTER(LEN=*),PARAMETER                        :: TAGNAME="electric_field" 
    TYPE ( berryPhaseOutput_type )                    :: bp_loc_obj
    TYPE ( dipoleOutput_type )                        :: dip_loc_obj
    TYPE ( finiteFieldOut_type )                      :: finiteField_obj
    LOGICAL                                           :: bp_is = .FALSE. , finfield_is = .FALSE. , &
                                                         dipo_is = .FALSE.
    ! 
    
    IF (lberry .AND. PRESENT ( bp_obj))  THEN
       bp_is = .TRUE. 
       bp_loc_obj = bp_obj
    END IF 
    IF ( lelfield .AND. PRESENT(el_pol) .AND. PRESENT (ion_pol ) ) THEN 
       finfield_is=.TRUE.
       CALL qes_init_finiteFieldOut (finiteField_obj, "finiteElectricFieldInfo", el_pol, ion_pol)
    END IF 
    IF ( ldipole .AND. PRESENT( dipole_obj ) ) THEN
       dipo_is = .TRUE.
       dip_loc_obj=dipole_obj
    END IF 
    CALL  qes_init_outputElectricField(obj, TAGNAME, BerryPhase = bp_obj, &
                                    finiteElectricFieldInfo  = finiteField_obj, &
                                    dipoleInfo = dipole_obj, &
                                    GATEINFO =  gateInfo   )
    IF ( finfield_is) CALL qes_reset_finiteFieldOut( finiteField_obj) 
    !
    END SUBROUTINE qexsd_init_outputElectricField
    ! 
    !----------------------------------------------------------------------------------------
    SUBROUTINE qexsd_step_addstep(i_step, max_steps, ntyp, atm, ityp, nat, tau, alat, a1, a2, a3, &
                                  etot, eband, ehart, vtxc, etxc, ewald, degauss, demet, forces,  &
                                  stress, scf_has_converged, n_scf_steps, scf_error, efieldcorr, potstat_contr,      &
                                  fcp_force, fcp_tot_charge, gatefield_en)
    !-----------------------------------------------------------------------------------------
    !! This routing initializes le steps array containing up to max_steps elements of the step_type
    !! data structure. Each element contains structural and energetic info for m.d. trajectories and 
    !! structural minimization paths. All quantities must be provided directly in Hartree atomic units. 
    !! @Note updated on April 10th 2018 by Pietro Delugas
    IMPLICIT NONE 
    ! 
    INTEGER ,INTENT(IN)             :: i_step, max_steps, ntyp, nat, n_scf_steps, ityp(:)
    REAL(DP),INTENT(IN)             :: tau(3,nat), alat, a1(3), a2(3), a3(3), etot, eband, ehart, vtxc, &
                                       etxc, ewald, scf_error, forces(3,nat), stress(3,3) 
    LOGICAL,INTENT(IN)              :: scf_has_converged 
    REAL(DP),OPTIONAL,INTENT(IN)    :: degauss, demet, gatefield_en, efieldcorr
    REAL(DP),OPTIONAL,INTENT (IN)   :: potstat_contr, fcp_force, fcp_tot_charge       
    CHARACTER(LEN=*),INTENT(IN)     :: atm(:)
    TYPE (step_type)                :: step_obj
    TYPE ( scf_conv_type )          :: scf_conv_obj
    TYPE ( atomic_structure_type )  :: atomic_struct_obj
    TYPE ( total_energy_type )      :: tot_en_obj
    TYPE ( matrix_type )            :: mat_forces, mat_stress  
    !    
    IF ( i_step .EQ. 1 ) THEN 
       ALLOCATE (steps(max_steps))
       step_counter = 0
    END IF 
    step_counter = step_counter+1
    !
    step_obj%tagname="step"
    step_obj%n_step = step_counter
    !
    CALL qes_init_scf_conv( scf_conv_obj,"scf_conv", scf_has_converged, n_scf_steps, scf_error )
    !
    step_obj%scf_conv = scf_conv_obj 
    CALL qes_reset_scf_conv(scf_conv_obj)
    ! 
    CALL qexsd_init_atomic_structure(atomic_struct_obj, ntyp, atm, ityp, nat, tau, &
                                     alat, a1, a2, a3, 0)
    step_obj%atomic_structure=atomic_struct_obj
    CALL qes_reset_atomic_structure( atomic_struct_obj )
    ! 
    CALL qexsd_init_total_energy (tot_en_obj, etot, eband, ehart, &
          vtxc, etxc, ewald, degauss, demet, efieldcorr, potstat_contr, gatefield_en)  
    step_obj%total_energy=tot_en_obj
    CALL qes_reset_total_energy( tot_en_obj )
    ! 
    CALL  qes_init_matrix( mat_forces, "forces", [3, nat], forces ) 
    step_obj%forces=mat_forces
    CALL qes_reset_matrix ( mat_forces )
    ! 
    CALL qes_init_matrix( mat_stress, "stress", [3, 3], stress ) 
    step_obj%stress = mat_stress
    CALL qes_reset_matrix ( mat_stress ) 
    IF ( PRESENT ( fcp_force ) ) THEN 
       step_obj%FCP_force = fcp_force
       step_obj%FCP_force_ispresent = .TRUE.
    END IF 
    IF (PRESENT( fcp_tot_charge)) THEN 
       step_obj%FCP_tot_charge = fcp_tot_charge
       step_obj%FCP_tot_charge_ispresent = .TRUE. 
    END IF 
    !  
    ! 
    steps(step_counter) = step_obj
    steps(step_counter)%lwrite  = .TRUE.
    steps(step_counter)%lread   = .TRUE. 
    call qes_reset_step(step_obj)
    END SUBROUTINE qexsd_step_addstep 
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE qexsd_reset_steps()
       IMPLICIT NONE
       INTEGER  :: i_step
       IF (ALLOCATED(steps)) THEN
          DO i_step =1, SIZE(steps) 
            CALL qes_reset_step(steps(i_step))
          END DO
          DEALLOCATE (steps)
      END IF
   END SUBROUTINE
    !-------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_berryPhaseOutput( obj, gpar, gvec, nppstr, nkort, xk, pdl_ion,  &    
                                            mod_ion, pdl_ion_tot, mod_ion_tot, nstring, pdl_elec,  &
                                          mod_elec, wstring, pdl_elec_up, mod_elec_up, pdl_elec_dw,& 
                                          mod_elec_dw, pdl_elec_tot,mod_elec_tot, pdl_tot, mod_tot,&
                                          upol, rmod)
    !---------------------------------------------------------------------------------------------------
    !
    USE ions_base,            ONLY: nat, tau, atm, zv, ityp
    USE cell_base,            ONLY: omega
    USE noncollin_module,     ONLY : noncolin, nspin_lsda
    IMPLICIT NONE 
    ! 
    TYPE (berryPhaseOutput_type)                      :: obj
    REAL(DP),INTENT(IN)                               :: gpar(3), gvec, pdl_ion(nat), pdl_ion_tot, xk(3,*) 

    REAL(DP),INTENT(IN)                               :: pdl_elec(:), pdl_elec_up, pdl_elec_dw, pdl_elec_tot,    & 
                                                         pdl_tot, upol(3), rmod 
    !  
    INTEGER,INTENT(IN)                                :: mod_ion(nat), mod_ion_tot, mod_elec(:), mod_elec_up,    &
                                                         mod_elec_dw, mod_elec_tot, mod_tot, nppstr, nkort, nstring  
    !  
    REAL(DP),INTENT(IN)                               :: wstring(nstring)      
    ! 
    CHARACTER(LEN=*),PARAMETER                        :: TAGNAME = "BerryPhase"
    TYPE ( polarization_type)                         :: tot_pol_obj
    ! 
    TYPE ( electronicPolarization_type),ALLOCATABLE   :: str_pol_obj(:)
    TYPE ( ionicPolarization_type ),    ALLOCATABLE   :: ion_pol_obj(:)
    TYPE ( k_point_type )                             :: kp_obj
    TYPE ( phase_type)                                :: el_phase, ion_phase, tot_phase
    TYPE ( atom_type )                                :: atom_obj
    TYPE ( scalarQuantity_type )                      :: pol_val
    INTEGER                                           :: iat, istring, indstring, ispin 
    CHARACTER(10)                                     :: mod_string
    LOGICAL                                           :: spin_is = .FALSE. 
    !
    ALLOCATE (ion_pol_obj(nat), str_pol_obj(nat)) 
    DO iat =1, nat 
       WRITE(mod_string,'("(mod" ,I1,")")') mod_ion(iat) 
       CALL qes_init_phase(ion_phase,"phase", 0.d0,.FALSE.,0.d0,.FALSE.,TRIM(mod_string),.TRUE., pdl_ion(iat) )
       CALL qes_init_atom(atom_obj,"ion",name=TRIM(atm(ityp(iat))),position_ispresent=.FALSE.,atom = tau(:,iat), &
                          index_ispresent = .FALSE.)
       CALL qes_init_ionicPolarization(ion_pol_obj(iat), "ionicPolarization", atom_obj, zv(ityp(iat)), ion_phase )       
       CALL qes_reset_phase(ion_phase)
       CALL qes_reset_atom(atom_obj)
    END DO
    ! 
    IF ( nspin_lsda .EQ. 2 ) spin_is  = .TRUE.
    DO  istring= 1, nstring
        indstring = 1+(istring-1)*nppstr
        WRITE(mod_string,'("(mod ",I1,")")') mod_elec(istring)
        CALL qes_init_phase(el_phase, "phase", 0.d0, .FALSE., 0.d0, .FALSE., TRIM (mod_string), .TRUE., &
                            pdl_elec(istring) )
        IF (istring .LE. nstring/nspin_lsda) THEN 
           ispin = 1 
        ELSE 
           ispin = 2 
        END IF
        CALL qes_init_k_point(kp_obj, "firstKeyPoint", wstring(istring), .TRUE., "",.FALSE., xk(:,indstring))
        CALL qes_init_electronicPolarization(str_pol_obj(istring),"electronicPolarization", kp_obj, spin_is, ispin, &
                                             el_phase )
        CALL qes_reset_phase ( el_phase ) 
        CALL qes_reset_k_point(kp_obj)
    END DO
    ! 
    WRITE(mod_string,'("(mod ",I1,")")') mod_tot
    CALL qes_init_phase(tot_phase, "totalPhase", pdl_ion_tot, .TRUE. , pdl_elec_tot, .TRUE., TRIM(mod_string), &
                        .TRUE., pdl_tot)
    ! 
    CALL qes_init_scalarQuantity ( pol_val, "polarization", Units="e/bohr^2", scalarQuantity=(rmod/omega)*pdl_tot )
    !
    CALL qes_init_polarization(tot_pol_obj, "totalPolarization", pol_val, modulus = (rmod/omega)*dble(mod_tot), &
                               direction = upol )  
    ! 
    CALL qes_init_berryPhaseOutput( obj, TAGNAME, tot_pol_obj, tot_phase, nat, ion_pol_obj, nstring, str_pol_obj )
    ! 
    DO istring=1,nstring     
       CALL  qes_reset_electronicPolarization(str_pol_obj(istring))
    END DO 
    DEALLOCATE (str_pol_obj)
    DO iat=1, nat
       CALL qes_reset_ionicPolarization(ion_pol_obj(iat))
    END DO
    DEALLOCATE (ion_pol_obj)
    CALL qes_reset_polarization(tot_pol_obj)
    CALL qes_reset_scalarQuantity(pol_val)
    CALL qes_reset_phase(tot_phase) 
    !
    END SUBROUTINE qexsd_init_berryPhaseOutput
    !
    !-------------------------------------------------------------------------------------------------         
    SUBROUTINE qexsd_set_status(status_int)
    !-------------------------------------------------------------------------------------------------
    IMPLICIT NONE 
    !
    INTEGER      :: status_int
    END SUBROUTINE qexsd_set_status 
    !
    !--------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_set_closed() 
    ! 
    IMPLICIT NONE 
    CHARACTER(LEN=9)                  :: cdate, time_string
    CHARACTER(LEN=12)                 :: date_string
    !
    CALL date_and_tim( cdate, time_string ) 
    date_string = cdate(1:2) // ' ' // cdate(3:5) // ' ' // cdate (6:9)
    CALL qes_init_closed (qexsd_closed_element, "closed", date_string, time_string,&
                          "")
    END SUBROUTINE qexsd_set_closed 
    
!-----------------------------------------------------------------------------------
SUBROUTINE qexsd_init_gate_info(obj, tagname, gatefield_en, zgate_, nelec_, alat_, at_, bg_, zv_, ityp_) 
   !--------------------------------------------------------------------------------
   USE kinds,         ONLY : DP
   USE constants,     ONLY : tpi
   !
   IMPLICIT NONE
   TYPE (gateInfo_type),INTENT(INOUT)      :: obj;
   CHARACTER(LEN=*)                        :: tagname
   REAL(DP), INTENT(IN)                    :: gatefield_en, zgate_, alat_, at_(3,3), bg_(3,3), zv_(:), nelec_ 
   INTEGER,INTENT(IN)                      :: ityp_(:) 
   ! 
   REAL(DP)                                :: bmod, area, ionic_charge, gateamp, gate_gate_term
   ! 
   bmod=SQRT(bg_(1,3)**2+bg_(2,3)**2+bg_(3,3)**2)
   ionic_charge = SUM( zv_(ityp_(:)) )
   area = ABS((at_(1,1)*at_(2,2)-at_(2,1)*at_(1,2))*alat_**2) 
   gateamp = (-(nelec_-ionic_charge)/area*tpi)
   gate_gate_term =  (- (nelec_-ionic_charge) * gateamp * (alat_/bmod) / 6.0)
   obj = gateInfo_type( TAGNAME = TRIM(tagname), lwrite = .TRUE., lread = .FALSE., POT_PREFACTOR = gateamp, &
                        GATE_ZPOS = zgate_,  GATE_GATE_TERM = gate_gate_term, GATEFIELDENERGY = gatefield_en) 
   ! 
END SUBROUTINE qexsd_init_gate_info 



END MODULE qexsd_module

!
!----------------
! ... dummy defs
!----------------
!

! 
! 

