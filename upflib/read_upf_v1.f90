
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
      MODULE read_upf_v1_module
!=----------------------------------------------------------------------------=!

!  this module handles the reading of pseudopotential data

! ...   declare modules
        USE upf_kinds,    ONLY: DP
        USE upf_io,       ONLY: stdout
        USE pseudo_types, ONLY: pseudo_upf
        USE upf_utils,    ONLY: matches
!
        IMPLICIT NONE
        SAVE
        PRIVATE
        PUBLIC :: read_upf_v1, scan_begin, scan_end
      CONTAINS
!
!---------------------------------------------------------------------
SUBROUTINE read_upf_v1 ( file_pseudo, upf, ierr )
  !---------------------------------------------------------------------
  !
  !   read pseudopotential "upf" in the Unified Pseudopotential Format
  !   from file "file_pseudo" - File is closed on output.
  !   return error code in "ierr" (success: ierr=0)
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: file_pseudo
  INTEGER, INTENT(OUT) :: ierr 
  TYPE (pseudo_upf), INTENT(INOUT)       :: upf
  !
  !     Local variables
  !
  INTEGER :: ios, iunps
  CHARACTER (len=80) :: dummy  
  !
  ! Open the file
  !
  OPEN ( NEWUNIT = iunps, FILE = file_pseudo, STATUS = 'old', &
       FORM = 'formatted', iostat=ios )
  IF ( ios /= 0 ) GO TO 200
  !
  ! First check if this pseudo-potential has spin-orbit or GIPAW information
  !
  upf%q_with_l=.false.
  upf%has_so=.false.
  upf%has_gipaw = .false.
  addinfo_loop: do while (ios == 0)  
     read (iunps, *, iostat = ios, err = 200) dummy  
     if (matches ("<PP_ADDINFO>", dummy) ) then
        upf%has_so=.true.
     endif
     if ( matches ( "<PP_GIPAW_RECONSTRUCTION_DATA>", dummy ) ) then
        upf%has_gipaw = .true.
     endif
     if (matches ("<PP_QIJ_WITH_L>", dummy) ) then
        upf%q_with_l=.true. 
     endif
  enddo addinfo_loop
  
  !------->Search for Header
  !     This version doesn't use the new routine scan_begin
  !     because this search must set extra flags for
  !     compatibility with other pp format reading

  ios = 0
  rewind(iunps)
  header_loop: do while (ios == 0)  
     read (iunps, *, iostat = ios, err = 200) dummy  
     if (matches ("<PP_HEADER>", dummy) ) then  
        call read_pseudo_header (iunps, upf, ierr)  
        exit header_loop
     endif
  enddo header_loop
  if (ierr > 0) GO TO 200
  !
  ! this should be read from the PP_INFO section
  !
  upf%generated='Generated by new atomic code, or converted to UPF format'
  
  call scan_end (iunps, "HEADER",ierr)  
  if (ierr > 0) GO TO 200

  ! Compatibility with later formats:
  upf%has_wfc = .false.

  !-------->Search for mesh information
  call scan_begin (iunps, "MESH", .true., ierr)  
  if ( ierr > 0 ) go to 200
  call read_pseudo_mesh (iunps, upf, ierr)  
  if ( ierr > 0 ) go to 200
  call scan_end (iunps, "MESH",ierr)  
  if (ierr > 0) GO TO 200
  !-------->If  present, search for nlcc
  if ( upf%nlcc ) then  
     call scan_begin (iunps, "NLCC", .true.,ierr)  
     if ( ierr > 0 ) go to 200
     call read_pseudo_nlcc (iunps, upf, ierr)  
     if ( ierr > 0 ) go to 200
     call scan_end (iunps, "NLCC",ierr)  
     if ( ierr > 0 ) go to 200
  else
     ALLOCATE( upf%rho_atc( upf%mesh ) )
     upf%rho_atc = 0.0_DP
  endif
  !-------->Fake 1/r potential: do not read PP
  if (.not. matches ("1/r", upf%typ) ) then
  !-------->Search for Local potential
     call scan_begin (iunps, "LOCAL", .true.,ierr)  
     if ( ierr > 0 ) go to 200
     call read_pseudo_local (iunps, upf, ierr)  
     if ( ierr > 0 ) goto 200
     call scan_end (iunps, "LOCAL",ierr)  
     if ( ierr > 0 ) goto 200
  !-------->Search for Nonlocal potential
     call scan_begin (iunps, "NONLOCAL", .true., ierr)  
     if ( ierr > 0 ) goto 200
     call read_pseudo_nl (iunps, upf, ierr)  
     if ( ierr > 0 ) go to 200
     call scan_end (iunps, "NONLOCAL", ierr)
     if ( ierr > 0 ) go to 200
  !--------
  else
     call set_coulomb_nonlocal(upf)
  end if
  !-------->Search for atomic wavefunctions
  call scan_begin (iunps, "PSWFC", .true.,ierr)  
  if ( ierr > 0 ) go to 200
  call read_pseudo_pswfc (iunps, upf, ierr)
  if ( ierr > 0 ) go to 200
  call scan_end (iunps, "PSWFC",ierr)  
  if ( ierr > 0 ) go to 200
  !-------->Search for atomic charge
  call scan_begin (iunps, "RHOATOM", .true.,ierr)  
  if ( ierr > 0 ) go to 200
  call read_pseudo_rhoatom (iunps, upf, ierr)
  if ( ierr > 0 ) go to 200
  call scan_end (iunps, "RHOATOM",ierr)  
  if ( ierr > 0 ) go to 200
  !-------->Search for add_info
  if (upf%has_so) then
     call scan_begin (iunps, "ADDINFO", .true.,ierr)  
     if ( ierr > 0 ) go to 200
     call read_pseudo_addinfo (iunps, upf, ierr)  
     if ( ierr > 0 ) go to 200
     call scan_end (iunps, "ADDINFO",ierr)  
     if ( ierr > 0 ) go to 200
  endif
  !-------->GIPAW data
  IF ( upf%has_gipaw ) then
     CALL scan_begin ( iunps, "GIPAW_RECONSTRUCTION_DATA", .false.,ierr )
     if ( ierr > 0 ) go to 200
     CALL read_pseudo_gipaw ( iunps, upf, ierr )
     if ( ierr > 0 ) go to 200
     CALL scan_end ( iunps, "GIPAW_RECONSTRUCTION_DATA",ierr )
     if ( ierr > 0 ) go to 200
  END IF
  !--- Try to get the core radius if not present. Needed by the 
  !    atomic code for old pseudo files
  IF (upf%nbeta>0) THEN ! rcutus may be unallocated if nbeta=0
     IF(upf%rcutus(1)<1.e-9_DP) THEN 
        call scan_begin (iunps, "INFO", .true.,ierr)  
        if ( ierr > 0 ) go to 200
        call read_pseudo_ppinfo (iunps, upf, ierr)  
        if ( ierr > 0 ) go to 200
        call scan_end (iunps, "INFO",ierr)
        if ( ierr > 0 ) go to 200
     ENDIF
  ENDIF

200 CLOSE (UNIT=iunps)

end subroutine read_upf_v1
!---------------------------------------------------------------------
subroutine scan_begin (iunps, string, rew, ierr )  
  !---------------------------------------------------------------------
  !
  implicit none
  ! Unit of the input file
  integer, intent(in) :: iunps  
  ! Label to be matched
  character (len=*), intent(in) :: string  
  ! Flag if .true. rewind the file
  logical, intent(in) :: rew  
  ! Return error code
  integer, intent(out), optional :: ierr
  ! String read from file
  character (len=75) :: rstring  
  integer :: ios

  ios = 0
  if (rew) rewind (iunps)  
  do while (ios==0)  
     read (iunps, *, iostat = ios, err = 300) rstring  
     if (matches ("<PP_"//string//">", rstring) ) then
        if ( present(ierr) ) ierr = ios
        return  
     end if
  enddo
  return
300 write(stdout, '("scan_begin: No ",a," block")') trim(string) 
  if ( present(ierr) ) ierr = 1
  !
end subroutine scan_begin
!
!---------------------------------------------------------------------
subroutine scan_end (iunps, string, ios)  
  !---------------------------------------------------------------------
  implicit none
  ! Unit of the input file
  integer :: iunps
  ! Label to be matched
  character (len=*) :: string  
  ! String read from file
  character (len=75) :: rstring
  ! Return error code
  integer, intent(out), optional:: ios

  if (present(ios)) ios = 0 
  read (iunps, '(a)', end = 300, err = 300) rstring  
  if (matches ("</PP_"//string//">", rstring) ) return  
  return
  300 if (present(ios)) ios = 1
  write(stdout, '("scan_end: No ",a," end statement, corrupted file?")') trim(string) 
  !
end subroutine scan_end
!
!---------------------------------------------------------------------
subroutine read_pseudo_header (iunps, upf, ierr)
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  integer :: nw,pqn(0:4),l
  character (len=80) :: dummy
  !
  ierr = 1
  ! GTH analytical format: obviously not true in this case
  upf%is_gth=.false.
  ! PP is assumed to be multi-projector
  upf%is_multiproj = .true. 
  ! Version number (presently ignored)
  read (iunps, *, err = 100, end = 100) upf%nv , dummy  
  ! Element label
  read (iunps, *, err = 100, end = 100) upf%psd , dummy  
  ! Type of pseudo (1/r cannot be read with default format!!!)
  read (iunps, '(a80)', err = 100, end = 100) dummy
  upf%typ=trim(adjustl(dummy))
  !
  if (matches ('US', upf%typ) ) then
     upf%tvanp = .true.  
     upf%tpawp = .false.  
     upf%tcoulombp = .false.
  else if (matches ('PAW', upf%typ) ) then
     ! Note: if tvanp is set to false the results are wrong!
     upf%tvanp = .true.  
     upf%tpawp = .true.  
     upf%tcoulombp = .false.
  else if (matches ('NC', upf%typ) ) then
     upf%tvanp = .false.  
     upf%tpawp = .false.  
     upf%tcoulombp = .false.
  else if (matches ('1/r', upf%typ) ) then
     upf%tvanp = .false.  
     upf%tpawp = .false.
     upf%tcoulombp = .true.
  else
     write(stdout,'("read_pseudo_header: unknown pseudo type")')
     return
  endif

  read (iunps, *, err = 100, end = 100) upf%nlcc , dummy  

  read (iunps, '(a20,t24,a)', err = 100, end = 100) upf%dft, dummy  

  read (iunps, * ) upf%zp , dummy  
  read (iunps, * ) upf%etotps, dummy  
  read (iunps, * ) upf%ecutwfc, upf%ecutrho
  read (iunps, * ) upf%lmax , dummy
  read (iunps, *, err = 100, end = 100) upf%mesh , dummy
  read (iunps, *, err = 100, end = 100) upf%nwfc, upf%nbeta , dummy
  read (iunps, '(a)', err = 100, end = 100) dummy
  ALLOCATE( upf%els( upf%nwfc ), upf%lchi( upf%nwfc ), upf%oc( upf%nwfc ),&
       upf%nchi( upf%nwfc ) )
  ! set default n
  do l = 0, 4
    pqn(l) = l + 1
  enddo
  do nw = 1, upf%nwfc  
     read (iunps, * ) upf%els (nw), upf%lchi (nw), upf%oc (nw)  
     dummy = trim(adjustl(upf%els (nw)))
     if (dummy(1:1) == '1') pqn(upf%lchi (nw)) = 1
     if (dummy(1:1) == '2') pqn(upf%lchi (nw)) = 2
     if (dummy(1:1) == '3') pqn(upf%lchi (nw)) = 3
     if (dummy(1:1) == '4') pqn(upf%lchi (nw)) = 4
     if (dummy(1:1) == '5') pqn(upf%lchi (nw)) = 5
     if (dummy(1:1) == '6') pqn(upf%lchi (nw)) = 6
     if (dummy(1:1) == '7') pqn(upf%lchi (nw)) = 7
     upf%nchi (nw) = pqn(upf%lchi (nw))
     pqn(upf%lchi (nw)) = pqn(upf%lchi (nw)) + 1
     !upf%nchi (nw) = upf%lchi(nw)+1
  enddo
  ! next lines for compatibility with upf v.2
  ALLOCATE( upf%rcut_chi( upf%nwfc ), upf%rcutus_chi( upf%nwfc ), &
       upf%epseu( upf%nwfc ) )
  upf%rcut_chi = 0.0_dp; upf%rcutus_chi=0.0_dp; upf%epseu = 0.0_dp

  ierr = 0 
  return  
  !
100 write(stdout,'("read_pseudo_header: error reading pseudo file")')
  !
end subroutine read_pseudo_header
!
!---------------------------------------------------------------------
subroutine read_pseudo_mesh (iunps, upf, ierr)  
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  integer :: ir
  !
  ierr = 1
  ALLOCATE( upf%r( upf%mesh ), upf%rab( upf%mesh ) )
  upf%r   = 0.0_DP
  upf%rab = 0.0_DP

  call scan_begin (iunps, "R", .false.)  
  read (iunps, *, err = 100, end = 100) (upf%r(ir), ir=1,upf%mesh )
  call scan_end (iunps, "R")  
  call scan_begin (iunps, "RAB", .false.)  
  read (iunps, *, err = 101, end = 101) (upf%rab(ir), ir=1,upf%mesh )
  call scan_end (iunps, "RAB")  

  ierr = 0 
  return  

100 write(stdout,'("read_pseudo_mesh: error reading file (R) for ",a)') upf%psd
  return
101 write(stdout,'("read_pseudo_mesh: error reading file (RAB) for ",a)') upf%psd
  !
end subroutine read_pseudo_mesh
!---------------------------------------------------------------------
subroutine read_pseudo_nlcc (iunps, upf, ierr)
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  integer :: ir
  !
  ierr = 1
  ALLOCATE( upf%rho_atc( upf%mesh ) )
  upf%rho_atc = 0.0_DP
  read (iunps, *, err = 100, end = 100) (upf%rho_atc(ir), ir=1,upf%mesh )
  ierr = 0
  return
100 write(stdout,'("read_pseudo_nlcc: error reading pseudo file")')
  !
end subroutine read_pseudo_nlcc

!---------------------------------------------------------------------
subroutine read_pseudo_local (iunps, upf, ierr)
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  integer :: ir
  !
  ALLOCATE( upf%vloc( upf%mesh ) )
  upf%vloc = 0.0_DP
  ierr = 1
  read (iunps, *, err=100, end=100) (upf%vloc(ir) , ir=1,upf%mesh )
  ierr = 0
  return
100 write(stdout,'("read_pseudo_local: error reading pseudo file")')
end subroutine read_pseudo_local
!
!---------------------------------------------------------------------
subroutine read_pseudo_nl (iunps, upf, ierr)
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  integer :: nb, mb, ijv, n, ir, ios, idum, ldum, icon, lp, i, ikk, l, l1,l2, nd
  ! counters
  character (len=75) :: dummy  
  !
  ierr = 1
  ! Threshold for qfunc to be considered zero (inserted in version UPF v2)
  upf%qqq_eps = -1._dp
  !
  if ( upf%nbeta == 0) then
     upf%nqf = 0
     upf%nqlc= 0
     upf%kkbeta = 0
     ALLOCATE( upf%kbeta( 1 ) )
     ALLOCATE( upf%lll( 1 ) )
     ALLOCATE( upf%beta( upf%mesh, 1 ) )
     ALLOCATE( upf%dion( 1, 1 ) )
     ALLOCATE( upf%rinner( 1 ) )
     ALLOCATE( upf%qqq   ( 1, 1 ) )
     ALLOCATE( upf%qfunc ( upf%mesh, 1 ) )
     ALLOCATE( upf%qfcoef( 1, 1, 1, 1 ) )
     ALLOCATE( upf%rcut( 1 ) )
     ALLOCATE( upf%rcutus( 1 ) )
     ALLOCATE( upf%els_beta( 1 ) )
     ierr = 0
     return
  end if
  ALLOCATE( upf%kbeta( upf%nbeta ) )
  ALLOCATE( upf%lll( upf%nbeta ) )
  ALLOCATE( upf%beta( upf%mesh, upf%nbeta ) )
  ALLOCATE( upf%dion( upf%nbeta, upf%nbeta ) )
  ALLOCATE( upf%rcut( upf%nbeta ) )
  ALLOCATE( upf%rcutus( upf%nbeta ) )
  ALLOCATE( upf%els_beta( upf%nbeta ) )

  upf%kkbeta = 0  
  upf%lll    = 0  
  upf%beta   = 0.0_DP
  upf%dion   = 0.0_DP
  upf%rcut   = 0.0_DP
  upf%rcutus = 0.0_DP
  upf%els_beta = '  '

  do nb = 1, upf%nbeta 
     call scan_begin (iunps, "BETA", .false.)  
     read (iunps, *, err = 100, end = 100) idum, upf%lll(nb), dummy
     read (iunps, *, err = 100, end = 100) ikk  
     upf%kbeta(nb) = ikk
     upf%kkbeta = MAX ( upf%kkbeta, upf%kbeta(nb) )  
     read (iunps, *, err = 100, end = 100) (upf%beta(ir,nb), ir=1,ikk)

     read (iunps, *, err=200,iostat=ios) upf%rcut(nb), upf%rcutus(nb)
     read (iunps, *, err=200,iostat=ios) upf%els_beta(nb)
     call scan_end (iunps, "BETA")  
200  continue
  enddo

  call scan_begin (iunps, "DIJ", .false.)  
  read (iunps, *, err = 101, end = 101) nd, dummy  
  do icon = 1, nd
     !! FIXME: dangerous syntax, are we sure mb has the expected value?
     read (iunps, *, err = 101, end = 101) nb, mb, upf%dion(nb,mb)
     upf%dion (mb,nb) = upf%dion (nb,mb)  
  enddo
  call scan_end (iunps, "DIJ")  

  if ( upf%tvanp .or. upf%tpawp) then
     call scan_begin (iunps, "QIJ", .false.)  
     read (iunps, *, err = 102, end = 102) upf%nqf
     upf%nqlc = 2 * upf%lmax  + 1
     ALLOCATE( upf%rinner( upf%nqlc ) )
     ALLOCATE( upf%qqq   ( upf%nbeta, upf%nbeta ) )
     IF (upf%q_with_l .or. upf%tpawp) then
        ALLOCATE( upf%qfuncl ( upf%mesh, upf%nbeta*(upf%nbeta+1)/2, 0:2*upf%lmax ) )
        upf%qfuncl  = 0.0_DP
     ELSE
        ALLOCATE( upf%qfunc ( upf%mesh, upf%nbeta*(upf%nbeta+1)/2 ) )
        upf%qfunc  = 0.0_DP
     ENDIF
     IF ( upf%nqf > 0 ) THEN
        ALLOCATE( upf%qfcoef(upf%nqf, upf%nqlc, upf%nbeta, upf%nbeta ) )
     ELSE
        ALLOCATE( upf%qfcoef(1,1,1,1) )
     ENDIF
     upf%rinner = 0.0_DP
     upf%qqq    = 0.0_DP
     upf%qfcoef = 0.0_DP
     if ( upf%nqf > 0) then
        call scan_begin (iunps, "RINNER", .false.)  
        read (iunps,*,err=103,end=103) ( idum, upf%rinner(i), i=1,upf%nqlc )
        call scan_end (iunps, "RINNER")  
     end if
     do nb = 1, upf%nbeta
        do mb = nb, upf%nbeta

           read (iunps,*,err=102,end=102) idum, idum, ldum, dummy
           !"  i    j   (l)"
           if (ldum /= upf%lll(mb) ) go to 102

           read (iunps,*,err=104,end=104) upf%qqq(nb,mb), dummy
           ! "Q_int"
           upf%qqq(mb,nb) = upf%qqq(nb,mb)  
           ! ijv is the combined (nb,mb) index
           ijv = mb * (mb-1) / 2 + nb
           IF (upf%q_with_l .or. upf%tpawp) THEN
              l1=upf%lll(nb)
              l2=upf%lll(mb)
              DO l=abs(l1-l2),l1+l2
                 read (iunps, *, err=105, end=105) (upf%qfuncl(n,ijv,l), &
                                                    n=1,upf%mesh)
              END DO
           ELSE
              read (iunps, *, err=105, end=105) (upf%qfunc(n,ijv), n=1,upf%mesh)
           ENDIF

           if ( upf%nqf > 0 ) then
              call scan_begin (iunps, "QFCOEF", .false.)  
              read (iunps,*,err=106,end=106) &
                        ( ( upf%qfcoef(i,lp,nb,mb), i=1,upf%nqf ), lp=1,upf%nqlc )
              do i = 1, upf%nqf
                 do lp = 1, upf%nqlc
                    upf%qfcoef(i,lp,mb,nb) = upf%qfcoef(i,lp,nb,mb)
                 end do
              end do
              call scan_end (iunps, "QFCOEF")  
           end if

        enddo
     enddo
     call scan_end (iunps, "QIJ")  
  else
     upf%nqf  = 1
     upf%nqlc = 2 * upf%lmax  + 1
     ALLOCATE( upf%rinner( upf%nqlc ) )
     ALLOCATE( upf%qqq   ( upf%nbeta, upf%nbeta ) )
     ALLOCATE( upf%qfunc ( upf%mesh, upf%nbeta*(upf%nbeta+1)/2 ) )
     ALLOCATE( upf%qfcoef( upf%nqf, upf%nqlc, upf%nbeta, upf%nbeta ) )
     upf%rinner = 0.0_DP
     upf%qqq    = 0.0_DP
     upf%qfunc  = 0.0_DP
     upf%qfcoef = 0.0_DP
  endif
  ierr = 0
  return  

100 write(stdout,'("read_pseudo_nl: error reading pseudo file (BETA)")'  )
  return
101 write(stdout,'("read_pseudo_nl: error reading pseudo file (DIJ) ")'  )  
  return
102 write(stdout,'("read_pseudo_nl: error reading pseudo file (QIJ) ")'  )
  return
103 write(stdout,'("read_pseudo_nl: error reading pseudo file (RINNER)")')
  return
104 write(stdout,'("read_pseudo_nl: error reading pseudo file (qqq)")'   )
  return
105 write(stdout,'("read_pseudo_nl: error reading pseudo file (qfunc)")' )
  return
106 write(stdout,'("read_pseudo_nl: error reading pseudo file (qfcoef)")')
  !
end subroutine read_pseudo_nl
!---------------------------------------------------------------------
subroutine read_pseudo_pswfc (iunps, upf, ierr)
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  character (len=75) :: dummy  
  integer :: nb, ir

  ierr = 1
  ALLOCATE( upf%chi( upf%mesh, upf%nwfc ) )
  upf%chi = 0.0_DP
  do nb = 1, upf%nwfc  
     read (iunps, *, err=100, end=100) dummy  !Wavefunction labels
     read (iunps, *, err=100, end=100) ( upf%chi(ir,nb), ir=1,upf%mesh )
  enddo
  ierr = 0
  return  
100 write(stdout,'("read_pseudo_pswfc: error reading pseudo file")')
  !
end subroutine read_pseudo_pswfc

!---------------------------------------------------------------------
subroutine read_pseudo_rhoatom (iunps, upf, ierr)
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  integer :: ir
  !
  ierr =  1
  ALLOCATE( upf%rho_at( upf%mesh ) )
  upf%rho_at = 0.0_DP
  read (iunps,*,err=100,end=100) ( upf%rho_at(ir), ir=1,upf%mesh )
  ierr = 0
  return  
100 write(stdout,'("read_pseudo_rhoatom: error reading pseudo file")')
  !
end subroutine read_pseudo_rhoatom
!
!---------------------------------------------------------------------
subroutine read_pseudo_addinfo (iunps, upf, ierr)
!---------------------------------------------------------------------
!
!     This routine reads from the new UPF file,
!     and the total angular momentum jjj of the beta and jchi of the
!     wave-functions.
!
  implicit none
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  integer :: nb
  
  ALLOCATE( upf%jchi(upf%nwfc) )
  ALLOCATE( upf%jjj(upf%nbeta) )

  ierr = 1
  upf%jchi=0.0_DP
  do nb = 1, upf%nwfc
     read (iunps, *,err=100,end=100) upf%els(nb),  &
          upf%nchi(nb), upf%lchi(nb), upf%jchi(nb), upf%oc(nb)
     if ( abs ( upf%jchi(nb)-upf%lchi(nb)-0.5_dp ) > 1.0d-7 .and. &
          abs ( upf%jchi(nb)-upf%lchi(nb)+0.5_dp ) > 1.0d-7      ) then
        write(stdout,'(5x,"read_pseudo_upf: obsolete ADDINFO section ignored")')
        upf%has_so = .false.
        return
     end if
  enddo
  
  upf%jjj=0.0_DP
  do nb = 1, upf%nbeta
     read (iunps, *, err=100,end=100) upf%lll(nb), upf%jjj(nb)
     if ( abs ( upf%lll(nb)-upf%jjj(nb)-0.5_dp) > 1.0d-7 .and. &
          abs ( upf%lll(nb)-upf%jjj(nb)+0.5_dp) > 1.0d-7       ) then
        write(stdout,'(5x,"read_pseudo_upf: obsolete ADDINFO section ignored")')
        upf%has_so = .false.
        return
     end if
  enddo
  
  read(iunps, *) upf%xmin, upf%rmax, upf%zmesh, upf%dx
  ierr = 0
  return
100 write(stdout,'("read_pseudo_addinfo: error reading pseudo file")')
  !
end subroutine read_pseudo_addinfo
!
!---------------------------------------------------------------------
SUBROUTINE read_pseudo_gipaw ( iunps, upf, ierr )
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE ( pseudo_upf ), INTENT ( INOUT ) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  REAL (dp) :: version
  !
  ierr = 1
  CALL scan_begin ( iunps, "GIPAW_FORMAT_VERSION", .false. )
  READ ( iunps, *, err=100, end=100 ) version
  upf%gipaw_data_format = INT(version)
  CALL scan_end ( iunps, "GIPAW_FORMAT_VERSION" )
  
  IF ( upf%gipaw_data_format == 1 .or. upf%gipaw_data_format == 0 ) THEN
     CALL read_pseudo_gipaw_core_orbitals ( iunps, upf, ierr )
     CALL read_pseudo_gipaw_local ( iunps, upf, ierr )
     CALL read_pseudo_gipaw_orbitals ( iunps, upf, ierr )
  ELSE
     write(stdout,'("read_pseudo_gipaw: UPF/GIPAW in unknown format")')
     return
  END IF
  ierr = 0
  RETURN
100 write(stdout,'("read_pseudo_gipaw: error  reading pseudo file")')
  !
END SUBROUTINE read_pseudo_gipaw

!---------------------------------------------------------------------
SUBROUTINE read_pseudo_gipaw_core_orbitals ( iunps, upf, ierr )
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  TYPE ( pseudo_upf ), INTENT ( INOUT ) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  CHARACTER ( LEN = 75 ) :: dummy1, dummy2
  INTEGER :: nb, ir

  ierr = 1
  CALL scan_begin ( iunps, "GIPAW_CORE_ORBITALS", .false. )
  READ ( iunps, *, err=100, end=100 ) upf%gipaw_ncore_orbitals
  
  ALLOCATE ( upf%gipaw_core_orbital_n(upf%gipaw_ncore_orbitals) )
  ALLOCATE ( upf%gipaw_core_orbital_l(upf%gipaw_ncore_orbitals) )
  ALLOCATE ( upf%gipaw_core_orbital_el(upf%gipaw_ncore_orbitals) )
  ALLOCATE ( upf%gipaw_core_orbital(upf%mesh,upf%gipaw_ncore_orbitals) )
  upf%gipaw_core_orbital = 0.0_dp
  
  DO nb = 1, upf%gipaw_ncore_orbitals
     CALL scan_begin ( iunps, "GIPAW_CORE_ORBITAL", .false. )
     READ (iunps, *, err=100, end=100) &
          upf%gipaw_core_orbital_n(nb), upf%gipaw_core_orbital_l(nb), &
          dummy1, dummy2, upf%gipaw_core_orbital_el(nb)
     READ ( iunps, *, err=100, end=100 ) &
          ( upf%gipaw_core_orbital(ir,nb), ir = 1, upf%mesh )
     CALL scan_end ( iunps, "GIPAW_CORE_ORBITAL" )
  END DO
  
  CALL scan_end ( iunps, "GIPAW_CORE_ORBITALS" )
  ierr = 0
  RETURN
  !
100 write(stdout,'("read_pseudo_gipaw_core_orbital: error reading file")')
  !
END SUBROUTINE read_pseudo_gipaw_core_orbitals

!---------------------------------------------------------------------
SUBROUTINE read_pseudo_gipaw_local ( iunps, upf, ierr )
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  TYPE ( pseudo_upf ), INTENT ( INOUT ) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  INTEGER :: ir
  
  ierr = 1
  CALL scan_begin ( iunps, "GIPAW_LOCAL_DATA", .false. )
  
  ALLOCATE ( upf%gipaw_vlocal_ae(upf%mesh) )
  ALLOCATE ( upf%gipaw_vlocal_ps(upf%mesh) )
  
  CALL scan_begin ( iunps, "GIPAW_VLOCAL_AE", .false. )
  
  READ ( iunps, *, err=100, end=100 ) &
       ( upf%gipaw_vlocal_ae(ir), ir = 1, upf%mesh )
  
  CALL scan_end ( iunps, "GIPAW_VLOCAL_AE" )
  
  CALL scan_begin ( iunps, "GIPAW_VLOCAL_PS", .false. )
  
  READ ( iunps, *, err=100, end=100 ) &
       ( upf%gipaw_vlocal_ps(ir), ir = 1, upf%mesh )
  
  CALL scan_end ( iunps, "GIPAW_VLOCAL_PS" )
  
  CALL scan_end ( iunps, "GIPAW_LOCAL_DATA" )
  ierr = 0
  RETURN
  !
100 write(stdout,'("read_pseudo_gipaw_local: error reading pseudo file")')
  !
END SUBROUTINE read_pseudo_gipaw_local

!---------------------------------------------------------------------
SUBROUTINE read_pseudo_gipaw_orbitals ( iunps, upf, ierr )
  !---------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  TYPE ( pseudo_upf ), INTENT ( INOUT ) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  !
  CHARACTER ( LEN = 75 ) :: dummy
  INTEGER :: nb, ir
  
  ierr = 1
  CALL scan_begin ( iunps, "GIPAW_ORBITALS", .false. )
  READ ( iunps, *, err=100, end=100 ) upf%gipaw_wfs_nchannels
  
  ALLOCATE ( upf%gipaw_wfs_el(upf%gipaw_wfs_nchannels) )
  ALLOCATE ( upf%gipaw_wfs_ll(upf%gipaw_wfs_nchannels) )
  ALLOCATE ( upf%gipaw_wfs_rcut(upf%gipaw_wfs_nchannels) )
  ALLOCATE ( upf%gipaw_wfs_rcutus(upf%gipaw_wfs_nchannels) )
  ALLOCATE ( upf%gipaw_wfs_ae(upf%mesh,upf%gipaw_wfs_nchannels) )
  ALLOCATE ( upf%gipaw_wfs_ps(upf%mesh,upf%gipaw_wfs_nchannels) )
  
  inquire ( unit = iunps, name = dummy )
  DO nb = 1, upf%gipaw_wfs_nchannels
     CALL scan_begin ( iunps, "GIPAW_AE_ORBITAL", .false. )
     READ (iunps, *, err=100, end=100) &
          upf%gipaw_wfs_el(nb), upf%gipaw_wfs_ll(nb)
     READ ( iunps, *, err=100, end=100 ) &
          ( upf%gipaw_wfs_ae(ir,nb), ir = 1, upf%mesh )
     CALL scan_end ( iunps, "GIPAW_AE_ORBITAL" )
     
     CALL scan_begin ( iunps, "GIPAW_PS_ORBITAL", .false. )
     READ (iunps, *, err=100, end=100) &
          upf%gipaw_wfs_rcut(nb), upf%gipaw_wfs_rcutus(nb)
     READ ( iunps, *, err=100, end=100 ) &
          ( upf%gipaw_wfs_ps(ir,nb), ir = 1, upf%mesh )
     CALL scan_end ( iunps, "GIPAW_PS_ORBITAL" )
  END DO
  CALL scan_end ( iunps, "GIPAW_ORBITALS" )
  ierr = 0
  RETURN
  !
100 write(stdout,'("read_pseudo_gipaw_orbitals: error reading pseudo file")')
  !
END SUBROUTINE read_pseudo_gipaw_orbitals
!</apsi>

subroutine read_pseudo_ppinfo (iunps, upf, ierr)
  !---------------------------------------------------------------------
  !
  implicit none
  !
  TYPE (pseudo_upf), INTENT(INOUT) :: upf
  integer, intent(in) :: iunps  
  integer, intent(out):: ierr
  character (len=80) :: dummy  
  real(dp) :: rdummy
  integer :: idummy, nb, ios

  ios=0
  DO while (ios==0) 
     READ (iunps, '(a)', err = 100, end = 100, iostat=ios) dummy  
     IF (matches ("Rcut", dummy) ) THEN
        DO nb=1,upf%nbeta
           READ (iunps, '(a2,2i3,f6.2,3f19.11)',err=100, end=100,iostat=ios) &
               upf%els_beta(nb), idummy, &
               idummy, rdummy, upf%rcut(nb), upf%rcutus (nb), rdummy
        ENDDO
        ios=100
     ELSE 
        nb = 1 
     ENDIF
  ENDDO
100 CONTINUE 
  ! PP_INFO reports values only for bound valence states,  
  IF ( nb .LE.  upf%nbeta ) THEN 
     IF (upf%els_beta(nb) == "</") THEN 
        upf%els_beta(nb:upf%nbeta)="__"
        upf%rcut(nb:upf%nbeta)=-1.0
        upf%rcutus(nb:upf%nbeta) = -1.0 
     END IF 
  END IF
  ierr = 0
  RETURN

  END SUBROUTINE read_pseudo_ppinfo

  SUBROUTINE set_coulomb_nonlocal(upf)
  IMPLICIT NONE
  TYPE(pseudo_upf) :: upf

  upf%nqf = 0
  upf%nqlc= 0
  upf%qqq_eps= -1._dp
  upf%kkbeta = 0
  ALLOCATE( upf%kbeta(1),         &
            upf%lll(1),           &
            upf%beta(upf%mesh,1), &
            upf%dion(1,1),        &
            upf%rinner(1),        &
            upf%qqq(1,1),         &
            upf%qfunc(upf%mesh,1),&
            upf%qfcoef(1,1,1,1),  &
            upf%rcut(1),          &
            upf%rcutus(1),        &
            upf%els_beta(1) )
   RETURN
   END SUBROUTINE set_coulomb_nonlocal
!
!=----------------------------------------------------------------------------=!
      END MODULE read_upf_v1_module
!=----------------------------------------------------------------------------=!
