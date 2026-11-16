!
! Copyright (C) 2016-2026 Quantum ESPRESSO Foundation 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_base
  !----------------------------------------------------------------------------
  !
  !! Subroutines used to read and write binary data produced by QE.  
  !! Author: Paolo Giannozzi, based on previous work by Carlo Cavazzoni
  !
  USE kinds,     ONLY : dp
  !
  IMPLICIT NONE
  !
  INTEGER, ALLOCATABLE::  index_file(:,:,:) 
  !! contains indices of G-components as a function of the Miller indices
  !! Note: defined here to avoid trouble with negative indices
  INTEGER, ALLOCATABLE::  index_l2f(:) 
  !! contains indices of local G-components in the file
  PRIVATE
  PUBLIC :: write_wfc, read_wfc, write_rhog, read_rhog
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc( iuni, filename, root_in_group, intra_group_comm, &
         ik, xk, ispin, nspin, wfc, ngw, gamma_only, nbnd, igl, ngwl, &
         b1,b2,b3, mill_k, scalef )
      !------------------------------------------------------------------------
      !
      !! Collects wfc, distributed on "intra_group_comm", writes them
      !! together with related information to file "filename.*"
      !! (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Only processor "root_in_group" collects data and writes to file
      !!
      USE mp_wave,    ONLY : mergewf, mergekg
      USE mp,         ONLY : mp_size, mp_rank, mp_max
      !
#if defined(__HDF5)
      USE qeh5_base_module,  ONLY  : qeh5_file, qeh5_dataset, qeh5_openfile, qeh5_open_dataset, &
                             qeh5_add_attribute, qeh5_write_dataset, qeh5_close, qeh5_set_space, &
                             qeh5_set_file_hyperslab              
#endif

      INTEGER,            INTENT(IN) :: iuni
      CHARACTER(LEN=*),   INTENT(IN) :: filename
      INTEGER,            INTENT(IN) :: ik, ispin, nspin
      REAL(DP),           INTENT(IN) :: xk(:)
      COMPLEX(DP),        INTENT(IN) :: wfc(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      INTEGER,            INTENT(IN) :: mill_k(:,:)
      REAL(DP),           INTENT(IN) :: b1(3), b2(3), b3(3)    
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      INTEGER,            INTENT(IN) :: root_in_group, intra_group_comm
      !
      LOGICAL                  :: ionode_in_group
      INTEGER                  :: igwx, npwx, npol, j
      INTEGER                  :: me_in_group, nproc_in_group, my_group
      INTEGER, ALLOCATABLE     :: itmp(:,:)
      COMPLEX(DP), ALLOCATABLE, TARGET :: wtmp(:)
      COMPLEX(DP), POINTER             :: wtmp2(:)
      !
#if defined(__HDF5) 
      TYPE (qeh5_file)         :: h5file
      TYPE (qeh5_dataset)      :: evc_dset, igw_dset
      !  
#endif 
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      npol = 1
      IF ( nspin == 4 ) npol = 2
      npwx = SIZE( wfc, 1 ) / npol
      !
      IF ( ionode_in_group ) THEN
#if defined  __HDF5
         CALL qeh5_openfile(h5file, TRIM(filename)//'.hdf5',action = 'write') 
         CALL qeh5_add_attribute( h5file%id, "ik", ik )
         CALL qeh5_add_attribute( h5file%id, "xk", xk, 1, [3]) 
         CALL qeh5_add_attribute( h5file%id, "ispin", ispin ) 
         IF (gamma_only) THEN 
            CALL qeh5_add_attribute(h5file%id, "gamma_only", ".TRUE.") 
         ELSE 
            CALL qeh5_add_attribute( h5file%id, "gamma_only", ".FALSE." )
         END IF 
         CALL qeh5_add_attribute( h5file%id, "scale_factor", scalef ) 
         CALL qeh5_add_attribute( h5file%id, "ngw", ngw )
         CALL qeh5_add_attribute( h5file%id, "igwx", igwx )
         CALL qeh5_add_attribute( h5file%id, "npol", npol )
         CALL qeh5_add_attribute( h5file%id, "nbnd", nbnd )
#else
         OPEN ( UNIT = iuni, FILE = TRIM(filename)//'.dat', &
              FORM='unformatted', STATUS = 'unknown' )
         WRITE(iuni) ik, xk, ispin, gamma_only, scalef
         WRITE(iuni) ngw, igwx, npol, nbnd
#endif
         !
      END IF
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( itmp( 3, MAX (igwx,1) ) )
      ELSE
         ! not used: some compiler do not like passing unallocated arrays
         ALLOCATE( itmp( 3, 1 ) )
      END IF
      itmp (:,:) = 0
      CALL mergekg( mill_k, itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         igw_dset%name = "MillerIndices"
         CALL qeh5_set_space( igw_dset, itmp(1,1), RANK = 2, DIMENSIONS = [3,igwx])
         CALL qeh5_open_dataset (h5file, igw_dset, ACTION = 'write') 
         CALL qeh5_add_attribute( igw_dset%id, "bg1", b1, RANK =1, DIMS = [3]) 
         CALL qeh5_add_attribute( igw_dset%id, "bg2", b2, RANK =1, DIMS = [3])
         CALL qeh5_add_attribute( igw_dset%id, "bg3", b3, RANK =1, DIMS = [3]) 
         CALL qeh5_add_attribute( igw_dset%id, "doc","Miller Indices of the wave-vectors, &
                                   same ordering as wave-function components") 
         CALL qeh5_write_dataset(itmp, igw_dset)
         CALL qeh5_close(igw_dset)   
       
#else
         WRITE(iuni) b1, b2, b3
         WRITE(iuni) itmp(1:3,1:igwx)
#endif
      END IF
      DEALLOCATE( itmp )
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( wtmp( MAX( npol*igwx, 1 ) ) )
         IF ( npol == 2 ) wtmp2 => wtmp( igwx+1:2*igwx )
      ELSE
         ALLOCATE( wtmp( 1 ) )
         IF ( npol == 2 ) wtmp2 => wtmp( 1:1 )
      ENDIF
      wtmp = 0.0_DP
      !
#if defined(__HDF5)
      IF ( ionode_in_group) THEN 
         CALL qeh5_set_space ( evc_dset, wtmp(1), 2, [npol*igwx, nbnd], MODE = 'f')
         CALL qeh5_set_space ( evc_dset, wtmp(1), 1, [npol*igwx], MODE = 'm')
         CALL qeh5_open_dataset (h5file, evc_dset, ACTION = 'write', NAME = 'evc' )       
         CALL qeh5_add_attribute( evc_dset%id, "doc:","Wave Functions, (npwx,nbnd), &
                               each contiguous line represents a wave function,  &
                               each complex coefficient is given by a couple of contiguous floats") 
      END IF
#endif 
      DO j = 1, nbnd
         !
         IF ( npol == 2 ) THEN
            !
            ! Quick-and-dirty noncolinear case - mergewf should be modified
            ! Collect into wtmp(1:igwx) the first set of plane waves components
            !
            CALL mergewf( wfc(1:npwx,       j), wtmp , ngwl, igl,&
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            !
            ! Collect into wtmp(igwx+1:2*igwx) the second set of plane waves
            ! components - pointer wtmp2 is used instead of wtmp(igwx+1:2*igwx)
            ! in order to avoid a bogus out-of-bound error
            !
            CALL mergewf( wfc(npwx+1:2*npwx,j), wtmp2, ngwl, igl,&
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            !
         ELSE
            !
            CALL mergewf( wfc(:,j), wtmp, ngwl, igl, me_in_group, &
                 nproc_in_group, root_in_group, intra_group_comm )
            !
         END IF
         !
         IF ( ionode_in_group ) THEN
#if defined(__HDF5)
            CALL qeh5_set_file_hyperslab ( evc_dset,  OFFSET = [0,j-1], COUNT = [2*npol*igwx,1] ) 
            CALL qeh5_write_dataset ( wtmp, evc_dset)   
#else
            WRITE(iuni) wtmp(1:npol*igwx)
#endif
         END IF
         !
      END DO
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         CALL qeh5_close ( evc_dset) 
         CALL qeh5_close (h5file)   
#else 
         CLOSE (UNIT = iuni, STATUS = 'keep' )
#endif
      END IF
      !
      IF ( npol == 2 ) NULLIFY ( wtmp2 )
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE write_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, filename, root_in_group, intra_group_comm,  &
         ik, xk, ispin, npol, wfc, ngw, gamma_only, nbnd, igl, ngwl, &
         b1, b2, b3, mill_k, scalef, ierr )
      !
      !! Processor "root_in_group" reads wfc and related information from file 
      !! "filename.*" (* = dat if fortran binary, * = hdf5 if HDF5),
      !! distributes wfc on "intra_group_comm"
      !! if ierr is present, return 0 if everything is ok, /= 0 if not
      !------------------------------------------------------------------------
      !
      USE mp,          ONLY : mp_bcast, mp_size, mp_rank, mp_max
      !
#if defined (__HDF5)
      USE  qeh5_base_module
#endif

      INTEGER,            INTENT(IN)    :: iuni
      CHARACTER(LEN=*),   INTENT(IN)    :: filename
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm
      INTEGER,            INTENT(IN)    :: ik
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, npol
      COMPLEX(DP),        INTENT(OUT)   :: wfc(:,:)
      INTEGER,            INTENT(IN)    :: igl(:)
      REAL(DP),           INTENT(OUT)   :: scalef
      REAL(DP),           INTENT(OUT)   :: xk(3)
      REAL(DP),           INTENT(OUT)   :: b1(3), b2(3), b3(3)
      INTEGER,            INTENT(INOUT) :: mill_k(:,:)
      LOGICAL,            INTENT(OUT)   :: gamma_only
      INTEGER, OPTIONAL,  INTENT(OUT)   :: ierr
      !
      INTEGER                           :: j
      INTEGER, ALLOCATABLE              :: mill_g(:,:)
      COMPLEX(DP), ALLOCATABLE          :: wtmp_g(:), wtmp(:)
      INTEGER                           :: ierr_
      INTEGER                           :: igwx_, npwx, ik_, nbnd_
      INTEGER                           :: me_in_group, nproc_in_group
      LOGICAL                           :: ionode_in_group
#if defined(__HDF5)
      TYPE (qeh5_file)    ::   h5file
      TYPE (qeh5_dataset) ::   h5dset_wfc, h5dset_mill
      CHARACTER(LEN=8)    ::   char_buf 
#endif  
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      !
      IF ( ionode_in_group ) THEN
#if !defined __HDF5
         OPEN ( UNIT = iuni, FILE=TRIM(filename)//'.dat', &
                FORM='unformatted', STATUS = 'old', IOSTAT = ierr_)
#else
         CALL qeh5_openfile( h5file, TRIM(filename)//'.hdf5', ACTION = 'read', ERROR = ierr_)
#endif
      END IF
      CALL mp_bcast( ierr_, root_in_group, intra_group_comm )
      IF ( PRESENT(ierr) ) THEN
         ierr = ierr_
         IF ( ierr /= 0 ) RETURN
      ELSE
         CALL errore( 'read_wfc ', &
              'cannot open restart file ' // TRIM(filename) //' for reading', ierr_ )
      END IF
      !
      IF ( ionode_in_group ) THEN
#if defined  __HDF5
         CALL qeh5_read_attribute (h5file%id, "ik", ik_)
         CALL qeh5_read_attribute (h5file%id, "xk",xk, RANK =1, DIMENSIONS = [3])
         CALL qeh5_read_attribute (h5file%id, "ispin", ispin)
         CALL qeh5_read_attribute (h5file%id, "gamma_only", char_buf, MAXLEN = len(char_buf) )
         IF (TRIM(char_buf) =='.TRUE.' .OR. TRIM(char_buf)=='.true.') THEN 
            gamma_only = .TRUE. 
         ELSE 
            gamma_only = .FALSE.
         END IF
         CALL qeh5_read_attribute (h5file%id, "scale_factor",scalef)
         CALL qeh5_read_attribute (h5file%id, "ngw", ngw)
         CALL qeh5_read_attribute (h5file%id, "nbnd", nbnd_)
         CALL qeh5_read_attribute (h5file%id, "npol",npol)
         CALL qeh5_read_attribute (h5file%id, "igwx",igwx_)
#else
         READ (iuni) ik_, xk, ispin, gamma_only, scalef
         READ (iuni) ngw, igwx_, npol, nbnd_
#endif
      END IF
      !
      CALL mp_bcast( ik_,    root_in_group, intra_group_comm )
      CALL mp_bcast( xk,     root_in_group, intra_group_comm )
      CALL mp_bcast( ispin,  root_in_group, intra_group_comm )
      CALL mp_bcast( gamma_only, root_in_group, intra_group_comm )
      CALL mp_bcast( scalef, root_in_group, intra_group_comm )
      CALL mp_bcast( ngw,    root_in_group, intra_group_comm )
      CALL mp_bcast( igwx_,  root_in_group, intra_group_comm )
      CALL mp_bcast( npol,   root_in_group, intra_group_comm )
      CALL mp_bcast( nbnd_,   root_in_group, intra_group_comm )
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( mill_g(3,igwx_) )
#if defined(__HDF5)
       CALL qeh5_open_dataset(h5file, h5dset_mill, ACTION = 'read', NAME = 'MillerIndices')
       IF ( h5dset_mill%filespace%dims(2) > igwx_ ) &
          CALL errore ( 'read_wfc', 'real dimensions of Miller indices dataset do not match igwx attribute', 8) 
       ! no reading of b1, b2, and b3 from file. They should be already set. 
       CALL qeh5_read_dataset ( mill_g, h5dset_mill) 
       CALL qeh5_close ( h5dset_mill) 
#else
         READ (iuni) b1, b2, b3
         READ (iuni) mill_g(1:3,1:igwx_)
#endif
      END IF
      !
      !! fill array "index_l2f" containing the correspondence between
      !! the local and on-file indices of G-components 
      !! requires Miller indices for both file and local G-components
      !
      ALLOCATE ( index_l2f(ngwl) )
      CALL fill_index_l2f ( root_in_group, intra_group_comm, &
           igwx_, mill_g, gamma_only, ngwl, mill_k, index_l2f )
      IF (ALLOCATED(mill_g))  DEALLOCATE( mill_g )
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( wtmp_g( npol*igwx_ ) )
#if defined (__HDF5) 
         CALL qeh5_open_dataset( h5file, h5dset_wfc, ACTION = 'read', NAME = 'evc')
         CALL qeh5_set_space ( h5dset_wfc, wtmp_g(1), RANK = 1, DIMENSIONS = [npol*igwx_], MODE = 'm')
#endif
      ELSE
         ALLOCATE( wtmp_g(1) )
      ENDIF
      !
      ALLOCATE( wtmp(ngwl*npol) )
      nbnd = nbnd_ 
      npwx = SIZE( wfc, 1 ) / npol
      !
      DO j = 1, nbnd_ 
         !
         IF ( j <= SIZE( wfc, 2 ) ) THEN
            !
            IF ( ionode_in_group ) THEN 
#if defined __HDF5
               
               CALL qeh5_set_file_hyperslab (h5dset_wfc, OFFSET = [0,j-1], COUNT = [2*npol*igwx_,1] )
               CALL qeh5_read_dataset (wtmp_g, h5dset_wfc )
#else
               READ (iuni) wtmp_g(1:npol*igwx_) 
#endif
            END IF
            !
            !! first dimension of wfc (npol*npwx) may exceed the number of
            !! components (ngwl): ensure that no garbage is left for n>ngwl
            !
            wfc(:,j) = (0.0_dp, 0.0_dp)
            !
            !! Distribute each wave-function read from file across processors,
            !! ensuring consistency with the ordering of the current calculation
            !! also convert from half to full G-vector sphere, or vice versa
            !
            CALL distribute_vect ( root_in_group, intra_group_comm, &
                 igwx_, wtmp_g, gamma_only, ngwl, index_l2f, npol, wtmp ) 
            !
            !! copy data to wfc
            wfc(1:ngwl,j) = wtmp(1:ngwl)
            !! noncolinear case: spin-up and spin-down PWs are not stored 
            !! contiguously (there can be a gap in the middle and at the end)
            IF ( npol == 2 ) THEN
               wfc(npwx+1:npwx+ngwl,j) = wtmp(ngwl+1:npol*ngwl)
            END IF
            !
         END IF
         !
      END DO
      !
      IF ( ionode_in_group ) THEN
#if defined (__HDF5)
         CALL qeh5_close(h5dset_wfc) 
         CALL qeh5_close(h5file)
#else
         CLOSE ( UNIT = iuni, STATUS = 'keep' )
#endif
      END IF
      !
      DEALLOCATE ( index_l2f )
      DEALLOCATE( wtmp )
      DEALLOCATE( wtmp_g )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rhog ( filename, root_in_group, intra_group_comm, &
         b1, b2, b3, gamma_only, mill, ig_l2g, rho )
      !------------------------------------------------------------------------
      !! Collects rho(G), distributed on "intra_group_comm", writes it
      !! together with related information to file "filename".*
      !! (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Processor "root_in_group" collects data and writes to file
      !
      USE mp,                   ONLY : mp_sum, mp_bcast, mp_size, mp_rank
      USE mp_wave,              ONLY : mergewf, mergekg
#if defined (__HDF5)
      USE qeh5_base_module
#endif
      !
      CHARACTER(LEN=*), INTENT(IN) :: filename
      !! name of file written (to which a suffix is added)
      INTEGER,            INTENT(IN) :: root_in_group
      !! root processor that collects and writes
      INTEGER,            INTENT(IN) :: intra_group_comm
      !! rho(G) is distributed over this group of processors
      REAL(dp),         INTENT(IN) :: b1(3), b2(3), b3(3)
      !!  b1, b2, b3 are the three primitive vectors in a.u.
      INTEGER,          INTENT(IN) :: mill(:,:)
      !! Miller indices for local G-vectors
      !! G = mill(1)*b1 + mill(2)*b2 + mill(3)*b3
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! local-to-global indices, for machine- and mpi-independent ordering
      !! on this processor, G(ig) maps to G(ig_l2g(ig)) in global ordering
      LOGICAL,          INTENT(IN) :: gamma_only
      !! if true, only the upper half of G-vectors (z >=0) is present
      COMPLEX(dp),      INTENT(IN) :: rho(:,:)
      !! rho(G) on this processor
      !
      COMPLEX(dp), ALLOCATABLE :: rhoaux(:)
      !! Local rho(G), with LSDA workaround
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      !! Global rho(G) collected on root proc
      INTEGER, ALLOCATABLE     :: mill_g(:,:)
      !! Global Miller indices collected on root proc
      INTEGER                  :: me_in_group, nproc_in_group
      LOGICAL                  :: ionode_in_group
      INTEGER                  :: ngm, nspin, ngm_g, igwx
      INTEGER                  :: iun, ns, ig, ierr
      !
#if defined __HDF5
      TYPE (qeh5_file)          ::  h5file
      TYPE (qeh5_dataset)       ::  h5dset_mill, h5dset_rho_g
      CHARACTER(LEN=10)          :: bool_char = ".FALSE.", datasets(4)        
      !
#endif
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      ngm  = SIZE (rho, 1)
      IF (ngm /= SIZE (mill, 2) .OR. ngm /= SIZE (ig_l2g, 1) ) &
         CALL errore('write_rhog', 'inconsistent input dimensions', 1)
      nspin= SIZE (rho, 2)
#if defined(__HDF5) 
      IF ( nspin <=2) THEN 
         datasets(1:2) = ["rhotot_g  ", "rhodiff_g "]
      ELSE
         datasets(1) =      "rhotot_g"
         datasets(2) =      "m_x" 
         datasets(3) =      "m_y" 
         datasets(4) =      "m_z"
      END IF  
#endif
      iun  = 4
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      CALL mp_sum( ngm_g, intra_group_comm )
      !
      ierr = 0
#if defined (__HDF5)
      IF ( ionode_in_group ) CALL qeh5_openfile(h5file, FILE = &
           TRIM(filename)//'.hdf5', ACTION = 'write', ERROR = ierr) 
#else
      IF ( ionode_in_group ) OPEN ( UNIT = iun, FILE = TRIM(filename)//'.dat', &
                FORM = 'unformatted', STATUS = 'unknown', iostat = ierr )
#endif
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error opening file ' &
           & // TRIM( filename ), 1 )
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
          IF ( gamma_only) bool_char = '.TRUE.'
          CALL qeh5_add_attribute (h5file%id, NAME = "gamma_only", TEXT = TRIM(bool_char) )
          CALL qeh5_add_attribute (h5file%id, "ngm_g",  ngm_g ) 
          CALL qeh5_add_attribute (h5file%id, "nspin",  nspin )       
#else
          WRITE (iun, iostat=ierr) gamma_only, ngm_g, nspin
          WRITE (iun, iostat=ierr) b1, b2, b3
#endif
      END IF
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
           & // TRIM( filename ), 1 )
      !
      ! ... collect all G-vectors across processors within the band group
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( mill_g( 3, ngm_g ) )
      ELSE
         ! not used: some compiler do not like passing unallocated arrays
         ALLOCATE( mill_g( 3, 1 ) )
      END IF
      !
      ! ... mergekg collects distributed array mill(1:3,ig) where ig is the
      ! ... local index, into array mill_g(1:3,ig_g), where ig_g=ig_l2g(ig)
      ! ... is the global index. mill_g is collected on root_bgrp only
      !
      CALL mergekg( mill, mill_g, ngm, ig_l2g, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      !
      ! ... write G-vectors
      !
      IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         CALL qeh5_set_space ( h5dset_mill, mill_g(1,1), RANK = 2, DIMENSIONS = [3,ngm_g] ) 
         CALL qeh5_open_dataset ( h5file, h5dset_mill, NAME = "MillerIndices" , ACTION = 'write')
         !
         CALL qeh5_add_attribute(h5dset_mill%id, NAME = 'bg1', VALUE = b1(1), RANK = 1, DIMS = [3]) 
         CALL qeh5_add_attribute(h5dset_mill%id, NAME = 'bg2', VALUE = b2(1), RANK = 1, DIMS = [3])
         CALL qeh5_add_attribute(h5dset_mill%id, NAME = 'bg3', VALUE = b3(1), RANK = 1, DIMS = [3]) 
         !
         CALL qeh5_write_dataset( mill_g, h5dset_mill )
         ! 
         CALL qeh5_close( h5dset_mill)    
#else
         WRITE (iun, iostat=ierr) mill_g(1:3,1:ngm_g)
#endif
      END IF
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... deallocate to save memory
      !
      DEALLOCATE( mill_g )
      !
      ! ... now collect all G-vector components of the charge density
      ! ... (one spin at the time to save memory) using the same logic
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( rho_g( ngm_g ) )
      ELSE
         ALLOCATE( rho_g( 1 ) )
      END IF
      ALLOCATE (rhoaux(ngm))
      !
      DO ns = 1, nspin
         !
         DO ig = 1, ngm
               rhoaux(ig) = rho(ig,ns)
         END DO
         !
         rho_g = 0
         CALL mergewf( rhoaux, rho_g, ngm, ig_l2g, me_in_group, &
              nproc_in_group, root_in_group, intra_group_comm )
         !
         IF ( ionode_in_group ) THEN
#if defined(__HDF5)
         CALL qeh5_set_space ( h5dset_rho_g, rho_g(1), RANK = 1 , DIMENSIONS = [ngm_g] ) 
         CALL qeh5_open_dataset( h5file, h5dset_rho_g, NAME = TRIM(datasets(ns)) , ACTION = 'write', ERROR = ierr )
         if (ierr /= 0 ) CALL infomsg('write_rho:', 'error while opening h5 dataset in charge_density.hdf5') 
         CALL qeh5_write_dataset(rho_g, h5dset_rho_g) 
         CALL qeh5_close( h5dset_rho_g)     
#else
            WRITE (iun, iostat=ierr) rho_g(1:ngm_g)
#endif
         END IF
         CALL mp_bcast( ierr, root_in_group, intra_group_comm )
         IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
              & // TRIM( filename ), 2+ns )
         !
      END DO
      !
#if defined(__HDF5) 
      IF (ionode_in_group) CALL qeh5_close(h5file) 
#else
      IF (ionode_in_group) CLOSE (UNIT = iun, status ='keep' )
#endif
      !
      DEALLOCATE( rhoaux )
      DEALLOCATE( rho_g )
      !
      RETURN
      !
    END SUBROUTINE write_rhog
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rhog ( filename, root_in_group, intra_group_comm, &
         ig_l2g, nspin, rho, gamma_only, ier_ )
      !------------------------------------------------------------------------
      !! Read and distribute rho(G) from file  "filename".* 
      !! (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Processor "root_in_group" reads from file, distributes to
      !! all processors in the intra_group_comm communicator 
      !
      USE mp,         ONLY : mp_size, mp_rank, mp_bcast
      USE gvect,      ONLY : ngm_g, ngm, mill
      !
#if defined (__HDF5) 
      USE qeh5_base_module
#endif
      CHARACTER(LEN=*), INTENT(IN) :: filename
      !! name of file read (to which a suffix is added)
      INTEGER,          INTENT(IN) :: root_in_group
      !! root processor that reads and distributes
      INTEGER,          INTENT(IN) :: intra_group_comm
      !! rho(G) is distributed over this group of processors
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! unused, to be removed from call
      INTEGER,          INTENT(IN) :: nspin
      !! read up to nspin components
      COMPLEX(dp), INTENT(INOUT)   :: rho(:,:)
      !! must be allocated on input, is filled with G components on output
      LOGICAL, INTENT(IN), OPTIONAL:: gamma_only
      !! true if half G-vectors are used - may differ from value read from file
      !! used only for informative message, not in actual reading/distributing
      INTEGER, OPTIONAL, INTENT(OUT):: ier_
      !
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      REAL(dp)                 :: b1(3), b2(3), b3(3)
      INTEGER                  :: nspin_, isup, isdw
      INTEGER                  :: iun, ns, ig, ierr, ierr_open
      INTEGER                  :: me_in_group, nproc_in_group
      LOGICAL                  :: ionode_in_group, gamma_only_
      INTEGER                  :: ngm_g_
      INTEGER, ALLOCATABLE     :: mill_g(:,:)
      !
#if defined __HDF5
      TYPE ( qeh5_file)       :: h5file
      TYPE ( qeh5_dataset)    :: h5dset_mill, h5dset_rho_g
      CHARACTER(LEN=10)       :: tempchar, datasets(4)
#endif
      !
      IF ( ngm /= SIZE (rho, 1) )  &
           CALL errore('read_rhog', 'inconsistent input dimensions', 1)
      !
      iun  = 4
      ierr = 0
      ierr_open = 0
      IF ( PRESENT(ier_) ) ier_ = 0
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      !
      IF ( ionode_in_group ) THEN
#if defined (__HDF5) 
         CALL qeh5_openfile(h5file, TRIM(filename)//'.hdf5', ACTION = 'read', &
              error = ierr)
         IF ( ierr /= 0 .AND. PRESENT(ier_) ) THEN
            ierr_open = ierr
            GO TO 10
         END IF
         CALL qeh5_read_attribute (h5file%id, "gamma_only", tempchar, MAXLEN = len(tempchar)  )
         CALL qeh5_read_attribute (h5file%id, "ngm_g", ngm_g_ ) 
         CALL qeh5_read_attribute (h5file%id, "nspin", nspin_)  
         SELECT CASE (TRIM(tempchar) )  
         CASE ('.true.', '.TRUE.' ) 
            gamma_only_ = .TRUE.
         CASE DEFAULT
            gamma_only_ = .FALSE.
         END SELECT
#else
         OPEN ( UNIT = iun, FILE = TRIM( filename ) // '.dat', &
              FORM = 'unformatted', STATUS = 'old', iostat = ierr )
         IF ( ierr /= 0 ) THEN
            IF ( PRESENT(ier_) ) THEN
               ierr_open = ierr
               GO TO 10
            END IF
            ierr = 1
            GO TO 10
         END IF
         READ (iun, iostat=ierr) gamma_only_, ngm_g_, nspin_
         IF ( ierr /= 0 ) THEN
            ierr = 2
            GO TO 10
         END IF
         READ (iun, iostat=ierr) b1, b2, b3
         IF ( ierr /= 0 ) ierr = 3
#endif
10       CONTINUE 
      END IF
      !
      ! open failure with ier_ present: all ranks must return, not just the reader
      CALL mp_bcast( ierr_open, root_in_group, intra_group_comm )
      IF ( ierr_open /= 0 ) THEN
         ier_ = ierr_open
         RETURN
      END IF
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), ierr )
      CALL mp_bcast( ngm_g_, root_in_group, intra_group_comm )
      CALL mp_bcast( nspin_, root_in_group, intra_group_comm )
      CALL mp_bcast( gamma_only_, root_in_group, intra_group_comm )
      !
      IF ( nspin > nspin_ ) &
         CALL infomsg('read_rhog', 'some spin components not found on file')
      !
      ! ... allocate and read Miller indices - useful when there is a mismatch
      ! ... between input gamma tricks and gamma tricks read from file, but
      ! ... also to avoid mismatch in ordering that may occur in rare cases
      !
      IF (ionode_in_group) THEN
         ALLOCATE (mill_g(3,ngm_g_))
#if defined (__HDF5)
         CALL qeh5_open_dataset( h5file, h5dset_mill, &
              NAME = "MillerIndices", ACTION = 'read', ERROR = ierr)
         CALL qeh5_read_dataset ( mill_g , h5dset_mill )
         CALL qeh5_close ( h5dset_mill )
#else
         READ (iun, iostat=ierr) mill_g(1:3,1:ngm_g_)
#endif
      ELSE
         ALLOCATE (mill_g(1,1))
      END IF
      !
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... now read, broadcast and re-order G-vector components
      ! ... of the charge density (one spin at the time to save memory)
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( rho_g(ngm_g_) )
      ELSE
         ALLOCATE( rho_g( 1 ) )
      END IF
#if defined(__HDF5)
      IF (nspin_ <= 2) THEN 
        datasets(1:2) =["rhotot_g  ", "rhodiff_g "]
      ELSE
        datasets(1)  = "rhotot_g"; datasets(2:4) = ["m_x", "m_y", "m_z"] 
      END IF
#endif 
      !
      !! fill array "index_l2f" containing the correspondence between
      !! the local and on-file indices of G-components 
      !! requires Miller indices for both file and local G-components
      !
      ALLOCATE ( index_l2f(ngm) )
      CALL fill_index_l2f ( root_in_group, intra_group_comm, &
           ngm_g_, mill_g, gamma_only_, ngm, mill, index_l2f )
      IF (ALLOCATED(mill_g))  DEALLOCATE( mill_g )
      !
      IF ( PRESENT(gamma_only) ) THEN
         IF ( gamma_only .and. .NOT. gamma_only_ ) THEN 
            call infomsg('read_rhog','Conversion: K charge to Gamma charge') 
         ELSE IF ( gamma_only_ .and. .NOT. gamma_only ) THEN 
            call infomsg ('read_rhog', 'Conversion: Gamma charge to K charge')
         !ELSE
         !call infomsg ('read_rhog', 'Reordering of the charge')
         ENDIF 
      ENDIF 
      !
      DO ns = 1, nspin_
         !
         IF ( ionode_in_group ) THEN
#if defined(__HDF5)
           CALL qeh5_open_dataset( h5file, h5dset_rho_g, NAME = TRIM(datasets(ns)), ACTION = 'read', ERROR = ierr) 
           CALL qeh5_read_dataset ( rho_g , h5dset_rho_g )
           CALL qeh5_close ( h5dset_rho_g )  
#else 
           READ (iun, iostat=ierr) rho_g(1:ngm_g_)
#endif
         END IF
         CALL mp_bcast( ierr, root_in_group, intra_group_comm )
         IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
              & // TRIM( filename ), 2+ns )
         !
         !! Distribute the charge density read from file across processors,
         !! ensuring consistency with the ordering of the current calculation
         !! also convert from half to full G-vector sphere, or vice versa
         !
         CALL distribute_vect ( root_in_group, intra_group_comm, &
           ngm_g_, rho_g, gamma_only_, ngm, index_l2f, 1, rho(:,ns) )
         !
      END DO
      !
#if defined(__HDF5)
      IF ( ionode_in_group ) CALL qeh5_close( h5file)
#else
      IF ( ionode_in_group ) CLOSE (UNIT = iun, status ='keep' )
#endif
      !
      DEALLOCATE ( index_l2f )
      DEALLOCATE( rho_g )
      !
      RETURN
      !
    END SUBROUTINE read_rhog
    !
    !-----------------------------------------------------------------
    SUBROUTINE distribute_vect ( root_in_group, intra_group_comm, &
         ngm_g_file, vect_g_file, gamma_only_file, ngm, index_l2f, npol, &
         vect )
    !-----------------------------------------------------------------
      !
      !! Distribute array "vect_g_file", containing G-vector components
      !! as read from file, from process root_in_group to all the MPI
      !! processes in intra_group_comm. Output in array "vect".
      !! index_l2f contains the correspondence between local and file
      !! ordering of G components: vect(n) = vect_g_file(index_l2f(n))
      !! npol consecutive sets of data, of input size ngm_g_file and
      !! output file ngm, are distributed and ordered
      !
      USE io_global,     ONLY : stdout
      USE mp,            ONLY : mp_size, mp_rank, mp_get
      !
      INTEGER, INTENT(in) :: root_in_group
      INTEGER, INTENT(in) :: intra_group_comm
      INTEGER, INTENT(in) :: ngm_g_file  
      !! total number of G components found in file 
      INTEGER, INTENT(in) :: ngm
      !! required number of g vectors (for this processor)
      INTEGER, INTENT(in) :: index_l2f(:)
      !! n-th local component G(n) in file at position m=index_l2f(n)
      INTEGER, INTENT(in) :: npol
      !! Order npol consecutive sets of data 
      COMPLEX(kind=DP), INTENT(in) :: vect_g_file(:)
      !! vect(G) read from file 
      COMPLEX(kind=DP), INTENT(out) :: vect(:)
      !! output vect(G) distributed according to indices index_l2f
      LOGICAL, INTENT(in) :: gamma_only_file 
      !
      !! local variables for G components distribution in parallel case 
      !
      COMPLEX(kind=DP), ALLOCATABLE :: vect_l(:)
      INTEGER, ALLOCATABLE::  index_l2f_l(:) 
      INTEGER             :: ngm_l
      INTEGER             :: me_in_group, nproc_in_group, npr
      INTEGER             :: ig, ig_file, ig_read, ig_miss, ierr
      LOGICAL             :: ionode_in_group
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      !
      ig_read = 0
      ig_miss = 0
      DO npr = 0, nproc_in_group-1
         IF ( npr == root_in_group ) THEN
            IF ( ionode_in_group ) THEN
               !! here fill vect on ionode_in_group
               CALL fill_vect (ngm, npol, index_l2f, &
                    ngm_g_file, vect_g_file, vect, ierr )
               ig_read = ig_read + ngm
               ig_miss = ig_miss + ierr
            END IF
         ELSE 
            !! here send ngm, index_l2f, from processor "npr"
            !!      receive ngm_l, index_l2f_l, on processor "root_in_group"
            CALL mp_get (ngm_l, ngm, me_in_group, root_in_group,&
                    npr, 1, intra_group_comm)
            IF ( ionode_in_group ) THEN
               ALLOCATE ( vect_l(npol*ngm_l), index_l2f_l(ngm_l) )
            ELSE
               ALLOCATE ( vect_l(1), index_l2f_l(1) )
            END IF
            CALL mp_get (index_l2f_l, index_l2f, me_in_group, root_in_group, npr, 2, intra_group_comm)
            !! here fill vect_l (on ionode_in_group only)
            IF ( ionode_in_group ) THEN
               CALL fill_vect (ngm_l, npol, index_l2f_l, &
                    ngm_g_file, vect_g_file, vect_l, ierr )
               ig_read = ig_read + ngm_l
               ig_miss = ig_miss + ierr
            END IF
            !! here send vect_l back to processor "npr" from "root_in_group"
            CALL mp_get( vect, vect_l, me_in_group, npr, root_in_group, 3, intra_group_comm)
            DEALLOCATE ( vect_l, index_l2f_l )
         END IF
      END DO
      !
      IF ( ig_miss > 0 ) THEN
         WRITE (stdout,'(5X,"read_rhog: ",I6," G-components read, ",I6," not found (set to 0)")') ig_read,ig_miss
         WRITE (stdout,'(5X,"Zero-padding will be used")')
      END IF
      !
    END SUBROUTINE distribute_vect
    !
    !
    SUBROUTINE fill_vect (ngm, npol, index_l2f, ngm_g_file, &
                    vect_g_file, vect, ig_miss )
      !
      !! Fill vect with G-vector components vect_g_file read from file
      !! Negative indices in index_l2f mean "take the conjugate (-G)"
      !
      INTEGER, INTENT(in) :: ngm
      !! number of output G-components
      INTEGER, INTENT(in) :: ngm_g_file
      !! number of G-component found in file 
      INTEGER, INTENT(in) :: npol
      !! number of dataset to be reordered
      INTEGER, INTENT(in) :: index_l2f(:)
      !! n-th local component G(n) in file at position m=index_l2f(n)
      COMPLEX(kind=DP), INTENT(in) :: vect_g_file(:)
      !! vect(G) read from file 
      COMPLEX(kind=DP), INTENT(out) :: vect(:)
      !! output vect(G) distributed according to indices index_l2f
      INTEGER, INTENT(out) :: ig_miss
      !! count missing components, if any
      INTEGER :: ig, ig_file, ig_max, ipol
      !
      ig_miss= 0
      ig_max = SIZE(vect_g_file,1)
      vect (:) = (0.0_dp, 0.0_dp)
      DO ig = 1, ngm
         ig_file = index_l2f(ig)
         IF ( ig_file > 0 .AND. ig_file <= ig_max ) THEN
            do ipol = 1, npol
               vect(ig+(ipol-1)*ngm) = vect_g_file (ig_file+(ipol-1)*ngm_g_file)
            end do
         ELSE IF ( ig_file < 0 .AND. -ig_file <= ig_max) THEN
            do ipol = 1, npol
               vect(ig+(ipol-1)*ngm) = &
                    CONJG(vect_g_file (-ig_file+(ipol-1)* ngm_g_file))
            end do
         ELSE
            ig_miss = ig_miss + 1
         END IF
      END DO
      !
    END SUBROUTINE fill_vect
    !
    !-----------------------------------------------------------------
    SUBROUTINE fill_index_l2f ( root_in_group, intra_group_comm, &
         ngm_g_file, mill_g_file, gamma_only_file, ngm, mill, index_l2f )
    !-----------------------------------------------------------------
      !
      !! Compute index_l2f giving the correspondence from data read from
      !! file and local date to be distributed across MPI processors
      !! Requires Miller indices for both the file and the local data
      !
      USE io_global,     ONLY : stdout
      USE mp,            ONLY : mp_size, mp_rank, mp_get
      !
      INTEGER, INTENT(in) :: root_in_group
      INTEGER, INTENT(in) :: intra_group_comm
      INTEGER, INTENT(in) :: ngm_g_file  
      !! total number of G components found in file 
      INTEGER, INTENT(in) :: mill_g_file(:,:)
      !! Miller indices of G components read from file
      INTEGER, INTENT(in) :: ngm
      !! required number of g vectors (for this processor)
      INTEGER, INTENT(in) :: mill(:,:)
      !! Miller indices of required G components (for this processor)
      INTEGER, INTENT(out) :: index_l2f(:)
      !! output local-to-file indices
      LOGICAL, INTENT(in) :: gamma_only_file 
      !
      !! local variables for G components distribution in parallel case 
      !
      INTEGER, ALLOCATABLE::  mill_l(:,:) 
      INTEGER, ALLOCATABLE::  index_l2f_l(:) 
      INTEGER             :: ngm_l
      INTEGER             :: me_in_group, nproc_in_group, npr
      INTEGER             :: ig, ig_file, ig_read, ig_miss, ierr
      LOGICAL             :: ionode_in_group
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      !
      !! next call allocates and fills the table "index_file",
      !! on ionode ONLY! it is a big table, for all G components
      IF(ionode_in_group) CALL fill_index_file (ngm_g_file, mill_g_file)
      !
      ig_read = 0
      ig_miss = 0
      DO npr = 0, nproc_in_group-1
         IF ( npr == root_in_group ) THEN
            IF ( ionode_in_group ) THEN

               !! here fill index_l2f on ionode_in_group
               CALL fill_index_l (ngm, mill, gamma_only_file, &
                    ngm_g_file, index_l2f, ierr )

               ig_read = ig_read + ngm
               ig_miss = ig_miss + ierr
            END IF
         ELSE 
            !! here send ngm, mill, from processor "npr"
            !!      receive ngm_l, mill_l, on processor "root_in_group"
            CALL mp_get (ngm_l, ngm, me_in_group, root_in_group,&
                    npr, 1, intra_group_comm)
            IF ( ionode_in_group ) THEN
               ALLOCATE ( index_l2f_l(ngm_l), mill_l(3,ngm_l) )
            ELSE
               ALLOCATE ( index_l2f_l(1), mill_l(3,1) )
            END IF
            CALL mp_get (mill_l, mill, me_in_group, root_in_group, npr, 2, intra_group_comm)
            !! here fill index_l2f_l (on ionode_in_group only)
            IF ( ionode_in_group ) THEN
               CALL fill_index_l (ngm_l, mill_l, gamma_only_file, &
                    ngm_g_file, index_l2f_l, ierr )
               ig_read = ig_read + ngm_l
               ig_miss = ig_miss + ierr
            END IF
            !! here send index_l2f_l back to processor "npr" from "root_in_group"
            CALL mp_get( index_l2f, index_l2f_l, me_in_group, npr, root_in_group, 3, intra_group_comm)
            DEALLOCATE ( index_l2f_l, mill_l )
         END IF
      END DO
      !
      IF ( ig_miss > 0 ) THEN
         WRITE (stdout,'(5X,"read_rhog: ",I6," G-components read, ",I6," not found (set to 0)")') ig_read,ig_miss
         WRITE (stdout,'(5X,"Zero-padding will be used")')
      END IF
      !
      IF(ionode_in_group) DEALLOCATE(index_file)
      !
    END SUBROUTINE fill_index_l2f
    !-----------------------------------------------------------------
    SUBROUTINE fill_index_file (ngm_g_file, mill_g_file)
    !-----------------------------------------------------------------
    !
    INTEGER, INTENT(IN) :: ngm_g_file, mill_g_file(:,:)
    INTEGER :: n1l,n2l,n3l, n1u,n2u,n3u, n(3), ig
    !
    !! find lower and upper bounds for file "index_file"
    n1l = MINVAL(mill_g_file(1,:))
    n2l = MINVAL(mill_g_file(2,:))
    n3l = MINVAL(mill_g_file(3,:))
    n1u = MAXVAL(mill_g_file(1,:))
    n2u = MAXVAL(mill_g_file(2,:))
    n3u = MAXVAL(mill_g_file(3,:))
    ALLOCATE(index_file(n1l:n1u,n2l:n2u,n3l:n3u)) 
    index_file = 0
    !
    DO ig = 1, ngm_g_file
       n(:) = mill_g_file(:,ig)
       index_file ( n(1), n(2), n(3) ) = ig
    END DO
    !
    END SUBROUTINE fill_index_file
    !
    !-----------------------------------------------------------------
    SUBROUTINE fill_index_l (ngm, mill, gamma_only_file, ngm_g_file, &
                    index_l2f, ig_miss )
      !-----------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: ngm
      !! number of output G-vectors 
      INTEGER, INTENT(in) :: ngm_g_file
      !! number of G-vectors found in file 
      INTEGER, INTENT(in) :: mill(:,:)
      !! Miller indices of G components read from file 
      LOGICAL, INTENT(in) :: gamma_only_file 
      !! true if half G-vector sphere read from file 
      INTEGER, INTENT(out) :: index_l2f(:)
      !! output local-to-file indices
      INTEGER, INTENT(out) :: ig_miss
      !! count missing components, if any
      INTEGER :: ig, ig_file, ipol, n(3), nlo(3), nup(3)
      !! bounds
      LOGICAL :: is_gm
      !! is_gm=.true. : G-component is not there, look for -G instead
      !
      nlo = LBOUND (index_file)
      nup = UBOUND (index_file)
      ig_miss = 0
      DO ig = 1, ngm
         n(:) = mill(:,ig)
         !! BEWARE: gamma_only must follow here the same logic as in "ggen"
         is_gm= ( (n(1) < 0) .OR. (n(1) == 0 .AND. n(2) < 0 ) .OR. &
                  (n(1) == 0 .AND. n(2) == 0 .AND. n(3) < 0 ) ) .AND. &
                  gamma_only_file 
         ig_file = 0 
         IF ( is_gm ) THEN
            IF ( ALL( -n >= nlo ) .AND. ALL ( -n <= nup ) ) &
                 ig_file =-index_file(-n(1),-n(2),-n(3))
         ELSE IF ( ALL( n >= nlo ) .AND. ALL( n <= nup ) ) THEN
            ig_file = index_file( n(1), n(2), n(3))
         END IF
         index_l2f(ig) = ig_file
         IF ( ig_file == 0 .AND. ig_file <= ngm_g_file) THEN
            ig_miss = ig_miss + 1
         END IF
      END DO
      !
    END SUBROUTINE fill_index_l
    
  END MODULE io_base

