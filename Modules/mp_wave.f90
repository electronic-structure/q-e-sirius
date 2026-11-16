!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
    MODULE mp_wave
      !
      !! MPI management of wave function related arrays.
      !
      IMPLICIT NONE
      SAVE

    CONTAINS

      SUBROUTINE mergewf ( pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )

      !! This subroutine merges the pieces of a wave functions (pw) splitted across 
      !! processors into a total wave function (pwt) containing al the components
      !! in a pre-defined order (the same as if only one processor is used).

      USE kinds
      USE parallel_include

      IMPLICIT NONE

      COMPLEX(DP), intent(in) :: PW(:)
      !! piece of wave function
      COMPLEX(DP), intent(out) :: PWT(:)
      !! total wave function
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor ( the one that should receive the data )
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)

      INTEGER :: ierr, i, ip, ngw_ip, ngw_lmax, itmp, igwx, gid

#if defined __MPI
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE( ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( igwx, itmp, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif

      IF ( mpime == root .AND. igwx > SIZE( pwt ) ) &
        CALL errore(' mergewf ',' wrong size for pwt ',SIZE(pwt) )

#if defined __MPI

      DO ip = 1, nproc

        IF( (ip-1) /= root ) THEN

! ...     In turn each processors send to root the wave components and their indexes in the 
! ...     global array
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid, IERR )
            CALL MPI_SEND( pw(1), ngwl, MPI_DOUBLE_COMPLEX, ROOT, IP+NPROC, gid, IERR )
          END IF
          IF ( mpime == root) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(pw_ip(ngw_lmax))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_RECV( pw_ip, ngw_lmax, MPI_DOUBLE_COMPLEX, (ip-1), IP+NPROC, gid, istatus, IERR )
            CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_COMPLEX, ngw_ip, ierr ) 
            DO I = 1, ngw_ip
              PWT(ig_ip(i)) = pw_ip(i)
            END DO
            DEALLOCATE(ig_ip)
            DEALLOCATE(pw_ip)
          END IF

        ELSE

          IF(mpime == root) THEN
            DO I = 1, ngwl
              PWT(ig_l2g(i)) = pw(i)
            END DO
          END IF

        END IF

        CALL MPI_BARRIER( gid, IERR )

      END DO

#elif ! defined __MPI

      DO I = 1, ngwl
        ! WRITE( stdout,*) 'MW ', ig_l2g(i), i
        PWT( ig_l2g(i) ) = pw(i)
      END DO

#else

      CALL errore(' MERGEWF ',' no communication protocol ',0)

#endif

      RETURN
      END SUBROUTINE mergewf

!=----------------------------------------------------------------------------=!
      
      SUBROUTINE mergekg ( mill, millt, ngwl, ig_l2g, mpime, nproc, root, comm )

      !! Same logic as for \(\texttt{mergewf}\), for Miller indices.

      USE kinds
      USE parallel_include

      IMPLICIT NONE

      INTEGER, intent(in) :: mill(:,:)
      !! Miller indices: distributed input
      INTEGER, intent(out):: millt(:,:)
      !! Miller indices: collected output
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      INTEGER, ALLOCATABLE :: mill_ip(:,:)

      INTEGER :: ierr, i, ip, ngw_ip, ngw_lmax, itmp, igwx, gid

#if defined __MPI
      INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE( ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE( igwx, itmp, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif
      IF ( mpime == root .AND. igwx > SIZE( millt, 2 ) ) &
        CALL errore(' mergekg',' wrong size for millt ',SIZE(millt,2) )

#if defined __MPI

      DO ip = 1, nproc

        IF( (ip-1) /= root ) THEN

! ...     In turn each processors send to root the wave components and their indexes in the 
! ...     global array
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid, IERR )
            CALL MPI_SEND( mill,3*ngwl, MPI_INTEGER, ROOT, IP+NPROC, gid, IERR )
          END IF
          IF ( mpime == root) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(mill_ip(3,ngw_lmax))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT( istatus, MPI_INTEGER, ngw_ip, ierr ) 
            CALL MPI_RECV( mill_ip,3*ngw_lmax, MPI_INTEGER, (ip-1), IP+NPROC, gid, istatus, IERR )
            DO I = 1,ngw_ip
              millt(:,ig_ip(i)) = mill_ip(:,i)
            END DO
            DEALLOCATE(ig_ip)
            DEALLOCATE(mill_ip)
          END IF

        ELSE

          IF(mpime == root) THEN
            DO I = 1, ngwl
              millt(:,ig_l2g(i)) = mill(:,i)
            END DO
          END IF

        END IF

        CALL MPI_BARRIER( gid, IERR )

      END DO

#elif ! defined __MPI

      DO I = 1, ngwl
        ! WRITE( stdout,*) 'MW ', ig_l2g(i), i
         millt(:,ig_l2g(i) ) = mill(:,i)
      END DO

#else

      CALL errore(' mergekg ',' no communication protocol ',0)

#endif

      RETURN
    END SUBROUTINE mergekg

!=----------------------------------------------------------------------------=!

      SUBROUTINE splitwf ( pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )

      !! This subroutine splits a total wave function (PWT) containing al the components
      !! in a pre-defined order (the same as if only one processor is used), across 
      !! processors (PW).

      USE kinds
      USE parallel_include
      IMPLICIT NONE

      COMPLEX(DP), INTENT(OUT) :: PW(:)
      !! piece of wave function
      COMPLEX(DP), INTENT(IN) :: PWT(:)
      !! total wave function
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl

      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)

      INTEGER ierr, i, ngw_ip, ip, ngw_lmax, gid, igwx, itmp, size_pwt

#if defined __MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

!
! ... Subroutine Body
!

      igwx = MAXVAL( ig_l2g(1:ngwl) )

#if defined __MPI

      gid = comm

! ... Get local and global wavefunction dimensions
      CALL MPI_ALLREDUCE(ngwl, ngw_lmax, 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      CALL MPI_ALLREDUCE(igwx, itmp    , 1, MPI_INTEGER, MPI_MAX, gid, IERR )
      igwx = itmp

#endif

      IF ( mpime == root .AND. igwx > SIZE(pwt )) &
        CALL errore (' splitwf ',' wrong size for pwt', SIZE(pwt) )

#if defined __MPI

      DO ip = 1, nproc
! ...   In turn each processor send to root the the indexes of its wavefunction conponents
! ...   Root receive the indexes and send the componens of the wavefunction read from the disk (pwt)
        IF ( (ip-1) /= root ) THEN
          IF ( mpime == (ip-1) ) THEN
            CALL MPI_SEND( ig_l2g, ngwl, MPI_INTEGER, ROOT, IP, gid,IERR)
            CALL MPI_RECV( pw(1), ngwl, MPI_DOUBLE_COMPLEX, ROOT, IP+NPROC, gid, istatus, IERR )
          END IF
          IF ( mpime == root ) THEN
            ALLOCATE(ig_ip(ngw_lmax))
            ALLOCATE(pw_ip(ngw_lmax))
            CALL MPI_RECV( ig_ip, ngw_lmax, MPI_INTEGER, (ip-1), IP, gid, istatus, IERR )
            CALL MPI_GET_COUNT(istatus, MPI_INTEGER, ngw_ip, ierr)
            DO i = 1, ngw_ip
                  pw_ip(i) = PWT(ig_ip(i)) 
            END DO
            CALL MPI_SEND( pw_ip, ngw_ip, MPI_DOUBLE_COMPLEX, (ip-1), IP+NPROC, gid, IERR )
            DEALLOCATE(ig_ip)
            DEALLOCATE(pw_ip)
          END IF
        ELSE
          IF ( mpime == root ) THEN
            DO i = 1, ngwl
                 pw(i) = PWT(ig_l2g(i))  
            END DO
          END IF
        END IF
        CALL MPI_BARRIER(gid, IERR)
      END DO

#elif ! defined __MPI

      DO I = 1, ngwl
           pw(i) = pwt( ig_l2g(i) ) 
      END DO

#else

      CALL errore(' SPLITWF ',' no communication protocol ',0)

#endif

      RETURN
      END SUBROUTINE splitwf

!=----------------------------------------------------------------------------=!

SUBROUTINE redistwf( c_dist_pw, c_dist_st, npw_p, nst_p, comm, idir )
   !
   !! Redistribute wave function.
   !
   USE kinds
   USE parallel_include

   implicit none

   COMPLEX(DP) :: c_dist_pw(:,:)
   !! the wave functions with plane waves distributed over processors 
   COMPLEX(DP) :: c_dist_st(:,:)
   !! the wave functions with electronic states distributed over processors 
   INTEGER, INTENT(IN) :: npw_p(:)
   !! the number of plane wave on each processor
   INTEGER, INTENT(IN) :: nst_p(:)
   !! the number of states on each processor
   INTEGER, INTENT(IN) :: comm
   !! group communicator
   INTEGER, INTENT(IN) :: idir
   !! direction of the redistribution:  
   !! \(\text{idir}>0\):  \(\text{c_dist_pw}\rightarrow\text{c_dist_st}\)  
   !! \(\text{idir}<0\):  \(\text{c_dist_pw}\leftarrow\text{c_dist_st}\)

   INTEGER :: mpime, nproc, ierr, npw_t, nst_t, proc, i, j, ngpww, ii
   INTEGER, ALLOCATABLE :: rdispls(:),  recvcount(:)
   INTEGER, ALLOCATABLE :: sendcount(:),  sdispls(:)
   COMPLEX(DP), ALLOCATABLE :: ctmp( : )

#if defined(__MPI)
   CALL mpi_comm_rank( comm, mpime, ierr )
   IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_comm_rank ', ierr )
   CALL mpi_comm_size( comm, nproc, ierr )
   IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_comm_size ', ierr )

   ALLOCATE( rdispls( nproc ), recvcount( nproc ), sendcount( nproc ), sdispls( nproc ) )

   npw_t = 0
   nst_t = 0
   DO proc=1,nproc
      sendcount(proc) = npw_p(mpime+1) * nst_p(proc)
      recvcount(proc) = npw_p(proc) * nst_p(mpime+1)
      npw_t = npw_t + npw_p(proc)
      nst_t = nst_t + nst_p(proc)
   END DO
   sdispls(1)=0
   rdispls(1)=0
   DO proc=2,nproc
      sdispls(proc) = sdispls(proc-1) + sendcount(proc-1)
      rdispls(proc) = rdispls(proc-1) + recvcount(proc-1)
   END DO

   ALLOCATE( ctmp( npw_t * nst_p( mpime + 1 ) ) )

   IF( idir > 0 ) THEN
      !
      ! ... Step 1. Communicate to all Procs so that each proc has all
      ! ... G-vectors and some states instead of all states and some
      ! ... G-vectors. This information is stored in the 1-d array ctmp.
      !
      CALL MPI_BARRIER( comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_barrier ', ierr )
      !
      CALL MPI_ALLTOALLV( c_dist_pw, sendcount, sdispls, MPI_DOUBLE_COMPLEX,             &
           &             ctmp, recvcount, rdispls, MPI_DOUBLE_COMPLEX, comm, ierr)
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_alltoallv ', ierr )
      !
      !   Step 2. Convert the 1-d array ctmp into a 2-d array consistent with the
      !   original notation c(ngw,nbsp). Psitot contains ntot = SUM_Procs(ngw) G-vecs
      !   and nstat states instead of all nbsp states
      !
      ngpww = 0
      DO proc = 1, nproc
         DO i = 1, nst_p(mpime+1)
            ii = (i-1) * npw_p(proc) 
            DO j = 1, npw_p(proc)
               c_dist_st( j + ngpww, i ) = ctmp( rdispls(proc) + j + ii )
            END DO
         END DO
         ngpww = ngpww + npw_p(proc)
      END DO

   ELSE
      !
      !   Step 4. Convert the 2-d array c_dist_st into 1-d array
      !
      ngpww = 0
      DO proc = 1, nproc
         DO i = 1, nst_p(mpime+1) 
            ii = (i-1) * npw_p(proc)
            DO j = 1, npw_p(proc)
               ctmp( rdispls(proc) + j + ii ) = c_dist_st( j + ngpww, i )
            END DO
         END DO
         ngpww = ngpww + npw_p(proc)
      END DO
      !        
      !   Step 5. Redistribute among processors. The result is stored in 2-d
      !   array c_dist_pw consistent with the notation c(ngw,nbsp)
      !
      CALL MPI_BARRIER( comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_barrier ', ierr )

      CALL MPI_ALLTOALLV( ctmp, recvcount, rdispls, MPI_DOUBLE_COMPLEX,          &
          &               c_dist_pw, sendcount , sdispls, MPI_DOUBLE_COMPLEX, comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_alltoallv ', ierr )


   END IF

   DEALLOCATE( ctmp )
   DEALLOCATE( rdispls, recvcount, sendcount, sdispls )
#endif
   RETURN
END SUBROUTINE redistwf

!=----------------------------------------------------------------------------=!

SUBROUTINE redistwfr( c_dist_pw, c_dist_st, npw_p, nst_p, comm, idir )
   !
   !!  Redistribute wave function.
   !
   USE kinds
   USE parallel_include

   implicit none

   REAL(DP) :: c_dist_pw(:,:)
   !! the wave functions with plane waves distributed over processors 
   REAL(DP) :: c_dist_st(:,:)
   !! the wave functions with electronic states distributed over processors 
   INTEGER, INTENT(IN) :: npw_p(:)
   !! the number of plane wave on each processor
   INTEGER, INTENT(IN) :: nst_p(:)
   !! the number of states on each processor
   INTEGER, INTENT(IN) :: comm
   !! group communicator
   INTEGER, INTENT(IN) :: idir
   !! direction of the redistribution:  
   !! \(\text{idir}>0\):  \(\text{c_dist_pw}\rightarrow\text{c_dist_st}\)  
   !! \(\text{idir}<0\):  \(\text{c_dist_pw}\leftarrow\text{c_dist_st}\)

   INTEGER :: mpime, nproc, ierr, npw_t, nst_t, proc, i, j, ngpww
   INTEGER, ALLOCATABLE :: rdispls(:),  recvcount(:)
   INTEGER, ALLOCATABLE :: sendcount(:),  sdispls(:)
   REAL(DP), ALLOCATABLE :: ctmp( : )

#if defined(__MPI)
   CALL mpi_comm_rank( comm, mpime, ierr )
   IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_comm_rank ', ierr )
   CALL mpi_comm_size( comm, nproc, ierr )
   IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_comm_size ', ierr )

   ALLOCATE( rdispls( nproc ), recvcount( nproc ), sendcount( nproc ), sdispls( nproc ) )

   npw_t = 0
   nst_t = 0
   DO proc=1,nproc
      sendcount(proc) = npw_p(mpime+1) * nst_p(proc)
      recvcount(proc) = npw_p(proc) * nst_p(mpime+1)
      npw_t = npw_t + npw_p(proc)
      nst_t = nst_t + nst_p(proc)
   END DO
   sdispls(1)=0
   rdispls(1)=0
   DO proc=2,nproc
      sdispls(proc) = sdispls(proc-1) + sendcount(proc-1)
      rdispls(proc) = rdispls(proc-1) + recvcount(proc-1)
   END DO

   ALLOCATE( ctmp( npw_t * nst_p( mpime + 1 ) ) )

   IF( idir > 0 ) THEN
      !
      ! ... Step 1. Communicate to all Procs so that each proc has all
      ! ... G-vectors and some states instead of all states and some
      ! ... G-vectors. This information is stored in the 1-d array ctmp.
      !
      CALL MPI_BARRIER( comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_barrier ', ierr )
      !
      CALL MPI_ALLTOALLV( c_dist_pw, sendcount, sdispls, MPI_DOUBLE_PRECISION,             &
           &             ctmp, recvcount, rdispls, MPI_DOUBLE_PRECISION, comm, ierr)
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_alltoallv ', ierr )
      !
      !   Step 2. Convert the 1-d array ctmp into a 2-d array consistent with the
      !   original notation c(ngw,nbsp). Psitot contains ntot = SUM_Procs(ngw) G-vecs
      !   and nstat states instead of all nbsp states
      !
      ngpww = 0
      DO proc = 1, nproc
         DO i = 1, nst_p(mpime+1)
            DO j = 1, npw_p(proc)
               c_dist_st( j + ngpww, i ) = ctmp( rdispls(proc) + j + (i-1) * npw_p(proc) )
            END DO
         END DO
         ngpww = ngpww + npw_p(proc)
      END DO

   ELSE
      !
      !   Step 4. Convert the 2-d array c_dist_st into 1-d array
      !
      ngpww = 0
      DO proc = 1, nproc
         DO i = 1, nst_p(mpime+1) 
            DO j = 1, npw_p(proc)
               ctmp( rdispls(proc) + j + (i-1) * npw_p(proc) ) = c_dist_st( j + ngpww, i )
            END DO
         END DO
         ngpww = ngpww + npw_p(proc)
      END DO
      !        
      !   Step 5. Redistribute among processors. The result is stored in 2-d
      !   array c_dist_pw consistent with the notation c(ngw,nbsp)
      !
      CALL MPI_BARRIER( comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_barrier ', ierr )

      CALL MPI_ALLTOALLV( ctmp, recvcount, rdispls, MPI_DOUBLE_PRECISION,          &
          &               c_dist_pw, sendcount , sdispls, MPI_DOUBLE_PRECISION, comm, ierr )
      IF( ierr /= 0 ) CALL errore( ' wf_redist ', ' mpi_alltoallv ', ierr )


   END IF

   DEALLOCATE( ctmp )
   DEALLOCATE( rdispls, recvcount, sendcount, sdispls )
#endif
   RETURN
END SUBROUTINE redistwfr

!=----------------------------------------------------------------------------=!

    END MODULE mp_wave

