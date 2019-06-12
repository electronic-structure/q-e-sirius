  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------
  SUBROUTINE readwfc( ipool, recn, evc0 )
  !--------------------------------------------------------
  !
  !  open wfc files as direct access, read, and close again
  !
  ! RM - Nov/Dec 2014
  ! Imported the noncolinear case implemented by xlzhang
  !
  !-------------------------------------------------------------
  !
  USE kinds,    ONLY : DP
  USE io_files, ONLY : prefix, tmp_dir
  USE units_lr, ONLY : lrwfc, iuwfc
  USE wvfct,    ONLY : npwx
  USE pwcom,    ONLY : nbnd
  USE noncollin_module, ONLY : npol
  USE mp_global,        ONLY : nproc_pool, me_pool, npool
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: recn
  !! kpoint number
  INTEGER, INTENT(in) :: ipool
  !! poolfile number to be read (not used in serial case)
  !
  COMPLEX(DP), INTENT(out) :: evc0(npwx*npol,nbnd)
  !! wavefunction is read from file
  !
  ! Local variables
  !
  INTEGER :: unf_recl, ios
  REAL(DP) :: dummy 
  CHARACTER(len=256) :: tempfile
  CHARACTER(len=3) :: nd_nmbr0
  ! file number for shuffle
  !
  !  open the wfc file, read and close
  !
  CALL set_ndnmbr( ipool, me_pool, nproc_pool, npool, nd_nmbr0 )
  !
#if defined(__MPI)
  tempfile = trim(tmp_dir) // trim(prefix) // '.wfc' // nd_nmbr0
# else
  tempfile = trim(tmp_dir) // trim(prefix) // '.wfc'
#endif
  INQUIRE (IOLENGTH = unf_recl) dummy 
  unf_recl = unf_recl * lrwfc
  !
  OPEN(iuwfc, file = tempfile, form = 'unformatted', &
       access = 'direct', iostat = ios, recl = unf_recl)
  IF (ios /= 0) CALL errore('readwfc', 'error opening wfc file', iuwfc)
  READ (iuwfc, rec = recn) evc0
  CLOSE(iuwfc, status = 'keep')
  !
  RETURN
  !
  END SUBROUTINE readwfc
