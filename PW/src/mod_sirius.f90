module mod_sirius
use input_parameters, only : use_sirius, sirius_cfg
use sirius
implicit none

logical :: use_sirius_radial_integration_beta = .true.
logical :: use_sirius_beta_projectors         = .true.
logical :: use_sirius_q_operator              = .false.
logical :: use_sirius_ks_solver               = .true.
logical :: use_sirius_density                 = .true.
logical :: use_sirius_density_matrix          = .true.
! initialize G-vectors once or at eeach step of ionic relaxation
logical :: init_gvec_once                     = .false.

! inverse of the reciprocal lattice vectors matrix
real(8) bg_inv(3,3)
! id of the k-point set for ground state calculations
integer kset_id
! total number of k-points
integer num_kpoints
real(8), allocatable :: kpoints(:,:)
real(8), allocatable :: wkpoints(:)

type atom_type_t
  ! atom label
  character(len=1, kind=C_CHAR) :: label(100) !c_string(len_trim(f_string) + 1)
  ! nh(iat) in the QE notation
  integer                 :: num_beta_projectors
  ! plane-wave coefficients of Q-operator
  complex(8), allocatable :: qpw(:, :)
end type atom_type_t

type(atom_type_t), allocatable :: atom_type(:)

contains

! this is not used at the moment as QE generates the effective potential itself
subroutine put_xc_functional_to_sirius
implicit none

  !if (get_meta().ne.0.or.get_inlc().ne.0) then
  !  write(*,*)get_igcx()
  !  write(*,*)get_igcc()
  !  write(*,*)get_meta()
  !  write(*,*)get_inlc()
  !  stop ("interface for this XC functional is not implemented")
  !endif

  !!== write(*,*)"xc_funtionals:", get_iexch(), get_icorr()

  !if (get_iexch().ne.0.and.get_igcx().eq.0) then
  !  select case(get_iexch())
  !  case(0)
  !  case(1)
  !    call sirius_add_xc_functional(c_str("XC_LDA_X"))
  !  case default
  !    stop ("interface for this exchange functional is not implemented")
  !  end select
  !endif

  !if (get_iexch().ne.0.and.get_igcx().ne.0) then
  !  select case(get_igcx())
  !  case(0)
  !  case(2)
  !    call sirius_add_xc_functional(c_str("XC_GGA_X_PW91"))
  !  case(3)
  !    call sirius_add_xc_functional(c_str("XC_GGA_X_PBE"))
  !  case default
  !    write(*,*)get_igcx()
  !    stop ("interface for this gradient exchange functional is not implemented")
  !  end select
  !endif

  !if (get_icorr().ne.0.and.get_igcc().eq.0) then
  !  select case(get_icorr())
  !  case(0)
  !  case(1)
  !    call sirius_add_xc_functional(c_str("XC_LDA_C_PZ"))
  !  case(4)
  !    call sirius_add_xc_functional(c_str("XC_LDA_C_PW"))
  !  case default
  !    stop ("interface for this correlation functional is not implemented")
  !  end select
  !endif

  !if (get_icorr().ne.0.and.get_igcc().ne.0) then
  !  select case(get_igcc())
  !  case(0)
  !  case(2)
  !    call sirius_add_xc_functional(c_str("XC_GGA_C_PW91"))
  !  case(4)
  !    call sirius_add_xc_functional(c_str("XC_GGA_C_PBE"))
  !  case default
  !    stop ("interface for this gradient correlation functional is not implemented")
  !  end select
  !endif
end subroutine put_xc_functional_to_sirius

subroutine setup_sirius()
use cell_base, only : alat, at, bg
use funct, only : get_iexch, get_icorr, get_inlc, get_meta, get_igcc, get_igcx
use ions_base, only : tau, nsp, atm, zv, amass, ityp, nat
use uspp_param, only : upf, nhm, nh
use atom, only : rgrid, msh
use fft_base, only :  dfftp
use klist, only : nks, xk, nkstot, wk
use gvect, only : ngm_g, ecutrho, ngm, mill
use gvecw, only : ecutwfc
use control_flags, only : gamma_only
use mp_pools, only : inter_pool_comm, npool
use mp_images,        only : nproc_image, intra_image_comm
use mp, only : mp_sum, mp_bcast
use wvfct, only : nbnd
use parallel_include
use sirius
use input_parameters, only : sirius_cfg
use noncollin_module, only : noncolin, npol, angle1, angle2
use lsda_mod, only : lsda, nspin, starting_magnetization
use cell_base, only : omega
use symm_base, only : nosym
use spin_orb,  only : lspinorb
use esm,       only : esm_local, esm_bc, do_comp_esm
use mp_diag, only : nproc_ortho
implicit none
!
integer :: dims(3), i, ia, iat, rank, ierr, ijv, ik, li, lj, mb, nb, j, l,&
     ilast, ir, num_gvec, num_ranks_k, vt(3), iwf, num_kp
real(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), tmp
real(8), allocatable :: dion(:, :), qij(:,:,:), vloc(:), wk_tmp(:), xk_tmp(:,:)
integer, allocatable :: nk_loc(:)
integer :: ih, jh, ijh, lmax_beta
logical(C_BOOL) bool_var

if (allocated(atom_type)) then
  do iat = 1, nsp
    if (allocated(atom_type(iat)%qpw)) deallocate(atom_type(iat)%qpw)
  enddo
  deallocate(atom_type)
endif

allocate(atom_type(nsp))

do iat = 1, nsp
  atom_type(iat)%label = c_str(atm(iat))
enddo

! create context of simulation
call sirius_create_simulation_context(c_str(trim(adjustl(sirius_cfg))), c_str("pseudopotential"), intra_image_comm)

! set number of first-variational states
call sirius_set_num_bands(nbnd)

bool_var = gamma_only
call sirius_set_gamma_point(bool_var)

num_ranks_k = nproc_image / npool
i = sqrt(dble(num_ranks_k) + 1d-10)
if (i * i .ne. num_ranks_k) then
  !stop ("not a square MPI grid")
  dims(1) = num_ranks_k
  dims(2) = 1
else 
  dims(1) = i
  dims(2) = i
endif
call sirius_set_mpi_grid_dims(2, dims(1))

!write(*,*)'i=',i
!write(*,*)'num_ranks_k=',num_ranks_k
!write(*,*)'nproc_ortho=',nproc_ortho

! set |G| cutoff of the dense FFT grid
! convert from G^2/2 Rydbergs to |G| in [a.u.^-1]
call sirius_set_pw_cutoff(sqrt(ecutrho))

! set |G+k| cutoff for the wave-functions
! convert from |G+k|^2/2 Rydbergs to |G+k| in [a.u.^-1]
call sirius_set_gk_cutoff(sqrt(ecutwfc))

if (lspinorb) then
   call sirius_set_num_mag_dims(3)
   call sirius_set_so_correction(.true.)
else
   if (noncolin) then
      call sirius_set_num_mag_dims(3)
   else
      if (nspin.eq.2) then
         call sirius_set_num_mag_dims(1)
      else
         call sirius_set_num_mag_dims(0)
      endif
   endif
endif

! set lattice vectors of the unit cell (length is in [a.u.])
a1(:) = at(:, 1) * alat
a2(:) = at(:, 2) * alat
a3(:) = at(:, 3) * alat
call sirius_set_lattice_vectors(a1(1), a2(1), a3(1))

vlat(:, 1) = a1(:)
vlat(:, 2) = a2(:)
vlat(:, 3) = a3(:)
! get the inverse of Bravais lattice vectors
call invert_mtrx(vlat, vlat_inv)
! get the inverse of reciprocal lattice vectors
call invert_mtrx(bg, bg_inv)
! get MPI rank associated with the distribution of k-points
call mpi_comm_rank(inter_pool_comm, rank, ierr)

! initialize atom types
do iat = 1, nsp

  ! get lmax_beta for this atom type
  lmax_beta = -1
  do i = 1, upf(iat)%nbeta
    lmax_beta = max(lmax_beta, upf(iat)%lll(i))
  enddo
  !if (upf(iat)%lmax .ne. lmax_beta) then
  !  write(*,*)
  !  write(*,'("Mismatch between lmax_beta and upf%lmax for atom type", I2)')iat
  !  write(*,'("  lmax =", I2)')upf(iat)%lmax
  !  write(*,'("  lmax_beta =", I2)')lmax_beta
  !endif

  ! add new atom type
  bool_var = upf(iat)%has_so
  call sirius_add_atom_type(atom_type(iat)%label, zn=nint(zv(iat)+0.001d0), mass=amass(iat), spin_orbit=bool_var)

  ! set radial grid
  call sirius_set_atom_type_radial_grid(atom_type(iat)%label, upf(iat)%mesh, upf(iat)%r(1))

  ! set beta-projectors
  do i = 1, upf(iat)%nbeta
    call sirius_add_atom_type_beta_radial_function(atom_type(iat)%label, upf(iat)%lll(i),&
                                                  &upf(iat)%beta(1, i), upf(iat)%kbeta(i))
  enddo

  ! set the atomic radial functions
  do iwf = 1, upf(iat)%nwfc
    l = upf(iat)%lchi(iwf)
    call sirius_add_atom_type_ps_atomic_wf(c_str(atm(iat)), l, upf(iat)%chi(1, iwf), 0.d0, msh(iat))
  enddo

  allocate(dion(upf(iat)%nbeta, upf(iat)%nbeta))
  ! convert to hartree
  do i = 1, upf(iat)%nbeta
    do j = 1, upf(iat)%nbeta
      dion(i, j) = upf(iat)%dion(i, j) / 2.d0
    end do
  end do
  ! sed d^{ion}_{i,j}
  call sirius_set_atom_type_dion(c_str(atm(iat)), upf(iat)%nbeta, dion(1, 1))
  deallocate(dion)

  ! set radial function of augmentation charge
  if (upf(iat)%tvanp) then
    !do l = 0, upf(iat)%nqlc - 1
    do l = 0, 2 * lmax_beta
      do i = 0, upf(iat)%nbeta - 1
        do j = i, upf(iat)%nbeta - 1
          ijv = j * (j + 1) / 2 + i + 1
          call sirius_add_atom_type_q_radial_function(c_str(atm(iat)), l, i, j,&
                                                     &upf(iat)%qfuncl(1, ijv, l), upf(iat)%kkbeta)
        enddo
      enddo
    enddo
  endif

  if (upf(iat)%tpawp) then ! TODO: cleaup this
    call sirius_set_atom_type_paw_data(c_str(atm(iat)), upf(iat)%aewfc(1,1), upf(iat)%pswfc(1,1),&
         &upf(iat)%nbeta, upf(iat)%mesh, upf(iat)%paw%iraug,&
         &upf(iat)%paw%core_energy, upf(iat)%paw%ae_rho_atc(1),&
         &upf(iat)%mesh, upf(iat)%paw%oc(1), upf(iat)%nbeta )
  endif
  
  ! TODO: pass PW coefficients of this functions

  !! set non-linear core correction
  !if (associated(upf(iat)%rho_atc)) then
  !  call sirius_set_atom_type_rho_core(c_str(atm(iat)), upf(iat)%mesh, upf(iat)%rho_atc(1))
  !else
  !  allocate(vloc(upf(iat)%mesh))
  !  vloc = 0.d0
  !  call sirius_set_atom_type_rho_core(c_str(atm(iat)), upf(iat)%mesh, vloc(1))
  !  deallocate(vloc)
  !endif

  !! set total charge density of a free atom (to compute initial rho(r))
  !call sirius_set_atom_type_rho_tot(c_str(atm(iat)), upf(iat)%mesh, upf(iat)%rho_at(1))

  !allocate(vloc(upf(iat)%mesh)) ! TODO: cut vloc to 10 a.u. here, not in the SIRIUS code
  !! convert to Hartree                ! issue a warning in SIRIUS if the tail of vloc is not -z/r
  !do i = 1, upf(iat)%mesh
  !  vloc(i) = upf(iat)%vloc(i) / 2.d0
  !end do
  !! set local part of pseudo-potential
  !call sirius_set_atom_type_vloc(c_str(atm(iat)), upf(iat)%mesh, vloc(1))
  !deallocate(vloc)
enddo
  
! add atoms to the unit cell
! WARNING: sirius accepts only fractional coordinates;
!          if QE stores coordinates in a different way, the conversion must be made here
do ia = 1, nat
  iat = ityp(ia)
  ! Cartesian coordinates
  v1(:) = tau(:, ia) * alat
  ! fractional coordinates
  v1(:) = matmul(vlat_inv, v1)
  ! reduce coordinates to [0, 1) interval
  call sirius_reduce_coordinates(v1(1), v2(1), vt(1))
  if (noncolin) then
    v1(1) = zv(iat) * starting_magnetization(iat) * sin(angle1(iat)) * cos(angle2(iat))
    v1(2) = zv(iat) * starting_magnetization(iat) * sin(angle1(iat)) * sin(angle2(iat))
    v1(3) = zv(iat) * starting_magnetization(iat) * cos(angle1(iat))
  else
    v1 = 0
    v1(3) = zv(iat) * starting_magnetization(iat)
  endif
  call sirius_add_atom(c_str(atm(iat)), v2(1), v1(1))
enddo

! QE is taking care of symmetry
!if (nosym) then
  call sirius_set_use_symmetry(0)
!endif

bool_var = do_comp_esm
call sirius_set_esm(bool_var, esm_bc)

! initialize global variables/indices/arrays/etc. of the simulation
call sirius_initialize_simulation_context()

call sirius_create_potential
  
!! get number of g-vectors of the dense fft grid
!call sirius_get_num_gvec(num_gvec)
!
!! TODO: number of G-vectors can be different; adapt the code wo work in this situation
!if (.not.((num_gvec .eq. ngm_g) .or. (num_gvec * 2 - 1 .eq. ngm_g))) then
!  write(*,*)"wrong number of g-vectors"
!  write(*,*)"num_gvec=",num_gvec
!  write(*,*)"ngm_g=",ngm_g
!endif

!call sirius_get_fft_grid_size(dims(1)) ! TODO: size of FFT box is not very relevant and in principle can be slightly different
!if (dims(1).ne.dfftp%nr1.or.dims(2).ne.dfftp%nr2.or.dims(3).ne.dfftp%nr3) then
!  write(*,*)"wrong fft grid dimensions"
!  write(*,*)"qe: ", dfftp%nr1,  dfftp%nr2,  dfftp%nr3
!  write(*,*)"sirius: ", dims
!  stop 111
!endif

!!== i = 1
!!== if (nosym) i = 0
!!== kmesh(:) = (/nk1, nk2, nk3/)
!!== kshift(:) = (/k1, k2, k3/)
!!== call sirius_create_irreducible_kset(kmesh, kshift, i, kset_id)

if (.true.) then
  do i = 1, num_kpoints
    write(*,*)'ik=',i,' kpoint=',matmul(bg_inv,kpoints(:,i))
  enddo
endif

allocate(wk_tmp(nkstot))
allocate(xk_tmp(3, nkstot))
! weights of k-points in SIRIUS must sum to one
do i = 1, nkstot
  if (nspin.eq.1) then
    wk_tmp(i) = wk(i) / 2.d0
  else
    wk_tmp(i) = wk(i)
  endif
  xk_tmp(:,i) = xk(:,i)
end do

call mpi_bcast(wk_tmp(1),        nkstot, mpi_double, 0, inter_pool_comm, ierr)
call mpi_bcast(xk_tmp(1, 1), 3 * nkstot, mpi_double, 0, inter_pool_comm, ierr)

! convert to fractional coordinates
do ik = 1, nkstot
  xk_tmp(:, ik) = matmul(bg_inv, xk_tmp(:, ik))
end do

allocate(nk_loc(0:npool-1))
nk_loc = 0
nk_loc(rank) = nks
call mp_sum(nk_loc, inter_pool_comm)
if (nspin.eq.2) then
  nk_loc(:) = nk_loc(:)
endif

if (nspin.eq.2) then
  num_kp = nkstot / 2
else
  num_kp = nkstot
endif
call sirius_create_kset(num_kp, xk_tmp(1, 1), wk_tmp(1), 1, kset_id)

! create Density class
call sirius_create_density()
 
! create Potential class
call sirius_create_potential()

! create ground-state class
call sirius_create_ground_state(kset_id)

! create ground-state class
! create a set of k-points
! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice
!if (nspin.eq.2) then
!  call sirius_create_kset(nkstot / 2, xk_tmp(1, 1), wk_tmp(1), 1, kset_id, nk_loc(0))
!else
!  call sirius_create_kset(nkstot, xk_tmp(1, 1), wk_tmp(1), 1, kset_id, nk_loc(0))
!endif
deallocate(wk_tmp)
deallocate(xk_tmp)
deallocate(nk_loc)


end subroutine setup_sirius


subroutine get_q_operator_from_sirius
use uspp_param, only : upf, nh, nhm
use ions_base,  only : nsp
use gvect,      only : ngm, mill
use uspp,       only : qq_nt
implicit none
integer iat, ih, jh, ijh, i

do iat = 1, nsp
  call sirius_get_num_beta_projectors(atom_type(iat)%label, atom_type(iat)%num_beta_projectors)
  if (nh(iat).ne.atom_type(iat)%num_beta_projectors) then
    stop 'wrong number of beta projectors'
  endif
  if (upf(iat)%tvanp) then
    i = atom_type(iat)%num_beta_projectors
    allocate(atom_type(iat)%qpw(ngm, i * (i + 1) / 2))
    ijh = 0
    do ih = 1, atom_type(iat)%num_beta_projectors
      do jh = ih, atom_type(iat)%num_beta_projectors
        ijh = ijh + 1
        call sirius_get_q_operator(atom_type(iat)%label, ih, jh, ngm, mill(1, 1), atom_type(iat)%qpw(1, ijh))
        call sirius_get_q_operator_matrix(atom_type(iat)%label, qq_nt(1, 1, iat), nhm)
      enddo
    enddo
  endif
enddo

end subroutine get_q_operator_from_sirius

subroutine invert_mtrx(vlat, vlat_inv)
  implicit none
  real(8), intent(in) :: vlat(3,3)
  real(8), intent(out) :: vlat_inv(3, 3)
  real(8) d1

  d1 = vlat(1,2)*vlat(2,3)*vlat(3,1)-vlat(1,3)*vlat(2,2)*vlat(3,1)+vlat(1,3)*vlat(2,1)*vlat(3,2) &
  &   -vlat(1,1)*vlat(2,3)*vlat(3,2)+vlat(1,1)*vlat(2,2)*vlat(3,3)-vlat(1,2)*vlat(2,1)*vlat(3,3)
  d1 = 1.d0 / d1
  vlat_inv(1,1)=(vlat(2,2)*vlat(3,3)-vlat(2,3)*vlat(3,2))*d1
  vlat_inv(1,2)=(vlat(1,3)*vlat(3,2)-vlat(1,2)*vlat(3,3))*d1
  vlat_inv(1,3)=(vlat(1,2)*vlat(2,3)-vlat(1,3)*vlat(2,2))*d1
  vlat_inv(2,1)=(vlat(2,3)*vlat(3,1)-vlat(2,1)*vlat(3,3))*d1
  vlat_inv(2,2)=(vlat(1,1)*vlat(3,3)-vlat(1,3)*vlat(3,1))*d1
  vlat_inv(2,3)=(vlat(1,3)*vlat(2,1)-vlat(1,1)*vlat(2,3))*d1
  vlat_inv(3,1)=(vlat(2,1)*vlat(3,2)-vlat(2,2)*vlat(3,1))*d1
  vlat_inv(3,2)=(vlat(1,2)*vlat(3,1)-vlat(1,1)*vlat(3,2))*d1
  vlat_inv(3,3)=(vlat(1,1)*vlat(2,2)-vlat(1,2)*vlat(2,1))*d1
end subroutine

subroutine get_band_energies_from_sirius
  !
  use wvfct,    only : nbnd, et
  use klist,    only : nkstot, nks
  use lsda_mod, only : nspin
  use sirius
  !
  implicit none
  !
  integer, external :: global_kpoint_index
  !
  real(8), allocatable :: band_e(:,:)
  integer :: ik, nk, nb, nfv

  allocate(band_e(nbnd, nkstot))
  
  ! get band energies
  if (nspin.ne.2) then
    ! non-magnetic or non-collinear case
    do ik = 1, nkstot
      call sirius_get_band_energies(kset_id, ik, 0, band_e(1, ik))
    end do
  else
    ! collinear magnetic case
    nk = nkstot / 2
    ! get band energies
    do ik = 1, nk
      call sirius_get_band_energies(kset_id, ik, 0, band_e(1, ik))
      call sirius_get_band_energies(kset_id, ik, 1, band_e(1, nk + ik))
    end do
  
  endif
  
  ! convert to Ry
  do ik = 1, nks
    et(:, ik) = 2.d0 * band_e(:, global_kpoint_index(nkstot, ik))
  enddo
  
  deallocate(band_e)

end subroutine get_band_energies_from_sirius

subroutine put_band_occupancies_to_sirius
  !
  use wvfct,    only : nbnd, wg
  use klist,    only : nkstot, nks, wk
  use lsda_mod, only : nspin
  use mp_pools, only : inter_pool_comm
  use parallel_include
  use sirius
  !
  implicit none
  !
  integer, external :: global_kpoint_index
  !
  real(8), allocatable :: bnd_occ(:, :)
  real(8) :: maxocc
  integer :: ik, ierr, nk, nb
  
  ! compute occupancies
  allocate(bnd_occ(nbnd, nkstot))
  bnd_occ = 0.d0
  ! define a maximum band occupancy (2 in case of spin-unpolarized, 1 in case of spin-polarized)
  maxocc = 2.d0
  if (nspin.gt.1) then
    maxocc = 1.d0
  endif
  do ik = 1, nks
    bnd_occ(:, global_kpoint_index(nkstot, ik)) = maxocc * wg(:, ik) / wk(ik)
  enddo
  call mpi_allreduce(MPI_IN_PLACE, bnd_occ(1, 1), nbnd * nkstot, MPI_DOUBLE, MPI_SUM, inter_pool_comm, ierr)
  
  if (nspin.ne.2) then
    ! set band occupancies
    do ik = 1, nkstot
      call sirius_set_band_occupancies(kset_id, ik, 0, bnd_occ(1, ik))
    enddo
  else 
    nk = nkstot / 2
    do ik = 1, nk
      call sirius_set_band_occupancies(kset_id, ik, 0, bnd_occ(1, ik))
      call sirius_set_band_occupancies(kset_id, ik, 1, bnd_occ(1, ik + nk))
    enddo
  endif
  
  deallocate(bnd_occ)

end subroutine put_band_occupancies_to_sirius

subroutine get_density_matrix_from_sirius
use scf,        only : rho
use ions_base,  only : nat, nsp, ityp
use uspp_param, only : nhm, nh
use lsda_mod,   only : nspin
implicit none
complex(8), allocatable :: dens_mtrx(:,:,:)
integer iat, na, ijh, ih, jh, ispn
! complex density matrix in SIRIUS has at maximum three components
allocate(dens_mtrx(nhm, nhm, 3))
do iat = 1, nsp
  do na = 1, nat
    if (ityp(na).eq.iat.and.allocated(rho%bec)) then
      rho%bec(:, na, :) = 0.d0
      call sirius_get_density_matrix(na, dens_mtrx(1, 1, 1), nhm)

      ijh = 0
      do ih = 1, nh(iat)
        do jh = ih, nh(iat)
          ijh = ijh + 1
          if (nspin.le.2) then
            do ispn = 1, nspin
              rho%bec(ijh, na, ispn) = dreal(dens_mtrx(ih, jh, ispn))
            enddo
          endif
          if (nspin.eq.4) then
            rho%bec(ijh, na, 1) = dreal(dens_mtrx(ih, jh, 1) + dens_mtrx(ih, jh, 2))
            rho%bec(ijh, na, 4) = dreal(dens_mtrx(ih, jh, 1) - dens_mtrx(ih, jh, 2))
            rho%bec(ijh, na, 2) = 2.d0 * dreal(dens_mtrx(ih, jh, 3))
            rho%bec(ijh, na, 3) = -2.d0 * dimag(dens_mtrx(ih, jh, 3))
          endif
          ! off-diagonal elements have a weight of 2
          if (ih.ne.jh) then
            do ispn = 1, nspin
              rho%bec(ijh, na, ispn) = rho%bec(ijh, na, ispn) * 2.d0
            enddo
          endif
        enddo
      enddo
    endif
  enddo
enddo
deallocate(dens_mtrx)
end subroutine get_density_matrix_from_sirius

subroutine get_density_from_sirius
  !
  use scf,        only : rho
  use gvect,      only : mill, ngm
  use mp_bands,   only : intra_bgrp_comm
  use lsda_mod,   only : nspin
  use ions_base,  only : nat, nsp, ityp
  use uspp_param, only : nhm, nh
  use sirius
  !
  implicit none
  !
  integer iat, ig, ih, jh, ijh, na, ispn
  complex(8) z1, z2
  
  ! get rho(G)
  call sirius_get_pw_coeffs(c_str("rho"), rho%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
  
  if (nspin.eq.2) then
    call sirius_get_pw_coeffs(c_str("magz"), rho%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
    ! convert to rho_{up}, rho_{dn}
    do ig = 1, ngm
      z1 = rho%of_g(ig, 1)
      z2 = rho%of_g(ig, 2)
      rho%of_g(ig, 1) = 0.5 * (z1 + z2)
      rho%of_g(ig, 2) = 0.5 * (z1 - z2)
    enddo
  endif
  if (nspin.eq.4) then
    call sirius_get_pw_coeffs(c_str("magx"), rho%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_get_pw_coeffs(c_str("magy"), rho%of_g(1, 3), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_get_pw_coeffs(c_str("magz"), rho%of_g(1, 4), ngm, mill(1, 1), intra_bgrp_comm)
  endif
  ! get density matrix
  call get_density_matrix_from_sirius
end subroutine get_density_from_sirius

!subroutine put_density_to_sirius
!  !
!  use scf,        only : rho
!  use gvect,      only : mill, ngm
!  use mp_bands,   only : intra_bgrp_comm
!  use lsda_mod,   only : nspin
!  use ions_base,  only : nat, nsp, ityp
!  use uspp_param, only : nhm, nh
!  use sirius
!  implicit none
!  !
!  complex(8), allocatable :: rho_tot(:), mag(:)
!  complex(8), allocatable :: dens_mtrx(:,:,:)
!  integer iat, ig, ih, jh, ijh, na, ispn
!  real(8) :: fact
!  !
!  if (nspin.eq.1.or.nspin.eq.4) then
!    call sirius_set_pw_coeffs(c_str("rho"), rho%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
!  endif
!
!  if (nspin.eq.2) then
!    allocate(rho_tot(ngm))
!    allocate(mag(ngm))
!    do ig = 1, ngm
!      rho_tot(ig) = rho%of_g(ig, 1) + rho%of_g(ig, 2)
!      mag(ig) = rho%of_g(ig, 1) - rho%of_g(ig, 2)
!    enddo
!    call sirius_set_pw_coeffs(c_str("rho"), rho_tot(1), ngm, mill(1, 1), intra_bgrp_comm)
!    call sirius_set_pw_coeffs(c_str("magz"), mag(1), ngm, mill(1, 1), intra_bgrp_comm)
!    deallocate(rho_tot)
!    deallocate(mag)
!  endif
!
!  if (nspin.eq.4) then
!    call sirius_set_pw_coeffs(c_str("magx"), rho%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
!    call sirius_set_pw_coeffs(c_str("magy"), rho%of_g(1, 3), ngm, mill(1, 1), intra_bgrp_comm)
!    call sirius_set_pw_coeffs(c_str("magz"), rho%of_g(1, 4), ngm, mill(1, 1), intra_bgrp_comm)
!  endif
!
!  !!== ! set density matrix
!  !!== ! complex density matrix in SIRIUS has at maximum three components
!  !!== allocate(dens_mtrx(nhm, nhm, 3))
!  !!== do iat = 1, nsp
!  !!==   do na = 1, nat
!  !!==     if (ityp(na).eq.iat.and.allocated(rho%bec)) then
!  !!==       dens_mtrx = (0.d0, 0.d0)
!  !!==       ijh = 0
!  !!==       do ih = 1, nh(iat)
!  !!==         do jh = ih, nh(iat)
!  !!==           ijh = ijh + 1
!  !!==           ! off-diagonal elements have a weight of 2
!  !!==           if (ih.ne.jh) then
!  !!==             fact = 0.5d0
!  !!==           else
!  !!==             fact = 1.d0
!  !!==           endif
!  !!==           if (nspin.le.2) then
!  !!==             do ispn = 1, nspin
!  !!==               dens_mtrx(ih, jh, ispn) = fact * rho%bec(ijh, na, ispn)
!  !!==               dens_mtrx(jh, ih, ispn) = fact * rho%bec(ijh, na, ispn)
!  !!==             enddo
!  !!==           endif
!  !!==           if (nspin.eq.4) then
!  !!==             ! 0.5 * (rho + mz)
!  !!==             dens_mtrx(ih, jh, 1) = fact * 0.5 * (rho%bec(ijh, na, 1) + rho%bec(ijh, na, 4))
!  !!==             dens_mtrx(jh, ih, 1) = fact * 0.5 * (rho%bec(ijh, na, 1) + rho%bec(ijh, na, 4))
!  !!==             ! 0.5 * (rho - mz)
!  !!==             dens_mtrx(ih, jh, 2) = fact * 0.5 * (rho%bec(ijh, na, 1) - rho%bec(ijh, na, 4))
!  !!==             dens_mtrx(jh, ih, 2) = fact * 0.5 * (rho%bec(ijh, na, 1) - rho%bec(ijh, na, 4))
!  !!==             ! 0.5 * (mx - I * my)
!  !!==             dens_mtrx(ih, jh, 3) = fact * 0.5 * dcmplx(rho%bec(ijh, na, 2), -rho%bec(ijh, na, 3))
!  !!==             dens_mtrx(jh, ih, 3) = fact * 0.5 * dcmplx(rho%bec(ijh, na, 2), -rho%bec(ijh, na, 3))
!  !!==           endif
!  !!==         enddo
!  !!==       enddo
!  !!==       call sirius_set_density_matrix(na, dens_mtrx(1, 1, 1), nhm)
!  !!==     endif
!  !!==   enddo
!  !!== enddo
!  !!== deallocate(dens_mtrx)
!
!end subroutine put_density_to_sirius

subroutine put_potential_to_sirius
  use scf,                  only : v, vltot
  use gvect,                only : mill, ngm
  use mp_bands,             only : intra_bgrp_comm
  use lsda_mod,             only : nspin
  use noncollin_module,     only : nspin_mag
  USE wavefunctions_module,  ONLY: psic

  use fft_base,             only : dfftp
  USE fft_interfaces,       ONLY : fwfft
  use uspp,                 only : deeq
  use uspp_param,           only : nhm
  use paw_variables,        only : okpaw
  use ions_base,            only : nat
  use sirius
  !
  implicit none
  !
  complex(8), allocatable :: vxcg(:)
  complex(8) :: z1, z2
  integer ig, is, ir, i, ia, j
  character(10) label
  real(8), allocatable :: deeq_tmp(:,:)
  real(8) :: d1,d2
  !
  if (nspin.eq.1.or.nspin.eq.4) then
    ! add local part of the potential and transform to PW domain
    psic(:) = v%of_r(:, 1) + vltot(:)
    call fwfft('Rho', psic, dfftp)
    ! convert to Hartree
    do ig = 1, ngm
      v%of_g(ig, 1) = psic(dfftp%nl(ig)) * 0.5d0
    enddo
    ! set effective potential
    call sirius_set_pw_coeffs(c_str("veff"), v%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
  endif

  if (nspin.eq.2) then
    do is = 1, 2
      ! add local part of the potential and transform to PW domain
      psic(:) = v%of_r(:, is) + vltot(:)
      call fwfft('Rho', psic, dfftp)
      ! convert to Hartree
      do ig = 1, ngm
         v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
      enddo
    enddo

    do ig = 1, ngm
      z1 = v%of_g(ig, 1) 
      z2 = v%of_g(ig, 2)
      v%of_g(ig, 1) = 0.5 * (z1 + z2)
      v%of_g(ig, 2) = 0.5 * (z1 - z2)
    enddo
    ! set effective potential and magnetization
    call sirius_set_pw_coeffs(c_str("veff"), v%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_set_pw_coeffs(c_str("bz"),   v%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
  endif

  if (nspin.eq.4) then
    do is = 2, nspin_mag
      psic(:) = v%of_r(:, is)
      call fwfft('Rho', psic, dfftp)
      ! convert to Hartree
      do ig = 1, ngm
        v%of_g(ig, is) = psic(dfftp%nl(ig)) * 0.5d0
      enddo
      if (is.eq.2) label="bx"
      if (is.eq.3) label="by"
      if (is.eq.4) label="bz"

      call sirius_set_pw_coeffs(c_str(label), v%of_g(1, is), ngm, mill(1, 1), intra_bgrp_comm)
    enddo
  endif

!  ! convert Vxc to plane-wave domain
!  if (nspin.eq.1.or.nspin.eq.4) then
!    do ir = 1, dfftp%nnr
!      psic(ir) = vxc(ir, 1)
!    enddo
!  else
!    do ir = 1, dfftp%nnr
!      psic(ir) = 0.5d0 * (vxc(ir, 1) + vxc(ir, 2))
!    enddo
!  endif
!  call fwfft('Rho', psic, dfftp)
!  allocate(vxcg(ngm))
!  ! convert to Hartree
!  do ig = 1, ngm
!     vxcg(ig) = psic(dfftp%nl(ig)) * 0.5d0
!  end do
!  ! set XC potential
!  call sirius_set_pw_coeffs(c_str("vxc"), vxcg(1), ngm, mill(1, 1), intra_bgrp_comm)
!  deallocate(vxcg)
!
  ! update D-operator matrix
  !call sirius_generate_d_operator_matrix()
  !if (okpaw) then
    allocate(deeq_tmp(nhm, nhm))
  !  !! get D-operator matrix
  !  !do ia = 1, nat
  !  !  do is = 1, nspin
  !  !    call sirius_get_d_operator_matrix(ia, is, deeq(1, 1, ia, is), nhm)
  !  !  enddo
  !  !  if (nspin.eq.2) then
  !  !    do i = 1, nhm
  !  !      do j = 1, nhm
  !  !        d1 = deeq(i, j, ia, 1)
  !  !        d2 = deeq(i, j, ia, 2)
  !  !        deeq(i, j, ia, 1) = d1 + d2
  !  !        deeq(i, j, ia, 2) = d1 - d2
  !  !      enddo
  !  !    enddo
  !  !  endif
  !  !  ! convert to Ry
  !  !  deeq(:, :, ia, :) = deeq(:, :, ia, :) * 2
  !  !enddo
  !  !call add_paw_to_deeq(deeq)
    do ia = 1, nat
      do is = 1, nspin
        if (nspin.eq.2.and.is.eq.1) then
          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) + deeq(:, :, ia, 2)) / 2 ! convert to Ha
        endif
        if (nspin.eq.2.and.is.eq.2) then
          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) - deeq(:, :, ia, 2)) / 2 ! convert to Ha
        endif
        if (nspin.eq.1.or.nspin.eq.4) then
          deeq_tmp(:, :) = deeq(:, :, ia, is) / 2 ! convert to Ha
        endif
        call sirius_set_d_operator_matrix(ia, is, deeq_tmp(1, 1), nhm)
      enddo
    enddo
    deallocate(deeq_tmp)
  !endif

end subroutine put_potential_to_sirius

subroutine put_q_operator_matrix_to_sirius
use uspp,       only : qq_nt
use uspp_param, only : upf, nhm
use ions_base,  only : nsp, atm
implicit none
integer iat

do iat = 1, nsp
  if (upf(iat)%tvanp) then
    call sirius_set_q_operator_matrix(c_str(atm(iat)), qq_nt(1, 1, iat), nhm)
  endif
enddo

end subroutine put_q_operator_matrix_to_sirius

subroutine get_wave_functions_from_sirius
use klist, only : nkstot, nks, ngk, igk_k
use gvect, only : mill
use buffers, only : save_buffer
use io_files, only : iunwfc, nwordwfc
use bp, only : lelfield
use noncollin_module, only : npol
use wvfct, only : npwx, nbnd
use wavefunctions_module, only : evc
implicit none
integer, external :: global_kpoint_index
integer, allocatable :: gvl(:,:)
integer ig, ik, ik_, i, j
complex(8) z1

allocate(gvl(3, npwx))
do ik = 1, nks
  do ig = 1, ngk(ik)
    gvl(:,ig) = mill(:, igk_k(ig, ik))
  enddo
  !
  ik_ = global_kpoint_index(nkstot, ik)
  call sirius_get_wave_functions(kset_id, ik_, ngk(ik), gvl(1, 1), evc(1, 1), npwx * npol) 
  !write(*,*)'checking wfs for k-point ', ik_
  !do i = 1, nbnd
  !  do j = 1, nbnd
  !    z1 = 0.d0
  !    do ig = 1, ngk(ik)
  !      z1 = z1 + conjg(evc(ig, i)) * evc(ig, j)
  !    enddo
  !    if (i.eq.j) z1 = z1 - 1.d0
  !    if (abs(z1).gt.1e-10) then
  !      write(*,*)'not orthogonal ',i,j,' diff: ',z1
  !    endif
  !  enddo
  !enddo
  !
  IF ( nks > 1 .OR. lelfield ) &
    CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
enddo
deallocate(gvl)

end subroutine get_wave_functions_from_sirius

!subroutine put_vltot_to_sirius
!  use scf,       only : vltot
!  use gvect, only : mill, ngm
!  use wavefunctions_module, only : psic
!  use fft_interfaces,       only : fwfft, invfft
!  use fft_base,             only : dfftp
!  use mp_bands, only : intra_bgrp_comm
!  use sirius
!  !
!  implicit none
!  !
!  complex(8), allocatable :: vg(:)
!  integer ig
!  !
!  allocate(vg(ngm))
!  psic(:) = vltot(:)
!  call fwfft('Rho', psic, dfftp)
!  ! convert to Hartree
!  do ig = 1, ngm
!    vg(ig) = psic(dfftp%nl(ig)) * 0.5d0 ! convert to Ha
!  enddo
!  ! set local potential
!  call sirius_set_pw_coeffs(c_str("vloc"), vg(1), ngm, mill(1, 1), intra_bgrp_comm)
!  deallocate(vg)
!end subroutine put_vltot_to_sirius

!subroutine get_rhoc_from_sirius
!  use uspp_param,only : upf
!  use ener,      only : etxcc
!  use scf,       only : rho_core, rhog_core
!  use control_flags, only : gamma_only
!  use wavefunctions_module, only : psic
!  use gvect, only : mill, ngm
!  use scf, only : rho_core, rhog_core
!  use mp_bands, only : intra_bgrp_comm
!  use ions_base, only : ntyp => nsp
!  use fft_interfaces,only : invfft
!  use fft_base,  only : dfftp
!  use sirius
!  !
!  implicit none
!  
!  etxcc = 0.0d0
!  if (any(upf(1:ntyp)%nlcc)) then
!    call sirius_get_pw_coeffs(c_str("rhoc"), rhog_core(1), ngm, mill(1, 1), intra_bgrp_comm)
!    psic(:) = (0.d0, 0.d0)
!    psic(dfftp%nl(:)) = rhog_core(:)
!    if (gamma_only) psic(dfftp%nlm(:)) = conjg(rhog_core(:))
!    call invfft ('Rho', psic, dfftp)
!    rho_core(:) = psic(:)
!  else 
!    rhog_core(:) = 0.0d0
!    rho_core(:)  = 0.0d0
!  endif
!
!end subroutine get_rhoc_from_sirius

!subroutine set_vloc_sirius
!use sirius
!use gvect, only : ngm, mill, igtongl, ngl
!use vlocal, only : vloc
!use mp_bands, only : intra_bgrp_comm
!use ions_base, only : atm
!integer nt,i
!real(8), allocatable :: tmp(:)
!
!allocate(tmp(ngm))
!vloc(:,:) = 0.d0
!do nt = 1, ntyp
!  call sirius_get_pw_coeffs_real(c_str(atm(nt)), c_str("vloc"), tmp(1), ngm, mill(1, 1), intra_bgrp_comm)
!  do i = 1, ngm
!    vloc(igtongl(i), nt) = tmp(i) * 2 ! convert to Ry
!  enddo
!enddo
!
!deallocate(tmp, tmp1)
!
!call set_vloc_sirius
!CALL setlocal()
!end subroutine

!subroutine get_vloc_from_sirius
!  use wavefunctions_module, only : psic
!  use gvect, only : mill, ngm, gg
!  use scf, only: vltot, v_of_0
!  use fft_interfaces, only : fwfft, invfft
!  use fft_base, only : dfftp
!  use constants, only : eps8
!  use control_flags, only : gamma_only
!  use mp_bands, only : intra_bgrp_comm
!  use mp, only : mp_bcast, mp_sum
!  use sirius
!  !
!  implicit none
!  !
!  complex(8), allocatable :: vpw(:)
!  allocate(vpw(ngm))
!  call sirius_get_pw_coeffs(c_str("vloc"), vpw(1), ngm, mill(1, 1), intra_bgrp_comm)
!  psic(:) = 0.d0
!  psic(dfftp%nl(:)) = vpw(:)
!  if (gamma_only) psic(dfftp%nlm(:)) = conjg(vpw(:))
!  call invfft('Rho', psic, dfftp)
!  vltot(:) = dble(psic(:)) * 2 ! convert to Ry
!  v_of_0=0.d0
!  IF (gg(1) < eps8) v_of_0 = dble(vpw(1))
!  !
!  call mp_sum(v_of_0, intra_bgrp_comm)
!  deallocate(vpw)
!
!end subroutine get_vloc_from_sirius

!subroutine get_density_matrix_from_sirius
!implicit none
!real(8), allocatable :: dens_mtrx(:, :, :)
!integer iat, na
!
!allocate(dens_mtrx(nhm, nhm, 3))
!do iat = 1, nsp
!  do na = 1, nat
!    if (ityp(na).eq.iat.and.allocated(rho%bec)) then
!      rho%bec(:, na, :) = 0.d0
!      call sirius_get_density_matrix(na, dens_mtrx(1, 1, 1), nhm)
!
!      ijh = 0
!      do ih = 1, nh(iat)
!        do jh = ih, nh(iat)
!          ijh = ijh + 1
!          if (nspin.le.2) then
!            do ispn = 1, nspin
!              rho%bec(ijh, na, ispn) = dreal(dens_mtrx(ih, jh, ispn))
!            enddo
!          endif
!          if (nspin.eq.4) then
!            rho%bec(ijh, na, 1) = dreal(dens_mtrx(ih, jh, 1) + dens_mtrx(ih, jh, 2))
!            rho%bec(ijh, na, 4) = dreal(dens_mtrx(ih, jh, 1) - dens_mtrx(ih, jh, 2))
!            rho%bec(ijh, na, 2) = 2.d0 * dreal(dens_mtrx(ih, jh, 3))
!            rho%bec(ijh, na, 3) = -2.d0 * dimag(dens_mtrx(ih, jh, 3))
!          endif
!          ! off-diagonal elements have a weight of 2
!          if (ih.ne.jh) then
!            do ispn = 1, nspin
!              rho%bec(ijh, na, ispn) = rho%bec(ijh, na, ispn) * 2.d0
!            enddo
!          endif
!        enddo
!      enddo
!    endif
!  enddo
!enddo
!deallocate(dens_mtrx)
!
!end subroutine get_density_matrix_from_sirius


end module
