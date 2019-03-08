module mod_sirius
use input_parameters, only : use_sirius, sirius_cfg
use sirius
implicit none

! use SIRIUS to compute radial integrals of beta projectors
logical :: use_sirius_radial_integration_beta = .true.
! use SIRIUS to compute Q-operator radial integrals
logical :: use_sirius_radial_integration_q    = .true.
! use SIRIUS to compute radial integrals of Vloc(r)
logical :: use_sirius_radial_integration_vloc = .true.
! use SIRIUS to compute radial integrals of rho_core(r)
logical :: use_sirius_radial_integration_rhoc = .true.
! use SIRIUS to get radial integrals of beta-projectors
logical :: use_sirius_radial_integrals_beta   = .true.
! use SIRIUS to get radial integrals of Q-operator
logical :: use_sirius_radial_integrals_q      = .true.
! use SIRIUS to compute beta projectors
logical :: use_sirius_beta_projectors         = .false.
! use SIRIUS to compute Q-operator
logical :: use_sirius_q_operator              = .false.
! use SIRIUS to solve KS equations
logical :: use_sirius_ks_solver               = .true.
! use SIRIUS to generate density
logical :: use_sirius_density                 = .true.
! use SIRIUS to generate effective potential; WARNING: currently must be always set to .false.
logical :: use_sirius_potential               = .false.
! use SIRIUS to generate density matrix ('bec' thing in QE) WARNING: currently must be set to the value of use_sirius_density
logical :: use_sirius_density_matrix          = .true.
! use SIRIUS to generate D-operator matrix (non-local part of pseudopotential)
logical :: use_sirius_d_operator_matrix       = .true.
! use SIRIUS to compute local part of pseudopotential
logical :: use_sirius_vloc                    = .true.
! use SIRIUS to compute core charge density
logical :: use_sirius_rho_core                = .true.
! use SIRIUS to compute plane-wave coefficients of atomic charge density
logical :: use_sirius_rho_atomic              = .true.
! use SIRIUS to compute forces
logical :: use_sirius_forces                  = .true.
! use SIRIUS to compute stress tensor
logical :: use_sirius_stress                  = .true.
! initialize G-vectors once or at each step of ionic relaxation
logical :: recompute_gvec                     = .false.

! inverse of the reciprocal lattice vectors matrix
real(8) bg_inv(3,3)
! id of the k-point set for ground state calculations
!integer kset_id
! total number of k-points
integer num_kpoints
real(8), allocatable :: kpoints(:,:)
real(8), allocatable :: wkpoints(:)
! phase factors
complex(8), allocatable ::eigts(:,:,:)

type atom_type_t
  ! atom label
  character(len=1, kind=C_CHAR) :: label(100)
  ! nh(iat) in the QE notation
  integer                 :: num_beta_projectors
  ! plane-wave coefficients of Q-operator
  complex(8), allocatable :: qpw(:, :)
end type atom_type_t

type(atom_type_t), allocatable :: atom_type(:)

type(C_PTR) :: sctx = C_NULL_PTR
type(C_PTR) :: gs_handler = C_NULL_PTR
type(C_PTR) :: ks_handler = C_NULL_PTR

contains

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
use control_flags, only : gamma_only, diago_full_acc
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
use ldaU, only : lda_plus_U, Hubbard_J, Hubbard_U, Hubbard_alpha, &
     & Hubbard_beta, is_Hubbard, lda_plus_u_kind, Hubbard_J0, U_projection, Hubbard_l
use esm,       only : esm_local, esm_bc, do_comp_esm
use mp_diag, only : nproc_ortho
use control_flags, only : iverbosity
implicit none
!
integer :: dims(3), i, ia, iat, rank, ierr, ijv, li, lj, mb, nb, j, l,&
     ilast, ir, num_gvec, num_ranks_k, vt(3), iwf, num_kp, nmagd
real(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), tmp
real(8), allocatable :: dion(:, :), qij(:,:,:), vloc(:), wk_tmp(:), xk_tmp(:,:)
integer, allocatable :: nk_loc(:)
integer :: ih, jh, ijh, lmax_beta
INTEGER,EXTERNAL           :: set_hubbard_l,set_hubbard_n
real(8), external :: hubbard_occ

allocate(atom_type(nsp))

do iat = 1, nsp
  atom_type(iat)%label = string(atm(iat))
enddo

! create context of simulation
sctx = sirius_create_context(intra_image_comm)
! set type of caclulation
call sirius_import_parameters(sctx, string('{"parameters" : {"electronic_structure_method" : &
                                             "pseudopotential"}}'))
! import config file
call sirius_import_parameters(sctx, string(trim(adjustl(sirius_cfg))))

! derive the number of magnetic dimensions
nmagd = 0
if (nspin.eq.2) then
  nmagd = 1
endif
if (lspinorb.or.noncolin) then
  nmagd = 3
endif
! set basic parameters
! set |G| cutoff of the dense FFT grid: convert from G^2/2 Rydbergs to |G| in [a.u.^-1]
! set |G+k| cutoff for the wave-functions: onvert from |G+k|^2/2 Rydbergs to |G+k| in [a.u.^-1]
! disable symmetry on SIRIUS side; QE is taking care of the symmetrization
call sirius_set_parameters(sctx, num_bands=nbnd, num_mag_dims=nmagd, gamma_point=bool(gamma_only),&
                          &use_symmetry=bool(.false.), so_correction=bool(lspinorb),&
                          &pw_cutoff=sqrt(ecutrho), gk_cutoff=sqrt(ecutwfc),&
                          &hubbard_correction=bool(lda_plus_U), hubbard_correction_kind=lda_plus_u_kind,&
                          &hubbard_orbitals=string(U_projection), verbosity=min(1, iverbosity))
if (do_comp_esm) then
  call sirius_set_parameters(sctx, esm_bc=esm_bc)
endif

num_ranks_k = nproc_image / npool
i = sqrt(dble(num_ranks_k) + 1d-10)
if (i * i .ne. num_ranks_k) then
  dims(1) = 1
  dims(2) = num_ranks_k
else
  dims(1) = i
  dims(2) = i
endif
call sirius_set_mpi_grid_dims(sctx, 2, dims(1))

if (diago_full_acc) then
  call sirius_set_parameters(sctx, iter_solver_tol_empty=0.d0)
endif

! set lattice vectors of the unit cell (length is in [a.u.])
a1(:) = at(:, 1) * alat
a2(:) = at(:, 2) * alat
a3(:) = at(:, 3) * alat
call sirius_set_lattice_vectors(sctx, a1(1), a2(1), a3(1))

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

  ! add new atom type
   call sirius_add_atom_type(sctx, atom_type(iat)%label, &
        & zn=nint(zv(iat)+0.001d0), &
        & mass=amass(iat), &
        & spin_orbit=bool(upf(iat)%has_so))


  ! set radial grid
  call sirius_set_atom_type_radial_grid(sctx, atom_type(iat)%label, upf(iat)%mesh, upf(iat)%r(1))

  ! set beta-projectors
  do i = 1, upf(iat)%nbeta
    l = upf(iat)%lll(i);
    if (upf(iat)%has_so) then
      if (upf(iat)%jjj(i) .le. upf(iat)%lll(i)) then
        l = - upf(iat)%lll(i)
      endif
    endif
    call sirius_add_atom_type_radial_function(sctx, &
         & atom_type(iat)%label, &
         & string("beta"), &
         & upf(iat)%beta(1, i), &
         & upf(iat)%kbeta(i), &
         & l=l)
  enddo

  ! set the atomic radial functions
  do iwf = 1, upf(iat)%nwfc
    l = upf(iat)%lchi(iwf)
    if (upf(iat)%has_so) then
      if (upf(iat)%jchi(iwf) < l) then
        l = -l
      endif
    endif
    if (associated(upf(iat)%nchi)) then
      i = upf(iat)%nchi(iwf)
    else
      i = -1
    endif
    call sirius_add_atom_type_radial_function(sctx, &
         & atom_type(iat)%label, &
         & string("ps_atomic_wf"), &
         & upf(iat)%chi(1, iwf), &
         & msh(iat), &
         & l=l, &
         & occ=upf(iat)%oc(iwf), &
         & n=i)
  enddo

  if (is_hubbard(iat)) then
     call sirius_set_atom_type_hubbard(sctx, &
          & atom_type(iat)%label, &
          & Hubbard_l(iat), &
          & set_hubbard_n(upf(iat)%psd), &
          & hubbard_occ ( upf(iat)%psd ), &
          & Hubbard_U(iat), &
          & Hubbard_J(1,iat), &
          & Hubbard_alpha(iat), &
          & Hubbard_beta(iat), &
          & Hubbard_J0(iat))
  endif

  allocate(dion(upf(iat)%nbeta, upf(iat)%nbeta))
  ! convert to hartree
  do i = 1, upf(iat)%nbeta
    do j = 1, upf(iat)%nbeta
      dion(i, j) = upf(iat)%dion(i, j) / 2.d0
    end do
  end do
  ! sed d^{ion}_{i,j}
  call sirius_set_atom_type_dion(sctx, atom_type(iat)%label, upf(iat)%nbeta, dion(1, 1))
  deallocate(dion)

  ! get lmax_beta for this atom type
  lmax_beta = -1
  do i = 1, upf(iat)%nbeta
    lmax_beta = max(lmax_beta, upf(iat)%lll(i))
  enddo

  ! set radial function of augmentation charge
  if (upf(iat)%tvanp) then
    !do l = 0, upf(iat)%nqlc - 1
    do l = 0, 2 * lmax_beta
      do i = 0, upf(iat)%nbeta - 1
        do j = i, upf(iat)%nbeta - 1
          ijv = j * (j + 1) / 2 + i + 1
          call sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, string("q_aug"),&
                                                   &upf(iat)%qfuncl(1, ijv, l), upf(iat)%kkbeta,&
                                                   &l=l, idxrf1=i, idxrf2=j)
        enddo
      enddo
    enddo
  endif

  if (upf(iat)%tpawp) then
    do i = 1, upf(iat)%nbeta
      call sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, string("ae_paw_wf"),&
                                               &upf(iat)%aewfc(1,i), upf(iat)%paw%iraug)
      call sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, string("ps_paw_wf"),&
                                               &upf(iat)%pswfc(1,i), upf(iat)%paw%iraug)
    enddo
    call sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, string("ae_paw_core"),&
                                             &upf(iat)%paw%ae_rho_atc(1),upf(iat)%mesh)

    call sirius_set_atom_type_paw(sctx, atom_type(iat)%label, upf(iat)%paw%core_energy / 2,&
                                 &upf(iat)%paw%oc(1), upf(iat)%nbeta)
  endif

  ! set non-linear core correction
  if (use_sirius_rho_core) then
    allocate(vloc(upf(iat)%mesh))
    vloc = 0.d0
    if (associated(upf(iat)%rho_atc)) then
      do i = 1, msh(iat)
        vloc(i) = upf(iat)%rho_atc(i)
      enddo
    endif
    call sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, string("ps_rho_core"),&
                                             &vloc(1), upf(iat)%mesh)
    deallocate(vloc)
  endif

  ! set total charge density of a free atom (to compute initial rho(r))
  call sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, string("ps_rho_total"),&
                                           &upf(iat)%rho_at(1), upf(iat)%mesh)

  ! the hack is done in Modules/readpp.f90
  if (use_sirius_vloc) then
    allocate(vloc(upf(iat)%mesh))
    do i = 1, msh(iat)
      vloc(i) = upf(iat)%vloc(i)
    enddo
    ! convert to Hartree
    vloc = vloc / 2.d0
    ! add a correct tail
    do i = msh(iat) + 1, upf(iat)%mesh
      vloc(i) = -zv(iat) / upf(iat)%r(i)
    enddo
    ! set local part of pseudo-potential
    call sirius_add_atom_type_radial_function(sctx, atom_type(iat)%label, string("vloc"), vloc(1), upf(iat)%mesh)
    deallocate(vloc)
  endif
  !call test_integration(msh(iat), upf(iat)%r, upf(iat)%rab)
  !call test_integration(upf(iat)%mesh, upf(iat)%r, upf(iat)%rab)
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
  !call sirius_reduce_coordinates(v1(1), v2(1), vt(1))
  v2 = v1
  if (noncolin) then
    v1(1) = zv(iat) * starting_magnetization(iat) * sin(angle1(iat)) * cos(angle2(iat))
    v1(2) = zv(iat) * starting_magnetization(iat) * sin(angle1(iat)) * sin(angle2(iat))
    v1(3) = zv(iat) * starting_magnetization(iat) * cos(angle1(iat))
  else
    v1 = 0
    v1(3) = zv(iat) * starting_magnetization(iat)
  endif
  call sirius_add_atom(sctx, atom_type(iat)%label, v2(1), v1(1))
enddo

! initialize global variables/indices/arrays/etc. of the simulation
call sirius_initialize_context(sctx)

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

!if (.true.) then
!  do i = 1, num_kpoints
!    write(*,*)'ik=',i,' kpoint=',matmul(bg_inv,kpoints(:,i))
!  enddo
!endif

!allocate(wk_tmp(nkstot))
!allocate(xk_tmp(3, nkstot))
!! weights of k-points in SIRIUS must sum to one
!do i = 1, nkstot
!  if (nspin.eq.1) then
!    wk_tmp(i) = wk(i) / 2.d0
!  else
!    wk_tmp(i) = wk(i)
!  endif
!  xk_tmp(:,i) = xk(:,i)
!end do
!
!call mpi_bcast(wk_tmp(1),        nkstot, mpi_double, 0, inter_pool_comm, ierr)
!call mpi_bcast(xk_tmp(1, 1), 3 * nkstot, mpi_double, 0, inter_pool_comm, ierr)
!
!! convert to fractional coordinates
!do ik = 1, nkstot
!  xk_tmp(:, ik) = matmul(bg_inv, xk_tmp(:, ik))
!end do

!allocate(nk_loc(0:npool-1))
!nk_loc = 0
!nk_loc(rank) = nks
!call mp_sum(nk_loc, inter_pool_comm)
!if (nspin.eq.2) then
!  nk_loc(:) = nk_loc(:)
!endif

!if (nspin.eq.2) then
!  num_kp = nkstot / 2
!else
!  num_kp = nkstot
!endif

!allocate(xk_tmp(3, num_kpoints))
!do i = 1, num_kpoints
!  xk_tmp(:, i) =  matmul(bg_inv, kpoints(:, i))
!enddo
!deallocate(xk_tmp)

! create k-point set
! WARNING: k-points must be provided in fractional coordinates of the reciprocal lattice
ks_handler = sirius_create_kset(sctx, num_kpoints, kpoints(1, 1), wkpoints(1), bool(.true.))

! create ground-state class
gs_handler = sirius_create_ground_state(ks_handler)

end subroutine setup_sirius


subroutine update_sirius
use cell_base, only : alat, at, bg
use ions_base, only : tau, nat
implicit none
real(8) :: a1(3), a2(3), a3(3), vlat(3, 3), vlat_inv(3, 3), v1(3), v2(3), tmp
integer ia
! set lattice vectors of the unit cell (length is in [a.u.])
a1(:) = at(:, 1) * alat
a2(:) = at(:, 2) * alat
a3(:) = at(:, 3) * alat
call sirius_set_lattice_vectors(sctx, a1(1), a2(1), a3(1))
!
vlat(:, 1) = a1(:)
vlat(:, 2) = a2(:)
vlat(:, 3) = a3(:)
! get the inverse of Bravais lattice vectors
call invert_mtrx(vlat, vlat_inv)
! get the inverse of reciprocal lattice vectors
call invert_mtrx(bg, bg_inv)

do ia = 1, nat
  v1(:) = tau(:, ia) * alat
  ! fractional coordinates
  v1(:) = matmul(vlat_inv, v1)
  call sirius_set_atom_position(sctx, ia, v1(1))
enddo
call sirius_update_ground_state(gs_handler)

end subroutine update_sirius


subroutine test_integration(nr, r, rab)
use mod_spline
implicit none
!
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: rab(nr)
!
real(8) x, x0, x1, exact_val, val1, val2, val3, pi
real(8), allocatable :: f(:)
integer ir
!
pi = 3.1415926535897932385d0
allocate(f(nr))
!-- start of generated code
write(*,*)"testing f(x)=1"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=1
enddo
exact_val=-x0 + x1
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=x"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=x
enddo
exact_val=-x0**2/2. + x1**2/2.
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=x**2"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=x**2
enddo
exact_val=-x0**3/3. + x1**3/3.
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=x**3"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=x**3
enddo
exact_val=-x0**4/4. + x1**4/4.
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=1/(1 + x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=1/(1 + x)
enddo
exact_val=Log((1 + x1)/(1 + x0))
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=(1 + x)**(-2)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=(1 + x)**(-2)
enddo
exact_val=1/(1 + x0) - 1/(1 + x1)
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=Sqrt(x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=Sqrt(x)
enddo
exact_val=(-2*(x0**1.5 - x1**1.5))/3.
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=exp(-x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=exp(-x)
enddo
exact_val=exp(-x0) - exp(-x1)
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=Sin(x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=Sin(x)
enddo
exact_val=cos(x0) - cos(x1)
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=Sin(2*x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=Sin(2*x)
enddo
exact_val=(cos(2*x0) - cos(2*x1))/2.
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=exp(x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=exp(x)
enddo
exact_val=-exp(x0) + exp(x1)
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=exp(-x**2)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=exp(-x**2)
enddo
exact_val=(Sqrt(Pi)*(-Erf(x0) + Erf(x1)))/2.
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=exp(x - x**2)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=exp(x - x**2)
enddo
exact_val=(exp(0.25)*Sqrt(Pi)*(Erf(0.5 - x0) - Erf(0.5 - x1)))/2.
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=x/exp(x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=x/exp(x)
enddo
exact_val=(1 + x0)/exp(x0) - (1 + x1)/exp(x1)
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=x**2/exp(x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=x**2/exp(x)
enddo
exact_val=(2 + x0*(2 + x0))/exp(x0) - (2 + x1*(2 + x1))/exp(x1)
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=Log(1 + x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=Log(1 + x)
enddo
exact_val=x0 - x1 - (1 + x0)*Log(1 + x0) + (1 + x1)*Log(1 + x1)
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
write(*,*)"testing f(x)=Sin(x)/exp(2*x)"
x0=r(1)
x1=r(nr)
do ir=1,nr
  x=r(ir)
  f(ir)=Sin(x)/exp(2*x)
enddo
exact_val=((cos(x0) + 2*Sin(x0))/exp(2*x0) - (cos(x1) + 2*Sin(x1))/exp(2*x1))/5.
call simpson(nr,f,rab,val1)
call integrate(nr,f,r,val2)
call sirius_integrate(0,nr,r(1),f(1),val3)
if (abs(exact_val-val3).lt.abs(exact_val-val2)) write(*,*)'sirius spline ok'
write(*,'("        simpson err: ", G18.10)')abs(exact_val-val1)
write(*,'(" fortran spline err: ", G18.10)')abs(exact_val-val2)
write(*,'("  sirius spline err: ", G18.10)')abs(exact_val-val3)
!-- end of generated code
deallocate(f)
end subroutine


subroutine clear_sirius
use ions_base, only : nsp
implicit none
integer iat

call sirius_free_handler(gs_handler)
call sirius_free_handler(ks_handler)
call sirius_free_handler(sctx)

if (allocated(atom_type)) then
  do iat = 1, nsp
    if (allocated(atom_type(iat)%qpw)) deallocate(atom_type(iat)%qpw)
  enddo
  deallocate(atom_type)
endif

end subroutine


subroutine get_q_operator_from_sirius
use uspp_param, only : upf, nh, nhm
use ions_base,  only : nsp
use gvect,      only : ngm, mill
use uspp,       only : qq_nt
implicit none
integer iat, ih, jh, ijh, i

do iat = 1, nsp
  atom_type(iat)%num_beta_projectors = sirius_get_num_beta_projectors(sctx, atom_type(iat)%label)
  if (nh(iat).ne.atom_type(iat)%num_beta_projectors) then
     write(*,*) nh(iat), atom_type(iat)%num_beta_projectors
     stop 'wrong number of beta projectors'
  endif
  if (upf(iat)%tvanp) then
    i = atom_type(iat)%num_beta_projectors
    if (allocated(atom_type(iat)%qpw)) deallocate(atom_type(iat)%qpw)
    allocate(atom_type(iat)%qpw(ngm, i * (i + 1) / 2))
    ijh = 0
    do ih = 1, atom_type(iat)%num_beta_projectors
      do jh = ih, atom_type(iat)%num_beta_projectors
        ijh = ijh + 1
        call sirius_get_q_operator(sctx, atom_type(iat)%label, ih, jh, ngm, mill(1, 1), atom_type(iat)%qpw(1, ijh))
      enddo
    enddo
  endif
enddo

call get_q_operator_matrix_from_sirius

end subroutine get_q_operator_from_sirius


subroutine get_q_operator_matrix_from_sirius
use uspp_param, only : upf, nh, nhm
use ions_base,  only : nsp, ityp, nat
use uspp,       only : qq_nt, qq_at
implicit none
integer iat, ih, jh, ijh, ia

qq_nt = 0
do iat = 1, nsp
  atom_type(iat)%num_beta_projectors = sirius_get_num_beta_projectors(sctx, atom_type(iat)%label)
  if (nh(iat).ne.atom_type(iat)%num_beta_projectors) then
    stop 'wrong number of beta projectors'
  endif
  if (upf(iat)%tvanp) then
    ijh = 0
    do ih = 1, atom_type(iat)%num_beta_projectors
      do jh = ih, atom_type(iat)%num_beta_projectors
        ijh = ijh + 1
        call sirius_get_q_operator_matrix(sctx, atom_type(iat)%label, qq_nt(1, 1, iat), nhm)
      enddo
    enddo
  endif
enddo
do ia = 1, nat
   qq_at(:, :, ia) = qq_nt(:, :, ityp(ia))
end do

end subroutine get_q_operator_matrix_from_sirius


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
      call sirius_get_band_energies(ks_handler, ik, 0, band_e(1, ik))
    end do
  else
    ! collinear magnetic case
    nk = nkstot / 2
    ! get band energies
    do ik = 1, nk
      call sirius_get_band_energies(ks_handler, ik, 0, band_e(1, ik))
      call sirius_get_band_energies(ks_handler, ik, 1, band_e(1, nk + ik))
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
      call sirius_set_band_occupancies(ks_handler, ik, 0, bnd_occ(1, ik))
    enddo
  else
    nk = nkstot / 2
    do ik = 1, nk
      call sirius_set_band_occupancies(ks_handler, ik, 0, bnd_occ(1, ik))
      call sirius_set_band_occupancies(ks_handler, ik, 1, bnd_occ(1, ik + nk))
    enddo
  endif

  deallocate(bnd_occ)

end subroutine put_band_occupancies_to_sirius


subroutine get_density_matrix_from_sirius
  !
  use scf,        only : rho
  use ions_base,  only : nat, nsp, ityp
  use uspp_param, only : nhm, nh
  use lsda_mod,   only : nspin
  !
  implicit none
  !
  complex(8), allocatable :: dens_mtrx(:,:,:)
  integer iat, na, ijh, ih, jh, ispn
  ! complex density matrix in SIRIUS has at maximum three components
  allocate(dens_mtrx(nhm, nhm, 3))
  do iat = 1, nsp
    do na = 1, nat
      if (ityp(na).eq.iat.and.allocated(rho%bec)) then
        rho%bec(:, na, :) = 0.d0
        call sirius_get_density_matrix(gs_handler, na, dens_mtrx(1, 1, 1), nhm)
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
  call sirius_get_pw_coeffs(gs_handler, string("rho"), rho%of_g(1, 1), ngm, mill(1, 1), intra_bgrp_comm)
  if (nspin.eq.2) then
    call sirius_get_pw_coeffs(gs_handler, string("magz"), rho%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
    !! convert to rho_{up}, rho_{dn}
    !do ig = 1, ngm
    !  z1 = rho%of_g(ig, 1)
    !  z2 = rho%of_g(ig, 2)
    !  rho%of_g(ig, 1) = 0.5 * (z1 + z2)
    !  rho%of_g(ig, 2) = 0.5 * (z1 - z2)
    !enddo
  endif
  if (nspin.eq.4) then
    call sirius_get_pw_coeffs(gs_handler, string("magx"), rho%of_g(1, 2), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_get_pw_coeffs(gs_handler, string("magy"), rho%of_g(1, 3), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_get_pw_coeffs(gs_handler, string("magz"), rho%of_g(1, 4), ngm, mill(1, 1), intra_bgrp_comm)
  endif
  ! get density matrix
  call get_density_matrix_from_sirius
end subroutine get_density_from_sirius


subroutine put_density_to_sirius
  !
  use scf,        only : rho
  use gvect,      only : mill, ngm
  use mp_bands,   only : intra_bgrp_comm
  use lsda_mod,   only : nspin
  use ions_base,  only : nat, nsp, ityp
  use uspp_param, only : nhm, nh
  use sirius
  implicit none
  !
  complex(8), allocatable :: rho_tot(:), mag(:)
  integer iat, ig, ih, jh, ijh, na, ispn
  real(8) :: fact
  !
  !if (nspin.eq.1.or.nspin.eq.4) then
    call sirius_set_pw_coeffs(gs_handler, string("rho"), rho%of_g(1, 1), bool(.true.),&
                             &ngm, mill(1, 1), intra_bgrp_comm)
  !endif

  if (nspin.eq.2) then
    !allocate(rho_tot(ngm))
    !allocate(mag(ngm))
    !do ig = 1, ngm
    !  rho_tot(ig) = rho%of_g(ig, 1) + rho%of_g(ig, 2)
    !  mag(ig) = rho%of_g(ig, 1) - rho%of_g(ig, 2)
    !enddo
    !call sirius_set_pw_coeffs(gs_handler, string("rho"), rho_tot(1), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    !call sirius_set_pw_coeffs(gs_handler, string("magz"), mag(1), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    !deallocate(rho_tot)
    !deallocate(mag)
    call sirius_set_pw_coeffs(gs_handler, string("magz"), rho%of_g(1, 2), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
  endif

  if (nspin.eq.4) then
    call sirius_set_pw_coeffs(gs_handler, string("magx"), rho%of_g(1, 2), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_set_pw_coeffs(gs_handler, string("magy"), rho%of_g(1, 3), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_set_pw_coeffs(gs_handler, string("magz"), rho%of_g(1, 4), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
  endif
end subroutine put_density_to_sirius


subroutine put_density_matrix_to_sirius
  !
  use scf,        only : rho
  use ions_base,  only : nat, nsp, ityp
  use lsda_mod,   only : nspin
  use uspp_param, only : nhm, nh
  use uspp,       only : becsum
  implicit none
  !
  integer iat, na, ih, jh, ijh, ispn
  complex(8), allocatable :: dens_mtrx(:,:,:)
  real(8), allocatable :: dens_mtrx_tmp(:, :, :)
  real(8) fact
  ! set density matrix
  ! complex density matrix in SIRIUS has at maximum three components
  allocate(dens_mtrx_tmp(nhm * (nhm + 1) / 2, nat, nspin))
  !if (allocated(rho%bec)) then
  !  dens_mtrx_tmp = rho%bec
  !else
    dens_mtrx_tmp = becsum
  !endif

  allocate(dens_mtrx(nhm, nhm, 3))
  do iat = 1, nsp
    do na = 1, nat
      if (ityp(na).eq.iat) then
        dens_mtrx = (0.d0, 0.d0)
        ijh = 0
        do ih = 1, nh(iat)
          do jh = ih, nh(iat)
            ijh = ijh + 1
            ! off-diagonal elements have a weight of 2
            if (ih.ne.jh) then
              fact = 0.5d0
            else
              fact = 1.d0
            endif
            if (nspin.le.2) then
              do ispn = 1, nspin
                dens_mtrx(ih, jh, ispn) = fact * dens_mtrx_tmp(ijh, na, ispn)
                dens_mtrx(jh, ih, ispn) = fact * dens_mtrx_tmp(ijh, na, ispn)
              enddo
            endif
            if (nspin.eq.4) then
              ! 0.5 * (rho + mz)
              dens_mtrx(ih, jh, 1) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) + dens_mtrx_tmp(ijh, na, 4))
              dens_mtrx(jh, ih, 1) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) + dens_mtrx_tmp(ijh, na, 4))
              ! 0.5 * (rho - mz)
              dens_mtrx(ih, jh, 2) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) - dens_mtrx_tmp(ijh, na, 4))
              dens_mtrx(jh, ih, 2) = fact * 0.5 * (dens_mtrx_tmp(ijh, na, 1) - dens_mtrx_tmp(ijh, na, 4))
              ! 0.5 * (mx - I * my)
              dens_mtrx(ih, jh, 3) = fact * 0.5 * dcmplx(dens_mtrx_tmp(ijh, na, 2), -dens_mtrx_tmp(ijh, na, 3))
              dens_mtrx(jh, ih, 3) = fact * 0.5 * dcmplx(dens_mtrx_tmp(ijh, na, 2), -dens_mtrx_tmp(ijh, na, 3))
            endif
          enddo
        enddo
        call sirius_set_density_matrix(gs_handler, na, dens_mtrx(1, 1, 1), nhm)
      endif
    enddo
  enddo
  deallocate(dens_mtrx)
  deallocate(dens_mtrx_tmp)
end subroutine put_density_matrix_to_sirius


subroutine put_d_matrix_to_sirius
use uspp_param,           only : nhm
use ions_base,            only : nat
use lsda_mod,             only : nspin
use uspp,                 only : deeq
implicit none
real(8), allocatable :: deeq_tmp(:,:)
integer ia, is
allocate(deeq_tmp(nhm, nhm))
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
    call sirius_set_d_operator_matrix(sctx, ia, is, deeq_tmp(1, 1), nhm)
  enddo
enddo
deallocate(deeq_tmp)
end subroutine put_d_matrix_to_sirius


subroutine get_d_matrix_from_sirius
use uspp_param,           only : nhm
use ions_base,            only : nat
use lsda_mod,             only : nspin
use uspp,                 only : deeq
implicit none
real(8) d1, d2
integer ia, is, i, j
! get D-operator matrix
do ia = 1, nat
  do is = 1, nspin
    call sirius_get_d_operator_matrix(sctx, ia, is, deeq(1, 1, ia, is), nhm)
  enddo
  if (nspin.eq.2) then
    do i = 1, nhm
      do j = 1, nhm
        d1 = deeq(i, j, ia, 1)
        d2 = deeq(i, j, ia, 2)
        deeq(i, j, ia, 1) = d1 + d2
        deeq(i, j, ia, 2) = d1 - d2
      enddo
    enddo
  endif
  ! convert to Ry
  deeq(:, :, ia, :) = deeq(:, :, ia, :) * 2
enddo
end subroutine get_d_matrix_from_sirius

subroutine put_potential_to_sirius
  use scf,                  only : v, vltot
  use gvect,                only : mill, ngm
  use mp_bands,             only : intra_bgrp_comm
  use lsda_mod,             only : nspin
  use noncollin_module,     only : nspin_mag
  USE wavefunctions,        ONLY : psic

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
    call sirius_set_pw_coeffs(gs_handler, string("veff"), v%of_g(1, 1), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
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
    call sirius_set_pw_coeffs(gs_handler, string("veff"), v%of_g(1, 1), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
    call sirius_set_pw_coeffs(gs_handler, string("bz"),   v%of_g(1, 2), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
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

      call sirius_set_pw_coeffs(gs_handler, string(label), v%of_g(1, is), bool(.true.), ngm, mill(1, 1), intra_bgrp_comm)
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
!  call sirius_set_pw_coeffs(string("vxc"), vxcg(1), ngm, mill(1, 1), intra_bgrp_comm)
!  deallocate(vxcg)
!
!  ! update D-operator matrix
!  !call sirius_generate_d_operator_matrix()
!  !if (okpaw) then
!    allocate(deeq_tmp(nhm, nhm))
!  !  !! get D-operator matrix
!  !  !do ia = 1, nat
!  !  !  do is = 1, nspin
!  !  !    call sirius_get_d_operator_matrix(ia, is, deeq(1, 1, ia, is), nhm)
!  !  !  enddo
!  !  !  if (nspin.eq.2) then
!  !  !    do i = 1, nhm
!  !  !      do j = 1, nhm
!  !  !        d1 = deeq(i, j, ia, 1)
!  !  !        d2 = deeq(i, j, ia, 2)
!  !  !        deeq(i, j, ia, 1) = d1 + d2
!  !  !        deeq(i, j, ia, 2) = d1 - d2
!  !  !      enddo
!  !  !    enddo
!  !  !  endif
!  !  !  ! convert to Ry
!  !  !  deeq(:, :, ia, :) = deeq(:, :, ia, :) * 2
!  !  !enddo
!  !  !call add_paw_to_deeq(deeq)
!    do ia = 1, nat
!      do is = 1, nspin
!        if (nspin.eq.2.and.is.eq.1) then
!          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) + deeq(:, :, ia, 2)) / 2 ! convert to Ha
!        endif
!        if (nspin.eq.2.and.is.eq.2) then
!          deeq_tmp(:, :) = 0.5 * (deeq(:, :, ia, 1) - deeq(:, :, ia, 2)) / 2 ! convert to Ha
!        endif
!        if (nspin.eq.1.or.nspin.eq.4) then
!          deeq_tmp(:, :) = deeq(:, :, ia, is) / 2 ! convert to Ha
!        endif
!        call sirius_set_d_operator_matrix(ia, is, deeq_tmp(1, 1), nhm)
!      enddo
!    enddo
!    deallocate(deeq_tmp)
!  !endif

end subroutine put_potential_to_sirius

subroutine put_q_operator_matrix_to_sirius
use uspp,       only : qq_nt
use uspp_param, only : upf, nhm
use ions_base,  only : nsp, atm
implicit none
integer iat

do iat = 1, nsp
  if (upf(iat)%tvanp) then
    call sirius_set_q_operator_matrix(sctx,string(atm(iat)), qq_nt(1, 1, iat), nhm)
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
use wavefunctions, only : evc
use lsda_mod, only : isk, lsda
use mp_pools, only : inter_pool_comm
use parallel_include
implicit none
integer, external :: global_kpoint_index
integer, allocatable :: gvl(:,:)
integer ig, ik, ik_, i, j, ispn, rank, ierr, nksmax
complex(8) z1

! rank of communicator that distributes k-points
call mpi_comm_rank(inter_pool_comm, rank, ierr)
call mpi_allreduce(nks, nksmax, 1, MPI_INTEGER, MPI_MAX, inter_pool_comm, ierr)

!CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )

allocate(gvl(3, npwx))
do ik = 1, nksmax
  if (ik.le.nks) then
    do ig = 1, ngk(ik)
      gvl(:,ig) = mill(:, igk_k(ig, ik))
    enddo
    !
    ik_ = global_kpoint_index(nkstot, ik)
    ispn = isk(ik)
    if (lsda.and.ispn.eq.2) then
      ik_ = ik_ - nkstot / 2
    endif
    call sirius_get_wave_functions(ks_handler, ik_, ispn, ngk(ik), gvl(1, 1), evc(1, 1), npwx, npol)
    if (nks > 1 .or. lelfield) then
      call save_buffer ( evc, nwordwfc, iunwfc, ik )
    endif
  else
    call sirius_get_wave_functions(ks_handler, -1, -1, -1, -1, z1, -1, -1)
  endif
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
!  call sirius_set_pw_coeffs(string("vloc"), vg(1), ngm, mill(1, 1), intra_bgrp_comm)
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
!    call sirius_get_pw_coeffs(string("rhoc"), rhog_core(1), ngm, mill(1, 1), intra_bgrp_comm)
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
!  call sirius_get_pw_coeffs_real(string(atm(nt)), string("vloc"), tmp(1), ngm, mill(1, 1), intra_bgrp_comm)
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
!  call sirius_get_pw_coeffs(string("vloc"), vpw(1), ngm, mill(1, 1), intra_bgrp_comm)
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
  !    call sirius_add_xc_functional(string("XC_LDA_X"))
  !  case default
  !    stop ("interface for this exchange functional is not implemented")
  !  end select
  !endif

  !if (get_iexch().ne.0.and.get_igcx().ne.0) then
  !  select case(get_igcx())
  !  case(0)
  !  case(2)
  !    call sirius_add_xc_functional(string("XC_GGA_X_PW91"))
  !  case(3)
  !    call sirius_add_xc_functional(string("XC_GGA_X_PBE"))
  !  case default
  !    write(*,*)get_igcx()
  !    stop ("interface for this gradient exchange functional is not implemented")
  !  end select
  !endif

  !if (get_icorr().ne.0.and.get_igcc().eq.0) then
  !  select case(get_icorr())
  !  case(0)
  !  case(1)
  !    call sirius_add_xc_functional(string("XC_LDA_C_PZ"))
  !  case(4)
  !    call sirius_add_xc_functional(string("XC_LDA_C_PW"))
  !  case default
  !    stop ("interface for this correlation functional is not implemented")
  !  end select
  !endif

  !if (get_icorr().ne.0.and.get_igcc().ne.0) then
  !  select case(get_igcc())
  !  case(0)
  !  case(2)
  !    call sirius_add_xc_functional(string("XC_GGA_C_PW91"))
  !  case(4)
  !    call sirius_add_xc_functional(string("XC_GGA_C_PBE"))
  !  case default
  !    stop ("interface for this gradient correlation functional is not implemented")
  !  end select
  !endif
end subroutine put_xc_functional_to_sirius

subroutine write_json()
use ener
use force_mod
use ions_base
implicit none
integer i

open(200, file="output.json", action="write", form="formatted")
write(200,'("{")')
write(200,'("  ""energy"": {")')
write(200,'("    ""total"": ", G18.10)')etot
write(200, '("  },")')
write(200,'("  ""stress"": {")')
write(200,'("    ""total"": [")')
write(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') sigma(1, 1), sigma(1, 2), sigma(1, 3)
write(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') sigma(2, 1), sigma(2, 2), sigma(2, 3)
write(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "]")')  sigma(3, 1), sigma(3, 2), sigma(3, 3)
write(200,'("    ]")')
write(200,'("  },")')
write(200,'("  ""force"": {")')
write(200,'("    ""total"": [")')
do i = 1, nat
  if (i.eq.nat) then
    write(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "]")') force(1, i), force(2, i), force(3, i)
  else
    write(200,'("        [", G18.10, ",", G18.10, ",", G18.10, "],")') force(1, i), force(2, i), force(3, i)
  endif
enddo
write(200,'("    ]")')
write(200,'("  }")')
write(200,'("}")')
close(200)

end subroutine

function idx_m_qe(m) result(m1)
  implicit none
  integer :: m
  integer :: m1

  if (m .gt. 0) then
     m1 = 2 * m - 1
  else
     m1 = -2 * m
  endif
end function idx_m_qe

subroutine qe_to_sirius_real(ns, ns_sirius)
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin

  complex(8), intent(out) :: ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  REAL(8), intent(in) :: ns(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  integer :: m1, m2, mm2, mm1, is, nt, na
  ns_sirius(:, :, :, :) = 0.0
  do na = 1, nat
     nt = ityp (na)
     if (is_hubbard(nt)) then
        do m1 = -Hubbard_l(nt), Hubbard_l(nt)
           mm1 = idx_m_qe(m1)
           do m2 = -Hubbard_l(nt), Hubbard_l(nt)
              mm2 = idx_m_qe(m2)
              DO is = 1, nspin
                 ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, is, na) = ns(mm1 + 1, mm2 + 1, is, na)
              enddo
           enddo
        enddo
     endif
  enddo
end subroutine qe_to_sirius_real

subroutine qe_to_sirius_complex(ns, ns_sirius)
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin

  complex(8), intent(out) :: ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  complex(8), intent(in) :: ns(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  integer :: m1, m2, mm2, mm1, is, nt, na

  ns_sirius(:, :, :, :) = 0.0
  do na = 1, nat
     nt = ityp (na)
     if (is_hubbard(nt)) then
        do m1 = -Hubbard_l(nt), Hubbard_l(nt)
           mm1 = idx_m_qe(m1)
           do m2 = -Hubbard_l(nt), Hubbard_l(nt)
              mm2 = idx_m_qe(m2)
              ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 1, na) = ns(mm1 + 1, mm2 + 1, 1, na)
              ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 2, na) = ns(mm1 + 1, mm2 + 1, 4, na)
              ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 3, na) = ns(mm1 + 1, mm2 + 1, 2, na)
              ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 4, na) = ns(mm1 + 1, mm2 + 1, 3, na)
           enddo
        enddo
     endif
  enddo
end subroutine qe_to_sirius_complex

subroutine sirius_to_qe_real(ns_sirius, ns)
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin
  integer :: m1, m2, mm2, mm1, is, nt, na

  complex(8), intent(in) :: ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  real(8), intent(out) :: ns(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)

  ns(:, :, :, :) = 0.0
  do na = 1, nat
     nt = ityp (na)
     if (is_hubbard(nt)) then
        do m1 = -Hubbard_l(nt), Hubbard_l(nt)
           mm1 = idx_m_qe(m1)
           do m2 = -Hubbard_l(nt), Hubbard_l(nt)
              mm2 = idx_m_qe(m2)
              DO is = 1, nspin
                 ns(mm1 + 1, mm2 + 1, is, na) = REAL(ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, is, na))
              enddo
           enddo
        enddo
     endif
  enddo
end subroutine sirius_to_qe_real

subroutine sirius_to_qe_complex(ns_sirius, ns)
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin

  complex(8), intent(in) :: ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat)
  complex(8), intent(out) :: ns(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat)
  integer :: m1, m2, mm2, mm1, is, nt, na

  ns(:, :, :, :) = 0.0
  do na = 1, nat
     nt = ityp (na)
     if (is_hubbard(nt)) then
        do m1 = -Hubbard_l(nt), Hubbard_l(nt)
           mm1 = idx_m_qe(m1)
           do m2 = -Hubbard_l(nt), Hubbard_l(nt)
              mm2 = idx_m_qe(m2)
              ns(mm1 + 1, mm2 + 1, 1, na) = ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 1, na)
              ns(mm1 + 1, mm2 + 1, 4, na) = ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 2, na)
              ns(mm1 + 1, mm2 + 1, 2, na) = ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 3, na)
              ns(mm1 + 1, mm2 + 1, 3, na) = ns_sirius(m1 + Hubbard_l(nt) + 1, m2 + Hubbard_l(nt) + 1, 4, na)
           enddo
        enddo
     endif
  enddo
end subroutine sirius_to_qe_complex

subroutine qe_sirius_set_hubbard_occupancy(rho)
  USE scf,              ONLY : scf_type
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE lsda_mod,             ONLY : nspin

  IMPLICIT NONE

  TYPE(scf_type), INTENT(IN) :: rho  ! the valence charge
  complex(8), allocatable :: ns_sirius(:, :, :, :)


  if (noncolin) then
     allocate(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat))
     call qe_to_sirius_complex(rho%ns_nc(1, 1, 1, 1), ns_sirius(1, 1, 1, 1))
  else
     allocate(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat))
     call qe_to_sirius_real(rho%ns(1, 1, 1, 1), ns_sirius(1, 1, 1, 1))
  endif
  call sirius_set_hubbard_occupancies(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)
  deallocate(ns_sirius)
end subroutine qe_sirius_set_hubbard_occupancy

subroutine qe_sirius_get_hubbard_occupancy(rho)
  USE scf,              ONLY : scf_type
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  USE lsda_mod,             ONLY : nspin
  IMPLICIT NONE

  TYPE(scf_type), INTENT(out) :: rho  ! the valence charge
  complex(8), allocatable :: ns_sirius(:, :, :, :)


  if (noncolin) then
     allocate(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat))
     call sirius_get_hubbard_occupancies(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)
     call sirius_to_qe_complex(ns_sirius(1, 1, 1, 1), rho%ns_nc(1, 1, 1, 1))
  else
     allocate(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 2, nat))
     call sirius_get_hubbard_occupancies(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)
     call sirius_to_qe_real(ns_sirius(1, 1, 1, 1), rho%ns(1, 1, 1, 1))
  endif

  deallocate(ns_sirius)
end subroutine qe_sirius_get_hubbard_occupancy

subroutine qe_sirius_set_hubbard_potential(v)
  USE scf,              ONLY : scf_type
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  IMPLICIT NONE
  TYPE(scf_type), INTENT(IN) :: v  ! the valence charge
  complex(8), allocatable :: ns_sirius(:, :, :, :)

  if (noncolin) then
     allocate(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat))
     call qe_to_sirius_complex(v%ns_nc(1, 1, 1, 1), ns_sirius(1, 1, 1, 1))
  else
     allocate(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat))
     call qe_to_sirius_real(v%ns(1, 1, 1, 1), ns_sirius(1, 1, 1, 1))
  endif
  ns_sirius(:,:,:,:) = 0.5 * ns_sirius(:,:,:,:)
  call sirius_set_hubbard_potential(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)

  deallocate(ns_sirius)
end subroutine qe_sirius_set_hubbard_potential

subroutine qe_sirius_get_hubbard_potential(v)
  USE scf,              ONLY : scf_type
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
       Hubbard_J, Hubbard_alpha, lda_plus_u_kind,&
       Hubbard_J0, Hubbard_beta,  is_hubbard
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module, ONLY : noncolin, nspin_lsda
  IMPLICIT NONE

  complex(8), allocatable :: ns_sirius(:, :, :, :)
  TYPE(scf_type), INTENT(out) :: v  ! the valence charge

  call sirius_get_hubbard_potential(gs_handler, ns_sirius(1, 1, 1, 1), 2 * hubbard_lmax + 1)

  if (noncolin) then
     allocate(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, 4, nat))
     ns_sirius(:,:,:,:) = 2.0 * ns_sirius(:,:,:,:)
     call sirius_to_qe_complex(ns_sirius(1, 1, 1, 1), v%ns_nc(1, 1, 1, 1))

  else
     allocate(ns_sirius(2 * Hubbard_lmax + 1, 2 * Hubbard_lmax + 1, nspin, nat))
     ns_sirius(:,:,:,:) = 2.0 * ns_sirius(:,:,:,:)
     call sirius_to_qe_real(ns_sirius(1, 1, 1, 1), v%ns(1, 1, 1, 1))
  endif

  deallocate(ns_sirius)
end subroutine qe_sirius_get_hubbard_potential


end module mod_sirius
