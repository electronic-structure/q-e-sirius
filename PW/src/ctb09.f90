!
! Copyright (C) 2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE ctb09
  !--------------------------------------------------------------------------
  !! Mixing parameter (c) of the Tran-Blaha (TB09) meta-GGA exchange
  !! potential [F. Tran and P. Blaha, PRL 102, 226401 (2009)], computed from
  !! the density at every evaluation rather than left at the Libxc default.
  !
  !! c is a function of \(g=\langle|\nabla\rho|/\rho\rangle\), an average
  !! dominated by the core region, where a pseudo density has no structure to
  !! average over. The all-electron correction to it is a property of each
  !! pseudopotential alone, so it is computed once, when they are read, and
  !! cached here for tb09_c to add to the plane-wave average. For each
  !! species this module builds the all-electron core charge of the neutral
  !! atom and stores the atom-centred difference
  !
  !! \[ A = \int 4\pi r^2 \left( \frac{|\nabla\rho_\text{AE}|}
  !!    {\rho_\text{AE}} - \frac{|\nabla\rho_\text{PS}|}
  !!    {\rho_\text{PS}} \right) dr \]
  !
  !! where \(\rho_\text{PS}\) is the atomic pseudo valence density plus the
  !! core charge of the nonlinear core correction, that is, exactly the
  !! density the plane-wave average is taken over, and \(\rho_\text{AE}\)
  !! replaces that core charge by the true core density. The valence density
  !! is common to the two and cancels pointwise, which is why the integrand,
  !! and not the two integrals separately, is what gets integrated: at the
  !! outer edge the two integrals are of order \(10^3\) and their difference
  !! of order \(10\).
  !
  !! The core charge comes from a spherical all-electron solve of the frozen
  !! ground-state atom, in LDA, on its own logarithmic grid.
  !
  USE kinds,        ONLY : DP
  USE radial_grids, ONLY : radial_grid_type
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: tb09_c, ae_correction, ae_deallocate
  !
  REAL(DP), ALLOCATABLE :: ae_grad_int(:)
  !! \(A\), the atom-centred correction to \(g\), per species (bohr^2)
  LOGICAL, ALLOCATABLE :: ae_valid(:)
  !! TRUE if \(A\) could be computed for that species
  CHARACTER(LEN=64), ALLOCATABLE :: ae_note(:)
  !! provenance, or the reason why the species was skipped
  !
  REAL(DP), PARAMETER :: rcut_ae = 10.0_DP
  !! the integral is truncated here, or at the end of the pseudopotential
  !! mesh if that comes first: outside, AE and PS densities coincide
  REAL(DP), PARAMETER :: rmax_ae = 100.0_DP
  REAL(DP), PARAMETER :: dx_ae = 0.008_DP
  !! grid parameters of the all-electron solve, the ld1.x defaults
  REAL(DP), PARAMETER :: density_floor = 1.0E-12_DP
  REAL(DP), PARAMETER :: charge_eps = 1.0E-6_DP
  !
  INTEGER, PARAMETER :: nshell = 19
  !! shells up to 7p, enough for every element with a pseudopotential
  !
  ! ... shells in spectroscopic order, by n and then by l. This is the order
  ! ... in which they are spatially nested, and therefore the order in which
  ! ... a pseudopotential freezes them into its core. It is NOT the order of
  ! ... increasing binding energy: in tungsten 4f lies 3.4 Ry ABOVE 5s, yet a
  ! ... 14-electron pseudopotential puts 4f in the core and 5s in the valence
  !
  INTEGER, PARAMETER :: n_spec(nshell) = &
     (/ 1, 2,2, 3,3,3, 4,4,4,4, 5,5,5,5, 6,6,6, 7,7 /)
  INTEGER, PARAMETER :: l_spec(nshell) = &
     (/ 0, 0,1, 0,1,2, 0,1,2,3, 0,1,2,3, 0,1,2, 0,1 /)
  !
  ! ... shells in Madelung order, by n+l and then by n: the filling order of
  ! ... the neutral ground state. The handful of elements that violate it,
  ! ... chromium and gold among them, do so by moving one electron between
  ! ... their outermost s and d or f shells, which no pseudopotential treats
  ! ... as core, so the core charge this module needs is unaffected
  !
  INTEGER, PARAMETER :: n_mad(nshell) = &
     (/ 1, 2, 2,3, 3,4, 3,4,5, 4,5,6, 4,5,6,7, 5,6,7 /)
  INTEGER, PARAMETER :: l_mad(nshell) = &
     (/ 0, 0, 1,0, 1,0, 2,1,0, 2,1,0, 3,2,1,0, 3,2,1 /)
  !
CONTAINS
!
!----------------------------------------------------------------------------
FUNCTION tb09_c( nnr, nspin, rho, rho_core, grho ) RESULT( c )
  !--------------------------------------------------------------------------
  !! Compute the mixing parameter as a functional of the density.
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE fft_base,  ONLY : dfftp
  USE mp,        ONLY : mp_sum
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nnr
  !! number of points of the dense FFT grid, local to this process
  INTEGER, INTENT(IN) :: nspin
  !! number of spin components
  REAL(DP), INTENT(IN) :: rho(nnr)
  !! total valence charge density, in real space
  REAL(DP), INTENT(IN) :: rho_core(nnr)
  !! core charge of the nonlinear core correction, in real space
  REAL(DP), INTENT(IN) :: grho(3,nnr,nspin)
  !! gradient of the spin densities, core charge included
  REAL(DP) :: c
  !! the mixing parameter
  !
  ! ... local variables
  !
  INTEGER  :: ia, k
  REAL(DP) :: gx, gy, gz, grho2, gbar, gcore, rhok
  REAL(DP), PARAMETER :: alpha = -0.012_DP, beta = 1.023_DP
  !! intercept and slope, in bohr^(1/2), of the c(g) fit
  REAL(DP), PARAMETER :: rho_min = 1.0E-10_DP
  !! below this the point is skipped rather than floored: in vacuum both
  !! the gradient and the density vanish, the ratio is numerical noise,
  !! and dividing by a floor turns that noise into a large contribution.
  !! Same threshold and same treatment as ABINIT, m_xc_tb09.F90
  !
  gbar = 0._DP
  !
  ! ... numerator and denominator must contain the same density the
  ! ... functional is evaluated with, that is, the valence density plus
  ! ... the core charge of the nonlinear core correction
  !
  IF ( nspin == 1 ) THEN
    !$acc parallel loop reduction(+:gbar) present(rho,rho_core,grho) &
    !$acc& private(gx,gy,gz,grho2,rhok)
    DO k = 1, nnr
      gx = grho(1,k,1)
      gy = grho(2,k,1)
      gz = grho(3,k,1)
      grho2 = gx*gx + gy*gy + gz*gz
      rhok = ABS( rho(k) + rho_core(k) )
      !
      IF ( rhok > rho_min ) gbar = gbar + SQRT(grho2) / rhok
      !
    ENDDO
  ELSE
    !$acc parallel loop reduction(+:gbar) present(rho,rho_core,grho) &
    !$acc& private(gx,gy,gz,grho2,rhok)
    DO k = 1, nnr
      gx = grho(1,k,1) + grho(1,k,2)
      gy = grho(2,k,1) + grho(2,k,2)
      gz = grho(3,k,1) + grho(3,k,2)
      grho2 = gx*gx + gy*gy + gz*gz
      rhok = ABS( rho(k) + rho_core(k) )
      !
      IF ( rhok > rho_min ) gbar = gbar + SQRT(grho2) / rhok
      !
    ENDDO
  ENDIF
  !
  CALL mp_sum( gbar, intra_bgrp_comm )
  !
  gbar = gbar / DBLE( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
  !
  ! ... alpha and beta were fitted to all-electron densities. The pseudo
  ! ... density is missing the core contribution to the average; add back
  ! ... the radial integral computed once per species by ae_correction,
  ! ... when the pseudopotentials were read. It is a
  ! ... volume integral, so unlike gbar itself it does not cancel omega,
  ! ... which is read here rather than saved so that vc-relax follows the
  ! ... cell. All-or-nothing: a partial correction would mix an
  ! ... all-electron core on some sites with a pseudo core on others
  !
  gcore = 0._DP
  IF ( ALLOCATED(ae_valid) ) THEN
    IF ( ALL(ae_valid(1:ntyp)) ) THEN
      DO ia = 1, nat
        gcore = gcore + ae_grad_int( ityp(ia) )
      ENDDO
      gcore = gcore / omega
    ENDIF
  ENDIF
  !
  c = alpha + beta * SQRT( MAX(gbar+gcore,0._DP) )
  !
  ! ... clamp to the range the functional was fitted in, as ABINIT does in
  ! ... m_xc_tb09.F90. c = 1 is the original Becke-Johnson potential, so the
  ! ... lower bound falls back to it rather than to something unphysical:
  ! ... the Becke-Johnson term carries a factor (3c-2), which changes sign
  ! ... below c = 2/3. The upper bound guards against a g inflated by a
  ! ... pathological density
  !
  c = MIN( MAX( c, 1._DP ), 2._DP )
  !
  IF ( gcore /= 0._DP ) THEN
    WRITE( stdout, '(/5X,"TB09: <|grad rho|/rho> = ",F12.8," + ",F12.8, &
                    &" (AE core), c = ",F12.8)' ) gbar, gcore, c
  ELSE
    WRITE( stdout, '(/5X,"TB09: <|grad rho|/rho> = ",F12.8,", c = ",F12.8)' ) &
           gbar, c
  ENDIF
  !
  RETURN
  !
END FUNCTION tb09_c
!
!----------------------------------------------------------------------------
SUBROUTINE ae_correction( )
  !--------------------------------------------------------------------------
  !! Fill the cache, once per run, and report it. Call it behind dft_is_ctb09:
  !! it is unconditional work, and no other functional needs the result.
  !
  USE io_global,    ONLY : stdout, ionode
  USE ions_base,    ONLY : ntyp => nsp
  USE uspp_param,   ONLY : upf
  !
  IMPLICIT NONE
  !
  INTEGER :: nt
  !
  ! ... meta-GGA is norm-conserving only, see setup.f90, which is where the
  ! ... error for this input belongs. Leave the cache unallocated and let it
  ! ... speak there rather than half-way through reading pseudopotentials
  !
  IF ( ANY(upf(1:ntyp)%tvanp) .OR. ANY(upf(1:ntyp)%tpawp) ) RETURN
  !
  CALL allocate_cache( ntyp )
  !
  DO nt = 1, ntyp
    CALL species( nt, upf(nt) )
  ENDDO
  !
  IF ( .NOT. ionode ) RETURN
  !
  WRITE( stdout, '(/5X,"CTB09: all-electron correction to <|grad rho|/rho>")' )
  DO nt = 1, ntyp
    IF ( ae_valid(nt) ) THEN
      WRITE( stdout, '(5X,I3,2X,A2,F14.8,"  bohr^2  ",A)' ) &
        nt, upf(nt)%psd, ae_grad_int(nt), TRIM(ae_note(nt))
    ELSE
      WRITE( stdout, '(5X,I3,2X,A2,"   unavailable    ",A)' ) &
        nt, upf(nt)%psd, TRIM(ae_note(nt))
    ENDIF
  ENDDO
  IF ( .NOT. ALL(ae_valid(1:ntyp)) ) WRITE( stdout, &
     '(5X,"CTB09: correction NOT applied - c is computed from the ", &
      &"pseudo density alone")' )
  !
END SUBROUTINE ae_correction
!
!----------------------------------------------------------------------------
SUBROUTINE allocate_cache( ntyp )
  !--------------------------------------------------------------------------
  !! Reset the per-species cache.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ntyp
  !
  CALL ae_deallocate()
  !
  ALLOCATE( ae_grad_int(ntyp), ae_valid(ntyp), ae_note(ntyp) )
  !
  ae_grad_int(:) = 0._DP
  ae_valid(:) = .FALSE.
  ae_note(:) = 'not computed'
  !
END SUBROUTINE allocate_cache
!
!----------------------------------------------------------------------------
SUBROUTINE ae_deallocate()
  !--------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  IF ( ALLOCATED(ae_grad_int) ) DEALLOCATE( ae_grad_int )
  IF ( ALLOCATED(ae_valid)    ) DEALLOCATE( ae_valid )
  IF ( ALLOCATED(ae_note)     ) DEALLOCATE( ae_note )
  !
END SUBROUTINE ae_deallocate
!
!----------------------------------------------------------------------------
SUBROUTINE species( nt, upf_in )
  !--------------------------------------------------------------------------
  !! Compute and cache \(A\) for one species.
  !
  USE pseudo_types, ONLY : pseudo_upf
  USE radial_grids, ONLY : do_mesh, deallocate_radial_grid, nullify_radial_grid
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nt
  TYPE(pseudo_upf), INTENT(IN) :: upf_in
  !
  INTEGER, EXTERNAL :: atomic_number
  !
  TYPE(radial_grid_type) :: grid
  CHARACTER(LEN=16) :: treatment
  LOGICAL  :: core(nshell)
  INTEGER  :: ierr, mesh, ncore, nwf, rel
  INTEGER  :: ll(nshell), nn(nshell)
  REAL(DP) :: enl(nshell), oc(nshell)
  REAL(DP) :: xmin, zcore, zed
  REAL(DP), ALLOCATABLE :: psi(:,:), rho_core_ae(:), rho_ps(:)
  !
  ae_grad_int(nt) = 0._DP
  ae_valid(nt) = .FALSE.
  !
  ! ... zmesh is meaningful only for pseudopotentials on a logarithmic mesh;
  ! ... those on a linear one, ONCVPSP among them, leave it at zero
  !
  zed = upf_in%zmesh
  IF ( zed < 1._DP ) zed = DBLE( atomic_number( upf_in%psd ) )
  IF ( zed < 1._DP .OR. zed > 103._DP ) THEN
    ae_note(nt) = 'element not recognized'
    RETURN
  ENDIF
  !
  zcore = zed - upf_in%zp
  IF ( zcore < 0.5_DP ) THEN
    ! ... no core to correct: an all-electron valence, so A vanishes exactly
    ae_valid(nt) = .TRUE.
    ae_note(nt) = 'no core electrons'
    RETURN
  ENDIF
  !
  ! ... which shells the pseudopotential froze into its core. This follows
  ! ... from the core charge alone and needs neither the shell labels of the
  ! ... file, which ONCVPSP does not write, nor an ordering by eigenvalue,
  ! ... which gets the 5d elements wrong
  !
  CALL core_shells( zcore, ncore, ierr )
  IF ( ierr /= 0 ) THEN
    WRITE( ae_note(nt), '("core charge ",F6.2," is not a closed shell")' ) zcore
    RETURN
  ENDIF
  !
  ! ... 'full' is treated scalar-relativistically: the Dirac solver of
  ! ... ld1.x is not available here, and the spin-orbit splitting of a
  ! ... closed core shell cancels in the spherical average anyway.
  ! ... The UPF v1 reader never fills upf%rel, so without a fallback every
  ! ... v1 pseudopotential would be solved nonrelativistically; for a heavy
  ! ... element that is a 10 percent error in the core average. Guess it the
  ! ... way upfconv does for the same files
  !
  IF ( LEN_TRIM(upf_in%rel) == 0 ) THEN
    rel = 0
    IF ( upf_in%has_so .OR. zed > 18._DP ) rel = 1
  ELSE
    rel = 0
    IF ( TRIM(upf_in%rel) == 'scalar' .OR. TRIM(upf_in%rel) == 'full' ) rel = 1
  ENDIF
  !
  treatment = 'nonrelativistic'
  IF ( rel > 0 ) treatment = 'scalar rel.'
  !
  CALL neutral_config( zed, ncore, nwf, nn, ll, oc, core, ierr )
  IF ( ierr /= 0 ) THEN
    ae_note(nt) = 'core shell is not full in the neutral atom'
    RETURN
  ENDIF
  !
  xmin = -7.0_DP
  IF ( rel > 0 ) xmin = -8.0_DP
  !
  ! ... radial_grid_type has no default initialization
  CALL nullify_radial_grid( grid )
  CALL do_mesh( rmax_ae, zed, xmin, dx_ae, 0, grid )
  mesh = grid%mesh
  !
  ALLOCATE( psi(mesh,nwf) )
  CALL scf( zed, rel, grid, nwf, nn, ll, oc, enl, psi, ierr )
  IF ( ierr /= 0 ) THEN
    ae_note(nt) = 'atomic scf not converged'
    DEALLOCATE( psi )
    CALL deallocate_radial_grid( grid )
    RETURN
  ENDIF
  !
  ALLOCATE( rho_core_ae(mesh), rho_ps(mesh) )
  CALL core_density( mesh, grid%r2, nwf, oc, core, psi, rho_core_ae )
  !
  ! ... the atomic pseudo valence density, interpolated onto the fine
  ! ... logarithmic grid the core density lives on
  !
  CALL pseudo_density( upf_in, mesh, grid%r, rho_ps )
  CALL integrate( nt, mesh, grid%r, grid%rab, rho_ps, rho_core_ae, &
                     upf_in%r(upf_in%mesh) )
  !
  DEALLOCATE( rho_core_ae, rho_ps, psi )
  CALL deallocate_radial_grid( grid )
  !
  IF ( ae_valid(nt) ) WRITE( ae_note(nt), '("Z=",I3,", ",I2," core shells, ",A)' ) &
    NINT(zed), ncore, TRIM(treatment)
  !
END SUBROUTINE species
!
!----------------------------------------------------------------------------
SUBROUTINE core_shells( zcore, ncore, ierr )
  !--------------------------------------------------------------------------
  !! How many shells, in spectroscopic order, the given core charge holds.
  !! The partition has to come out exact: a core that is not a set of closed
  !! shells means the assumption behind this whole correction is wrong for
  !! that pseudopotential.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: zcore
  INTEGER, INTENT(OUT) :: ncore, ierr
  !
  INTEGER :: i
  REAL(DP) :: charge
  !
  ierr = 1
  ncore = 0
  charge = 0._DP
  !
  DO i = 1, nshell
    charge = charge + DBLE( 2*(2*l_spec(i)+1) )
    IF ( ABS(charge-zcore) < charge_eps ) THEN
      ncore = i
      ierr = 0
      RETURN
    ENDIF
    IF ( charge > zcore ) EXIT
  ENDDO
  !
END SUBROUTINE core_shells
!
!----------------------------------------------------------------------------
SUBROUTINE neutral_config( zed, ncore, nwf, nn, ll, oc, core, ierr )
  !--------------------------------------------------------------------------
  !! Ground-state configuration of the neutral atom, filled in Madelung
  !! order, and the flag that marks the first \(\texttt{ncore}\) shells of
  !! the spectroscopic order within it. On exit the shells are ordered as
  !! they were filled.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: zed
  INTEGER, INTENT(IN) :: ncore
  INTEGER, INTENT(OUT) :: nwf, nn(nshell), ll(nshell), ierr
  REAL(DP), INTENT(OUT) :: oc(nshell)
  LOGICAL, INTENT(OUT) :: core(nshell)
  !
  INTEGER :: i, j
  REAL(DP) :: full, left
  !
  ierr = 0
  nwf = 0
  nn(:) = 0
  ll(:) = 0
  oc(:) = 0._DP
  core(:) = .FALSE.
  !
  left = zed
  DO i = 1, nshell
    IF ( left <= charge_eps ) EXIT
    full = DBLE( 2*(2*l_mad(i)+1) )
    nwf = nwf + 1
    nn(nwf) = n_mad(i)
    ll(nwf) = l_mad(i)
    oc(nwf) = MIN( full, left )
    left = left - oc(nwf)
  ENDDO
  !
  ! ... a core shell must be present and closed, or the core charge built
  ! ... from this configuration would not be the one the pseudopotential
  ! ... removed
  !
  DO i = 1, ncore
    DO j = 1, nwf
      IF ( nn(j) /= n_spec(i) .OR. ll(j) /= l_spec(i) ) CYCLE
      IF ( ABS( oc(j) - DBLE(2*(2*l_spec(i)+1)) ) > charge_eps ) EXIT
      core(j) = .TRUE.
      EXIT
    ENDDO
  ENDDO
  !
  IF ( COUNT( core(1:nwf) ) /= ncore ) ierr = 1
  !
END SUBROUTINE neutral_config
!
!----------------------------------------------------------------------------
SUBROUTINE core_density( mesh, r2, nwf, oc, core, psi, rho_core_ae )
  !--------------------------------------------------------------------------
  !! Spherical all-electron core charge from the core orbitals, which are
  !! radial wavefunctions multiplied by r. Not to be confused with the core
  !! charge of the nonlinear core correction, which \(\texttt{tb09_c}\)
  !! receives from \(\texttt{v_xc_meta}\) as \(\texttt{rho_core}\).
  !
  USE constants, ONLY : fpi
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: mesh, nwf
  REAL(DP), INTENT(IN) :: r2(mesh), oc(nwf), psi(mesh,nwf)
  LOGICAL, INTENT(IN) :: core(nwf)
  REAL(DP), INTENT(OUT) :: rho_core_ae(mesh)
  !
  INTEGER :: i
  !
  rho_core_ae(:) = 0._DP
  DO i = 1, nwf
    IF ( .NOT. core(i) ) CYCLE
    rho_core_ae(1:mesh) = rho_core_ae(1:mesh) + &
      oc(i) * psi(1:mesh,i)**2 / ( fpi * r2(1:mesh) )
  ENDDO
  !
END SUBROUTINE core_density
!
!----------------------------------------------------------------------------
SUBROUTINE integrate( nt, mesh, r, rab, rho_ps, rho_core_ae, rmax_pp )
  !--------------------------------------------------------------------------
  !! Integrate the difference of the two averages and store in the module
  !! variables. The all-electron density is the pseudo one with its core charge
  !! replaced by the true core, so everything else, the pseudization of the
  !! valence in particular, cancels pointwise.
  !
  USE constants, ONLY : fpi
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nt, mesh
  REAL(DP), INTENT(IN) :: r(mesh), rab(mesh), rho_ps(mesh), rho_core_ae(mesh)
  REAL(DP), INTENT(IN) :: rmax_pp
  !! end of the pseudopotential mesh: beyond it the pseudo density is not
  !! known and the integrand would be meaningless
  !
  INTEGER :: i, mesh_int
  REAL(DP), ALLOCATABLE :: grad_ae(:), grad_ps(:), integrand(:), rho_ae(:)
  !
  mesh_int = integration_mesh( mesh, r, rmax_pp )
  IF ( mesh_int < 5 ) THEN
    ae_note(nt) = 'radial grid too short'
    RETURN
  ENDIF
  !
  ALLOCATE( rho_ae(mesh), grad_ae(mesh), grad_ps(mesh), integrand(mesh) )
  !
  rho_ae(1:mesh) = rho_ps(1:mesh) + rho_core_ae(1:mesh)
  !
  CALL radial_gradient( rho_ae, grad_ae, r, mesh, 1 )
  CALL radial_gradient( rho_ps, grad_ps, r, mesh, 1 )
  !
  integrand(:) = 0._DP
  DO i = 1, mesh_int
    integrand(i) = fpi * r(i)**2 * &
      ( ABS(grad_ae(i)) / MAX( ABS(rho_ae(i)), density_floor ) - &
        ABS(grad_ps(i)) / MAX( ABS(rho_ps(i)), density_floor ) )
  ENDDO
  !
  CALL simpson( mesh_int, integrand, rab, ae_grad_int(nt) )
  ae_valid(nt) = .TRUE.
  !
  DEALLOCATE( rho_ae, grad_ae, grad_ps, integrand )
  !
END SUBROUTINE integrate
!
!----------------------------------------------------------------------------
FUNCTION integration_mesh( mesh, r, rmax_pp ) RESULT( mesh_int )
  !--------------------------------------------------------------------------
  !! Largest odd index with r below the truncation radius, so that Simpson
  !! integration is exact on the logarithmic grid. Outside the core region
  !! the two densities coincide and the integrand vanishes; the cut is also
  !! what keeps the integral away from the end of the pseudopotential mesh,
  !! where the pseudo density is only known by extrapolation.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: mesh
  REAL(DP), INTENT(IN) :: r(mesh), rmax_pp
  INTEGER :: mesh_int
  !
  INTEGER :: i
  REAL(DP) :: rcut
  !
  rcut = MIN( rcut_ae, rmax_pp )
  mesh_int = mesh
  DO i = 1, mesh
    IF ( r(i) > rcut ) THEN
      mesh_int = i - 1
      EXIT
    ENDIF
  ENDDO
  mesh_int = 2 * ( (mesh_int+1) / 2 ) - 1
  !
END FUNCTION integration_mesh
!
!----------------------------------------------------------------------------
SUBROUTINE pseudo_density( upf_in, mesh, r, rho_ps )
  !--------------------------------------------------------------------------
  !! The density the plane-wave average is taken over: the atomic pseudo
  !! valence density of the pseudopotential file, divided by \(4\pi r^2\),
  !! plus the core charge of the nonlinear core correction, interpolated
  !! onto the grid the all-electron core density lives on.
  !
  USE constants,    ONLY : fpi
  USE pseudo_types, ONLY : pseudo_upf
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf), INTENT(IN) :: upf_in
  INTEGER, INTENT(IN) :: mesh
  REAL(DP), INTENT(IN) :: r(mesh)
  REAL(DP), INTENT(OUT) :: rho_ps(mesh)
  !
  INTEGER :: i
  REAL(DP), ALLOCATABLE :: rho_pp(:)
  !
  ALLOCATE( rho_pp(upf_in%mesh) )
  !
  DO i = 1, upf_in%mesh
    IF ( upf_in%r(i) > 0._DP ) THEN
      rho_pp(i) = upf_in%rho_at(i) / ( fpi * upf_in%r(i)**2 )
    ELSE
      rho_pp(i) = 0._DP
    ENDIF
  ENDDO
  ! ... rho_at goes as r^2, so the ratio is finite but 0/0 at a grid
  ! ... starting exactly at the origin, as ONCVPSP grids do
  IF ( upf_in%r(1) <= 0._DP .AND. upf_in%mesh > 1 ) rho_pp(1) = rho_pp(2)
  !
  CALL interpolate( upf_in%mesh, upf_in%r, rho_pp, mesh, r, rho_ps, .FALSE. )
  !
  IF ( upf_in%nlcc ) CALL interpolate( upf_in%mesh, upf_in%r, &
                          upf_in%rho_atc, mesh, r, rho_ps, .TRUE. )
  !
  DEALLOCATE( rho_pp )
  !
END SUBROUTINE pseudo_density
!
!----------------------------------------------------------------------------
SUBROUTINE interpolate( nin, rin, fin, nout, rout, fout, accumulate )
  !--------------------------------------------------------------------------
  !! Cubic-spline interpolation of a radial function onto another grid.
  !! Outside the input range the end values are held constant.
  !
  USE splinelib, ONLY : spline, splint
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nin, nout
  REAL(DP), INTENT(IN) :: rin(nin), fin(nin), rout(nout)
  REAL(DP), INTENT(INOUT) :: fout(nout)
  LOGICAL, INTENT(IN) :: accumulate
  !! add to fout instead of overwriting it
  !
  INTEGER :: i
  REAL(DP) :: value
  REAL(DP), ALLOCATABLE :: d2y(:)
  !
  ALLOCATE( d2y(nin) )
  CALL spline( rin, fin, 0._DP, 0._DP, d2y )
  !
  DO i = 1, nout
    IF ( rout(i) < rin(1) ) THEN
      value = fin(1)
    ELSEIF ( rout(i) > rin(nin) ) THEN
      value = fin(nin)
    ELSE
      value = splint( rin, fin, d2y, rout(i) )
    ENDIF
    IF ( accumulate ) THEN
      fout(i) = fout(i) + value
    ELSE
      fout(i) = value
    ENDIF
  ENDDO
  !
  DEALLOCATE( d2y )
  !
END SUBROUTINE interpolate
!
!----------------------------------------------------------------------------
SUBROUTINE scf( zed, rel, grid, nwf, nn, ll, oc, enl, psi, ierr )
  !--------------------------------------------------------------------------
  !! Spin-restricted all-electron radial self-consistent field, in the local
  !! density approximation, with Anderson mixing of the screening potential.
  !
  !! LDA is used regardless of DFT functional from the PP.Only the
  !! core orbitals are kept, the valence density being taken from the
  !! pseudopotential file.
  !
  USE constants,     ONLY : e2, fpi, c_au
  USE radial_grids,  ONLY : hartree
  USE atomic_radial, ONLY : ascheq, lschps
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: zed
  INTEGER, INTENT(IN) :: rel, nwf
  TYPE(radial_grid_type), INTENT(IN) :: grid
  INTEGER, INTENT(IN) :: nn(nwf), ll(nwf)
  REAL(DP), INTENT(IN) :: oc(nwf)
  REAL(DP), INTENT(OUT) :: enl(nwf), psi(grid%mesh,nwf)
  INTEGER, INTENT(OUT) :: ierr
  !
  INTEGER, PARAMETER :: maxter = 200
  REAL(DP), PARAMETER :: beta = 0.2_DP, tr2 = 1.0E-14_DP, thresh = 1.0E-10_DP
  !
  LOGICAL  :: conv
  INTEGER  :: i, iter, mesh, n, nerr, nin, nstop
  REAL(DP) :: ze2
  REAL(DP), ALLOCATABLE :: rho(:), vh(:), vnuc(:), vpot(:)
  REAL(DP), ALLOCATABLE :: vin(:), vout(:), vin_old(:), res_old(:)
  !
  mesh = grid%mesh
  ze2 = - zed * e2
  ierr = 1
  conv = .FALSE.
  !
  ALLOCATE( rho(mesh), vh(mesh), vnuc(mesh), vpot(mesh) )
  ALLOCATE( vin(mesh), vout(mesh), vin_old(mesh), res_old(mesh) )
  psi(:,:) = 0._DP
  vin_old(:) = 0._DP
  res_old(:) = 0._DP
  !
  ! ... only the screening potential is mixed and tested for convergence.
  ! ... The nuclear term diverges at the origin, and the round-off of the
  ! ... sum there is larger than any sensible convergence threshold
  !
  vnuc(1:mesh) = ze2 / grid%r(1:mesh)
  CALL starting_potential( zed, mesh, grid%r, nwf, nn, ll, oc, enl, vpot )
  vin(1:mesh) = vpot(1:mesh) - vnuc(1:mesh)
  !
  DO iter = 1, maxter
    !
    vpot(1:mesh) = vnuc(1:mesh) + vin(1:mesh)
    !
    nerr = 0
    DO n = 1, nwf
      nin = mesh
      IF ( rel == 0 ) THEN
        CALL ascheq( nn(n), ll(n), enl(n), mesh, grid, vpot, ze2, thresh, &
                     psi(:,n), nstop )
      ELSE
        CALL lschps( 1, zed, thresh, grid, nin, nn(n), ll(n), enl(n), vpot, &
                     psi(:,n), nstop, c_au )
      ENDIF
      ! ... far from self-consistency the radial solver routinely misses a
      ! ... state; only the last iteration has to be clean
      IF ( nstop /= 0 ) THEN
        psi(:,n) = 0._DP
        nerr = nerr + 1
      ENDIF
    ENDDO
    !
    ! ... rho is the density multiplied by 4 pi r^2, as hartree expects
    !
    rho(:) = 0._DP
    DO n = 1, nwf
      IF ( oc(n) <= 0._DP ) CYCLE
      rho(1:mesh) = rho(1:mesh) + oc(n) * psi(1:mesh,n)**2
    ENDDO
    !
    CALL hartree( 0, 2, mesh, grid, rho, vh )
    !
    DO i = 1, mesh
      vout(i) = e2*vh(i) + vxc_lda( rho(i) / ( fpi * grid%r2(i) ) )
    ENDDO
    !
    CALL mix_anderson( mesh, iter, beta, tr2, vin, vout, vin_old, &
                          res_old, conv )
    !
    IF ( conv ) THEN
      IF ( nerr == 0 ) ierr = 0
      EXIT
    ENDIF
    !
  ENDDO
  !
  DEALLOCATE( rho, vh, vnuc, vpot )
  DEALLOCATE( vin, vout, vin_old, res_old )
  !
END SUBROUTINE scf
!
!----------------------------------------------------------------------------
SUBROUTINE mix_anderson( n, iter, beta, tr2, vin, vout, vin_old, res_old, conv )
  !--------------------------------------------------------------------------
  !! Anderson mixing with one history vector: atomic/src/dmixp.f90.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, iter
  REAL(DP), INTENT(IN) :: beta, tr2, vout(n)
  REAL(DP), INTENT(INOUT) :: vin(n), vin_old(n), res_old(n)
  LOGICAL, INTENT(OUT) :: conv
  !
  REAL(DP) :: den, num, theta
  REAL(DP), ALLOCATABLE :: res(:), vbar_in(:), vbar_out(:)
  !
  ALLOCATE( res(n), vbar_in(n), vbar_out(n) )
  !
  res(:) = vout(:) - vin(:)
  conv = ( DOT_PRODUCT( res, res ) / DBLE(n) < tr2 )
  !
  IF ( .NOT. conv ) THEN
    !
    vbar_in(:) = vin(:)
    vbar_out(:) = vout(:)
    !
    IF ( iter > 1 ) THEN
      den = DOT_PRODUCT( res-res_old, res-res_old )
      IF ( den > 0._DP ) THEN
        num = DOT_PRODUCT( res, res-res_old )
        theta = num / den
        ! ... an unbounded extrapolation is worse than none at all
        theta = MAX( MIN( theta, 2._DP ), -2._DP )
        vbar_in(:) = vin(:) - theta * ( vin(:) - vin_old(:) )
        vbar_out(:) = vout(:) - theta * ( res(:) - res_old(:) + &
                                          vin(:) - vin_old(:) )
      ENDIF
    ENDIF
    !
    vin_old(:) = vin(:)
    res_old(:) = res(:)
    vin(:) = vbar_in(:) + beta * ( vbar_out(:) - vbar_in(:) )
    !
  ENDIF
  !
  DEALLOCATE( res, vbar_in, vbar_out )
  !
END SUBROUTINE mix_anderson
!
!----------------------------------------------------------------------------
FUNCTION vxc_lda( rho ) RESULT( vxc )
  !--------------------------------------------------------------------------
  !! Slater exchange plus Perdew-Zunger correlation, in Ry.
  !
  USE constants, ONLY : e2, pi
  USE exch_lda,  ONLY : slater
  USE corr_lda,  ONLY : pz
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho
  REAL(DP) :: vxc
  !
  REAL(DP), PARAMETER :: rho_eps = 1.0E-30_DP
  REAL(DP) :: ec, ex, rs, vc, vx
  !
  vxc = 0._DP
  IF ( rho < rho_eps ) RETURN
  !
  rs = ( 3._DP / ( 4._DP * pi * rho ) )**( 1._DP/3._DP )
  !
  CALL slater( rs, ex, vx )
  CALL pz( rs, 1, ec, vc )
  !
  vxc = e2 * ( vx + vc )
  !
END FUNCTION vxc_lda
!
!----------------------------------------------------------------------------
SUBROUTINE starting_potential( zed, mesh, r, nwf, nn, ll, oc, enl, vpot )
  !--------------------------------------------------------------------------
  !! Generalized Thomas-Fermi atomic potential, as in starting_potential of
  !! ld1.x, with hydrogenic starting eigenvalues.
  !
  USE constants, ONLY : e2
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: mesh, nwf, nn(nwf), ll(nwf)
  REAL(DP), INTENT(IN) :: zed, r(mesh), oc(nwf)
  REAL(DP), INTENT(OUT) :: enl(nwf), vpot(mesh)
  !
  INTEGER :: i, n
  REAL(DP) :: enne, oce, t, x, zen
  !
  enne = 0._DP
  DO n = 1, nwf
    enne = enne + MAX( oc(n), 0._DP )
    zen = 0._DP
    DO i = 1, nwf
      oce = MAX( oc(i), 0._DP )
      IF ( nn(i) < nn(n) ) zen = zen + oce
      IF ( nn(i) == nn(n) .AND. ll(i) <= ll(n) ) zen = zen + oce
    ENDDO
    zen = MAX( zed - zen + 1._DP, 1._DP )
    enl(n) = - ( zen / nn(n) )**2
  ENDDO
  !
  DO i = 1, mesh
    x = r(i) * enne**(1._DP/3._DP) / 0.885_DP
    t = zed / ( 1._DP + SQRT(x)*(0.02747_DP-x*(0.1486_DP-0.007298_DP*x)) &
                + x*(1.243_DP+x*(0.2302_DP+0.006944_DP*x)) )
    t = MAX( 1._DP, t )
    vpot(i) = - e2 * t / r(i)
  ENDDO
  !
END SUBROUTINE starting_potential
!
!
END MODULE ctb09
