module ns_constraint

   USE upf_params,    ONLY : lqmax
   use kinds,     only : DP
   use io_global, only : stdout
   use ions_base, only : nat, ityp
   USE parameters,    ONLY : ntypx, natx 
   use ldaU,      only : Hubbard_l, Hubbard_lmax, starting_ns
   use lsda_mod,  only : nspin
   USE noncollin_module, ONLY : noncolin 
   
   SAVE
   
   REAL(DP), ALLOCATABLE :: Hubbard_constraints(:, :, :, :)
   REAL(DP), ALLOCATABLE :: Hubbard_occupations(:, :, :, :)
   LOGICAL :: Hubbard_constraining(natx)
   REAL(DP) :: input_Hubbard_occupations(lqmax, lqmax, 2, natx)
   REAL(DP) :: Hubbard_conv_thr
   REAL(DP) :: Hubbard_mixing_beta
   INTEGER :: Hubbard_maxstep
   REAL(DP) :: Hubbard_strength
   LOGICAL :: Hubbard_conv=.true.
   CHARACTER(len=30) :: Hubbard_constraint_type

contains
   ! Should be called after first ns_adj
   subroutine setup_constraints(ns)
      real(DP), intent(inout) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
      integer :: na, nt, is, m1, m2, ldim, i, j, l 
      logical :: full_mat
      
      Hubbard_constraining(:) = .false.
      Hubbard_conv = .true.
   
      
      full_mat = ANY(input_Hubbard_occupations /= -1.0_dp)
      
      if (full_mat .OR. (ANY(starting_ns /= -1.0_dp) .AND. Hubbard_constraint_type /= "occupations")) then
         ALLOCATE(Hubbard_constraints(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2, nat), source=0.0_dp)
         ALLOCATE(Hubbard_occupations(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2, nat), source=0.0_dp)
         Hubbard_conv = .false.
         DO na = 1, nat  
            nt = ityp(na)
            Hubbard_constraining(na) = ANY(input_Hubbard_occupations(:, :, :, na) /= -1.0_dp) .OR. &
              & (ANY(starting_ns(:, :, nt) >= 0) .AND. Hubbard_constraint_type /= "occupations")
            if (ANY(input_Hubbard_occupations(:, :, :, na) /= -1.0_dp)) then
               ldim = 2 * Hubbard_l(nt) + 1
               do is = 1, nspin
                  do m1 = 1, ldim
                     do m2 = 1, ldim
                        if (input_Hubbard_occupations(m1, m2, is, na) /= -1.0_dp) then
                           Hubbard_occupations(m1, m2, is, na) = input_Hubbard_occupations(m1, m2, is, na)
                        endif
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif
     
      ! If we're constraining, we fill the ns matrix with Hubbard_occupations, or
      ! if instead starting_ns is used, we fill Hubbard_occupations with what was computed
      ! in ns_adj
      if (ANY(Hubbard_constraining)) THEN
         write (stdout,*) "Modify starting ns matrices according to input values "
         do na = 1, nat
            if (Hubbard_constraining(na)) then
               nt = ityp(na)
               ldim = 2 * Hubbard_l(nt) + 1 
               do m1 = 1, ldim
             	   do m2 = 1, ldim
                     do is = 1, nspin
                        if (full_mat) then
                           ns(m1, m2, is, na) = Hubbard_occupations(m1, m2, is, na)
                        else
                           Hubbard_occupations(m1, m2, is, na) = ns(m1, m2, is, na)
                        endif
                     enddo
                  enddo
               enddo
            endif
         enddo
         ! This we don't do at the moment because it's not clear if we should have symmetries on the
         ! Occupation matrix when we are doing the constraining
         ! We now symmetrize rho%ns into Hubbard_occupations and then recopy it into rho%ns
         ! if (.NOT. noncolin) then
            ! Hubbard_occupations = 0.0_DP
         	 ! call symmetrize_ns(rho%ns, Hubbard_occupations, .false.)
         	 ! rho%ns(:, :, :, :) = Hubbard_occupations(:, :, :, :)
      	  ! endif
      endif
   end subroutine setup_constraints
   
   subroutine constrain_ns(nscur, hubbard_conv, dr2, tr2, iter)


      implicit none
      !
      logical, intent(inout) :: hubbard_conv
      real(DP), intent(in) ::dr2, tr2
      integer, intent(in) :: iter

      real(DP), intent(in) :: nscur(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
      integer, parameter :: ldmx = 7
      integer :: m1, m2, ldim, na, is, nt, m3
      real(DP) :: tmp, d, totdist

      totdist = 0.0_dp
      if (ALL(Hubbard_occupations == -1.0_DP) .OR. (Hubbard_constraint_type == 'ns') .OR. &
        & (Hubbard_constraint_type == 'occupations')) then
        hubbard_conv = .true.
        return
      endif
      if (.not. hubbard_conv) then
         WRITE( stdout,*) 'Constraints:'
         write(stdout, '(" Hubbard_constraint_type = ", A)') Hubbard_constraint_type
         write(stdout, '(" Hubbard_strength = ", f7.3)') Hubbard_strength
         do na = 1, nat
   		 if (Hubbard_constraining(na)) then
               nt = ityp(na)
               ldim = 2 * Hubbard_l(nt) + 1
               do is = 1, nspin
                  tmp = sum(abs(nscur(1:ldim, 1:ldim, is, na) - Hubbard_occupations(1:ldim, 1:ldim, is, na)))
                  write(stdout, '(" atom = ", i1, ", spin = ", i1, ", dist =", f7.3)') na, is, tmp
                     do m1 = 1, ldim
                        do m2 = m1, ldim
                           tmp = nscur(m1, m2, is, na) - Hubbard_occupations(m1, m2, is, na)
                           Hubbard_constraints(m1, m2, is, na) = Hubbard_constraints(m1, m2, is, na) + &
                           Hubbard_mixing_beta * tmp
                           Hubbard_constraints(m2, m1, is, na) = Hubbard_constraints(m1, m2, is, na)
                           totdist = totdist + abs(tmp)
                       enddo
                    enddo
                 WRITE( stdout,'(5x,"constraint matrix:")')
                 DO m1 = 1, ldim
                    WRITE( stdout,'(5x,7f7.3)') ( DBLE(Hubbard_constraints(m1,m2, is, na)), m2=1, ldim )
                 ENDDO
               enddo
            endif
         enddo
         
         if (totdist < Hubbard_conv_thr .AND. iter > 1) then
             hubbard_conv = .true.
             WRITE(stdout,'("Hubbard_conv true: Hubbard_conv_thr reached")')
         elseif (iter > Hubbard_maxstep) then
             hubbard_conv = .true.
             WRITE(stdout,'("Hubbard_conv true: Hubbard_maxstep reached")')
         endif
      endif
   end subroutine constrain_ns
end module ns_constraint
