# Makefile for CRPA

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
		  $(MOD_FLAG)../../PW/src \
		  $(MOD_FLAG)../../dft-d3 \
		  $(MOD_FLAG)../../LR_Modules \
		  $(MOD_FLAG)../../PHonon/PH \
		  $(MOD_FLAG)../wannier

CRPAOBJS = \
crpacom.o \
crpa_readin.o \
crpa_bcast_input.o \
crpa_solve_linear_system.o \
crpa_generate_grids.o \
crpa_q_points.o \
crpa_prepare_q.o \
crpa_R_points.o \
crpa_run_nscf.o \
crpa_load_q.o \
crpa_openfil_q.o \
crpa_init_q.o \
crpa_setup_q.o \
crpa_summary_q.o \
crpa_allocate_q.o \
crpa_dvpsi_pert.o \
crpa_stop_smoothly.o \
crpa_dealloc_1.o \
crpa_dealloc_2.o \
crpa_psym_dmag.o \
crpa_psymdvscf.o \
crpa_clean_q.o \
crpa_print_clock.o \
crpa_symdvscf.o \
crpa_sym_dmag.o \
crpa_dealloc_q.o \
crpa_close_q.o \
crpa_init.o \
crpa_compute_perturbation_q.o \
crpa_build_freq_grid.o \
crpa_coul_lanczos.o \
crpa_allocate_coul_lanczos.o \
crpa_restart_coul_lanczos.o \
crpa_restart_coul_cg.o \
crpa_average_w.o \
crpa_do_lanczos_coulomb.o \
crpa_calc_coulomb_lanczos.o \
crpa_extrapolate_coul.o \
crpa_diel_lanczos_chain_io.o \
crpa_coul_lanczos_chain_io.o \
crpa_run_spectrum_only.o \
crpa_symmetry.o \
crpa_summary.o \
crpa_setup_wannier.o \
crpa_prepare_perturbations_q.o \
crpa_postprocess_coulomb.o \
crpa_update_coulomb.o \
crpa_copy_original_xk.o \
lr_compute_matrix_elements.o \
crpa_q0correction.o \
crpa_phescf_w.o \
crpa_solve_e_w.o \
crpa_dielec.o \
crpa_plasma_frequency.o \
crpa_diel_lanczos.o \
crpa_allocate_diel_lanczos.o \
crpa_do_lanczos_epsilon.o \
crpa_psym_diel.o \
crpa_calc_chi_epsilon.o \
crpa_extrapolate_diel.o \
crpa_restart_diel_lanczos.o

# crpa_skip_nscf_by_rotate.o

PWOBJS = ../../PHonon/PH/libph.a ../../LR_Modules/liblrmod.a ../wannier/libw90mod.a \
	../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a \
	$(BASEMODS)

TLDEPS= crpalibs phlibs

all : tldeps libs-crpa crpa.x

libs-crpa : libcrpa.a

crpa.x : crpa_main.o libcrpa.a $(PWOBJS)
	$(LD) $(LDFLAGS) -o $@ crpa_main.o libcrpa.a \
		$(PWOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../CRPA/src/crpa.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

libcrpa.a : $(CRPAOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/crpa.x

include make.depend
# DO NOT DELETE
