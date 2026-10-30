# Makefile for HP

include ../../make.inc

VPATH = @srcdir@

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
	  $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../W90_Modules
HPOBJS = \
crpacom.o \
crpa_allocate_q.o \
crpa_compute_perturbation_q.o \
crpa_copy_original_xk.o \
crpa_clean_q.o \
crpa_close_q.o \
crpa_dealloc_1.o \
crpa_dealloc_2.o \
crpa_dealloc_q.o \
crpa_generate_grids.o \
crpa_init.o \
crpa_init_q.o \
crpa_load_q.o \
crpa_openfil_q.o \
crpa_postprocess_coulomb.o \
crpa_prepare_perturbations_q.o \
crpa_prepare_q.o \
crpa_print_clock.o \
crpa_q_points.o \
crpa_R_points.o \
crpa_run_nscf.o \
crpa_setup_q.o \
crpa_setup_wannier.o \
crpa_solve_linear_system.o \
crpa_stop_smoothly.o \
crpa_summary.o \
crpa_symmetry.o \
crpa_update_coulomb.o \
lr_compute_matrix_elements.o \
hp_allocate_q.o \
hp_bcast_input.o \
hpcom.o \
crpa_bcast_input.o \
crpa_readin.o \
crpa_run.o \
hp_clean_q.o \
hp_close_q.o \
hp_calc_chi.o \
hp_check_pert.o \
hp_check_type.o \
hp_dnstot_sum_q.o \
hp_dnsq.o \
hp_dealloc_q.o \
hp_dealloc_1.o \
hp_dealloc_2.o \
hp_dvpsi_pert.o \
hp_find_inequiv_sites.o \
hp_generate_grids.o \
hp_init.o \
hp_init_q.o \
hp_load_q.o \
hp_ns_trace.o \
hp_openfil_q.o \
hp_psym_dmag.o \
hp_psymdvscf.o \
hp_postproc.o \
hp_prepare_q.o \
hp_print_clock.o \
hp_q_points.o \
hp_R_points.o \
hp_rotate_dnsq.o \
hp_readin.o \
hp_read_chi.o \
hp_read_dnsq.o \
hp_run_nscf.o \
hp_solve_linear_system.o \
hp_symdnsq.o \
hp_symdvscf.o \
hp_sym_dmag.o \
hp_setup_q.o \
hp_summary.o \
hp_summary_q.o \
hp_stop_smoothly.o \
hp_write_chi.o \
hp_write_chi_full.o \
hp_write_dnsq.o 

LRMODS = ../../LR_Modules/liblrmod.a
W90MODS = ../../W90_Modules/libw90mod.a
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)

# see comments in file tldeps for explanations
TLDEPS=pwlibs lrmods w90mods
include ../../install/tldeps

all : hplibs hp.x

hplibs : libhp.a

hp.x : hp_main.o libhp.a $(PWOBJS) $(LRMODS) $(W90MODS)
	$(LD) $(LDFLAGS) -o $@ hp_main.o libhp.a \
		$(LRMODS) $(W90MODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../HP/src/hp.x . )

libhp.a : $(HPOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/hp.x 

ifneq ($(TOLERATE_MISSING_DEPEND),)
-include make.depend
else
ifeq ($(wildcard make.depend),)
$(error make.depend not found in $(CURDIR) -- run "make depend" first, or "make veryclean"/"make distclean" to clean without it)
endif
include make.depend
endif
