# Makefile for ELPH

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
		  $(MOD_FLAG)../../PW/src \
		  $(MOD_FLAG)../../dft-d3 \
		  $(MOD_FLAG)../../LR_Modules \
		  $(MOD_FLAG)../../PHonon/PH \
		  $(MOD_FLAG)../wannier

ELPHOBJS = \
elphcom.o \
elph_bcast_input.o \
elph_readin.o \
elph_compute_mel.o \
elph_prepare_q.o \
elph_q_points.o \
elph_load_q.o \


PWOBJS = ../../PHonon/PH/libph.a ../../LR_Modules/liblrmod.a ../wannier/libw90mod.a \
	../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a \
	$(BASEMODS)

TLDEPS= phlibs

all : tldeps libs-elph elph.x

libs-elph : libelph.a

elph.x : elph_main.o libelph.a $(PWOBJS)
	$(LD) $(LDFLAGS) -o $@ elph_main.o libelph.a \
		$(PWOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../CRPA/elph/elph.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

libelph.a : $(ELPHOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -rf ../../bin/elph.x

include make.depend
# DO NOT DELETE
