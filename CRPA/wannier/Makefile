# Makefile for PW

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
		  $(MOD_FLAG)../../dft-d3

W90OBJS = \
wannier2pw.o

W90LIBS = \
w90_interface.o \
w90_shift_center.o \
w90_multiply_iqr.o \
w90_wigner.o \
w90_real_space.o \
w90_readwrite_wan_Rr.o \
w90_wan_Rr_buffer.o \
w90_interpolate.o \


PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)

TLDEPS=bindir pwlibs mods

all : libw90mod.a wannier2pw.x

# wannier2pw.x : $(W90OBJS) libw90mod.a $(PWOBJS)
# 	$(LD) $(LDFLAGS) -o $@ \
# 	   libw90mod.a $(PWOBJS) $(QEMODS)
# 	- ( cd ../../bin; ln -fs ../CRPA/wannier/$@ . ; )
wannier2pw.x : wannier2pw.o libw90mod.a $(PWOBJS)
	$(LD) $(LDFLAGS) -o $@ wannier2pw.o libw90mod.a \
	    $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../CRPA/wannier/wannier2pw.x . )

libw90mod.a : $(W90LIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

#include make.depend

# DO NOT DELETE
