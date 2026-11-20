# Makefile for W90_Modules

include ../make.inc

VPATH = @srcdir@

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)$(BUILDDIR)/PW/src \
          $(MOD_FLAG)$(BUILDDIR)/dft-d3

W90OBJS = \
wan2pw.o

W90LIBS = \
w90_utils.o \
w90_interface.o \
w90_shift_center.o \
w90_multiply_iqr.o \
w90_wigner.o \
w90_real_space.o \
w90_readwrite_wan_Rr.o \
w90_wan_Rr_buffer.o \
w90_interpolate.o \
w90_hamiltonian.o

PWOBJS = ../PW/src/libpw.a ../KS_Solvers/libks_solvers.a ../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)

# see comments in file tldeps for explanations
TLDEPS=pwlibs mods
include ../install/tldeps

all : libw90mod.a wan2pw.x

wan2pw.x : wan2pw.o libw90mod.a $(PWOBJS)
	$(LD) $(LDFLAGS) -o $@ wan2pw.o libw90mod.a \
	    $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../bin ; ln -fs ../W90_Modules/wan2pw.x . )

libw90mod.a : $(W90LIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

include make.depend
