# Makefile for CRPA
sinclude ../make.inc

default: all

all: w90mod elph crpa

w90mod:	
	( cd wannier ; $(MAKE) all || exit 1 )

crpa-lib: w90mod
	( cd src ; $(MAKE) libs-crpa || exit 1 )

crpa: w90mod
	( cd src ; $(MAKE) all || exit 1 )

elph: w90mod
	( cd elph ; $(MAKE) all || exit 1 )


clean: w90mod_clean crpa_clean elph_clean

w90mod_clean:
	( cd wannier ; $(MAKE) clean )

crpa_clean:
	( cd src ; $(MAKE) clean )

elph_clean:
	( cd elph ; $(MAKE) clean )

# examples_clean:
#     if test -d examples ; then \
#     ( cd examples ; ./clean_all ) ; fi

# doc:
#     if test -d Doc ; then \
#     ( cd Doc ; $(MAKE) all || exit 1 ) ; fi

# doc_clean:
#     if test -d Doc ; then \
#     (cd Doc ; $(MAKE) clean ) ; fi

distclean: clean #doc_clean
