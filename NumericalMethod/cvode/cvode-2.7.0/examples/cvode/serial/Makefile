# -----------------------------------------------------------------
# $Revision: 1.12 $
# $Date: 2009/02/17 02:58:47 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# Makefile for CVODE serial examples
#
# examples/cvode/serial/Makefile.  Generated from Makefile.in by configure.
# -----------------------------------------------------------------

SHELL = /bin/bash

srcdir       = .
top_builddir = ../../../
builddir     = .
abs_builddir = /home/ericmjalbert/Documents/research/numerAnalysis/FindC/cvode-2.7.0/examples/cvode/serial
prefix       = /home/ericmjalbert/Documents/research/numerAnalysis/FindC/cvodeInstallDir
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib

INSTALL      = /usr/bin/install -c
INSTALL_PROG = ${INSTALL}
INSTALL_FILE = ${INSTALL} -m 644

LIBTOOL      = $(SHELL) $(top_builddir)/libtool
LIBTOOL_DEPS = config/ltmain.sh

SHARED_LIBS  = no

CPP      = cc -E
CPPFLAGS = 
CC       = cc
CFLAGS   = -g -O2
LDFLAGS  = 
LIBS     = -lm 

LAPACK_ENABLED = no
BLAS_LAPACK_LIBS = 

OBJ_EXT = .lo
EXE_EXT = 

EXS_INSTDIR = /home/ericmjalbert/Documents/research/numerAnalysis/FindC/cvodeInstallDir/examples

top_srcdir = $(srcdir)/../../..

SUNDIALS_INCS = -I$(top_srcdir)/include -I$(top_builddir)/include

SUNDIALS_LIBS = $(top_builddir)/src/cvode/libsundials_cvode.la            \
	        $(top_builddir)/src/nvec_ser/libsundials_nvecserial.la

mkinstalldirs = $(SHELL) $(top_srcdir)/config/mkinstalldirs
rminstalldirs = $(SHELL) $(top_srcdir)/config/rminstalldirs

EXAMPLES = cvAdvDiff_bnd \
	   cvDirectDemo_ls \
	   cvDiurnal_kry_bp \
	   cvDiurnal_kry \
	   cvKrylovDemo_ls \
	   cvKrylovDemo_prec \
	   cvRoberts_dns \
	   cvRoberts_dns_uw 

EXAMPLES_BL = cvAdvDiff_bndL \
	      cvRoberts_dnsL

OBJECTS = ${EXAMPLES:=${OBJ_EXT}}
OBJECTS_BL = ${EXAMPLES_BL:=${OBJ_EXT}}
EXECS = ${EXAMPLES:=${EXE_EXT}}
EXECS_BL = ${EXAMPLES_BL:=${EXE_EXT}}

# ----------------------------------------------------------------------------------------------------------------------

all:
	@for i in ${EXAMPLES} ; do \
	  $(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(SUNDIALS_INCS) $(CFLAGS) -c $(srcdir)/$${i}.c -o $${i}${OBJ_EXT} ; \
	  $(LIBTOOL) --mode=link $(CC) -o $(builddir)/$${i}$(EXE_EXT) $(builddir)/$${i}${OBJ_EXT} $(CFLAGS) $(LDFLAGS) $(SUNDIALS_LIBS) $(LIBS) $(BLAS_LAPACK_LIBS) ; \
	done
	@if test "X${LAPACK_ENABLED}" = "Xyes"; then \
	for i in ${EXAMPLES_BL} ; do \
	  $(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(SUNDIALS_INCS) $(CFLAGS) -c $(srcdir)/$${i}.c -o $${i}${OBJ_EXT} ; \
	  $(LIBTOOL) --mode=link $(CC) -o $(builddir)/$${i}$(EXE_EXT) $(builddir)/$${i}${OBJ_EXT} $(CFLAGS) $(LDFLAGS) $(SUNDIALS_LIBS) $(LIBS) $(BLAS_LAPACK_LIBS) ; \
	done ; \
	fi

install:
	$(mkinstalldirs) $(EXS_INSTDIR)/cvode/serial
	$(INSTALL_FILE) Makefile_ex $(EXS_INSTDIR)/cvode/serial/Makefile
	$(INSTALL_FILE) $(top_srcdir)/examples/cvode/serial/README $(EXS_INSTDIR)/cvode/serial/
	for i in ${EXAMPLES} ; do \
	  $(INSTALL_FILE) $(top_srcdir)/examples/cvode/serial/$${i}.c   $(EXS_INSTDIR)/cvode/serial/ ; \
	  $(INSTALL_FILE) $(top_srcdir)/examples/cvode/serial/$${i}.out $(EXS_INSTDIR)/cvode/serial/ ; \
	done
	@if test "X${LAPACK_ENABLED}" = "Xyes"; then \
	for i in ${EXAMPLES_BL} ; do \
	  $(INSTALL_FILE) $(top_srcdir)/examples/cvode/serial/$${i}.c   $(EXS_INSTDIR)/cvode/serial/ ; \
	  $(INSTALL_FILE) $(top_srcdir)/examples/cvode/serial/$${i}.out $(EXS_INSTDIR)/cvode/serial/ ; \
	done ; \
	fi

uninstall:
	rm -f $(EXS_INSTDIR)/cvode/serial/Makefile
	rm -f $(EXS_INSTDIR)/cvode/serial/README
	for i in ${EXAMPLES} ; do \
	  rm -f $(EXS_INSTDIR)/cvode/serial/$${i}.c ; \
	  rm -f $(EXS_INSTDIR)/cvode/serial/$${i}.out ; \
	done
	@if test "X${LAPACK_ENABLED}" = "Xyes"; then \
	for i in ${EXAMPLES_BL} ; do \
	  rm -f $(EXS_INSTDIR)/cvode/serial/$${i}.c ; \
	  rm -f $(EXS_INSTDIR)/cvode/serial/$${i}.out ; \
	done ; \
	fi
	$(rminstalldirs) $(EXS_INSTDIR)/cvode/serial
	$(rminstalldirs) $(EXS_INSTDIR)/cvode

clean:
	rm -rf .libs
	rm -f *.lo *.o
	rm -f ${OBJECTS} ${OBJECTS_BL}
	rm -f $(EXECS) $(EXECS_BL)

distclean: clean
	rm -f Makefile
	rm -f Makefile_ex

# ----------------------------------------------------------------------------------------------------------------------

libtool: $(top_builddir)/$(LIBTOOL_DEPS)
	@cd ${top_builddir} ;                 \
	${SHELL} ./config.status --recheck ; \
	cd ${abs_builddir}
