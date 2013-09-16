# -----------------------------------------------------------------
# $Revision: 1.3 $
# $Date: 2009/02/17 02:58:47 $
# -----------------------------------------------------------------
# Programmer: Radu Serban @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# Makefile for CVODE serial examples.
#
# This file is generated from a template using various variables
# set at configuration time. It can be used as a template for
# other user Makefiles.
#
# Note: if the solver was successfully configured with Blas/Lapack
# support, the Blas/Lapack libraries are specified through the 
# variable LIBRARIES_BL. Otherwise, this variable should contain
# an empty string. We include LIBRARIES_BL in the link line for
# all examples, whether they use the Lapack module or not, to
# address the case in which the SUNDIALS libraries are shared 
# objects. In that case, the solver library references Lapack 
# symbols which must be always resolved by linking against the
# Blas/Lapack libraries. If only static SUNDIALS libraries have 
# been built, it is not required to link the Blas/Lapack libraries
# for examples that do not use that module...
# -----------------------------------------------------------------

SHELL = /bin/bash

prefix       = /home/ejalbert/Documents/NumericalMethod/cvode/cvodeInstallDir
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib

CPP      = cc -E
CPPFLAGS = 
CC       = cc
CFLAGS   = -g -O2
LDFLAGS  = 
LIBS     = -lm 

INCLUDES = -I${includedir}
LIBRARIES = -lsundials_cvode -lsundials_nvecserial ${LIBS}
LIBRARIES_BL = 

#EXAMPLES = cvAdvDiff_bnd cvDirectDemo_ls cvDiurnal_kry_bp cvDiurnal_kry cvKrylovDemo_ls cvKrylovDemo_prec cvRoberts_dns cvRoberts_dns_uw solver
EXAMPLES = solver oldsolver

OBJECTS = ${EXAMPLES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CPPFLAGS} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------

all: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -o $${i} $${i}.o ${CFLAGS} ${LDFLAGS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL}" ; \
	  ${CC} -o $${i} $${i}.o ${CFLAGS} ${LDFLAGS} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL}; \
	done

clean:
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}

# -----------------------------------------------------------------------------------------

