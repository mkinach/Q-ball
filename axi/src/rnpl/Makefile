########################################################################
#
# Makefile for RNPL application 'qball'
#
# Requires following environment variables to be set appropriately:
#
# RNPL_RNPL
# RNPL_F77
# RNPL_F77LOAD
# RNPL_F77PP
# RNPL_FLIBS
#
########################################################################
.IGNORE:

SHELL := /bin/sh

# set 'APP' to application name stem (prefix) then execute 'make fix'
# to convert Makefile to use explicit targets
APP := qball

# if your application uses headers and/or libraries from non-system
# locations, define the following macros appropriately
USER_INC_PATHS :=
USER_LIB_PATHS :=

# if you want non-default flags for the compilers, set them here
RNPL     := $(RNPL_RNPL)
F77      := $(RNPL_F77)
F77_LOAD := $(RNPL_F77LOAD)
F77PP    := $(RNPL_F77PP)
FLIBS    := $(RNPL_FLIBS)

## uncomment this section to enable debug flags
##
## -ggdb: produce debugging info for GDB
## -Wall: verbose warnings at compile time
## -fcheck=bounds: runtime check array bounds
## -ffpe-trap=underflow,denormal,overflow,invalid: stop execution on FPE
## -finit-real=snan: check uninitialized variables
#DEBUGFLAGS  := -ggdb -fcheck=bounds -ffpe-trap=overflow,invalid,zero -finit-real=snan
##DEBUGFLAGS := $(DEBUGFLAGS) -ffpe-trap=underflow,denormal  # usually not a meaningful error
#DEBUGFLAGS  := $(DEBUGFLAGS) -Wall
#F77         := gfortran $(DEBUGFLAGS) -fno-second-underscore -Wno-argument-mismatch
##F77        := gfortran $(DEBUGFLAGS) -fno-second-underscore -fallow-argument-mismatch
##F77        := gfortran $(DEBUGFLAGS) -fno-second-underscore
#F77_LOAD    := $(F77) -L/usr/local/lib

########################################################################
# Main targets
########################################################################

# default targets
all: qball qball_init

# rename instances in this Makefile (see 'APP' above)
fix: Makefile
	sed "s@qball@$(APP)@g" < Makefile > .Makefile
	mv .Makefile Makefile

# define a default rule to automatically make .o files out of .f files
.f.o:
	$(F77) -c $*.f

qball: qball.o updates.o residuals.o interp.o $(wildcard *.inc)
	$(F77_LOAD) qball.o updates.o residuals.o interp.o \
	$(FLIBS) -lodepack -lvutil -llinpack -llapack -lblas \
	-o qball

qball.f: qball_rnpl
	$(RNPL) -l allf qball_rnpl
	rnpl_fix_f77 updates.f initializers.f residuals.f

updates.f:      qball_rnpl
residuals.f:    qball_rnpl
initializers.f: qball_rnpl
interp.f:       qball_rnpl

qball_init: qball_init.o updates.o initializers.o residuals.o interp.o
	$(F77_LOAD) qball_init.o updates.o residuals.o initializers.o interp.o \
	$(FLIBS) -lodepack -lvutil -llinpack -llapack -lblas \
	-o qball_init

qball_init.f: qball_rnpl

########################################################################
# Clean targets
########################################################################

clean:
	/bin/rm *.o > /dev/null 2>&1
	/bin/rm .rnpl.attributes > /dev/null 2>&1

vclean:
	make clean
	/bin/rm qball qball_init > /dev/null 2>&1
	/bin/rm residuals.f updates.f initializers.f qball.f qball_init.f > /dev/null 2>&1
	/bin/rm gfuni0.inc globals.inc other_glbs.inc sys_param.inc > /dev/null 2>&1
	/bin/rm allf.f allf_init.f > /dev/null 2>&1
	/bin/rm *.out > /dev/null 2>&1
