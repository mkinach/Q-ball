# Ignore errors in execution of recipes for all files
.IGNORE:

# Use the Bourne shell, 'sh'
SHELL = /bin/sh

# Set 'APP' to the RNPL application name stem (prefix) then execute
# 'make fix' to convert Makefile to use explicit targets
APP       =  qball

# Shorten the name of some of the environment variables (convenience).
# Also, if you want non-default flags for the compiler, set them here.
RNPL      = $(RNPL_RNPL)
F77       = $(RNPL_F77)
F77_LOAD  = $(RNPL_F77LOAD)
F77PP     = $(RNPL_F77PP)
FLIBS     = $(RNPL_FLIBS)
CC        = mpicc

#### Uncomment below to enable debug flags (Fortran)
## -ggdb: produce debugging info for GDB
## -Wall: verbose warnings at compile time
## -fcheck=bounds: runtime check array bounds
## -ffpe-trap=underflow,denormal,overflow,invalid: stop execution for floating-point exceptions
## -finit-real=snan: check uninitialized variables
#FDEBUGFLAGS  = -ggdb -fcheck=bounds -ffpe-trap=overflow,invalid,zero -finit-real=snan
##FDEBUGFLAGS := $(FDEBUGFLAGS) -ffpe-trap=underflow,denormal  # usually not a meaningful error
#FDEBUGFLAGS := $(FDEBUGFLAGS) -Wall
#F77          = mpifort $(FDEBUGFLAGS) -fno-second-underscore
#F77_LOAD     = mpifort $(FDEBUGFLAGS) -fno-second-underscore -L/usr/local/lib
#
#### Uncomment below to enable debug flags (C)
## -ggdb: produce debugging info for GDB
## -Wall: verbose warnings at compile time
## -Wfloat-equal: warning if floating point numbers are using in equality comparisons
#CDEBUGFLAGS  = -ggdb
#CDEBUGFLAGS := $(CDEBUGFLAGS) -Wall -Wfloat-equal
#CCFLAGS			 =
#CC           = mpicc $(CDEBUGFLAGS)

# Define some variables that are convenient for PAMR compilation
PAMR_LIBS = -lpamr -lamrd -lbbhutil -lm -lmpi \
            -lodepack -lvutil -llinpack -llapack -lblas -lgfortran
PAMR_OBJS = qball-pamr.o initializers.o updates.o fcn.o residuals.o \
            qtotcalc-pamr.o

# A general rule for building object files out of '.f' files
.f.o:
	$(F77) -c $*.f 

# A general rule for building object files out of '.c' files
.c.o:
	$(CC) $(CCCFLAGS) $(CCFLAGS) -c $*.c

all:
	@printf "\n    Please run 'make rnpl' and 'make pamr' separately\n\n"

# For renaming targets, as explained above
fix: Makefile
	  sed "s@qball@$(APP)@g" < Makefile > .Makefile 
	  mv .Makefile Makefile


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               R N P L
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#rnpl: qball qball_init
rnpl: qball

# Note: rnpl_fix_f77 is a Perl program to fix syntax errors when RNPL
# generates Fortran 77 output. See the RNPL docs
qball.f: qball_rnpl
	$(RNPL) -l allf  qball_rnpl
	rnpl_fix_f77 updates.f initializers.f residuals.f fcn.f

updates.f: qball_rnpl
residuals.f: qball_rnpl
initializers.f: qball_rnpl
qball_init.f: qball_rnpl
fcn.o: fcn.f fcn.inc

qball: qball.o updates.o residuals.o fcn.o
	$(F77_LOAD) qball.o updates.o residuals.o fcn.o $(FLIBS) -o qball


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               P A M R 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pamr: qball-pamr

qball-pamr.o: qball-pamr.h 

qball-pamr: $(PAMR_OBJS) 
	$(CC) $(CCFLAGS) $(LDFLAGS) $(PAMR_OBJS) $(PAMR_LIBS) -o qball-pamr


########################################################################
# Clean-up
########################################################################

clean:
	/bin/rm qball-pamr
	/bin/rm *_.c 
	/bin/rm *.o 
	/bin/rm *.param
	/bin/rm qball qball_init > /dev/null 2>&1
	/bin/rm residuals.f updates.f initializers.f qball.f qball_init.f > /dev/null 2>&1
	/bin/rm gfuni0.inc globals.inc other_glbs.inc sys_param.inc > /dev/null 2>&1
	/bin/rm *.sdf
	/bin/rm *.dat
	/bin/rm .rnpl_attributes
