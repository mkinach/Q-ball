########################################################################
#
# Makefile for PAMR application 'qball'
#
# Uses mpicc and mpifort compilers (run 'mpifort --showme' or 
# 'mpicc --showme' to see what they reduce to)
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
FD_DIR := ../fd

# if you want non-default flags for the compilers, set them here
CC     := mpicc -O3
CFLAGS := -I/usr/local/include -I/usr/include/X11
F77		 := mpifort -O6 -fno-second-underscore

## uncomment this section to enable debug flags
##
## Fortran:
## -ggdb: produce debugging info for GDB
## -Wall: verbose warnings at compile time
## -fcheck=bounds: runtime check array bounds
## -ffpe-trap=underflow,denormal,overflow,invalid: stop execution on FPE
## -finit-real=snan: check uninitialized variables
##
## C:
## -ggdb: produce debugging info for GDB
## -Wall: verbose warnings at compile time
## -Wfloat-equal: warning if floating point numbers are used in equality comparisons
#FDEBUGFLAGS  := -ggdb -fcheck=bounds -ffpe-trap=overflow,invalid,zero -finit-real=snan
##FDEBUGFLAGS := $(FDEBUGFLAGS) -ffpe-trap=underflow,denormal  # usually not a meaningful error
##FDEBUGFLAGS := $(FDEBUGFLAGS) -Wall
#F77         := mpifort $(FDEBUGFLAGS) -fno-second-underscore -Wno-argument-mismatch
##F77         := mpifort $(FDEBUGFLAGS) -fno-second-underscore -fallow-argument-mismatch
##F77         := mpifort $(FDEBUGFLAGS) -fno-second-underscore
#CDEBUGFLAGS := -ggdb 
##CDEBUGFLAGS := $(CDEBUGFLAGS) -Wall -Wfloat-equal
#CC          := mpicc $(CDEBUGFLAGS)

########################################################################
# Main targets
########################################################################

# default targets
all: qball

# rename instances in this makefile (see 'APP' above)
fix: Makefile
	  sed "s@qball@$(APP)@g" < Makefile > .Makefile
	  mv .Makefile Makefile

# define a default rule to automatically make .o files out of .f files
%.o: %.f
	$(F77) -I$(FD_DIR) -c $< -o $@

# define a default rule to automatically make .o files out of .c files
.c.o:
	$(CC) $(CFLAGS) -c $*.c

# define relevant variables for PAMR compilation
FD_OBJS   := $(patsubst $(FD_DIR)/%.f, $(FD_DIR)/%.o, $(wildcard $(FD_DIR)/*.f))
PAMR_LIBS := -lpamr -lamrd -lbbhutil -lm -lmpi \
					   -lodepack -lvutil -llinpack -llapack -lblas -lgfortran
PAMR_OBJS = qball.o interp.o init_qball.o Qtotcalc.o Etotcalc.o \
						num.o readwrite.o $(FD_OBJS)

qball: $(PAMR_OBJS) $(wildcard *.h) $(wildcard *.inc) $(wildcard *.mod)
	$(CC) $(PAMR_OBJS) $(LDFLAGS) $(PAMR_LIBS) -o qball

########################################################################
# Clean targets
########################################################################

clean:
	/bin/rm -f *.o *_.c *.mod > /dev/null 2>&1 
	/bin/rm -f *.sdf > /dev/null 2>&1 
	/bin/rm -f qball > /dev/null 2>&1
	/bin/rm -f core* > /dev/null 2>&1

vclean:
	make clean
