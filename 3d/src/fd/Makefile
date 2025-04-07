########################################################################
#
# Makefile for Maple language file 'gencode.mpl'
#
# Uses the FD toolkit (https://github.com/rmanak/FD/) and assumes that
# the binary 'maple' is available in your PATH
#
########################################################################
.IGNORE:

all:
	maple gencode.mpl

clean:
	/bin/rm *.h *_call > /dev/null 2>&1
	/bin/rm *.o > /dev/null 2>&1

vclean:
	make clean
	/bin/rm *.f > /dev/null 2>&1
