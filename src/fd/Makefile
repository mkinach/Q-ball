########################################################################
#
# Makefile for Maple language file 'gencode.mpl'
#
# Uses the FD toolkit (https://github.com/rmanak/FD/) and assumes that
# the binary 'maple' is available in your PATH
#
########################################################################

all:
	maple gencode.mpl

clean:
	/bin/rm -f *.h *_call > /dev/null 2>&1 
	/bin/rm -f *.o > /dev/null 2>&1 

vclean:
	/bin/rm -f *.f > /dev/null 2>&1 
	make clean
