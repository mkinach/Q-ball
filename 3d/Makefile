.IGNORE:

all:
	$(MAKE) -C src/fd   clean
	$(MAKE) -C src/pamr clean
	$(MAKE) -C src/fd
	$(MAKE) -C src/pamr

clean:
	$(MAKE) -C src/fd   clean
	$(MAKE) -C src/pamr clean
	/bin/rm run/output/*

vclean: clean
	make clean
	$(MAKE) -C src/fd   vclean
	$(MAKE) -C src/pamr vclean
