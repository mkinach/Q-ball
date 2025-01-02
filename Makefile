all:
	$(MAKE) -C src/fd   clean
	$(MAKE) -C src/pamr clean
	$(MAKE) -C src/fd
	$(MAKE) -C src/pamr

clean:
	$(MAKE) -C src/fd   clean
	$(MAKE) -C src/pamr clean
	/bin/rm -f run/output/*

vclean: clean
	$(MAKE) -C src/fd   vclean
	$(MAKE) -C src/pamr vclean
	/bin/rm -f run/output/*
	/bin/rm -f run/input/initdata*
