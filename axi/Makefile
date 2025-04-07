.IGNORE:

all:
	$(MAKE) -C src/rnpl vclean
	$(MAKE) -C src/rnpl
	$(MAKE) -C src/rnpl clean
	$(MAKE) -C src/pamr vclean
	$(MAKE) -C src/pamr

clean:
	$(MAKE) -C src/rnpl clean
	$(MAKE) -C src/pamr clean
	/bin/rm run/output/*

vclean: clean
	make clean
	$(MAKE) -C src/rnpl vclean
	$(MAKE) -C src/pamr vclean
