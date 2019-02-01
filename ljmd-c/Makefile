# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=Obj-serial/%.o)
############################################

default: split serial library

serial:
	$(MAKE) $(MFLAGS) -C Obj-$@

split:
	$(MAKE) $(MFLAGS) -C Obj-$@

library:
	$(MAKE) $(MFLAGS) -C $@

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) $(MFLAGS) -C Obj-split clean
	$(MAKE) $(MFLAGS) -C examples clean
	$(MAKE) $(MFLAGS) -C test clean
	$(MAKE) $(MFLAGS) -C library clean

check: split
	$(MAKE) $(MFLAGS) -C examples check


pretest: split
	$(MAKE) $(MFLAGS) -C test pretest

unittest: split
	$(MAKE) $(MFLAGS) -C test unittest

time: split
	$(MAKE) $(MFLAGS) -C examples time

python: library
	$(MAKE) $(MFLAGS) -C examples python

.PHONY: library
