# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=Obj-serial/%.o)
############################################

default: serial split

serial:
	$(MAKE) $(MFLAGS) -C Obj-$@

split:
	$(MAKE) $(MFLAGS) -C Obj-$@

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) $(MFLAGS) -C Obj-split clean
	$(MAKE) $(MFLAGS) -C examples clean
	$(MAKE) $(MFLAGS) -C test clean
	$(MAKE) $(MFLAGS) -C test_cell clean

check: serial split
	$(MAKE) $(MFLAGS) -C examples check

test_a:
	$(MAKE) $(MFLAGS) -C test test_a

cell:
	$(MAKE) $(MFLAGS) -C test_cell cell

time: split
	$(MAKE) $(MFLAGS) -C examples time
