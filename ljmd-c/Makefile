# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=Obj-serial/%.o)
############################################

default: split

serial:
	$(MAKE) $(MFLAGS) -C Obj-$@

split:
	$(MAKE) $(MFLAGS) -C Obj-$@

clean:
	$(MAKE) $(MFLAGS) -C Obj-serial clean
	$(MAKE) $(MFLAGS) -C Obj-split clean
	$(MAKE) $(MFLAGS) -C examples clean
	$(MAKE) $(MFLAGS) -C test clean

check: split
	$(MAKE) $(MFLAGS) -C examples check
	
pretest: split
	$(MAKE) $(MFLAGS) -C test pretest
