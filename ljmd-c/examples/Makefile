# -*- Makefile -*-

#EXE = ljmd-serial.x
EXE = ljmd-split.x

default: check

clean:
	rm -f *.xyz *.dat

check: ../$(EXE)
	../$(EXE) < argon_108.inp
	head -10 argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 1
