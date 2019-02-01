exe=../ljmd-serial.x
OBJDIR=../Obj-split
SRC=test_a.c force.c input.c utilities.c
VPATH=../src
OBJ=$(OBJDIR)/force.o $(OBJDIR)/input.o $(OBJDIR)/utilities.o $(OBJDIR)/verlet1.o $(OBJDIR)/verlet2.o
CC=gcc
CFLAGS=-I../include -fopenmp -Wall -O3 -ffast-math -fomit-frame-pointer -std=gnu11
LDLIBS=-lm
checklist=check_overall_energy.x check_force.x test_a.x test_b.x

default: pretest

### Run the pretests
pretest: main generate $(checklist) test_a.x test_b.x
	./run_test.sh

unittest: main generate $(checklist) test_a.x test_b.x
	python3 test_interface.py

### The program must be built in order to run tests
main:
	$(MAKE) $(MFLAGS) -C ..  


### Program that compares the output energy files 
check_overall_energy.x: check_overall_energy.c $(OBJ)
	$(CC) $^ -o $@ $(CFLAGS) $(LDLIBS)

### Chekc if force.c computes correctly.
check_force.x: check_force.c $(OBJ)
	$(CC) $^ -o $@ $(CFLAGS) $(LDLIBS)

# linker rule
test_a.x: test_a.c $(OBJ)
	$(CC) $^ -o $@ $(CFLAGS) $(LDLIBS) 

test_b.x: test_b.c $(OBJ)
	$(CC) $^ -o $@ $(CFLAGS) $(LDLIBS) -D THREE

test_b.x: test_b.h

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
%.o: ../include/%.h

clean:
	rm -f *.x *.test *.dat *.xyz *.o *.mod

### generate reference files
generate: main 
	$(exe) < argon_108.inp.gen > /dev/null

.PHONY: pretest clean generate main unittest



