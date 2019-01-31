exe=../ljmd-serial.x
OBJDIR=../Obj-split
SRC=test_a.c force.c input.c utilities.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)
CC=gcc
CFLAGS=-I../include -fopenmp -Wall -O3 -ffast-math -fomit-frame-pointer -std=gnu11
LDLIBS=-lm
checklist=check_overall_energy.x check_another_thing.x check_force.x

default: pretest

### Run the pretests
pretest: main generate $(checklist) test_a.x
	./run_test.sh

### The program must be built in order to run tests
main:
	$(MAKE) $(MFLAGS) -C ..  


### Program that compares the output energy files 
check_overall_energy.x: check_overall_energy.c $(OBJDIR)/input.o 
	$(CC) -o $@ $^ $(CFLAGS) $(LDLIBS)

###
check_another_thing.x: check_another_thing.c
	$(CC) -o $@ $^ $(CFLAGS)
	
### Chekc if force.c computes correctly.
check_force.x: check_force.c $(OBJDIR)/force.o $(OBJDIR)/input.o $(OBJDIR)/utilities.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDLIBS)

# linker rule
test_a.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
%.o: ../include/%.h
	
clean:
	rm *.x *.test *.dat *.xyz *.o *.mod

### generate reference files
generate: main 
	$(exe) < argon_108.inp.gen > /dev/null 

.PHONY: pretest clean generate main
	


