exe=../ljmd-serial.x
OBJDIR=../Obj-split
OBJS=$(OBJDIR)/input.o $(OBJDIR)/output.o
CFLAGS=-I../include -lm
checklist=check_overall_energy.x check_another_thing.x check_force.x

### Run the pretests
pretest: main generate $(checklist) 
	./run_test.sh

### The program must be built in order to run tests
main:
	$(MAKE) $(MFLAGS) -C ..  


### Program that compares the output energy files 
check_overall_energy.x: check_overall_energy.c $(OBJDIR)/input.o 
	gcc -o $@ $^ $(CFLAGS)

###
check_another_thing.x: check_another_thing.c
	gcc -o $@ $^ $(CFLAGS)
	
### Chekc if force.c computes correctly.
check_force.x: check_force.c $(OBJDIR)/force.o $(OBJDIR)/input.o $(OBJDIR)/utilities.o
	gcc -o $@ $^ $(CFLAGS)

clean:
	rm *.x *.test *.dat *.xyz

### generate reference files
generate: main 
	$(exe) < argon_108.inp.gen >/dev/null

.PHONY: pretest clean generate main
	

	
