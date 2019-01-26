exe=../ljmd-split.x
OBJDIR=../Obj-split
OBJS=$(OBJDIR)/input.o $(OBJDIR)/output.o
CFLAGS=-I../include
checklist=check_overall_energy.x check_another_thing.x

pretest: $(exe) $(checklist) 
	./run_test.sh

$(exe):
	$(MAKE) $(MFLAGS) -C ..  

check_overall_energy.x: check_overall_energy.c $(OBJS)
	gcc -o $@ $^ $(CFLAGS)
	
check_another_thing.x: check_another_thing.c
	gcc -o $@ $^ $(CFLAGS)

clean:
	rm *.x

.PHONY: pretest clean
	

	
