default:
	$(MAKE) $(MFLAGS) -C ljmd-c 
	
check:
	$(MAKE) $(MFLAGS) -C ljmd-c check

pretest:
	$(MAKE) $(MFLAGS) -C ljmd-c pretest

clean:
	$(MAKE) $(MFLAGS) -C ljmd-c clean

.PHONY: check pretest default clean

