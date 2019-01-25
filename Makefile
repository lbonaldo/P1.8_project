default:
	$(MAKE) $(MFLAGS) -C ljmd-c 
	
check:
	$(MAKE) $(MFLAGS) -C ljmd-c check

pretest:
	$(MAKE) $(MFLAGS) -C ljmd-c pretest

.PHONY: check pretest default

