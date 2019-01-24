default:
	$(MAKE) $(MFLAGS) -C ljmd-c 
	
check:
	$(MAKE) $(MFLAGS) -C ljmd-c check
	
.PHONY: test default

