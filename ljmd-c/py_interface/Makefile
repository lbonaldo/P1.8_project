# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-I../include -Wall -g -ggdb -O3 -ffast-math -fomit-frame-pointer -fPIC
LDLIBS=-lm -shared

# list of source files
SRC=force.c utilities.c verlet1.c verlet2.c input.c output.c
VPATH=.
OBJ=$(SRC:%.c=%.o)
MAIN=sysdlib.so

default: $(MAIN)

clean:
	rm -f *.mod *.o $(MAIN) *.xyz *.dat

# linker rule
$(MAIN): $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
%.o: ../include/%.h
