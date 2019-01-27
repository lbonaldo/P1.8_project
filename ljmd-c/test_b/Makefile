# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
ATOMS=THREE
CFLAGS=-I../include -Wall -g -O3 -D $(ATOMS)
LDLIBS=-lm

# list of source files
SRC=test_b.c input.c verlet1.c verlet2.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)
MAIN=test_b.x

default: $(MAIN)

clean:
	rm -f *.mod *.o $(MAIN)

# linker rule
$(MAIN): $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
%.o: ../include/%.h
