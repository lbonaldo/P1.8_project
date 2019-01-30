# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-I../include -Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

# list of source files
SRC=test_a.c force.c input.c utilities.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)
MAIN=test_a.x

default: $(MAIN)

clean:
	rm -f *.mod *.o $(MAIN)

test_a: test_a.x
	./test_a.x < test_a.inp 

# linker rule
$(MAIN): $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
%.o: ../include/%.h
