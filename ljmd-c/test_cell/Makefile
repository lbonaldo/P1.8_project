# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-I../include -Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

# list of source files
SRC=test_cell.c input.c utilities.c cell.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)
MAIN=test_cell.x

default: $(MAIN)

clean:
	rm -f *.mod *.o $(MAIN)

cell: test_cell.x
	./test_cell.x < test_cell.inp 

# linker rule
$(MAIN): $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
%.o: ../include/%.h
