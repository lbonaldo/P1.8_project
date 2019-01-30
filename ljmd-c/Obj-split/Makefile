# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-I../include -Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

# list of source files
SRC=main.c force.c input.c output.c utilities.c verlet1.c verlet2.c cell.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)
MAIN=../ljmd-split.x

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
