# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-I ../include -Wall -g -O3 -ffast-math -fomit-frame-pointer
LDLIBS=-lm

# list of source files
SRC=main.c force.c input.c output.c utilities.c verlet.c
VPATH=../src
OBJ=$(SRC:%.c=%.o)

default: ../main.x

clean:
	rm -f *.mod *.o ../main.x

# linker rule
../main.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compilation pattern rule for objects
%.o: %.c
	$(CC) -c $(CFLAGS) $<

# dependencies
main.o: main.c
%.o: ../include/%.h
