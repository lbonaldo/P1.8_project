#include <stdio.h>
#include <string.h>
#include <ctype.h>

/* generic file- or pathname buffer length */
#define BLEN 200

#ifndef INPUT_H
#define INPUT_H

/* helper function: read a line and then return
   the first string with whitespace stripped off */
static int get_a_line(FILE *fp, char *buf);

#endif
