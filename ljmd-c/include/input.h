#include <stdio.h>
#include <string.h>
#include <ctype.h>


#ifndef INPUT_H
#define INPUT_H

/* generic file- or pathname buffer length */
#define BLEN 200

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE *fp, char *buf);

#endif
