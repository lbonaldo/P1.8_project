#include <data_structures.h>

#ifndef CELL_H
#define CELL_H

/* to allocate the list of pairs */
void allocate_pairs(mdsys_t *sys);

/* to allocate atoms positions inside the cells*/
void allocate_cells(mdsys_t *sys);

/* to clean the cells */
void empty_cells(mdsys_t *sys);

#endif
