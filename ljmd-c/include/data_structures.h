#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

/* structure to store informations
 * cells */
struct _cell {
  int natoms;  //total number of atoms inside the cell
  int *idxlist; //atom indices belonging to the cell
};
typedef struct _cell cell_t;

/* structure to hold the complete information 
 * about the MD system */
struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;

    int ncells; //ncells along a SINGLE direction
    cell_t **clist;
    int *plist;
    double clen;
};
typedef struct _mdsys mdsys_t;

#endif
