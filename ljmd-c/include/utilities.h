#ifndef UTILITIES_H
#define UTILITIES_H

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */
/* helper function: zero out an array */
static void azzero(double *d, const int n);

/* helper function: apply minimum image convention */
static double pbc(double x, const double boxby2);

/* compute kinetic energy */
static void ekin(mdsys_t *sys);

#endif
