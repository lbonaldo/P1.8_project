#include <data_structures.h>

#ifndef UTILITIES_H
#define UTILITIES_H

/* a few physical constants */
static const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
static const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */
/* helper function: zero out an array */
void azzero(double *d, const int n);

/* helper function: apply minimum image convention */
inline double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* compute kinetic energy */
void ekin(mdsys_t *sys);

#endif
