#include <data_structures.h>

#ifndef UTILITIES_H
#define UTILITIES_H

/* a few physical constants */
static const double Na=6.022140857e23; /* Avogadro's constant, number of atoms in a mol  */
static const double jtokcal=0.00023900573614; /* 1 J/kcal  */
static const double kb_SI=1.38064852e-23; /* Boltzmann constant in international units J/K  */
static const double amu_kg=1.660539040e-27; /* 1 atomic mass in kg  */
static const double A_m=1e-10; /* 1 armstrong in m */

static const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K = kb_SI*Na*jtokcal */
/*
	m*v^2 	= m_{kg} * v_{m/s}^2 * kg * (m/s)^2
			= m_{kg} * v_{m/s}^2 * J
			= m_{amu} * v_{m/s}^2 * amu/kg * J
			= m_{amu} * v_{m/s}^2 * amu/kg * J/kcal * kcal
			= m_{amu} * v_{m/s}^2 * amu/kg * J/kcal * Na * kcal/mol
			= m_{amu} * v_{m/s}^2 * 2.390057e-7 * kcal/mol
			= m_{amu} * v_{A/T}^2 * ( (A/T)/(m/s) )^2 * 2.390057e-7 * kcal/mol
			= m_{amu} * v_{A/T}^2 * (s/T * 1e5 )^2 * 2390.057 * kcal/mol
			
	where A is armstrong and T is some unknown unit of time;
	then if Axel insists that mvsq2e=2390.057
	we must have that the velocities are given in units of A/T
	where T = 1e-5 seconds.
*/
static const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */


/* helper function: zero out an array */
void azzero(double *d, const int n);

/* helper function: apply minimum image convention */
inline double pbc(double x, const double boxby2,const double box)
{
    while (x >  boxby2) x -= box;
    while (x < -boxby2) x += box;
    return x;
}

/* compute kinetic energy */
void ekin(mdsys_t *sys);

/* 	bring all particles inside the origin box, 
	so pbc's while loop takes at
	most one step */
void normalize_positions(mdsys_t* sys);

#endif
