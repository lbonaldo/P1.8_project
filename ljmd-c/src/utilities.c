#include "utilities.h"

/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* 	helper function: apply minimum image convention 
	boxby2: half's box length
*/
double pbc(double x, const double boxby2,const double box)
{
    while (x >  boxby2) x -= box;
    while (x < -boxby2) x += box;
    return x;
}

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{   
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

/* 	bring all particles inside the origin box, 
	so pbc's while loop takes at
	most one step */
void normalize_positions(mdsys_t* sys){
	double box = sys->box, half_box = box * 0.5;
	for(int i=0;i<sys->natoms;++i){
		sys->rx[i] = half_box + pbc(sys->rx[i],half_box,box);
		sys->ry[i] = half_box + pbc(sys->ry[i],half_box,box);
		sys->rz[i] = half_box + pbc(sys->rz[i],half_box,box);
	}	
	return ;
}
