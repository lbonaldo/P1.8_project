#include "force.h"
#include <omp.h>

/* compute forces */
void force(mdsys_t *sys) 
{

	double c12,c6,rcsq;
    c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    c6=4.0*sys->epsilon*pow(sys->sigma,6.0);
    rcsq=(sys->rcut)*(sys->rcut);
      
    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

	#pragma omp parallel for 
    for(int i=0; i < sys->natoms; ++i) {
		
	//	printf("iter %d thread %d\n",i,omp_get_thread_num());
	
		double thread_epot = 0;
    	double rsq,ffac;
    	double rx,ry,rz;
		double r6,rinv;
        for(int j=0; j < sys->natoms; ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            rsq = rx*rx + ry*ry + rz*rz;
      
            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
	        rinv = 1.0/rsq; r6 = rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                thread_epot += 0.5*(c12*r6 - c6)*r6;

                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
            }
        }
		#pragma omp critical
		sys -> epot += thread_epot;
    }
}
