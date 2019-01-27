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
	
		double thread_epot = 0,fxi=0,fyi=0,fzi=0;
    	double rsq,ffac;
    	double rx,ry,rz,	
			rxi=sys->rx[i],
			ryi=sys->ry[i],
			rzi=sys->rz[i];
		double r6,rinv;
		double box=sys->box,boxby2=0.5*box;
		
        for(int j=0; j < sys->natoms; ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;
            
            /* get distance between particle i and j */
            rx=pbc(rxi - sys->rx[j], boxby2,box);
            ry=pbc(ryi - sys->ry[j], boxby2,box);
            rz=pbc(rzi - sys->rz[j], boxby2,box);
            rsq = rx*rx + ry*ry + rz*rz;
      
            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
	        	rinv = 1.0/rsq; r6 = rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                thread_epot += (c12*r6 - c6)*r6;

                fxi += rx*ffac;
                fyi += ry*ffac;
                fzi += rz*ffac;
            }
        }
		sys->fx[i]+=fxi;
		sys->fy[i]+=fyi;
		sys->fz[i]+=fzi;
		
		#pragma omp critical
		sys -> epot += thread_epot*0.5;
    }
}
