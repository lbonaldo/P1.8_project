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

    int i, j;

#pragma omp parallel for collapse(2) private(j)
    for(i=0; i < sys->Ncells; ++i){
      for(j=0; j < sys->natoms; ++j) {

	if(j>=sys->catoms[i]) continue;
	
	int start = 54*i;
	int c1 = sys->plist[start];
       
	int ii = sys->clist[c1][j];
	double thread_epot = 0,fxi=0,fyi=0,fzi=0;
    	double rsq,ffac;
    	double rx,ry,rz,	
	  rxi=sys->rx[ii],
	  ryi=sys->ry[ii],
	  rzi=sys->rz[ii];
	double r6,rinv;
	double box=sys->box,boxby2=0.5*box;
	
	for(int n=0; n<27; n++){ //close this loop
	  int c2 = sys->plist[start + (2*n + 1)];
       
	  for(int k=0; k<sys->catoms[c2]; k++){
	    int jj = sys->clist[c2][k];

            /* particles have no interactions with themselves */
            if (ii==jj) continue;
            
            /* get distance between particle i and j */
            rx=pbc(rxi - sys->rx[jj], boxby2,box);
            ry=pbc(ryi - sys->ry[jj], boxby2,box);
            rz=pbc(rzi - sys->rz[jj], boxby2,box);
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
	}
	sys->fx[ii]+=fxi;
	sys->fy[ii]+=fyi;
	sys->fz[ii]+=fzi;		
#pragma omp critical
	sys -> epot += thread_epot*0.5;
      }
    }
}
