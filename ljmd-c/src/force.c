#include "force.h"

/* compute forces */
void force(mdsys_t *sys) 
{
    double rsq,rcsq,ffac;
    double rx,ry,rz;
    double c12, c6;
    double r6, rinv;
    int i,j;

    //cell variables
    int c1, c2, ii, jj, k;

    c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    c6=4.0*sys->epsilon*pow(sys->sigma,6.0);
    rcsq=(sys->rcut)*(sys->rcut);
      
    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < sys->npair; ++i) {
      c1 = sys->plist[2*i];
      c2 = sys->plist[2*i+1];
      
      for(j=0; j < sys->catoms[c1]; ++j) {
	ii = sys->clist[c1][j];

	for(k=0; k<sys->catoms[c2]; k++){
	  jj = sys->clist[c2][k];

            /* particles have no interactions with themselves */
	    /* and Newton's 3rd law */
            if (ii>=jj) continue;
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[ii] - sys->rx[jj], 0.5*sys->box);
            ry=pbc(sys->ry[ii] - sys->ry[jj], 0.5*sys->box);
            rz=pbc(sys->rz[ii] - sys->rz[jj], 0.5*sys->box);
            rsq = rx*rx + ry*ry + rz*rz;
      
            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
	        rinv = 1.0/rsq; r6 = rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                sys->epot += (c12*r6 - c6)*r6;

                sys->fx[ii] += rx*ffac;
                sys->fy[ii] += ry*ffac;
                sys->fz[ii] += rz*ffac;

		sys->fx[jj] -= rx*ffac;
                sys->fy[jj] -= ry*ffac;
                sys->fz[jj] -= rz*ffac;
	    }
	}
      }
    }
}
