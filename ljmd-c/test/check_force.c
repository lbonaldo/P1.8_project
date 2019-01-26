#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include "utilities.h"
#include "data_structures.h"
#include "input.h"
#include "force.h"

const double eps = 1e-6,inf=1e-100;


int diff(double a,double b){
	if((fabs(a)+fabs(b)) < inf) return 0;
	return fabs(a-b)/(fabs(a)+fabs(b)) > eps;
}

/* compute forces benchmark */
void force_bench(mdsys_t *sys) 
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    for(i=0; i < (sys->natoms); ++i) {
        for(j=0; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            if (i==j) continue;
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);
      
            /* compute force and energy if within cutoff */
            if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);
                
                sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));

                sys->fx[i] += rx/r*ffac;
                sys->fy[i] += ry/r*ffac;
                sys->fz[i] += rz/r*ffac;
            }
        }
    }
}

void allocate(mdsys_t *sys){
    /* allocate memory */
    sys->rx=(double *)malloc(sys->natoms*sizeof(double));
    sys->ry=(double *)malloc(sys->natoms*sizeof(double));
    sys->rz=(double *)malloc(sys->natoms*sizeof(double));
    sys->vx=(double *)malloc(sys->natoms*sizeof(double));
    sys->vy=(double *)malloc(sys->natoms*sizeof(double));
    sys->vz=(double *)malloc(sys->natoms*sizeof(double));
    sys->fx=(double *)malloc(sys->natoms*sizeof(double));
    sys->fy=(double *)malloc(sys->natoms*sizeof(double));
    sys->fz=(double *)malloc(sys->natoms*sizeof(double));
	
}
void deallocate(mdsys_t *sys){
    free(sys->rx);
    free(sys->ry);
    free(sys->rz);
    free(sys->vx);
    free(sys->vy);
    free(sys->vz);
    free(sys->fx);
    free(sys->fy);
    free(sys->fz);
}

void copy_sys(mdsys_t* dest,const mdsys_t* orig){
	dest->natoms = orig->natoms;
	dest->nfi = orig->nfi;
	dest->nsteps = orig->nsteps;
	dest->dt = orig->dt;
	dest->mass = orig->mass;
	dest->epsilon = orig->epsilon;
	dest->sigma = orig->sigma;
	dest->box = orig->box;
	dest->rcut = orig->rcut;
	dest->temp = orig->temp;
	dest->epot = orig->epot;
	dest->ekin = orig->ekin;
	
	for(int i=0;i< orig->natoms;++i){
		dest->rx[i]=orig->rx[i];
		dest->ry[i]=orig->ry[i];
		dest->rz[i]=orig->rz[i];
		dest->vx[i]=orig->vx[i];
		dest->vy[i]=orig->vy[i];
		dest->vz[i]=orig->vz[i];
		dest->fx[i]=orig->fx[i];
		dest->fy[i]=orig->fy[i];
		dest->fz[i]=orig->fz[i];
	}
}

int cmp_sys(const mdsys_t* A,const mdsys_t* B){
	if( A->natoms 	!= 	B->natoms) 	return 1;
	if( A->nfi 		!= 	B->nfi ) 	return 2;
	if( A->nsteps	!= 	B->nsteps ) return 3;
	if( diff(A->dt	,  	B->dt))		return 4;
	if( diff(A->mass,  	B->mass)) 	return 5;
	if( diff(A->epsilon,B->epsilon))return 6;
	if( diff(A->sigma , B->sigma)) 	return 7;
	if( diff(A->box, 	B->box)	)	return 8;
	if( diff(A->rcut,	B->rcut) )	return 9;
	if( diff(A->temp,	B->temp))	return 10;
	if( diff(A->epot,	B->epot))	return 11;
	if( diff(A->ekin,	B->ekin))	return 12;
	
	for(int i=0;i< A->natoms;++i){
		if(diff(A->rx[i],B->rx[i] ))return 13;
		if(diff(A->ry[i],B->ry[i] ))return 13;
		if(diff(A->rz[i],B->rz[i] ))return 13;
		if(diff(A->vx[i],B->vx[i] ))return 14;
		if(diff(A->vy[i],B->vy[i] ))return 14;
		if(diff(A->vz[i],B->vz[i] ))return 14;
		if(diff(A->fx[i],B->fx[i] ))return 15;
		if(diff(A->fy[i],B->fy[i] ))return 15;
		if(diff(A->fz[i],B->fz[i] ))return 15;
	}
	
	return 0;
}


#define LEN 200
char line[LEN],restfile[LEN];
FILE *fp;

int main()
{
	mdsys_t sysA,sysB;
	
    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sysA.natoms=atoi(line);
    sysB.natoms=sysA.natoms;
	allocate(&sysA);
	allocate(&sysB);
	
    if(get_a_line(stdin,line)) return 1;
    sysA.mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sysA.epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sysA.sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sysA.rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sysA.box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
	
    
	/* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (int i=0; i<sysA.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sysA.rx+i, sysA.ry+i, sysA.rz+i);
        }
        for (int i=0; i<sysA.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sysA.vx+i, sysA.vy+i, sysA.vz+i);
        }
        fclose(fp);
    }else {
		fprintf(stderr,"Cannot read %s\n",restfile);
		return 1;
	}
	
	copy_sys(&sysB,&sysA);
	
	force_bench(&sysA);
	force(&sysB);
	
	int what;
	if(what=cmp_sys(&sysA,&sysB)){
		printf("exit with %d\n",what);
		return 1;
	}
	deallocate(&sysA);
	deallocate(&sysB);

	return 0;
}
