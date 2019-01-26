#include "output.h"
#include "input.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double eps = 1e-6;

int diff(double a,double b){
	return fabs(a-b)/(fabs(a)+fabs(b)) > eps;
}

int main(){
	
    int nprint,natoms,nsteps;
	double mass,epsilon,sigma,rcut,box,dt;
    char restfile[BLEN], test_trajfile[BLEN], test_ergfile[BLEN], line[BLEN], 
			ergfile[BLEN], trajfile[BLEN];
    FILE *fp,*traj,*erg,*test_traj,*test_erg;

    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,test_trajfile)) return 1;
    if(get_a_line(stdin,test_ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;

    erg=fopen(ergfile,"r");
    traj=fopen(trajfile,"r");
    test_erg=fopen(test_ergfile,"r");
    test_traj=fopen(test_trajfile,"r");
	
	for(int i=0;i<=nsteps;++i)if(i % nprint == 0){
		int ia,ib;
		double da,db;
		fscanf(erg,"%d",&ia);
		fscanf(test_erg,"%d",&ib);
		if(ia!=ib) return 1;
		
		for(int j=0;j<4;++j){
			fscanf(erg,"%lf",&da);
			fscanf(test_erg,"%lf",&db);
			if(diff(da,db)) return 1;
		}
	}
	

	fclose(erg);
	fclose(traj);
	fclose(test_erg);
	fclose(test_traj);

	return 0;
}
