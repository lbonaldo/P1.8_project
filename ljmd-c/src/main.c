/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <input.h>
#include <utilities.h>
#include <output.h>
#include <force.h>
#include <data_structures.h>
#include <verlet1.h>
#include <verlet2.h>
#include <cell.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stddef.h>

/* main */
int main(int argc, char **argv) 
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;

    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);
    
    if(get_a_line(stdin,line)) return 1;
    sys.clen=atof(line);

    assert(sys.clen >= sys.rcut);
    /* compute number of cells */
    sys.ncells = sys.box/sys.clen;
    assert(sys.ncells>2 && "Insufficient number of cells.");
    sys.Ncells = sys.ncells * sys.ncells * sys.ncells;

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));    
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    /* memory for cell operations */
    sys.clist=(int **)malloc(sys.Ncells*sizeof(int*)); //Pointers to cells
    sys.plist=(int *)malloc(sys.Ncells*54*sizeof(int)); //List of pairs
    sys.catoms=(int *)malloc(sys.Ncells*sizeof(int)); //atoms per cell

    /* initialize clist and catoms*/
    for(i=0; i<sys.Ncells; i++){
      sys.clist[i]=(int *)malloc(sizeof(int));
      sys.catoms[i]=0;
    }
    
    sys.npair = sys.Ncells * 27;
    /* create the pair list */
    allocate_pairs(&sys);
    
    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
        azzero(sys.fx, sys.natoms);
        azzero(sys.fy, sys.natoms);
        azzero(sys.fz, sys.natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* allocate particles inside the cells */
    allocate_cells(&sys);
    
    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);
    ekin(&sys);
    
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet1(&sys);

	/* allocate the atoms inside the cells */
	allocate_cells(&sys);

	force(&sys);
	velverlet2(&sys);
        ekin(&sys);
	
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    fclose(erg);
    fclose(traj);

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    free(sys.plist);
    free(sys.catoms);
    for(i=0; i<sys.Ncells; i++)
      free(sys.clist[i]);
    free(sys.clist);
    
    return 0;
}
