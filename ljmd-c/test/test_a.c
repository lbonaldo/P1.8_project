/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#include <input.h>
#include <utilities.h>
#include <force.h>
#include <data_structures.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/*What does the code test?*/
//if 2-particles are beyond cutoff, checks all the force components are null;
//if 2-particles are within cutoff, chechs:
//1. at least one force component is non-zero;
//2. forces are one opposed to the other;
//3. action is repulsive/attractive in accordance with (distance-potential minimum radius).
//Zero's are approximated by tolerance level 1e-8;

/* main */
int main(int argc, char **argv)
{
    int i, j;
    char line[BLEN];
    mdsys_t sys;
    double rx, ry, rz, r;
    double vmin, eps;
    double* fa[3]; //array of forces
    double sum[3], sp; //for vector sums and scalar products

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

    //set tolerance
    eps=1e-8;
    //minimum of the 2-partciles potential
    vmin=pow(2,1/6)*sys.sigma; 

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*sizeof(double));

    //Initialize particles of pbc{out} type
    sys.rx[0]=6.0; sys.ry[0]=-10.0; sys.rz[0]=12.0;
    sys.rx[1]=1.0; sys.ry[1]=-3.0; sys.rz[1]=-2.0;

    //check pbc
    rx=pbc(sys.rx[1] - sys.rx[0], 0.5*sys.box,sys.box);
    ry=pbc(sys.rx[1] - sys.rx[0], 0.5*sys.box,sys.box);
    rz=pbc(sys.rx[1] - sys.rx[0], 0.5*sys.box,sys.box);
    r=sqrt(rx*rx + ry*rx + rz*rz);
    
    //compute force
    force(&sys);
    
    //assert condition 1:
    //forces must be all zero
    fa[0]=sys.fx; fa[1]=sys.fy; fa[2]=sys.fz;
    for(i=0; i<3; i++)
      for(j=0; j<sys.natoms; j++)
	assert(fabs(fa[i][j])<eps);
	  
    //Initialize particles of pbc{in} type
    sys.rx[0]=0.0; sys.ry[0]=-10.0; sys.rz[0]=12.0;
    sys.rx[1]=1.0; sys.ry[1]=-3.0; sys.rz[1]=-2.0;

    //check pbc
    rx=pbc(sys.rx[1] - sys.rx[0], 0.5*sys.box,sys.box);
    ry=pbc(sys.rx[1] - sys.rx[0], 0.5*sys.box,sys.box);
    rz=pbc(sys.rx[1] - sys.rx[0], 0.5*sys.box,sys.box);
    r=sqrt(rx*rx + ry*rx + rz*rz);
    
    //compute force
    force(&sys);

    //assert condition 2:
    //at least one force component must be nonzero
    fa[0]=sys.fx; fa[1]=sys.fy; fa[2]=sys.fz;
    for(j=0; j<sys.natoms; j++)
    	assert(fabs(fa[0][j])>eps || fabs(fa[1][j])>eps || fabs(fa[2][j])>eps);

    //assert condition 3:
    //forces must have opposite sign
    for(i=0; i<3; i++){
      sum[i]=fabs(fa[i][0]+fa[i][1]);
      assert(sum[i]<eps);
    }
    //asser condition 4:
    //force must be attractive/repulsive
    //depending on the relative position between
    //mutual distance and potential minimum
    sp=sys.fx[1]*rx + sys.fy[1]*ry + sys.fz[1]*rz;
    if(r<vmin)
      assert(sp<0.0);
    else
      assert(sp>0.0);
      
    //exit
    printf("Test done succesfully.\n");

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
    
    return 0;
}
