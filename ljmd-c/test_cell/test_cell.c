#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stddef.h>

#include <input.h>
#include <data_structures.h>
#include <cell.h>

int main(){
    int i, j;
    char restfile[BLEN], line[BLEN];
    FILE *fp;
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

    //get cell lenght along a single direction
    if(get_a_line(stdin,line)) return 1;
    sys.clen=atof(line);

    //be sure that clen is not smaller than cutoff
    assert(sys.clen >= sys.rcut);
    /* compute number of cells */
    sys.ncells = sys.box/sys.clen;
    assert(sys.ncells>2 && "Insufficient number of cells.");
    sys.Ncells = sys.ncells * sys.ncells * sys.ncells;

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));

    /* memory for cell operations */
    sys.clist=(int **)malloc(sys.Ncells*sizeof(int*)); //Pointers to cells
    sys.plist=(int *)malloc(sys.Ncells*54*sizeof(int)); //List of pairs
    sys.catoms=(int *)malloc(sys.Ncells*sizeof(int)); //atoms per cell

    sys.npair = sys.Ncells * 27;
    printf("%ld\n", sizeof(sys.plist));
    printf("Old npair: %d\n",sys.npair);
    /* create the pair list */
    allocate_pairs(&sys);
    printf("New npair: %d\n",sys.npair);
    assert(sys.Ncells*27 >= sys.npair);
    printf("%ld\n", sizeof(sys.plist));
    assert(sizeof(sys.plist)*0.5/sizeof(int)==sys.npair);

    /* initialize clist and catoms*/
    for(i=0; i<sys.Ncells; i++){
      sys.clist[i]=(int *)malloc(sizeof(int));
      sys.catoms[i]=0;
    }
    
    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        fclose(fp);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* allocate the atoms inside the cells */
    allocate_cells(&sys);

    /* test that allocation works */
    /* check 1: cunt the correct total number of atoms */
    /* check 2: sum all atom indices and get (natoms-1)*natoms/2 */
    /* i.e. Gauss formula for integers sum. Recall that atom idx start from 0! */
    
    int sum1=0, sum2=0;
    for(i=0; i<sys.Ncells; i++){
      sum1 += sys.catoms[i];
      for(j=0; j<sys.catoms[i]; j++)
	sum2 += sys.clist[i][j];
    }
    assert(sum1 == sys.natoms && "wrong total number of atoms");
    assert(sum2 == (sys.natoms-1)*(sys.natoms)/2 && "wrong atom indices allocated");

    printf("Test completed successfully.\n");
    
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);

    free(sys.plist);
    free(sys.catoms);
    for(i=0; i<sys.Ncells; i++)
      free(sys.clist[i]);
    free(sys.clist);

    return 0;
}
