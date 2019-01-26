#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <data_structures.h>
#include "verlet1.h"
#include "verlet2.h"

void init_pos(mdsys_t, double *);
void init_vel(mdsys_t, double *);
void init_for(mdsys_t, double *);

void print_pos(const mdsys_t);
void print_vel(const mdsys_t);
void print_for(const mdsys_t);

int test(mdsys_t *sys, const int test_n);

void init_pos(mdsys_t *sys, double * pos_value){
  
  for (int i=0; i<sys.natoms; ++i) {
    sys.rx+i = pos_value[3*i];
    sys.ry+i = pos_value[3*i+1];
    sys.rz+i = pos_value[3*i+2];
  }
}

void init_vel(mdsys_t *sys, double * vel_value){
  
  for (int i=0; i<sys.natoms; ++i) {
    sys.vx+i = vel_value[3*i];
    sys.vy+i = vel_value[3*i+1];
    sys.vz+i = vel_value[3*i+2];
  }
}

void init_for(mdsys_t *sys, double * for_value){
  
  for (int i=0; i<sys.natoms; ++i) {
    sys.fx+i = for_value[3*i];
    sys.fy+i = for_value[3*i+1];
    sys.fz+i = for_value[3*i+2];
  }
}

void print_pos(const mdsys_t *sys)
{   
  for (int i=0; i<sys->natoms; ++i)
    fprintf(stdout, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
}

void print_vel(const mdsys_t *sys)
{   
  for (int i=0; i<sys->natoms; ++i)
    fprintf(stdout, "Ar  %20.8f %20.8f %20.8f\n", sys->vx[i], sys->vy[i], sys->vz[i]);
}

void print_for(const mdsys_t *sys)
{
  for (int i=0; i<sys->natoms; ++i)
    fprintf(stdout, "Ar  %20.8f %20.8f %20.8f\n", sys->fx[i], sys->fy[i], sys->fz[i]);
}

int test(mdsys_t *sys, const int test_n, const int eps){

  int status = 10;
  
  switch (test_n) {
    
    //test 1: v(t)=0, x(t)=0
  case 1:
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
      
      /* write output, if requested */
      if ((sys.nfi % nprint) == 0){
	print_pos(&sys);
	print_vel(&sys);
	print_for(&sys);
      }
      /* propagate system and recompute energies */
      velverlet1(&sys);
      
      if(
	 fabs(sys->vx[i] == 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass) < eps &&
	 fabs(sys->vy[i] == 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass) < eps &&
	 fabs(sys->vz[i] == 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass) < eps &&
	 fabs(sys->rx[i] == sys->dt*sys->vx[i]) < eps &&
	 fabs(sys->ry[i] == sys->dt*sys->vy[i]) < esp &&
	 fabs(sys->rz[i] == sys->dt*sys->vz[i]) < eps )
	printf("velverlet1 PASSED.\n");
      else
	printf("velverlet1 FAILED.\n");
      
      velverlet2(&sys);

      if(
	 fabs(sys->vx[i] - sys->dt / mvsq2e * sys->fx[i] / sys->mass) < eps &&
	 fabs(sys->vy[i] - sys->dt / mvsq2e * sys->fy[i] / sys->mass) < eps &&
	 fabs(sys->vz[i] - sys->dt / mvsq2e * sys->fz[i] / sys->mass) < eps )
	printf("velverlet2 PASSED.\n");
      else
	printf("velverlet2 FAILED.\n");
    }
    break;
    
    //test 2: v(t)=0, F(t)=0
  case 2:
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
      
      /* write output, if requested */
      if ((sys.nfi % nprint) == 0){
	print_pos(&sys);
	print_vel(&sys);
	print_for(&sys);
      }
      /* propagate system and recompute energies */
      velverlet1(&sys);

      if(
	 sys->vx[i] < eps &&
	 sys->vy[i] < eps &&
	 sys->vz[i] < eps &&
	 fabs(sys->rx[i] - pos_value[3*i]) < eps &&
	 fabs(sys->ry[i] - pos_value[3*i+1]) < eps &&
	 fabs(sys->rz[i] - pos_value[3*i+2]) < eps )
	printf("velverlet1 PASSED.\n");
      else
	printf("velverlet1 FAILED.\n");
      
      velverlet2(&sys);

      if(
	 sys->vx[i] < eps &&
	 sys->vy[i] < eps &&
	 sys->vz[i] < eps )
	printf("velverlet2 PASSED.\n");
      else
	printf("velverlet2 FAILED.\n");      
    }
    break;
    
    //test 3: a(t)=-0.5v(t)/dt
  case 3:
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
      
      /* write output, if requested */
      if ((sys.nfi % nprint) == 0){
	print_pos(&sys);
	print_vel(&sys);
	print_for(&sys);
      }
      /* propagate system and recompute energies */
      velverlet1(&sys);
      
      do_test_v1(&sys);

      velverlet2(&sys);

      do_test_v2(&sys);
    }
    break;
    
    //test 4: a(t)>-0.5v(t)/dt
  case 4:
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
      
      /* write output, if requested */
      if ((sys.nfi % nprint) == 0){
	print_pos(&sys);
	print_vel(&sys);
	print_for(&sys);
      }
      /* propagate system and recompute energies */
      velverlet1(&sys);
      do_test_v1(&sys);
      velverlet2(&sys);
      do_test_v2(&sys);
    }
    break;
    
    //test 5: a(t)<-0.5v(t)/dt
  case 5:
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
      
    /* write output, if requested */
      if ((sys.nfi % nprint) == 0){
	print_pos(&sys);
	print_vel(&sys);
	print_for(&sys);
      }
      /* propagate system and recompute energies */
      velverlet1(&sys);
      do_test_v1(&sys);
      velverlet2(&sys);
      do_test_v2(&sys);
    }
    break;
  }
}
