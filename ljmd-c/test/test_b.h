#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include <data_structures.h>
#include <verlet1.h>
#include <verlet2.h>
#include <input.h>
 
/*routines to initialize the structure with the testing values */
void init_pos(mdsys_t *, double *);
void init_vel(mdsys_t *, double *);
void init_for(mdsys_t *, double *);

/*routines to check the value inside the structure */
void print_pos(mdsys_t *);
void print_vel(mdsys_t *);
void print_for(mdsys_t *);

/* test function. */
/* inputs: */
/* - structure */
/* - test_n: 2 or 3 particles case */
/* - machine epsilon */
/* - pos_value, vel_value, for_value: testing values */

void test(mdsys_t *sys, const int test_n, const double eps, const double * pos_value, const double * vel_value, const double * for_value);

void init_pos(mdsys_t *sys, double * pos_value){
  
  for (int i=0; i<sys->natoms; ++i) {
    sys->rx[i] = pos_value[3*i];
    sys->ry[i] = pos_value[3*i+1];
    sys->rz[i] = pos_value[3*i+2];
  }
}

void init_vel(mdsys_t *sys, double * vel_value){
  
  for (int i=0; i<sys->natoms; ++i) {
    sys->vx[i] = vel_value[3*i];
    sys->vy[i] = vel_value[3*i+1];
    sys->vz[i] = vel_value[3*i+2];
  }
}

void init_for(mdsys_t *sys, double * for_value){
  
  for (int i=0; i<sys->natoms; ++i) {
    sys->fx[i] = for_value[3*i];
    sys->fy[i] = for_value[3*i+1];
    sys->fz[i] = for_value[3*i+2];
  }
}

void print_pos(mdsys_t *sys)
{   
  for (int i=0; i<sys->natoms; ++i)
    fprintf(stdout, "pos  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
}

void print_vel(mdsys_t *sys)
{   
  for (int i=0; i<sys->natoms; ++i)
    fprintf(stdout, "vel  %20.8f %20.8f %20.8f\n", sys->vx[i], sys->vy[i], sys->vz[i]);
}

void print_for(mdsys_t *sys)
{
  for (int i=0; i<sys->natoms; ++i)
    fprintf(stdout, "for  %20.8f %20.8f %20.8f\n", sys->fx[i], sys->fy[i], sys->fz[i]);
}

void test(mdsys_t *sys, const int test_n, const double eps, const double * pos_value, const double * vel_value, const double * for_value){
  
  
  /* switch function according to the type of test requested. See documentation calling the program without arguments.*/
  switch (test_n) {
    
    /* test 1: v(t)=0, x(t)=0 */
  case 1:
    
    /* print if needed as a first check */
    /* print_pos(sys); */
    /* print_vel(sys); */
    /* print_for(sys); */
    
    /* v(t+0.5dt)=v(t)+0.5a(t)dt */
    /* x(t+dt)=x(t)+v(t+0.5dt)dt */
    velverlet1(sys);
    
    /* do the test: check the results for the chosen values */
    for (int i=0; i<sys->natoms; ++i){
      if(
	 fabs(sys->vx[i] - 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass) < eps &&
	 fabs(sys->vy[i] - 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass) < eps &&
	 fabs(sys->vz[i] - 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass) < eps &&
	 fabs(sys->rx[i] - pos_value[3*i] - sys->dt*sys->vx[i]) < eps &&
	 fabs(sys->ry[i] - pos_value[3*i+1] - sys->dt*sys->vy[i]) < eps &&
	 fabs(sys->rz[i] - pos_value[3*i+2] - sys->dt*sys->vz[i]) < eps )
	printf("Verlet1-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
			printf("Verlet1-test (natoms: %d): FAILED.\n", sys->natoms);
			exit(1);
		}
    }
    
    /* v(t+dt)=v(t+0.5dt)+0.5a(t)dt */
    velverlet2(sys);
    
    /* do the second test: verlet2*/
    for (int i=0; i<sys->natoms; ++i){
      if(
	 fabs(sys->vx[i] - sys->dt / mvsq2e * sys->fx[i] / sys->mass) < eps &&
	 fabs(sys->vy[i] - sys->dt / mvsq2e * sys->fy[i] / sys->mass) < eps &&
	 fabs(sys->vz[i] - sys->dt / mvsq2e * sys->fz[i] / sys->mass) < eps )
	printf("Verlet2-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
	printf("Verlet2-test (natoms: %d): FAILED.\n", sys->natoms);
	exit(1);
	}
    }
    break;
    
    //test 2: v(t)=0, F(t)=0
  case 2:
    
    /* print_pos(sys); */
    /* print_vel(sys); */
    /* print_for(sys); */
    
    velverlet1(sys);
    
    for (int i=0; i<sys->natoms; ++i){
      if(
	 sys->vx[i] < eps &&
	 sys->vy[i] < eps &&
	 sys->vz[i] < eps &&
	 fabs(sys->rx[i] - pos_value[3*i]) < eps &&
	 fabs(sys->ry[i] - pos_value[3*i+1]) < eps &&
	 fabs(sys->rz[i] - pos_value[3*i+2]) < eps )
	printf("Verlet1-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
	printf("Verlet1-test (natoms: %d): FAILED.\n", sys->natoms);
	exit(1);
	}
    }
    
    velverlet2(sys);
    
   for (int i=0; i<sys->natoms; ++i){
     if(
	sys->vx[i] < eps &&
	sys->vy[i] < eps &&
	sys->vz[i] < eps )
       printf("Verlet2-test (natoms: %d): PASSED.\n", sys->natoms);
     else{
       printf("Verlet2-test (natoms: %d): FAILED.\n", sys->natoms);
	   exit(1);
	   }
   }
   break;
   
   //test 3: a(t)=-0.5v(t)/dt
  case 3:
    
    /* print_pos(sys); */
    /* print_vel(sys); */
    /* print_for(sys); */
    
    velverlet1(sys);
    
    for (int i=0; i<sys->natoms; ++i){
      if(
	 sys->vx[i] < eps &&
	 sys->vy[i] < eps &&
	 sys->vz[i] < eps &&
	 fabs(sys->rx[i] - pos_value[3*i]) < eps &&
	 fabs(sys->ry[i] - pos_value[3*i+1]) < eps &&
	 fabs(sys->rz[i] - pos_value[3*i+2]) < eps )
	printf("Verlet1-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
	printf("Verlet1-test (natoms: %d): FAILED.\n", sys->natoms);
	exit(1);
	}
    }

    velverlet2(sys);
    
    for (int i=0; i<sys->natoms; ++i){
      if(
	 fabs(sys->vx[i] + vel_value[3*i]) < eps &&
	 fabs(sys->vy[i] + vel_value[3*i+1]) < eps &&
	 fabs(sys->vz[i] + vel_value[3*i+2]) < eps )
	printf("Verlet2-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
	printf("Verlet2-test (natoms: %d): FAILED.\n", sys->natoms);
	exit(1);
	}
    }
    break;
    
    //test 4: a(t)>-0.5v(t)/dt
  case 4:
    
    /* print_pos(sys); */
    /* print_vel(sys); */
    /* print_for(sys); */
    
    velverlet1(sys);
    
    for (int i=0; i<sys->natoms; ++i){
      if(
	 sys->vx[i] > 0.0 &&
	 sys->vy[i] > 0.0 &&
	 sys->vz[i] > 0.0 &&
	 sys->rx[i] > pos_value[3*i] &&
	 sys->ry[i] > pos_value[3*i+1] &&
	 sys->rz[i] > pos_value[3*i+2] )
	printf("Verlet1-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
	printf("Verlet1-test (natoms: %d): FAILED.\n", sys->natoms);
	exit(1);
	}
      }
    
    velverlet2(sys);
    
    for (int i=0; i<sys->natoms; ++i){
      if(
	 (sys->vx[i] > 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass) &&
	 (sys->vy[i] > 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass) &&
	 (sys->vz[i] > 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass) )
	printf("Verlet2-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
	printf("Verlet2-test (natoms: %d): FAILED.\n", sys->natoms);
	exit(1);
	}
    }
    break;
    
    //test 5: a(t)<-0.5v(t)/dt
  case 5:
    
    /* print_pos(sys); */
    /* print_vel(sys); */
    /* print_for(sys); */
    
    velverlet1(sys);
    
    for (int i=0; i<sys->natoms; ++i){
      if(
	 sys->vx[i] < 0.0 &&
	 sys->vy[i] < 0.0 &&
	 sys->vz[i] < 0.0 &&
	 sys->rx[i] < pos_value[3*i] &&
	 sys->ry[i] < pos_value[3*i+1] &&
	 sys->rz[i] < pos_value[3*i+2] )
	printf("Verlet1-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
	printf("Verlet1-test (natoms: %d): FAILED.\n", sys->natoms);
	exit(1);
	}
    }
    
    velverlet2(sys);

    for (int i=0; i<sys->natoms; ++i){
      if(
	 (sys->vx[i] < 0.5*sys->dt / mvsq2e * sys->fx[i] / sys->mass) &&
	 (sys->vy[i] < 0.5*sys->dt / mvsq2e * sys->fy[i] / sys->mass) &&
	 (sys->vz[i] < 0.5*sys->dt / mvsq2e * sys->fz[i] / sys->mass) )
	printf("Verlet2-test (natoms: %d): PASSED.\n", sys->natoms);
      else{
		printf("Verlet2-test (natoms: %d): FAILED.\n", sys->natoms);
		exit(1);
		}
    }
    break;
  }
}
