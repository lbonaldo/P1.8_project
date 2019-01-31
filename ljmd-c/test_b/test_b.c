#include "test_b.h"

int main(int argc, char * argv[])
{

  /* machine epsilon */
  const double eps = 1e-8;
  char line[200];

  mdsys_t sys;
  
  if( argc < 2 ){
    printf("Error. The program runs as following: %s [TEST NUMBER]. \nCase 1: v(t)=0, x(t)=0\nCase 2: v(t)=0, F(t)=0\nCase 3: a(t)=-0.5v(t)/dt\nCase 4: a(t)>-0.5v(t)/dt\nCase 5: a(t)<-0.5v(t)/dt.\nProgram exit ...\n", argv[0]);
    exit(1);
  }

  /* type of tests */
  int test_n = atoi(argv[1]);

  printf("Test Number: %d\n", test_n);

  /* parameter specific for the test */
  if(get_a_line(stdin,line)) return 1;
  sys.mass=atof(line);
  if(get_a_line(stdin,line)) return 1;
  sys.rcut=atof(line);
  if(get_a_line(stdin,line)) return 1;
  sys.box=atof(line);
  sys.dt = 5.0;
  
#ifdef TWO
  sys.natoms = 2;
#endif
  
#ifdef THREE
  sys.natoms = 3;
#endif
  
  /* allocate memory for positions, velocity, force */
  sys.rx=(double *)malloc(sys.natoms*sizeof(double));
  sys.ry=(double *)malloc(sys.natoms*sizeof(double));
  sys.rz=(double *)malloc(sys.natoms*sizeof(double));
  sys.vx=(double *)malloc(sys.natoms*sizeof(double));
  sys.vy=(double *)malloc(sys.natoms*sizeof(double));
  sys.vz=(double *)malloc(sys.natoms*sizeof(double));
  sys.fx=(double *)malloc(sys.natoms*sizeof(double));
  sys.fy=(double *)malloc(sys.natoms*sizeof(double));
  sys.fz=(double *)malloc(sys.natoms*sizeof(double));
  
  /*allocate memory for the testing data */
  double * pos_value = (double *)malloc(sys.natoms*3*sizeof(double));
  double * vel_value = (double *)malloc(sys.natoms*3*sizeof(double));
  double * for_value = (double *)malloc(sys.natoms*3*sizeof(double));

  /* the swith is used to initialize the testing data according to the test chosen*/
  switch (test_n) {
    
    //test 1: v(t)=0, x(t)=0
  case 1:
#ifdef TWO

    pos_value[0] = 0.0;
    pos_value[1] = 0.0;
    pos_value[2] = 0.0;
    pos_value[3] = 10.0;
    pos_value[4] = 10.0;
    pos_value[5] = 10.0;

    vel_value[0] = 0.0;
    vel_value[1] = 0.0;
    vel_value[2] = 0.0;
    vel_value[3] = 0.0;
    vel_value[4] = 0.0;
    vel_value[5] = 0.0;

    for_value[0] = 10.0;
    for_value[1] = 10.0;
    for_value[2] = 10.0;
    for_value[3] = 10.0;
    for_value[4] = 10.0;
    for_value[5] = 10.0;

#endif
    
#ifdef THREE

    pos_value[0] = 0.0;
    pos_value[1] = 0.0;
    pos_value[2] = 0.0;
    pos_value[3] = 2.0;
    pos_value[4] = -2.0;
    pos_value[5] = 1.0;
    pos_value[6] = 4.0;
    pos_value[7] = -4.0;
    pos_value[8] = 2.0;

    vel_value[0] = 0.0;
    vel_value[1] = 0.0;
    vel_value[2] = 0.0;
    vel_value[3] = 0.0;
    vel_value[4] = 0.0;
    vel_value[5] = 0.0;
    vel_value[6] = 0.0;
    vel_value[7] = 0.0;
    vel_value[8] = 0.0;

    for_value[0] = 3.0;
    for_value[1] = 3.0;
    for_value[2] = 3.0;
    for_value[3] = 3.0;
    for_value[4] = 3.0;
    for_value[5] = 3.0;
    for_value[6] = 3.0;
    for_value[7] = 3.0;
    for_value[8] = 3.0;

#endif
    break;
    
    //test 2: v(t)=0, F(t)=0
  case 2:
#ifdef TWO

    pos_value[0] = 2.0;
    pos_value[1] = -2.0;
    pos_value[2] = 1.0;
    pos_value[3] = 4.0;
    pos_value[4] = -4.0;
    pos_value[5] = 2.0;

    vel_value[0] = 0.0;
    vel_value[1] = 0.0;
    vel_value[2] = 0.0;
    vel_value[3] = 0.0;
    vel_value[4] = 0.0;
    vel_value[5] = 0.0;

    for_value[0] = 0.0;
    for_value[1] = 0.0;
    for_value[2] = 0.0;
    for_value[3] = 0.0;
    for_value[4] = 0.0;
    for_value[5] = 0.0;
    
#endif
    
#ifdef THREE

    pos_value[0] = 2.0;
    pos_value[1] = -2.0;
    pos_value[2] = 1.0;
    pos_value[3] = 4.0;
    pos_value[4] = -4.0;
    pos_value[5] = 2.0;
    pos_value[6] = 6.0;
    pos_value[7] = -6.0;
    pos_value[8] = 3.0;

    vel_value[0] = 0.0;
    vel_value[1] = 0.0;
    vel_value[2] = 0.0;
    vel_value[3] = 0.0;
    vel_value[4] = 0.0;
    vel_value[5] = 0.0;
    vel_value[6] = 0.0;
    vel_value[7] = 0.0;
    vel_value[8] = 0.0;

    for_value[0] = 0.0;
    for_value[1] = 0.0;
    for_value[2] = 0.0;
    for_value[3] = 0.0;
    for_value[4] = 0.0;
    for_value[5] = 0.0;
    for_value[6] = 0.0;
    for_value[7] = 0.0;
    for_value[8] = 0.0;
    
#endif
    break;
    
    //test 3: a(t)=-0.5v(t)/dt
  case 3:
    #ifdef TWO
    
    pos_value[0] = 2.0;
    pos_value[1] = -2.0;
    pos_value[2] = 1.0;
    pos_value[3] = 4.0;
    pos_value[4] = -4.0;
    pos_value[5] = 2.0;

    vel_value[0] = -2.0e-04;
    vel_value[1] = 5.0e-5;
    vel_value[2] = -4.0e-4;
    vel_value[3] =  4.0e-5;
    vel_value[4] = 2.0e-6;
    vel_value[5] = -6.0e-5;
    
    for(int j = 0; j<6; ++j)
      for_value[j]=-2*sys.mass*mvsq2e/sys.dt*vel_value[j];
#endif
    
#ifdef THREE
    pos_value[0] = 2.0;
    pos_value[1] = -2.0;
    pos_value[2] = 1.0;
    pos_value[3] = 4.0;
    pos_value[4] = -4.0;
    pos_value[5] = 2.0;
    pos_value[6] = 6.0;
    pos_value[7] = -6.0;
    pos_value[8] = 3.0;
  
    vel_value[0] = -2.0e-04;
    vel_value[1] = 5.0e-5;
    vel_value[2] = -4.0e-4;
    vel_value[3] =  4.0e-5;
    vel_value[4] = 2.0e-6;
    vel_value[5] = -6.0e-5;
    vel_value[6] = -8.0e-5;
    vel_value[7] = 4.0e-5;
    vel_value[8] = -5.0e-5;

    for(int j = 0; j<9; ++j)
      for_value[j]=-2*sys.mass*mvsq2e/sys.dt*vel_value[j];
#endif
    break;
    
    //test 4: a(t)>-0.5v(t)/dt
  case 4:
#ifdef TWO

    pos_value[0] = 2.0;
    pos_value[1] = 1.0;
    pos_value[2] = 1.0;
    pos_value[3] = 0.0;
    pos_value[4] = 0.0;
    pos_value[5] = 0.0;

    vel_value[0] = -2.0e-3;
    vel_value[1] = 5.0e-4;
    vel_value[2] = -4.0e-4;
    vel_value[3] =  4.0e-4;
    vel_value[4] = 2.0e-5;
    vel_value[5] = -6.0e-4;

    for(int j = 0; j<6; ++j)
      for_value[j]=-2*sys.mass*mvsq2e/sys.dt*vel_value[j] + 2;
#endif
    
#ifdef THREE

    pos_value[0] = 4.0;
    pos_value[1] = 2.0;
    pos_value[2] = 2.0;
    pos_value[3] = 2.0;
    pos_value[4] = 1.0;
    pos_value[5] = 1.0;
    pos_value[6] = 0.0;
    pos_value[7] = 0.0;
    pos_value[8] = 0.0;

    vel_value[0] = -2.0e-3;
    vel_value[1] = 5.0e-4;
    vel_value[2] = -4.0e-4;
    vel_value[3] =  4.0e-4;
    vel_value[4] = 2.0e-5;
    vel_value[5] = -6.0e-4;
    vel_value[6] = -8.0e-4;
    vel_value[7] = 4.0e-4;
    vel_value[8] = -5.0e-4;

    for(int j = 0; j<9; ++j)
      for_value[j]=-2*sys.mass*mvsq2e/sys.dt*vel_value[j] + 2;
#endif
    break;
    
    //test 5: a(t)<-0.5v(t)/dt
  case 5:
#ifdef TWO

    pos_value[0] = 2.0;
    pos_value[1] = 1.0;
    pos_value[2] = 1.0;
    pos_value[3] = 0.0;
    pos_value[4] = 0.0;
    pos_value[5] = 0.0;

    vel_value[0] = -2.0e-3;
    vel_value[1] = 5.0e-4;
    vel_value[2] = -4.0e-4;
    vel_value[3] =  4.0e-4;
    vel_value[4] = 2.0e-5;
    vel_value[5] = -6.0e-4;

    for(int j = 0; j<6; ++j)
      for_value[j]=-2*sys.mass*mvsq2e/sys.dt*vel_value[j] - 2;
#endif
    
#ifdef THREE

    pos_value[0] = 4.0;
    pos_value[1] = 2.0;
    pos_value[2] = 2.0;
    pos_value[3] = 2.0;
    pos_value[4] = 1.0;
    pos_value[5] = 1.0;
    pos_value[6] = 0.0;
    pos_value[7] = 0.0;
    pos_value[8] = 0.0;

    vel_value[0] = -2.0e-3;
    vel_value[1] = 5.0e-4;
    vel_value[2] = -4.0e-4;
    vel_value[3] =  4.0e-4;
    vel_value[4] = 2.0e-5;
    vel_value[5] = -6.0e-4;
    vel_value[6] = -8.0e-4;
    vel_value[7] = 4.0e-4;
    vel_value[8] = -5.0e-4;

    for(int j = 0; j<9; ++j)
      for_value[j]=-2*sys.mass*mvsq2e/sys.dt*vel_value[j] - 2;
#endif
    break;
  }
  
  /* initialize positions, velocity and force of the structure*/
  
  init_pos(&sys, pos_value);
  init_vel(&sys, vel_value);
  init_for(&sys, for_value);
  
  printf("Starting test with %d atoms for %d steps and dt is %f.\n", sys.natoms, sys.nsteps, sys.dt);

  /* if needed as a check: it print the structure before starting the verlet algorithm*/
  /* print_pos(&sys); */
  /* print_vel(&sys); */
  /* print_for(&sys); */
  
  /**************************************************/
  
  test(&sys, test_n, eps, pos_value, vel_value, for_value);

  /**************************************************/
  
  /* clean up: close files, free memory */
  printf("Test Done.\n");
  
  free(sys.rx);
  free(sys.ry);
  free(sys.rz);
  free(sys.vx);
  free(sys.vy);
  free(sys.vz);
  free(sys.fx);
  free(sys.fy);
  free(sys.fz);
  
  free(pos_value);
  free(vel_value);
  free(for_value);
  
  return 0;
}
