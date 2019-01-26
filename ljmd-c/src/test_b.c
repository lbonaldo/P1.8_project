#include <test_b.h>

int main(int argc, char * argv[])
{
  int nprint, i;
  mdsys_t sys;

  double eps = 1e-10;
  
  int test_n = atoi(argv[1]);
  
  if( argc < 2 ){
    fprintf( stderr, "Error. The program runs as following: %s [TEST NUMBER]. \nCase 1:  v(t)=0, x(t)=0\nCase 2:  v(t)=0, F(t)=0\nCase 3: a(t)=-0.5v(t)/dt\nCase 4: a(t)>-0.5v(t)/dt\nCase 5: a(t)<-0.5v(t)/dt.\nProgram exit ...\n", argv[0]);
    exit(1);
  }
  
  /* read input file */
  /* if(get_a_line(stdin,line)) return 1; */
  /* sys.natoms=atoi(line); */
  if(get_a_line(stdin,line)) return 1;
  sys.mass=atof(line);
  /* if(get_a_line(stdin,line)) return 1; */
  /* sys.epsilon=atof(line); */
  /* if(get_a_line(stdin,line)) return 1; */
  /* sys.sigma=atof(line); */
  if(get_a_line(stdin,line)) return 1;
  sys.rcut=atof(line);
  if(get_a_line(stdin,line)) return 1;
  sys.box=atof(line);
  /* if(get_a_line(stdin,restfile)) return 1; */
  /* if(get_a_line(stdin,trajfile)) return 1; */
  /* if(get_a_line(stdin,ergfile)) return 1; */
  if(get_a_line(stdin,line)) return 1;
  sys.nsteps=atoi(line);
  /* if(get_a_line(stdin,line)) return 1; */
  /* sys.dt=atof(line); */
  if(get_a_line(stdin,line)) return 1;
  nprint=atoi(line);
  
  /* parameter specific for the test */
  sys.dt = 50.0;
  
#ifdef TWO
  int natoms = 2;
#end
  
#ifdef THREE
  int natoms = 3;
#end
  
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
  
  double * pos_value = (double *)malloc(sys.natoms*3*sizeof(double));
  double * vel_value = (double *)malloc(sys.natoms*3*sizeof(double));
  double * for_value = (double *)malloc(sys.natoms*3*sizeof(double));
  
  switch (test_n) {
    
    //test 1: v(t)=0, x(t)=0
  case 1:
#ifdef TWO
    pos_value = {0.0, 0.0, 0.0,
		 10.0, 10.0, 10.0};
    vel_value = {0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0};
    for_value = {10.0, 10.0, 10.0,
		 10.0, 10.0, 10.0};
#endif
    
#ifdef THREE
    pos_value = {0.0, 0.0, 0.0,
		 2.0, 2.0, 2.0,
		 4.0, 4.0, 4.0};
    vel_value = {0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0};
    for_value = {10.0, 10.0, 10.0,
		 10.0, 10.0, 10.0
                 10.0, 10.0, 10.0};
#endif
    break;
    
    //test 2: v(t)=0, F(t)=0
  case 2:
#ifdef TWO
    pos_value = {2.0, 2.0, 2.0,
		 4.0, 4.0, 4.0};
    vel_value = {0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0};
    for_value = {0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0};
#endif
    
#ifdef THREE
    pos_value = {2.0, 2.0, 2.0,
		 4.0, 4.0, 4.0,
                 6.0, 6.0, 6.0};
    vel_value = {0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0};
    for_value = {0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0,
		 0.0, 0.0, 0.0};
#endif
    break;
    
    //test 3: a(t)=-0.5v(t)/dt
  case 3:
    #ifdef TWO
    pos_value = {};
    vel_value = {};
    for_value = {};
#endif
    
#ifdef THREE
    pos_value = {};
    vel_value = {};
    for_value = {};
#endif
    break;
    
    //test 4: a(t)>-0.5v(t)/dt
  case 4:
#ifdef TWO
    pos_value = {};
    vel_value = {};
    for_value = {};
#endif
    
#ifdef THREE
    pos_value = {};
    vel_value = {};
    for_value = {};
#endif
    break;
    
    //test 5: a(t)<-0.5v(t)/dt
  case 5:
#ifdef TWO
    pos_value = {};
    vel_value = {};
    for_value = {};
#endif
    
#ifdef THREE
    pos_value = {};
    vel_value = {};
    for_value = {};
#endif
    break;
  }
  
  /* initialize positions, velocity and force of the structure*/
  
  init_pos(&sys, pos_value);
  init_vel(&sys, vel_value);
  init_for(&sys, for_value);
  
  /* initialize forces and energies.*/
  /* sys.nfi=0; */
  
  printf("Starting test with %d atoms for %d steps and dt is %f.\n", sys.natoms, sys.nsteps, sys.dt);
  
  print_pos(&sys);
  print_vel(&sys);
  print_for(&sys);
  
  /**************************************************/
  /* main MD loop */
  
  test(&sys, test_n, eps);
    
  /**************************************************/
  
  /* clean up: close files, free memory */
  printf("Simulation Done.\n");
  
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
