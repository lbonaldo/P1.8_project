import unittest
import subprocess
import sys

sys.path.append('../py_interface')

import interface

class Test(unittest.TestCase):

    @classmethod
    def SetUpClass(self):
        pass

    def test_check_input_arguments(self):
    
        input_file = '../examples/argon_108.inp'
        test_sys=interface.mdsys()
        restfile,trajfile,ergfile,nprint = interface.read_input_file(input_file,test_sys)
        arg_input = [test_sys.natoms,
                     test_sys.mass,
                     test_sys.epsilon,
                     test_sys.sigma,
                     test_sys.rcut,
                     test_sys.box,
                     restfile,
                     trajfile,
                     ergfile,
                     test_sys.nsteps,
                     test_sys.dt,
                     nprint]
        
        arg_test = [108,39.948,0.2379,3.405,8.5,17.1580,"argon_108.rest","argon_108.xyz","argon_108.dat",10000,5.0,100]
        
        self.assertListEqual(arg_input, arg_test)
        self.assertEqual(len(arg_input),12)

    def test_read_file(self):
        test_sys=interface.mdsys()
        
        rest='../examples/argon_108.rest'
        
        interface.read_input_file('../examples/argon_108.inp',test_sys)
        interface.allocate(test_sys)
        interface.read_file(rest,test_sys)
        
        f=open(rest,"r")
        
        for i in range(test_sys.natoms):
            x,y,z=f.readline().split()
            self.assertAlmostEqual(float(x),test_sys.rx[i]) 
            self.assertAlmostEqual(float(y),test_sys.ry[i]) 
            self.assertAlmostEqual(float(z),test_sys.rz[i]) 
        
        for i in range(test_sys.natoms):
            vx,vy,vz=f.readline().split()
            self.assertAlmostEqual(float(vx),test_sys.vx[i]) 
            self.assertAlmostEqual(float(vy),test_sys.vy[i]) 
            self.assertAlmostEqual(float(vz),test_sys.vz[i]) 
        f.close()
            
    def test_force(self):
        subprocess.check_output('./test_a.x < test_a.inp', shell=True)

    def test_verlet(self):
        subprocess.check_output('./test_b.x 1 < test_b.inp', shell=True)
        subprocess.check_output('./test_b.x 2 < test_b.inp', shell=True)
        subprocess.check_output('./test_b.x 3 < test_b.inp', shell=True)
        subprocess.check_output('./test_b.x 4 < test_b.inp', shell=True)
        subprocess.check_output('./test_b.x 5 < test_b.inp', shell=True)

    def test_check_force(self):
        subprocess.check_output('./check_force.x < argon_108.inp', shell=True)

    def test_check_energy(self):
        subprocess.check_output('../ljmd-split.x < argon_108.inp', shell=True)
        subprocess.check_output('./check_overall_energy.x < argon_108.inp ', shell=True)
        
if __name__ == '__main__':
    
    unittest.main()
