import unittest
import subprocess
import sys

sys.path.append('../py_interface')

import t_interface

class Test(unittest.TestCase):

    @classmethod
    def SetUpClass(self):
        pass

    def test_check_input_arguments(self):
        arg_input = [t_interface.sys_data.natoms,
                     t_interface.sys_data.mass,
                     t_interface.sys_data.epsilon,
                     t_interface.sys_data.sigma,
                     t_interface.sys_data.rcut,
                     t_interface.sys_data.box,
                     t_interface.restfile,
                     t_interface.trajfile,
                     t_interface.ergfile,
                     t_interface.sys_data.nsteps,
                     t_interface.sys_data.dt,
                     t_interface.nprint]
        
        arg_test = [108,39.948,0.2379,3.405,8.5,17.1580,"argon_108.rest","argon_108.xyz","argon_108.dat",10000,5.0,100]
        
        self.assertListEqual(arg_input, arg_test)
        self.assertEqual(len(arg_input),12)

    def test_read_file(self):
        rest = open('argon_108_rest.test', 'w+')
        t_interface.read_file(t_interface.restfile, t_interface.sys_data)
        for i in range(t_interface.sys_data.natoms):
            print("{:16.14f}\t{:16.14f}\t{:16.14f}".format(t_interface.sys_data.rx[i], t_interface.sys_data.ry[i], t_interface.sys_data.rz[i]), file=rest)
        for line1, line2 in zip(rest, t_interface.restfile):
            assertAlmostEqual(float(line1.split()[0]), float(line2.split()[0]))
            assertAlmostEqual(float(line1.split()[1]), float(line2.split()[1]))
            assertAlmostEqual(float(line1.split()[2]), float(line2.split()[2]))
        rest.close()
        if(not(rest.closed)):
            print("Error. Files not closed!")
            
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
