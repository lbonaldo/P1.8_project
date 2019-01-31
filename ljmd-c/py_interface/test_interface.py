import unittest
import interface
import argparse
import sys

class Test(unittest.TestCase):

    @classmethod
    def SetUpClass(self):
        pass

    def test_check_input_arguments(self):
        arg_input = [interface.sys_data.natoms,
                     interface.sys_data.mass,
                     interface.sys_data.epsilon,
                     interface.sys_data.sigma,
                     interface.sys_data.rcut,
                     interface.sys_data.box,
                     interface.restfile,
                     interface.trajfile,
                     interface.ergfile,
                     interface.sys_data.nsteps,
                     interface.sys_data.dt,
                     interface.nprint]
        
        arg_test = [108,39.948,0.2379,3.405,8.5,17.1580,"argon_108.rest","argon_108.xyz","argon_108.dat",1000,5.0,100]
        
        self.assertListEqual(arg_input, arg_test)
        self.assertEqual(len(arg_input),12)

    def test_read_file(self):
        rest = open('argon_108_rest.test', 'w+')
        #erg = open(argon_108_erg.test,"w")
        #traj = open(argon_108_traj.test,"w")
        interface.read_file(interface.restfile, interface.sys_data)
        for i in range(interface.sys_data.natoms):
            print("{:16.14f}\t{:16.14f}\t{:16.14f}".format(interface.sys_data.rx[i], interface.sys_data.ry[i], interface.sys_data.rz[i]), file=rest)
        for line1, line2 in zip(rest, interface.restfile):
            assertAlmostEqual(float(line1.split()[0]), float(line2.split()[0]))
            assertAlmostEqual(float(line1.split()[1]), float(line2.split()[1]))
            assertAlmostEqual(float(line1.split()[2]), float(line2.split()[2]))
        #output(interface.sys_data, erg, traj)
        #traj.close()
        #erg.close()
        rest.close()
        #if(!rest.closed and !traj.closed and !erg.closed):
        if(not(rest.closed)):
            print("Error. Files not closed!")

    def test_force(self):
        

if __name__ == '__main__':
    
    unittest.main()
