import unittest

import interface

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(self):
    	pass

    def test_check_input_arguments(self):
        arg_input = [sys_data.natoms, 
		     sys_data.mass,
		     sys_data.epsilon,
		     sys_data.sigma,
		     sys_data.rcut,
		     sys_data.box,
		     restfile,
		     trajfile,
		     ergfile,
		     sys_data.nsteps,
		     sys_data.dt,
		     nprint]

	arg_test = []
        self.assertListEqual(arg_input, arg_test)

if __name__ == '__main__':

    unittest.main()
