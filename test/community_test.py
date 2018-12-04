import unittest
from structs import community
import os
import geonomics as gnx
from structs import landscape



class CommunityTestCases(unittest.TestCase):
    """
    Unit tests for Community.py. 
    """
    def testMakeCommunity(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        land = landscape._make_landscape(params)
        com = community._make_community(land, params)
        self.assertEqual(com.n_pops, len(params.comm.pops))
        self.assertEqual(com.t, -1)
        self.assertIsInstance(com, community.Community)

if __name__ == '__main__':
    unittest.main()
