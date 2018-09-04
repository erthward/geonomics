import unittest
from structs import community
from structs import landscape
import params

class dataTestCases(unittest.TestCase):
    def testMakeCommunity(self):
        params = exec(open('./params.py', 'r').read())
        land = landscape.make_land(params)
        com = community.make_community(land, params)
        self.assertEqual(com.n_pops, len(params.comm.pops))
        self.assertEqual(com.t, -1)

if __name__ == '__main__':
    unittest.main()