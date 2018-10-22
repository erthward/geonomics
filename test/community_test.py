import unittest
from structs import community
from structs import landscape


class CommunityTestCases(unittest.TestCase):
    def testMakeCommunity(self):
        params = exec(open('./params.py', 'r').read())
        land = landscape._make_landscape(params)
        com = community._make_community(land, params)
        self.assertEqual(com.n_pops, len(params.comm.pops))
        self.assertEqual(com.t, -1)
        self.assertIsInstance(com, community.Community)

if __name__ == '__main__':
    unittest.main()