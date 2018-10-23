import unittest
from sim import model
from structs import landscape
from structs import population
import copy

class modelTestCases(unittest.TestCase):    
    def test_run(self):
        model.Model.run(self)
        return 

    def test_walk(self):
        model.Model.walk(self)
if __name__ == '__main__':
    unittest.main()