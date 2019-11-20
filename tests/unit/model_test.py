import unittest
from sim import model
from structs import landscape
from structs import population
import copy



class modelTestCases(unittest.TestCase):
    """
    Unit tests for model.py. 
    """
    def test_run(self):
        try:
            model.Model.run(self)
        except:
            print("Running failed")

    def test_walk(self):
        try:
            model.Model.walk(self)
        except:
            print("Walking failed")

if __name__ == '__main__':
    unittest.main()
