import unittest
import sim
from sim import burnin as burnin
from sim import model as model
import os
import geonomics as gnx
import copy

class BurinTestCases(unittest.TestCase):
    """
    Unit Tests for `burnin.py`.
    """

    def test_test_adf_threshold(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        pop = None #TODO: place holder for pop after testing model
        self.assertEquals(True, burnin._test_adf_threshold(pop, 1))


    def test_test_t_threshold(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        pop = None
        self.assertEquals(True, burnin._test_t_threshold(pop, 2))

    def test_test_t_threshold_2(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory =
        print("This is a test for the test!")

    def test_test_t_threshold(self):
        current_working_directory = os.getcwd()

if __name__ == '__main__':
    unittest.main()
