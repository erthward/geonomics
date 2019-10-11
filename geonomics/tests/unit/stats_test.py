import unittest
from structs import landscape
from structs import population
import geonomics as gnx
import os
from sim import stats


class StatstestCases(unittest.TestCase):
    """
    Unit tests for stats.py.
    """

    def test_write_stats(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        land = landscape._make_landscape(params)
        pop = population._make_population(land, params)
        stats._calc_ld(pop)

    def test_cal_het(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        land = landscape._make_landscape(params)
        pop = population._make_population(land, params)
        np_array = stats._calc_ld(pop)
        print(np_array)

if __name__ == '__main__':
    unittest.main()
