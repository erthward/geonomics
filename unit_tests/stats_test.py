import unittest
from structs import landscape
from structs import population
from sim import stats

class StatstestCases(unittest.TestCase):
    def test_write_stats(self):
        params = exec(open('./params.py', 'r').read())
        land = landscape.make_land(params)
        pop = population.make_population(land, params)
        stats.calc_ld(pop) 

    def test_cal_het(self):
        params = exec(open('./params.py', 'r').read())
        land = landscape.make_land(params)
        pop = population.make_population(land, params)
        np_array = stats.calc_ld(pop)
        print(np_array) 

if __name__ /== '__main__':
    unittest.main()