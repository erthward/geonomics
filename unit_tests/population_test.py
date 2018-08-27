import unittest
from sim import model
from structs import landscape
from structs import population
import copy

class populationTestCases(unittest.TestCase):
    
    def test_make_population(self): 
        params = exec(open('./params.py', 'r').read())
        land = landscape.make_land(params)
        pop = population.make_population(land, params)
        self.assertEqual(pop.get_age(), 0)
        self.assertEqual(pop.get_cells, None)
        self.assertEqual(pop.get_coords, None)
        self.assertEqual(pop.get_dom, None)
        self.assertEqual(pop.get_fitness, None)

if __name__ == '__main__':
    unittest.main()