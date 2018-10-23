import unittest
from sim import model
from structs import landscape
from structs import population
import copy

class PopulationTestCases(unittest.TestCase):
    
    def test_make_population(self): 
        params = exec(open('./params.py', 'r').read())
        land = landscape._make_landscape(params)
        pop = population._make_population(land, params)
        self.assertIsInstance(pop, population.Population)


if __name__ == '__main__':
    unittest.main()