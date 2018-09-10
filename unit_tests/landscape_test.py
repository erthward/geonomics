import unittest
from sim import burnin
from sim import model
from structs import landscape
import copy

class LandscapeTestCases(unittest.TestCase):
    def test_make_random_scape(self): 
        params = exec(open('./params.py', 'r').read())
        dim = params.land.main.dim 
        try:
            random_scape = landscape.make_random_scape(dim, 10)
        except:
            print("Can not make the random scape")

    def test_make_defined_scape(self):
        params = exec(open('./params.py', 'r').read())
        dim = params.land.main.dim
        res = params.land.main.res
        ulc = params.land.main.ulc
        #TODO: The make defined scape has wrong number of parameters as called in the function
        #defined_scape = landscape.make_defined_scape(dim,)

    def test_make_land(self): 
        params = exec(open('./params.py', 'r').read())
        dim = params.land.main.dim
        res = params.land.main.res
        ulc = params.land.main.ulc
        land = landscape.make_land(params, 2)
        self.assertIsInstance(land, landscape.Land)
    
        
    def test_get_gis_rasters(self):
        #TODO: get_gis_rasters has wrong number of parameters compared with the one called in function at line 344
   


if __name__ == '__main__':
    unittest.main()