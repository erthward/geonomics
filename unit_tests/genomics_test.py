import unittest
from sim import burnin
from sim import model
from structs import landscape
import geonomics as g 
import os 
import copy 


class GenomicsTestCases(unittest.TestCase):
    def make_params_file_test(self): 
        g.make_parameters_file()
    
    def read_params_test(self):
        params_filepath =  os.path.abspath("/home/xuboyue/Documents/geonomics/example_params.py") 
        params = g.read_params(params_filepath)
    def make_model_test(self):
        params_filepath =  os.path.abspath("/home/xuboyue/Documents/geonomics/example_params.py") 
        params = g.read_params(params_filepath)
        g.make_model(params)
        self.assertIsInstance(model.Model)
        return 
    def make_land_test(self):
        params = g.read_params()
        g.make_land()
        return 
    def make_genomic_architecture(self):
        return 
    def make_individual_test(self):
        return 
    def make_pop_test(self):
        return 
    def make_community(self):
        return 

if __name__ == '__main__':
    unittest.main()