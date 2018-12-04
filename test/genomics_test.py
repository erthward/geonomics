import unittest
from sim import model
from structs import landscape
import geonomics as g
import os
import copy



class GenomicsTestCases(unittest.TestCase):
    """
    Unit tests for genomics.py.
    """
    def make_params_file_test(self):
        try:
            g.make_parameters_file()
        except:
            print("System error, fail make parameter faile test.")

    def make_landscape_test(self):
        filenames = set(os.listdir('.'))
        g.make_parameters_file()
        new_filenames = set(os.listdir('.'))
        filename = [*new_filenames - filenames][0]
        try:
            g.make_landscape(filename)
        except:
            print("System error, fail make landscape test.")

    def make_population_test(self):
        filenames = set(os.listdir('.'))
        g.make_parameters_file()
        new_filenames = set(os.listdir('.'))
        filename = [*new_filenames - filenames][0]
        landscape = g.make_landscape(filename)
        try:
            g.make_population(landscape, filename)
        except:
            print("System error, fail make population test.")

    def make_community_test(self):
        filenames = set(os.listdir('.'))
        g.make_parameters_file()
        new_filenames = set(os.listdir('.'))
        filename = [*new_filenames - filenames][0]
        landscape = g.make_landscape(filename)
        try:
            g.make_community(landscape, filename)
        except:
            print("System error, fail make community test.")

    def make_genomic_architecture_test(self):
        filenames = set(os.listdir('.'))
        g.make_parameters_file()
        new_filenames = set(os.listdir('.'))
        filename = [*new_filenames - filenames][0]
        landscape = g.make_landscape(filename)
        try:
            g.make_genomic_architecture(filename, landscape)
        except:
            print("System error, fail make community test.")

    def make_model_test(self):
        try:
            g.make_model()
        except:
            print("System error, fail make model test without parameter.")

    def make_model_test_with_params(self):
        filenames = set(os.listdir('.'))
        g.make_parameters_file()
        new_filenames = set(os.listdir('.'))
        filename = [*new_filenames - filenames][0]
        try:
            g.make_model(filename)
        except:
            print("System error, fail make community test.")

    def make_individual_test(self):
        return


    def run_default_model_test(self):
        try:
            g.run_default_model()
        except:
            print("System error, fail running default model test")


    def read_parameters_file_test(self):
        return


if __name__ == '__main__':
    unittest.main()
