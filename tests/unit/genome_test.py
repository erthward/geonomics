import unittest
from structs import genome
import geonomics as gnx
import os
from ops import mutation


class GenomeTestCasese(unittest.TestCase):
    """
    Two ways of running single tests:
    python -m unittest discover -s test -p 'command_test.py'
    python test/command_test.py
    """
    def test_make_genomic_architecture(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        g_params = params.genome
        land = gnx.make_landscape(params)
        Genomic_Architecture = gnx.make_genomic_architecture(params, land)
        genome_arch = genome._make_genomic_architecture(params, land)
        self.assertIsInstance(genome_arch, Genomic_Architecture)

    def test_make_recomb_array(self):
        current_working_directory = os.getcwd()
        filepath = current_working_directory + "/GENOMICS_parameter.py"
        gnx.make_parameters_file(filepath)
        params = exec(open(filepath, 'r').read())
        g_params = params.genome
        try:
            genome._make_recomb_array(g_params, 1)
        except:
            print("Can't make the recomb array")

if __name__ == '__main__':
    unittest.main()
