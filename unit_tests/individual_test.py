import unittest
from sim import model
from structs import individual
from structs import genome

class IndividualTestCases(unittest.TestCase):
    def test_make_individual_simple_case(self):
        params = exec(open('./params.py', 'r').read())
        genome_arch = genome.make_genomic_architecture(params)
        ind = individual.make_individual(1, True, None, genome_arch)
        self.assertIsInstance(ind, individual.Individual)
        
    def test_make_individual_largeNumber(self):
        params = exec(open('./params.py', 'r').read())
        genome_arch = genome.make_genomic_architecture(params)
        ind = individual.make_individual(1000000, True, None, genome_arch)
        self.assertIsInstance(ind, individual.Individual)

if __name__ == '__main__':
    unittest.main()