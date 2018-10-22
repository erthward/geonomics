import unittest
import structs
from structs import genome
from sim import model
from structs import individual

class IndividualTestCases(unittest.TestCase):
    def test_make_individual_simple_case(self):
        params = exec(open('./params.py', 'r').read())
        genome_arch = genome.make_genomic_architecture(params)
        ind = individual._make_individual(1, True, None, genome_arch)
        self.assertIsInstance(ind, individual.Individual)
        
    def test_make_individual_largeNumber(self):
        params = exec(open('./params.py', 'r').read())
        genome_arch = genome.make_genomic_architecture(params)
        ind = individual._make_individual(1000000, True, None, genome_arch)
        self.assertIsInstance(ind, individual.Individual)

if __name__ == '__main__':
    from structs import genome
    unittest.main()