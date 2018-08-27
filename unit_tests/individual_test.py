import unittest
from sim import model
from structs import individual
from structs import genome

class individualTestCases(unittest.TestCase):
    def test_make_individual_simple_case(self):
        params = exec(open('./params.py', 'r').read())
        genome_arch = genome.make_genomic_architecture(params)
        ind = individual.make_individual(1, True, None, genome_arch)
        self.assertEqual(ind.idx, 1)
        self.assertEqual(ind.age, 0)
        # sex, x and y are both random number so no need to write cases
        self.assertEqual(ind.fitness, None)
        
        self.assertEqual(ind.genome, genome.draw_genome(genome_arch))
        self.assertEqual(ind.phenotype, None)
        self.assertEqual(ind.habitat, None)
        
    def test_make_individual_largeNumber(self):
