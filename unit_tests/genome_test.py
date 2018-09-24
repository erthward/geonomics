import unittest
from structs import genome 
from ops import mutation


class GenomeTestCasese(unittest.TestCase):        
    def test_make_genomic_architecture(self):
        params = exec(open('./params.py', 'r').read())
        g_params = params.genome
        genome_arch = genome.make_genomic_architecture(g_params)
        self.assertIsInstance(genome_arch, Genomic_Architecture)

    def test_make_recomb_array(self):
        params = exec(open('./params.py', 'r').read())
        g_params = params.genome
        try: 
            genome.make_recomb_array(g_params)
        except:
            print("Can't make the recomb array")

    

    
if __name__ == '__main__':
    unittest.main()