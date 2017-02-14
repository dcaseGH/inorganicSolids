import unittest
import numpy as np

class TestRDF(unittest.TestCase):
    def test_basicRDF(self):
        from coreClasses import Structure, UnitCell, Species
        from rdf import RDF

        structure = Structure(unitCell     =  UnitCell(vectors = np.diag((1.1, 2.22, 3.3))),
                              speciesList  = [Species(element = 'X',
                                                      fracCoord = np.zeros(3))]
                              )

        np.testing.assert_array_almost_equal(RDF.calculateRDF(structure,
                                                              [0],
                                                              returnUnbinnedData = True)[:7],
                                             np.array([1.1, 1.1, 2.2, 2.2, 2.22, 2.22, (1.1**2 + 2.22**2)**0.5]))

#        print RDF.calculateRDF(structure, [0], nBins = 42)
        self.assertTrue(RDF.calculateRDF(structure, [0], nBins = 42)[0].shape[0] == 42)

        

if __name__ == '__main__':
    unittest.main()
