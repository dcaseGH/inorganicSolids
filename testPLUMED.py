import unittest
import numpy as np

class TestPlumed(unittest.TestCase):
    def test_q6FromStructure(self):
        from plumedIO import PlumedInput
        from coreClasses import Structure, Species, UnitCell

        testStructure = Structure(unitCell    = UnitCell(vectors = 10000. * np.identity(3)),
                                  speciesList = [Species(element = 'C', cartCoord = np.zeros(3)),
                                                 Species(element = 'H', cartCoord = np.array([1., 1., 1.])),
                                                 Species(element = 'H', cartCoord = np.array([1., -1., -1.])),
                                                 Species(element = 'H', cartCoord = np.array([-1., 1., -1.])),
                                                 Species(element = 'H', cartCoord = np.array([-1., -1., 1.]))])

        q6Dict = PlumedInput.q6FromStructure(testStructure,
                                             Species(element='C'),
                                             Species(element='H'), cleanFiles = False)

        print "CLEAN THE FILES AFTER PLUMED"

if __name__ == '__main__':
    unittest.main()

