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

        q4Dict = PlumedInput.qlFromStructure(testStructure,
                                             Species(element='C'),
                                             Species(element='H'), 
                                             angMo = 4,
                                             r_0 = 1.5,
                                             nn = 4)

        print q4Dict

        q6Dict = PlumedInput.qlFromStructure(testStructure,
                                             Species(element='C'),
                                             Species(element='H'),
                                             angMo = 6,
                                             r_0 = 1.5,
                                             nn = 4)

        print q6Dict

        #test making the table

        np.testing.assert_array_almost_equal( PlumedInput.q4q6Array(testStructure,
                                                                    Species(element='C'),
                                                                    Species(element='H'),
                                                                    r_0 = 1.5,
                                                                    nn = 4),
                                              np.array([[0.430907, 0.639292]]))



if __name__ == '__main__':
    unittest.main()

