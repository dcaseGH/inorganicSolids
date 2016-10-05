import unittest
import numpy as np

class TestVASPIO(unittest.TestCase):
    def test_XML(self):
        from vaspIO import VaspXML
        testStructure = VaspXML().initialStructure(open('testFiles/AlTiO_0_K422.xml', 'r').read())
        self.assertTrue(testStructure.speciesList[2].element == 'Al')
        np.testing.assert_array_equal(testStructure.speciesList[2].fracCoord, 
                                      np.array([0.5000000000,   0.3627000000,   0.0645000000]))
        self.assertTrue(testStructure.speciesList[-2].element == 'Ti')
        np.testing.assert_array_equal(testStructure.speciesList[-2].fracCoord, 
                                      np.array([0.5000000000,   0.6883000000,   0.2500000000]))

        finalStructure = VaspXML().finalStructure(open('testFiles/AlTiO_0_K422.xml', 'r').read())
        self.assertTrue(finalStructure.speciesList[-1].element == 'Ti')
        finalStructure.setCartCoord()
        np.testing.assert_array_almost_equal(np.array([[3.6342859954725371,    0.0000153912678447,    0.0000029302011885],
                                                       [0.0000399628039769,    9.4602687719173009,    0.0003371485723345],
                                                       [0.0000079334022948,    0.0003471534509835,    9.7344558029469788]]),
                                             finalStructure.unitCell.vectors,
                                             decimal = 4)
        np.testing.assert_array_almost_equal(finalStructure.speciesList[-1].cartCoord, 
                                             np.array([        1.81717,      2.96276,      7.30102 ]),
                                             decimal = 5)
      



if __name__ == '__main__':
    unittest.main()

