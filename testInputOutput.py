import unittest
#from gulpBasics import *
#from inputOutput import *
import numpy as np

''' Tests only work on dhcs computer atm- make portable if others ever want them  '''

class TestCIFFileReading(unittest.TestCase):
    def test_tidyStringNumber(self):
        from cifIO import tidyLabelWithDigits
        self.assertTrue(tidyLabelWithDigits('Ti3') == 'Ti')
        self.assertTrue(tidyLabelWithDigits('Ti') == 'Ti')
        self.assertTrue(tidyLabelWithDigits('1111111Ti3') == '')

    def test_readCIF(self):
        self.assertTrue(True)

    def test_getUnitCellCIF(self):
        from cifIO import getUnitCellCIF
        cell = getUnitCellCIF('/Users/cmdc2-extra/Work/icsdFiles/EntryWithCollCode193041.cif')
        np.testing.assert_array_equal(cell.angles,  np.array([90., 90., 120.]))
        np.testing.assert_array_equal(cell.lengths, np.array([8.39187, 8.39187, 22.73297]))

    def test_getSpeciesListCIF(self):
        from cifIO import getSpeciesListCIF
        speciesList = getSpeciesListCIF('/Users/cmdc2-extra/Work/icsdFiles/EntryWithCollCode193041.cif')
        self.assertTrue(len(speciesList) == 11)
        self.assertTrue(speciesList[0].label   == 'Li1')
        self.assertTrue(speciesList[2].element == 'Ti')
        np.testing.assert_array_equal(speciesList[-1].fracCoord, np.array([0.9126, 0.1442, 0.2995]))

    def test_getSymmetryGroup(self):
        from cifIO import getSymmetryGroupCIF
        symmetryGroup = getSymmetryGroupCIF('/Users/cmdc2-extra/Work/icsdFiles/EntryWithCollCode193041.cif')
        self.assertTrue(symmetryGroup.labelHM == 'R -3 H')
        self.assertTrue(symmetryGroup.number  == 148)


if __name__ == '__main__':
    unittest.main()

