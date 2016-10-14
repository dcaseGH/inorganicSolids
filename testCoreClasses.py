import unittest
from coreClasses import *# Structure
from setTools import *

class TestHardcode(unittest.TestCase):
    def test_timeStringToFloat(self):
        from hardcode import timeStringToFloat
        self.assertAlmostEqual(1.e-9, timeStringToFloat('1.000n'))
        self.assertAlmostEqual(2.5 * 60, timeStringToFloat('2.5m'))

class TestSetTests(unittest.TestCase):

    def test_sameElementByAttributes(self):
        tempSpec1 = Species(element = 'Cl', fracCoord = [0.1,0.3,.4])
        tempSpec2 = Species(element = 'Cl', fracCoord = [0.10,0.3,.4])
        tempSpec3 = Species(element = 'Cl', core='shel', fracCoord = [0.1,0.3,.4])
        tempSpec4 = Species(element = 'XXXXX', core='shel', fracCoord = [0.1,0.3,.4])

        self.assertTrue( sameElementByAttributes(tempSpec1, tempSpec2, ['element']))
        self.assertTrue( sameElementByAttributes(tempSpec1, tempSpec2, ['element', 'fracCoord']))
        self.assertTrue( sameElementByAttributes(tempSpec1, tempSpec3, ['element', 'fracCoord']))
        self.assertFalse(sameElementByAttributes(tempSpec1, tempSpec3, ['element', 'fracCoord', 'core']))
        self.assertFalse(sameElementByAttributes(tempSpec1, tempSpec4, ['element', 'fracCoord', 'core']))

    def test_subsetByAttributes(self):
        ''' N.B. subset != proper subset '''
        tempSpec1 = Species(element = 'Cl', fracCoord = [0.1,0.3,.4])
        tempSpec2 = Species(element = 'Cl', fracCoord = [0.10,0.3,.4])
        tempSpec3 = Species(element = 'Cl', core='shel', fracCoord = [0.1,0.3,.4])
        tempSpec4 = Species(element = 'XXXXX', core='shel', fracCoord = [0.1,0.3,.4])

        self.assertTrue (subsetByAttributes([tempSpec1, tempSpec2, tempSpec3],
                                            [tempSpec1, tempSpec3],
                                            ['element', 'core']))
        self.assertFalse(subsetByAttributes([tempSpec4],
                                           [tempSpec1, tempSpec2, tempSpec3],
                                           ['element', 'core']))

    def test_sameByAttributes(self):
        tempSpec1 = Species(element = 'Cl', fracCoord = [0.1,0.3,.4])
        tempSpec2 = Species(element = 'Cl', fracCoord = [0.10,0.3,.4])
        tempSpec3 = Species(element = 'Cl', core='shel', fracCoord = [0.1,0.3,.4])
        tempSpec4 = Species(element = 'XXXXX', core='shel', fracCoord = [0.1,0.3,.4])

        self.assertFalse(sameSetByAttributes([tempSpec1, tempSpec2, tempSpec3],
                                             [tempSpec1, tempSpec3],
                                             ['element', 'core']))
        self.assertTrue(sameSetByAttributes([tempSpec1, tempSpec2, tempSpec3],
                                            [tempSpec3, tempSpec2, tempSpec1],
                                            ['element', 'core']))




#class TestSpeciesMethods(unittest.TestCase):
#    def test_sameAttributes(self):
#        now use element in set tests above

class TestStructureMethods(unittest.TestCase):
    def test_UniqueSpecies(self):
        ''' Test that tempStructure 1 can be reduced to list with same attributes as tempStructure2
            also, there is a test which compares to an incomplete list '''
        tempStructure1 = Structure(speciesList = [Species(element = 'Cl', fracCoord = [0.1,0.3,.4]),
                                                  Species(element = 'Cl', core = 'shel', fracCoord = [0.1,0.3,.4]),
                                                  Species(element = 'Na', fracCoord = [0.3,0.3,.1]),
                                                  Species(element = 'Na', fracCoord = [0.3,0.3,.10]),
                                                  Species(element = 'S', fracCoord = [0.2,0.3,.1], charge = 0.2),
                                                  Species(element = 'S', fracCoord = [0.2,0.3,.1], charge = 0)
                                                 ])

        tempStructure2 = Structure(speciesList = [Species(element = 'Cl', fracCoord = [0.1,0.3,.4]),
                                                  Species(element = 'Cl', core = 'shel', fracCoord = [0.1,0.3,.4]),
                                                  Species(element = 'Na', fracCoord = [0.3,0.3,.1]),
                                                  Species(element = 'S', fracCoord = [0.2,0.3,.1], charge = 0.2)
                                                 ])

        self.assertTrue (len(tempStructure1.uniqueSpecies(['element', 'core'])) == 4)
        self.assertTrue (sameSetByAttributes(tempStructure1.uniqueSpecies(['element', 'core']),
                                             tempStructure2.speciesList,
                                             ['element', 'core']))
        self.assertTrue ( subsetByAttributes(tempStructure1.uniqueSpecies(['element', 'core']),
                                             tempStructure2.speciesList,
                                             ['element', 'core']))
        self.assertFalse( subsetByAttributes(tempStructure1.uniqueSpecies(['element', 'core']),
                                             tempStructure2.speciesList[:2],
                                             ['element', 'core']))

    def test_UnitCell(self):
        ''' This had better work, as stolen from pymatget and test is crap'''
        u = UnitCell(lengths = np.array([1., 2., 3.]), angles = np.array([45., 46., 47.]))
        print 
        np.testing.assert_array_almost_equal(u.calculateVectors(np.array([1., 2., 3.]), np.array([45., 46., 47.])),
                                             np.array([[ 0.7193398,   0.,          0.69465837],
                                                       [ 0.53048842,  1.310947,    1.41421356],
                                                       [ 0.,          0.,          3.        ]]))

    def test_numberValenceElectrons(self):
        from dlpolyIO import dlpolyInput
        s = dlpolyInput.initFromCONFIG('testFiles/exampleCONFIG').parentStructure
        print len([x for x in s.speciesList if x.element == 'Al'])
        print len([x for x in s.speciesList if x.element == 'Li'])
        print len([x for x in s.speciesList if x.element == 'O'])
        print len([x for x in s.speciesList if x.element == 'P'])
        print len([x for x in s.speciesList if x.element == 'Ti'])


class TestSpeciesMethods(unittest.TestCase):
    def test_setThermalVelocity(self):
        from coreClasses import Species
        oCore  = Species(element = 'O', core = 'core')
        oShell = Species(element = 'O', core = 'shel')
        self.assertEqual(0., np.linalg.norm(oShell.setThermalVelocity(100.)))
        self.assertNotEqual(0., np.linalg.norm(oCore.setThermalVelocity(100.)))

        import fortranformat as ff
        line = ff.FortranRecordWriter('(3f20.0)')
        print oCore.cartVelocity, line.write(oCore.cartVelocity.tolist())

    def test_atomicValenceElectrons(self):
        from coreClasses import Species
        s = Species(element = 'Li')
        self.assertTrue(s.atomicValenceElectrons() == 1)

if __name__ == '__main__':
    unittest.main()
