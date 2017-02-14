import unittest
from coreClasses import *# Structure
from setTools import *

class TestXYZFiles(unittest.TestCase):
    def test_basicFunctions(self):
        self.assertTrue(XYZFile.standardSpeciesLine('Li2 23.4 1.e-10 100.0'))
        self.assertTrue(XYZFile.standardSpeciesLine('Li2 23.4 1.e-10 100.0\n'))
        self.assertFalse(XYZFile.standardSpeciesLine('Atoms 3 so there\n'))


class TestSymmetryGroup(unittest.TestCase):
    def test_symmetryGroup(self):
        sg = SymmetryGroup.fromCif('testFiles/example.cif')
        self.assertEqual(sg.number, 167)
        self.assertEqual(sg.labelHM, 'R-3c')
        self.assertEqual(len(sg.elementList), 36)

        testPt = np.array([2.3, -1.2, 3.4])
        np.testing.assert_array_almost_equal( sg.elementList[1].operate(testPt),
                                              np.array([-2.3,  1.2, -3.4]) )

        self.assertEqual(len( sg.generateUniquePoints([testPt])), 36 )
        self.assertEqual(len(  sg.generateUniquePoints([testPt] + [np.array([0.2903, 0., 0.25])], boundUnitCell = True)),
                               36 + 18)


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


class TestPotentials(unittest.TestCase):
    def test_buckingham(self):
        from savedPotentials import currentPotentialSet311016
        self.assertTrue('Li' in  [x.species1.element for x in currentPotentialSet311016 if x.__class__.__name__ == 'VBuckingham'])
        self.assertTrue(len(currentPotentialSet311016) == 7)
        currentPotentialSet311016 = VBuckingham.removeSpeciesFromListBucks(currentPotentialSet311016, Species(element='Li'))
        self.assertFalse('Li' in [x.species1.element for x in currentPotentialSet311016 if x.__class__.__name__ == 'VBuckingham'])
        self.assertTrue(len(currentPotentialSet311016) == 6)

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

    def test_removeListSpecies(self):
        tempStructure = Structure(speciesList = [Species(element = 'Cl', fracCoord = [0.1,0.3,.4]),
                                                 Species(element = 'Cl', core = 'shel', fracCoord = [0.1,0.3,.4]),
                                                 Species(element = 'Na', fracCoord = [0.3,0.3,.1]),
                                                 Species(element = 'Na', fracCoord = [0.3,0.3,.10]),
                                                 Species(element = 'S', fracCoord = [0.2,0.3,.1], charge = 0.2),
                                                 Species(element = 'S', fracCoord = [0.2,0.3,.1], charge = 0)
                                                 ])

        self.assertTrue(len(tempStructure.speciesList) ==6)
        removedList = tempStructure.removeListSpecies(Species(element = 'S'))
        self.assertTrue(len(tempStructure.speciesList) ==4 )
        self.assertTrue(len(removedList) ==2)
        np.testing.assert_array_almost_equal(removedList[1].fracCoord, 
                                             [0.2, 0.3, 0.1])

    def test_UnitCell(self):
        ''' This had better work, as stolen from pymatget and test is crap'''
        u = UnitCell(lengths = np.array([1., 2., 3.]), angles = np.array([45., 46., 47.]))
        np.testing.assert_array_almost_equal(u.calculateVectors(np.array([1., 2., 3.]), np.array([45., 46., 47.])),
                                             np.array([[ 0.7193398,   0.,          0.69465837],
                                                       [ 0.53048842,  1.310947,    1.41421356],
                                                       [ 0.,          0.,          3.        ]]))

        grid = u.createGrid()
        self.assertTrue(grid.shape == (24,3))
        np.testing.assert_array_almost_equal(grid[-1,:], np.dot([1,1,1], u.vectors))

    def test_numberValenceElectrons(self):
        from dlpolyIO import dlpolyInput
        s = dlpolyInput.initFromCONFIG('testFiles/exampleCONFIG').parentStructure
        print len([x for x in s.speciesList if x.element == 'Al'])
        print len([x for x in s.speciesList if x.element == 'Li'])
        print len([x for x in s.speciesList if x.element == 'O'])
        print len([x for x in s.speciesList if x.element == 'P'])
        print len([x for x in s.speciesList if x.element == 'Ti'])
        print 'finish test_numberValenceElectrons'
 
    def test_nearestNeighbourDistance(self):
        from coreClasses import Structure, Species
        structure = Structure().fromCIF('testFiles/exampleNeighbours.cif')
        np.testing.assert_almost_equal(structure.nearestNeighbourDistance(speciesListIndex = -1),
                                       1.9509368212716025)
        np.testing.assert_almost_equal(structure.nearestNeighbourDistance(speciesListIndex = -1,
                                                                          targetSpecies = Species(element = 'Al')),
                                       2.3585408060499997)

    def test_fromCIF(self):
        ''' Note two methods reorder things and dont always return cartesian coords  '''

        from coreClasses import Structure

        structure1 = Structure().fromCIF('testFiles/exampleNeighbours.cif')
        np.testing.assert_almost_equal(structure1.unitCell.lengths[2],
                                       20.974129)
        tempSpecies  = structure1.speciesList[0]
        self.assertTrue(tempSpecies.element == 'P')
        np.testing.assert_array_almost_equal(tempSpecies.fracCoord,
                                             np.array([0.2905, 0.0, 0.25]), decimal = 4)
        self.assertTrue(tempSpecies.cartCoord is None)

        structure2 = Structure().fromCIF('testFiles/exampleNeighbours.cif', expandFullCell = True)
        np.testing.assert_almost_equal(structure2.unitCell.lengths[2],
                                       20.974129)
        tempSpecies  = structure2.speciesList[0]
        self.assertTrue(tempSpecies.element == 'Li')
        np.testing.assert_array_almost_equal(tempSpecies.fracCoord,
                                             np.array([0.9891, 0.0002, 0.5028]), decimal = 4)
        np.testing.assert_array_almost_equal(tempSpecies.cartCoord,
                                             np.array([8.1429, 0.0015, 10.5466]), decimal = 4)


    def test_fromPMGStructure(self):
        from pymatgen.core import Structure as PMGS
        from coreClasses   import Structure
        pmgStructure = PMGS.from_file('testFiles/exampleNeighbours.cif')
        myStructure  = Structure().fromPMGStructure(pmgStructure)

        np.testing.assert_almost_equal(myStructure.unitCell.lengths[2],
                                       20.974129)

        tempSpecies  = myStructure.speciesList[0]
        self.assertTrue(tempSpecies.element == 'Li')
        np.testing.assert_array_almost_equal(tempSpecies.fracCoord,
                                             np.array([0.9891, 0.0002, 0.5028]), decimal = 4)
        np.testing.assert_array_almost_equal(tempSpecies.cartCoord,
                                             np.array([8.1429, 0.0015, 10.5466]), decimal = 4)


    def test_changeUnitCell(self):
        from coreClasses import Structure
        structure = Structure.fromCIF('testFiles/exampleNeighbours.cif')
        oldNAtoms = len(structure.speciesList)
        newLimits = np.array([[-0.5,0.5], [-0.5, 0.5], [1.5, 2.3]])
        structure.changeUnitCell(newLimits) 
        self.assertTrue(all([x.fracCoord[0] >= newLimits[0, 0] and x.fracCoord[0] < newLimits[0, 1] and 
                             x.fracCoord[1] >= newLimits[1, 0] and x.fracCoord[1] < newLimits[1, 1] and 
                             x.fracCoord[2] >= newLimits[2, 0] and x.fracCoord[2] < newLimits[2, 1] for x in structure.speciesList]))


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

    def test_initFromPMGAtom(self):
        from pymatgen.core import Structure as PMGS
        pmgStructure = PMGS.from_file('testFiles/exampleNeighbours.cif')

        from coreClasses import Species
        tempSpecies = Species().initFromPMGSite(pmgStructure._sites[0])
        self.assertTrue(tempSpecies.element == 'Li')
        np.testing.assert_array_almost_equal(tempSpecies.fracCoord,
                                             np.array([0.9891, 0.0002, 0.5028]), decimal = 4)
        np.testing.assert_array_almost_equal(tempSpecies.cartCoord,
                                             np.array([8.1429, 0.0015, 10.5466]), decimal = 4)

if __name__ == '__main__':
    unittest.main()
