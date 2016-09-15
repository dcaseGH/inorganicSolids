import unittest
#from coreClasses import *# Structure
import numpy as np

class TestDefects(unittest.TestCase):
    ''' Contains all defect Gulp/IO stuff  '''

    def testDefectsInputString(self):
        from gulpIO import InputFileGULP
        from coreClasses import Defect, Species

        a = InputFileGULP(defectSpecies = [Defect(species = Species(element = 'Fe',  fracCoord = [0.1,0.2,0.3]),
                                                  defectType    = 'Vacancy'),
                                           Defect(species = Species(label   = 'Fe2', fracCoord = [0.4,-0.2,0.3]),
                                                  defectType    = 'interstitial')],
                          defectDetails = {'centre':[0.2,0.3,0.0],
                                           'innerRadius': 10.,
                                           'outerRadius':20.}
                         )

        print "Currently doesnt do label vs element--- dunno how GULP works enough atm"
        self.assertEqual(a.defectsInputString(),
                        'Vacancy\nFe 0.1 0.2 0.3\ninterstitial\nNone 0.4 -0.2 0.3\ncentre\n0.2 0.3 0.0\nsize\n10.0 20.0\n')

    def testGULPOutputDefectEnergy(self):
        islamParams = {'LiOA': 632.1018,
                       'LiOrho': 0.2906,
                       'OOA': 22764.30,
                       'OOrho': 0.1490,
                       'OOC6': 27.89}

        currentParam = islamParams#dfFreq.loc[8]

        from coreClasses import Structure, VBuckingham, VSpring, VThreeBody, Species, Defect
        from gulpIO import InputFileGULP, OutputFileGULP
        from hardcode import GULP_EXE, CIF_FILE_DIRECTORY
        from subprocessHandling import RunCommand
        
        originalCIF = CIF_FILE_DIRECTORY + 'EntryWithCollCode155635.cif'
        xstal = Structure()
        xstal.inputCIF(originalCIF)

        qFeCore = 0#2.997                                                                                                                                       
        qOCore  = -2.86
        
        xstal.addShells(Species(element = 'O',  charge = qOCore))
        xstal.setCharges([Species(element = 'Fe', core='core', charge = 2. - qFeCore),
                          Species(element = 'P',  core='core', charge = 5.),
                          Species(element = 'Li', core='core', charge = 1.),
                          Species(element = 'O',  core='core', charge = -2 - qOCore)])
        
        inFile = InputFileGULP(fileName = 'LiFeP04',
                               title    = 'use to test stuff',
                               parentStructure = xstal,
                               keywords   = ['defect', 'regi', 'bulk', 'nodsym'],                                       
                               defectSpecies = [Defect(species = xstal.speciesList[0],
                                                      defectType    = 'vacancy')],
                               defectDetails = {'centre': xstal.speciesList[0].fracCoord,
                                                'innerRadius': 10.,
                                                'outerRadius': 20.},
                               potentials = [
                VBuckingham(species1 = Species(element = 'O', core='shel'),
                            species2 = Species(element = 'O', core='shel'),
                            A        = currentParam['OOA'],
                            rho      = currentParam['OOrho'],
                            C6       = currentParam['OOC6']),
                VBuckingham(species1 = Species(element = 'P'),
                            species2 = Species(element = 'O', core='shel'),
                            A        = 897.2648,
                            rho      = 0.3577,
                            C6       = 0.),
                VBuckingham(species1 = Species(element = 'Li'),
                            species2 = Species(element = 'O', core='shel'),
                            A        = currentParam['LiOA'],
                            rho      = currentParam['LiOrho'],
                            C6       = 0.),
                VBuckingham(species1 = Species(element = 'Fe'),
                            species2 = Species(element = 'O',  core='shel'),
                            A        = 1105.2409,
                            rho      = 0.3106,
                            C6       = 0.0),
                VThreeBody( species1 = Species(element = 'P'),
                            species2 = Species(element = 'O', core='shel'),
                            species3 = Species(element = 'O', core='shel'),
                            K        = 1.322626,
                            theta0   = 109.47,
                            cut12    = 2.,
                            cut13    = 2.,
                            cut23    = 4.),
                VSpring(    species1 = Species(element = 'O'),
                            K        = 65.)
                ])

        runner = RunCommand(GULP_EXE, inFile.stringForm())
        runner.run(timeout=10)
        np.testing.assert_almost_equal(OutputFileGULP.defectEnergy(runner.output),
                                       11.67965725)

    def testCheshireCells(self):
        ''' Test the matrix returns the same energy for the defect  '''
        pass

if __name__ == '__main__':
    unittest.main()
