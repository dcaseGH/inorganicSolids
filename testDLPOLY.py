import unittest
#from gulpBasics import *
#from inputOutput import *
import numpy as np

class TestDLPOLYOutput(unittest.TestCase):
    def testDataFrameOutput(self):
        ''' Assume have pandas '''
        from dlpolyIO import dlpolyOutput
        import pandas as pd
        egOutput  = dlpolyOutput(outputString = open('testFiles/OUTPUT310816', 'r').read())
        dataFrame = egOutput.dataFrameOutput('temp_tot', 'temp_shl')
#        print dataFrame.loc[:10]
        self.assertEqual(dataFrame.loc[0]['step'], 1)
        self.assertEqual(dataFrame.loc[0]['cpu-time'], 0.36)
        self.assertEqual(dataFrame.loc[0]['temp_shlAverage'], 2661.0)

class TestDLPOLYInput(unittest.TestCase):
    def testInitFromConfig(self):
        '''  '''
        from dlpolyIO import dlpolyInput
        newInput = dlpolyInput().initFromCONFIG('testFiles/exampleCONFIG', levcfg = 0)
        self.assertEqual(newInput.levcfg, 0)
        self.assertEqual(len(newInput.parentStructure.speciesList), 9000)
        self.assertTrue( newInput.speciesList[0].cartVelocity is None )
        np.testing.assert_array_almost_equal(newInput.speciesList[0].cartCoord,
                                             np.array([-0.2470917286E-01,   -0.2350931424,         10.40699720]))


class TestDLPOLYHistory(unittest.TestCase):
    def testMakePMGStructureList(self):
        from dlpolyIO import DLPOLYHistory
        pmgStructures = DLPOLYHistory.makePMGStructureList('testFiles/exampleLATPHISTORY')
        self.assertEqual(len(pmgStructures), 1)
        self.assertTrue( pmgStructures[0]._sites[0]._species._data.keys()[0].symbol == 'Al' )
        np.testing.assert_array_almost_equal(pmgStructures[0]._sites[0]._coords,
                                             np.array([-4.2284E+00, -2.4074E+00, -4.1802E+00]))
        self.assertTrue( pmgStructures[0]._sites[-1]._species._data.keys()[0].symbol == 'Ti' )
        np.testing.assert_array_almost_equal(pmgStructures[0]._sites[-1]._coords,
                                             np.array([-4.2747E+00, -2.3676E+00, -1.5118E+01]))

    def testInitPMGDiffusionAnalyzer(self):
        from dlpolyIO import DLPOLYHistory
        pmgDA = DLPOLYHistory.initPMGDiffusionAnalyzer(historyFileName = 'testFiles/exampleLATPHISTORY',
                                                       temperature     = 0,
                                                       specie          = 'Li',
                                                       time_step       = 5.e-9,
                                                       step_skip       = None)

    def testSelectOnlySpecie(self):
        from dlpolyIO import DLPOLYHistory
        pmgStructures = DLPOLYHistory.makePMGStructureList('testFiles/exampleLATPHISTORY',
                                                           selectOnlySpecie = 'Li')

        # it equals 423, because an extra point is included at origin
        self.assertEqual(len(pmgStructures[0]), 423)

    def testMakeListPositions(self):
        from dlpolyIO import DLPOLYHistory
        firstCartPosition = DLPOLYHistory.makeListPositions('testFiles/exampleLATPHISTORY', fracCoord = False, cartCoord = True)[0]
        np.testing.assert_array_almost_equal(np.array([-4.2284E+00, -2.4074E+00, -4.1802E+00]),
                                             firstCartPosition)

        cellVectors = np.array([[  42.18,      0.1496E-01, -0.1179E-01],
                                [ -21.08,       36.51,      0.7080E-02],
                                [-0.1341E-01,  0.2270E-02,   43.56]])

        firstFracPosition = DLPOLYHistory.makeListPositions('testFiles/exampleLATPHISTORY')[0]
        np.testing.assert_array_almost_equal(np.array([-4.2284E+00, -2.4074E+00, -4.1802E+00]),
                                             np.dot(firstFracPosition, cellVectors))


#        pmgDA
        #print pmgDA.__dict__

if __name__ == '__main__':
    unittest.main()

