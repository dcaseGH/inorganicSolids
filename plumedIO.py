''' Run plumed (on DHC mac) as ~/Work/Plumed/plumed2-master/src/lib/plumed driver --plumed plumed.dat --ixyz *.xyz 
    N.B. atom indexing begins at 1, not 0
    Also note- if stuff doesn't work, may need to go to src/ then touch [module].on then make again'''

import numpy as np
import os
from hardcode import PLUMED_EXE

class PlumedInput():
    '''  '''

    def __init__(self):
        pass

    @staticmethod
#    def q4q6Array(structure,
#                  centralSpecies,
#                  coordinateSpecies,
#                  kwargs):
    def q4q6Array(*args,
                  **kwargs):
#                  r_0 = 1.0,
#                  d_0 = 0.):
        ''' Make table of (q4, q6) 
            Maintain order so that points can be identified within structure '''

        kwargs['angMo'] = 4
        q4Dict = PlumedInput.qlFromStructure(*args, **kwargs)

        kwargs['angMo'] = 6
        q6Dict = PlumedInput.qlFromStructure(*args, **kwargs)

        orderedKeys = sorted( map(int, ([x.replace('.mean', '').replace('q6_', '') for x in q6Dict.keys()])) )

        return np.array([[q4Dict['q4_' + i + '.mean'], q6Dict['q6_' + i + '.mean']] for i in map(str, orderedKeys)])


    @staticmethod
    def qlFromStructure(structure,
                        centralSpecies,
                        coordinateSpecies,
                        angMo      = None,
                        neatenFile = True,
                        cleanFiles = True,
                        r_0        = 1.0,
                        nn         = None,
                        xyzFile    = 'plumed.xyz',
                        datFile    = 'plumed.dat',
                        colvarFile = 'COLVAR'):

        ''' Assume switching function is along the lines of (1 - ((r - d0)/r0)**n) / (1-((r-d0)/r0)**m) 
            Leave n, m as defaults (NN, MM in PLUMED) - d0 should aim to include 1st solvation shell
            r_0 -> makes it very sharp '''

        from subprocessHandling import RunCommandNew

        outDict = {}
#        print 'assert structure is P1'

        datString = ''
        # Plumed uses numbers from 1, not 0
        coordinateIndices = [1 + x for x in structure.speciesMatchIndices(targetSpecies = coordinateSpecies)]
        centralIndices    = [1 + x for x in structure.speciesMatchIndices(targetSpecies = centralSpecies)]

        if neatenFile:
            from hardcode import groupIntegers
            coordString = groupIntegers(coordinateIndices, returnString=True)
        else:
            coordString = ','.join(map(str, coordinateIndices))

        labelPrefix = "q%s_" % angMo
        for c in centralIndices:
            if nn is not None:
                datString += "Q%s SPECIESA=%s SPECIESB=%s R_0=%s NN=%s MEAN LABEL=%s\n" %(angMo, c, coordString, r_0, nn, labelPrefix + str(c))
            else:
                datString += "Q%s SPECIESA=%s SPECIESB=%s R_0=%s MEAN LABEL=%s\n" %(angMo, c, coordString, r_0, labelPrefix + str(c))

        datString += 'PRINT ARG=' + ','.join([labelPrefix+str(c)+'.mean' for c in centralIndices]) + " FILE=%s\n" % colvarFile

        with open(xyzFile, 'w') as outf:
            outf.write(structure.xyzString(includeCellVectors=True))

        with open(datFile, 'w') as outf:
            outf.write(datString)
        
        runner = RunCommandNew(PLUMED_EXE + " driver --plumed %s  --ixyz %s"%(datFile, xyzFile))
        runner.run(timeout=100.)

        with open(colvarFile, 'r') as inFile:            
            #change this if reading many lines
            colLines = inFile.readlines()
            assert(len(colLines) == 2) 
            keyLine = colLines[0].split('time')[1].split()
            datLine = colLines[1].split()[1:]
            assert(len(keyLine) == len(datLine))
            for i in xrange(len(keyLine)):
                outDict[keyLine[i]] = float(datLine[i])
#        print runner.output

        #remove files (from the current directory)
#        for f in [datFile, xyzFile]:
        for f in [datFile, xyzFile, colvarFile]:
            if cleanFiles and os.path.exists(f):
                os.remove(f)

        return outDict

#    @staticmethod
#    def xyzFromCif(cif):
#        ''' One cif -> xyz with appropriate dressing '''
#        pass

    


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Write input and output files mainly')
    parser.add_argument('-xc', '--xyzFromCif', default=None, help='write xyz from cif')
    parser.add_argument('-q6e', '--q6CentralElement', default=None, help='print Q6 for species of this element')
    parser.add_argument('-q6c', '--q6CoordElement',   default='O',  help='Q6 coordination species')
    parser.add_argument('-cif', '--cifFile', default=None, help='input cif to make P1 structure')


    args = parser.parse_args()
    outputBaseString = 'outputPlumed'

#    if args.cifFile:
#        structure = 

    if args.q6CentralElement:
        from coreClasses import Species, Structure
        q6Dict = PlumedInput.q6FromStructure(Structure.fromCIF(args.cifFile, expandFullCell=True),
                                             Species(element=args.q6CentralElement),
                                             Species(element=args.q6CoordElement))

        print q6Dict

    if args.xyzFromCif:
        # change this one day perhaps- conflicts with explicitly giving a -cif
        from coreClasses import Structure
        with open(outputBaseString + '.xyz', 'w') as outf:
            outf.write(Structure().fromCIF(args.xyzFromCif, expandFullCell=True).xyzString(includeCellVectors=True))
