''' Run plumed (on DHC mac) as ~/Work/Plumed/plumed2-master/src/lib/plumed driver --plumed plumed.dat --ixyz *.xyz '''

import numpy as np

class PlumedInput():
    '''  '''

    def __init__(self):
        pass

#    @staticmethod
#    def xyzFromCif(cif):
#        ''' One cif -> xyz with appropriate dressing '''
#        pass

    


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Write input and output files mainly')
    parser.add_argument('-xc', '--xyzFromCif', default=None, help='write xyz from cif')

    args = parser.parse_args()

    outputBaseString = 'outputPlumed'

    if args.xyzFromCif:
        from coreClasses import Structure
        with open(outputBaseString + '.xyz', 'w') as outf:
            outf.write(Structure().fromCIF(args.xyzFromCif, expandFullCell=True).xyzString(includeCellVectors=True))
