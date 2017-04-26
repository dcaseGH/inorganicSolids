''' run as python matrixToLengthsPrinter.py lx ly lz xy xz yz 
    returns a matrix as a string xx xy xz yx etc '''

import numpy as np
import sys
from coreClasses import UnitCell
#from lammpsIO import 

lx, ly, lz, xy, xz, yz = map(float, sys.argv[1:])

matrix = np.array([[lx, xy, xz],
                   [0., ly, yz],
                   [0., 0., lz]])

print UnitCell.generalLengthsAnglesMatrix(matrix)
