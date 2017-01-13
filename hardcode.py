''' Write your links and so on here '''

GULP_EXE           = '/Users/cmdc2-extra/Work/gulp-4.4/Src/gulp'
PLUMED_EXE         = '/Users/cmdc2-extra/Work/Plumed/plumed2-master/src/lib/plumed'
CIF_FILE_DIRECTORY = '/Users/cmdc2-extra/Work/icsdFiles/'

hartrees2eV   = 27.211396132
bohr2angstrom = 0.52917724900001

def timeStringToFloat(s):
    ''' Won't do anything that isn't number + string '''

    seconds = {'f': 1.e-15,
               'p': 1.e-12,
               'n': 1.e-9,
               's': 1.,
               'm': 60.,
               'h': 60.**2.,
               'd': 24.* 60.**2.,
               'w': 7. * 24. * 60.**2.}

    return float(s[:-1]) * seconds[s[-1].lower()]

def groupIntegers(listIn, returnString = True):
    ''' E.g. if you want to return the numbers of the atoms in some list '''
    from itertools import groupby
    from operator import itemgetter
    
    for k, g in groupby(enumerate(listIn), lambda (i, x): i-x):
        if returnString:
            #possibly need to save this as another object if using more than once???
            myList = map(itemgetter(1),g)
            return ",".join([str(min(myList)) + "-" + str(max(myList))])
        else:
            return map(itemgetter(1), g)

#### DO THIS ANOTHER TIME IF NEEDED
#def fractionalStoichiometryByWeight(components, returnOrderedDict = True):
#    ''' input is dictionary, e.g. {['LiV']: 20., ['BO2']: 80.}
#        output is dictionary of each element and also its stoichiometry '''
    
#    for i,j in components.items():
        
    
    #sum fractions to 1.
#    assert(abs(sum([x[1] for x in speciesFraction.items()]) - 1.) < 0.0001)
#    if returnOrderedDict:
#        from collections import OrderedDict
#        return
#    else:
#        return stoichiometryDict


atomicValenceElectrons = {
'H' : 1,
'Li': 1,
'O' : 6,
'P' : 5,
'Ti': 4
}

# N.B. only use default charges for very quick tests etc- specify for MD and GULP
defaultCharges = {
'O': -2.,
'Li': 1.,
'B':  3.,
'Ti': 4.,
'V':  5.,
'P':  5.,

}
