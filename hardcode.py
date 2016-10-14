''' Write your links and so on here '''

GULP_EXE = '/Users/cmdc2-extra/Work/gulp-4.4/Src/gulp'
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

atomicValenceElectrons = {
'H' : 1,
'Li': 1,
'O' : 6,
'P' : 5,
'Ti': 4
}
