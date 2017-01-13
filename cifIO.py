import numpy as np

def tidyStringNumber(stringIn):
    return float(stringIn.split('(')[0])


def tidyLabelWithDigits(stringIn):
    ''' If there is a digit in the string, chop off up to this point
        This is because Label = Number + element (ish) '''

    import re

    if bool(re.search(r'\d', stringIn)):
        return stringIn[:re.search("\d", stringIn).start()]
    else:
        return stringIn

def getSpeciesListCIF(inputFile):
    ''' inputFile is the name of an input file 
        use ASE at this point to quickly parse file 
        DOES NOTHING ABOUT FRACTIONAL OCCUPANCY '''

#    from pymatgen.io.cif import CifFile
    from ase.io.cif import parse_cif
    from coreClasses import Species

    aseParser = parse_cif(inputFile)
    labels    = aseParser[0][1]['_atom_site_label']
    elements  = map(tidyLabelWithDigits, aseParser[0][1]['_atom_site_label'])
    fracX     = aseParser[0][1]['_atom_site_fract_x']
    fracY     = aseParser[0][1]['_atom_site_fract_y']
    fracZ     = aseParser[0][1]['_atom_site_fract_z']

    return [Species(element   = elements[i],
                    label     = labels[i],
                    fracCoord = np.array([fracX[i],
                                          fracY[i],
                                          fracZ[i]], dtype = 'float64'))
            for i in xrange(len(labels))]


def getSpeciesListCIF_DEPRECATED(inputFile):
    ''' CIF files tend to be loop_ then _attributes then data, but for now assume 
        the third loop contains the atoms etc 
        Should consider that many things can be read, but for now just get label and fractional coordinates
        REPLACE THIS ONE DAY WITH SOMETHING FROM PMG OR ASE '''
    from coreClasses import Species
    import re
    speciesList   = []
    infoList          = []
#    speciesAttributes = [k.lower() for k in Species().__dict__]
#    loopCounter       = 0
    # should turn each loop into dictionary
    with open(inputFile, 'r') as iFile:
        for eachLoop in iFile.read().split("loop_\n"):
#            print eachLoop, 'eachloop'
            if not '_atom_site_label' in eachLoop.lower():
#                print '_atom_site_label not in', eachLoop.lower()
                continue
#            else:
#                print '_atom_site_label in', eachLoop.lower()
            for l in eachLoop.split("\n"):
                if "#end" in l.lower():
                   break
            #edited by DHC 110816 to disregard the _atom_site_aniso_ information
#            if "_atom_site_" in l and not "_atom_site_aniso_" in l:
                if "_atom_site_" in l:
                    infoList.append(l.replace("\n", "").split("_atom_site_")[1].lower())
                if not "_atom_site_" in l:
                    splitL = l.split()

                    # if there is a number in the label (should be) then cut off before to get element                    
                    if bool(re.search(r'\d', splitL[infoList.index('label')])):
                        element = splitL[infoList.index('label')][:re.search("\d", splitL[infoList.index('label')]).start()],
                    else:
                        element = splitL[infoList.index('label')]

                    speciesList.append(Species(label     = splitL[infoList.index('label')],
                                               element   = element,
                                               fracCoord = np.array(map(tidyStringNumber,
                                                                    [splitL[infoList.index('fract_x')],
                                                                     splitL[infoList.index('fract_y')],
                                                                     splitL[infoList.index('fract_z')]])))
                                       )

    return speciesList

def getSymmetryGroupCIF(inputFile):
    from coreClasses import SymmetryGroup
    with open(inputFile, 'r') as iFile:
        for l in iFile.readlines():
            if '_symmetry_space_group_name_h-m' in l.lower():
                labelHM = l.split("'")[1]
            if '_symmetry_int_tables_number' in l.lower():
                number = int(l.split()[1].replace("\n",""))

    return SymmetryGroup(labelHM = labelHM,
                         number  = number)

def loopToDict(loop):
    ''' Takes a loop_ from a cif, and returns a dictionary '''

    keys = [line for line in loop.split() if line.startswith('_')] 

    data = []
    for line in loop.split():
        if 'loop_' in line or any([k in line for k in keys]):
            continue
        else:
            data.append(line.split())
    print keys, data
    outDict = {}
    for i,x in enumerate(keys):
        outDict[x] = [d[i] for d in data]
    return outDict


def getSymmetryOperationsCIFstring(inputString):
    for l in inputString.split('loop_'):
        if '_symmetry_equiv_pos_site_id' not in l:
            continue
        return loopToDict(l)['_symmetry_equiv_pos_as_xyz']

def getUnitCellCIF(inputFile):
    ''' atm only reads things in precise order and format dhc 190716
        should ignore decimals in parentheses  
        include redundency/failure at some point '''
    from coreClasses import UnitCell
    with open(inputFile, 'r') as iFile:
        for l in iFile.readlines():
            if '_cell_length_a' in l.lower():
                lengths = [float(l.split()[1].split('(')[0])]
            if '_cell_length_b' in l.lower():
                lengths.append(float(l.split()[1].split('(')[0]))
            if '_cell_length_c' in l.lower():
                lengths.append(float(l.split()[1].split('(')[0]))

            if '_cell_angle_alpha' in l.lower():
                angles = [float(l.split()[1].split('(')[0])]
            if '_cell_angle_beta'  in l.lower():
                angles.append(float(l.split()[1].split('(')[0]))
            if '_cell_angle_gamma' in l.lower():
                angles.append(float(l.split()[1].split('(')[0]))

    if len(lengths) != 3:
        raise ValueError("Number of cell lengths in file %s is %s" %(inputFile, len(lengths)))

    if len(angles) != 3:
        raise ValueError("Number of cell angles in file %s is %s" %(inputFile, len(angles)))

    return UnitCell(lengths = np.array(lengths), angles = np.array(angles))

def quickCifString(structure):
    ''' P1 cell, no labels '''

    outString  = "data_\n" 
    outString += "_cell_length_a " + str(structure.unitCell.lengths[0]) + "\n"
    outString += "_cell_length_b " + str(structure.unitCell.lengths[1]) + "\n"
    outString += "_cell_length_c " + str(structure.unitCell.lengths[2]) + "\n"
    outString += "_cell_angle_alpha " + str(structure.unitCell.angles[0]) + "\n"
    outString += "_cell_angle_beta "  + str(structure.unitCell.angles[1]) + "\n"
    outString += "_cell_angle_gamma " + str(structure.unitCell.angles[2]) + "\n"
    outString += "loop_\n"
    outString += "_symmetry_equiv_pos_site_id\n"
    outString += "_symmetry_equiv_pos_as_xyz\n"
    outString += "1 'x, y, z'\n"
    outString += "loop_\n"
    outString += "_atom_site_label\n"
    outString += "_atom_site_fract_x\n"
    outString += "_atom_site_fract_y\n"
    outString += "_atom_site_fract_z\n"

    if structure.speciesList[0].fracCoord is None:
        structure.setFracCoord()

    for x in structure.speciesList:
        if x.core != 'core':
            continue
        outString += "%s %s %s %s \n"%(x.element, x.fracCoord[0], x.fracCoord[1], x.fracCoord[2])

    return outString + "#end\n"


def writeCIF(fileName, structure, quickCif = True):
    ''' At the moment just makes a quick P1 cell '''
    if quickCif:
        cifString = quickCifString(structure)
    with open(fileName, 'w') as outf:
        outf.write(cifString)

def writeCIFFromSpeciesList(fileName, 
                            speciesList, 
                            lengths      = np.array([50., 50.,50.]),
                            displacement = np.array([0.5, 0.5, 0.5])):
    ''' Puts the species in the list into a box and prints it 
        put into the middle with displacement = 0.5, 0.5, 0.5 '''
    from coreClasses import Structure, UnitCell

    tempStructure = Structure(unitCell = UnitCell(lengths = lengths,
                                                  angles  = np.array([90., 90., 90.])),
                              speciesList = speciesList)

    # moving to centre of cell is good idea for clarity
    tempStructure.displaceAllSpecies(fracDisplacement = displacement)

    with open(fileName, 'w') as outf:
        outf.write(quickCifString(tempStructure))
