import numpy as np

def tidyStringNumber(stringIn):
    return float(stringIn.split('(')[0])

def getSpeciesListCIF(inputFile):
    ''' CIF files tend to be loop_ then _attributes then data, but for now assume 
        the third loop contains the atoms etc 
        Should consider that many things can be read, but for now just get label and fractional coordinates '''
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
                    speciesList.append(Species(label     = splitL[infoList.index('label')],
                                               element   = splitL[infoList.index('label')][:re.search("\d", splitL[infoList.index('label')]).start()],
                                               fracCoord = np.array(map(tidyStringNumber,
                                                                    [splitL[infoList.index('fract_x')],
                                                                     splitL[infoList.index('fract_y')],
                                                                     splitL[infoList.index('fract_z')]])))
                                       )
#    with open(inputFile, 'r') as iFile:
#        for l in iFile.read().split("loop_\n")[-1].split("\n"):
#            if "#end" in l.lower():
#                break
#            #edited by DHC 110816 to disregard the _atom_site_aniso_ information
##            if "_atom_site_" in l and not "_atom_site_aniso_" in l:
#            if "_atom_site_" in l:
#                infoList.append(l.replace("\n", "").split("_atom_site_")[1].lower())
#            if not "_atom_site_" in l:
#                splitL = l.split()
#                speciesList.append(Species(label     = splitL[infoList.index('label')],
#                                           element   = splitL[infoList.index('label')][:re.search("\d", splitL[infoList.index('label')]).start()],
#                                           fracCoord = np.array(map(tidyStringNumber,
#                                                                    [splitL[infoList.index('fract_x')],
#                                                                     splitL[infoList.index('fract_y')],
#                                                                     splitL[infoList.index('fract_z')]])))
#                                  )
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

    outString  = "_cell_length_a " + str(structure.unitCell.lengths[0]) + "\n"
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

    return outString


def writeCIF(fileName, structure, quickCif = True):
    ''' At the moment just makes a quick P1 cell '''
    if quickCif:
        cifString = quickCifString(structure)
    with open(fileName, 'w') as outf:
        outf.write(cifString)

