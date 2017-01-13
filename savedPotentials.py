from coreClasses import VBuckingham, VSpring, VThreeBody, Species

''' List of saved potentials, with hopefully some notes about them '''

savedPotentials = {
'islamLFP':[
VBuckingham(species1 = Species(element = 'O', core='shel'),
                              species2 = Species(element = 'O', core='shel'),
                              A        = 22764.30,
                              rho      = 0.1490,
                              C6       = 27.89),
                  VBuckingham(species1 = Species(element = 'P'),
                              species2 = Species(element = 'O', core='shel'),
                              A        = 897.2648,
                              rho      = 0.3577,
                              C6       = 0.),
                  VBuckingham(species1 = Species(element = 'Li'),
                              species2 = Species(element = 'O', core='shel'),
                              A        = 632.1018,
                              rho      = 0.2906,
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
                              K        = 65.),
                  VSpring(    species1 = Species(element = 'Fe'),
                              K        = 19.26)
],
'lowestLTPSearch':[
                   VBuckingham(species1 = Species(element = 'O', core='shel'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 12759.206746,
                               rho      = 0.217390,
                               C6       = 29.975111),
                   VBuckingham(species1 = Species(element = 'P'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 897.2648,
                               rho      = 0.3577,
                               C6       = 0.),
                   VBuckingham(species1 = Species(element = 'Li'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 970.818704,
                               rho      = 0.266331,
                               C6       = 0.),
                   VBuckingham(species1 = Species(element = 'Ti'),
                               species2 = Species(element = 'O',  core='shel'),
                               A        = 6646.426774,
                               rho      = 0.243418,
                               C6       = 0.0),
                   VSpring(    species1 = Species(element = 'O'),
                               K        = 65.)
                  ],
'secondLTPSearch':[
                   VBuckingham(species1 = Species(element = 'O', core='shel'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 35733.031675,
                               rho      = 0.160849,
                               C6       = 25.744628),
                   VBuckingham(species1 = Species(element = 'P'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 897.2648,
                               rho      = 0.3577,
                               C6       = 0.),
                   VBuckingham(species1 = Species(element = 'Li'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 998.179946,
                               rho      = 0.257881,
                               C6       = 0.),
                   VBuckingham(species1 = Species(element = 'Ti'),
                               species2 = Species(element = 'O',  core='shel'),
                               A        = 3250.316127,
                               rho      = 0.281599,
                               C6       = 0.0),
                   VSpring(    species1 = Species(element = 'O'),
                               K        = 65.)
                  ],
'lowestLTPSearch_TiOModified':[
                   VBuckingham(species1 = Species(element = 'O', core='shel'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 12759.206746,
                               rho      = 0.217390,
                               C6       = 29.975111),
                   VBuckingham(species1 = Species(element = 'P'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 897.2648,
                               rho      = 0.3577,
                               C6       = 0.),
                   VBuckingham(species1 = Species(element = 'Li'),
                               species2 = Species(element = 'O', core='shel'),
                               A        = 970.818704,
                               rho      = 0.266331,
                               C6       = 0.),
                   VBuckingham(species1 = Species(element = 'Ti'),
                               species2 = Species(element = 'O',  core='shel'),
                               A        = 6.62236710e+03,
                               rho      = 2.54455157e-01,
                               C6       = 0.0),
                   VSpring(    species1 = Species(element = 'O'),
                               K        = 65.)
                  ],

'lowestAlTiOSearch':[
                     VBuckingham(species1 = Species(element = 'Al'),
                                 species2 = Species(element = 'O',  core='shel'),
                                 A        = 1761.207581 ,
                                 rho      = 0.292213 ,
                                 C6       = 0.0)],

'tempBVParameters120117':[
                          VBuckingham(species1 = Species(element = 'V'),
                                      species2 = Species(element = 'O', core='shel'),
                                      A        = 1.00535801e+03,
                                      rho      = 3.55733223e-01,
                                      C6       = 0.),
                          VBuckingham(species1 = Species(element = 'B'),
                                      species2 = Species(element = 'O', core='shel'),
                                      A        = 715.917969,
                                      rho      = 0.298132,
                                      C6       = 0.)
                          ]}
                          


currentPotentialSet141016 = savedPotentials['lowestLTPSearch'] + savedPotentials['lowestAlTiOSearch']
currentPotentialSet271016 = savedPotentials['lowestLTPSearch_TiOModified'] + savedPotentials['lowestAlTiOSearch']
currentPotentialSet311016 = savedPotentials['lowestLTPSearch_TiOModified'] + savedPotentials['lowestAlTiOSearch'] + [savedPotentials['islamLFP'][3]]
tempBVPotential120117     = currentPotentialSet311016 + savedPotentials['tempBVParameters120117']


''' islamLFP is modification of http://www.rsc.org/suppdata/c5/ta/c5ta09418f/c5ta09418f1.pdf
    lowestLTPSearch is the lowest one in in the Sobol search that passed the negative eigenvalue checks
    secondLTPSearch is 2nd lowest-- did this to check results for Li diffusion, and thermostats etc   
    lowestAlTiOSearch -> see Work/AlTiO/AlTiODataAnalysis.ipynb 
    modified TiO -> from script /Users/cmdc2-extra/Work/LATP/defectEnergies/OptimiseTiOParams.py (check details... Phonons (Rutile/anatase) and vacancies) 
    currentPotentialSet311016 is currentPotentialSet271016 + Islams Fe-O '''




