from coreClasses import VBuckingham, VSpring, VThreeBody, Species

''' List of saved potentials, with hopefully some notes about them '''

savedPotentials = {
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

}

''' lowestLTPSearch is the lowest one in in the Sobol search that passed the negative eigenvalue checks
    secondLTPSearch is 2nd lowest-- did this to check results for Li diffusion, and thermostats etc    '''



