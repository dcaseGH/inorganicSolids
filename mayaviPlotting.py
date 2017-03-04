import numpy as np
from mayavi import mlab

class atomGraphics():
    def __init__(self,
                 color        = None,
                 scale_factor = None):
        self.color       = color
        self.scale_factor = scale_factor

#    @staticmethod

defaultAtomObjects = {
'H':  atomGraphics(color        = (1, 1, 0),
                   scale_factor = 0.6),
'Li': atomGraphics(color        = (0.4, 0.5, 0.2),
                   scale_factor = 1.0),
'B':  atomGraphics(color        = (0., 0.5, 0.2),
                   scale_factor = 1.0),
'C':  atomGraphics(color        = (0.4, 0.5, 0.2),
                   scale_factor = 1.0),
'O':  atomGraphics(color        = (0.4, 0.5, 0.2),
                   scale_factor = 1.0),
'P':  atomGraphics(color        = (0.4, 0.5, 0.2),
                   scale_factor = 1.0),
'Al': atomGraphics(color        = (0.4, 0.5, 0.2),
                   scale_factor = 1.0),
'Ti': atomGraphics(color        = (0.2, 0.5, 0.2),
                   scale_factor = 1.0),
'V':  atomGraphics(color        = (0.2, 0.5, 0.2),
                   scale_factor = 1.0)
}

class coreObjectPlotting():

    def __init__(self,
                 structure = None):
        self.structure = structure

    def plotUnitCell(self,
                     color=(0,0,0),
                     tube_radius=0.1):

        cellEdges = np.array([[[0., 0., 0.], [1., 0., 0.]],
                              [[0., 0., 0.], [0., 1., 0.]],
                              [[0., 0., 0.], [0., 0., 1.]],
                              [[0., 0., 1.], [0., 1., 1.]],
                              [[0., 0., 1.], [1., 0., 1.]],
                              [[0., 1., 0.], [1., 1., 0.]],
                              [[0., 1., 0.], [0., 1., 1.]],
                              [[1., 0., 0.], [1., 1., 0.]],
                              [[1., 0., 0.], [1., 0., 1.]],
                              [[1., 1., 0.], [1., 1., 1.]],
                              [[0., 1., 1.], [1., 1., 1.]],
                              [[1., 0., 1.], [1., 1., 1.]]])
        for e in cellEdges:
            p = np.dot(e, self.structure.unitCell.vectors)
            mlab.plot3d(p[:, 0], p[:, 1], p[:, 2], color=color, tube_radius=tube_radius)

#        return

    def plotSpecies(self, **kwargs):
        ''' Possibly put kwargs later '''

        if self.structure.speciesList[0].cartCoord is None:
            self.structure.setCartCoord()
        
        for s in self.structure.speciesList:
#            print "plotting", s.cartCoord
            mlab.points3d(s.cartCoord[0],s.cartCoord[1],s.cartCoord[2],
                          color        = defaultAtomObjects[s.element].color,
                          scale_factor = defaultAtomObjects[s.element].scale_factor)

#        return 

    def plotPolyhedra(self, **kwargs):
        ''' Adjust radius - add this at some point 
            kwargs must include centralSpecies - or have default of not O '''

        from setTools      import sameElementByAttributes
        from scipy.spatial import ConvexHull
        from coreClasses   import Species

        if self.structure.speciesList[0].cartCoord is None:
            self.structure.setCartCoord()

        if 'centralSpecies' not in kwargs.keys():
            # default is non oxygen atoms
            kwargs['centralSpecies'] = [Species(element = x.element) for x in self.structure.uniqueSpecies(['element']) if x.element != 'O']

        if not isinstance(kwargs['centralSpecies'], list):
            kwargs['centralSpecies'] = list([kwargs['centralSpecies']])

        if 'cutoffRadius' not in kwargs.keys():
            kwargs['cutoffRadius'] = 2.

#        print kwargs['centralSpecies']
        for c in kwargs['centralSpecies']:
            for s in self.structure.speciesMatchIndices(targetSpecies   = c,
                                                        matchAttributes = ['element', 'core']):
                try:
                    # need 4 points for convexhull in R^{3}-- hack this for B for e.g. later
                    ch = ConvexHull(self.structure.radialDistributionFunction(s,
                                                                              cutoffRadius  = kwargs['cutoffRadius'],
                                                                              returnVectors = True)
                                    + self.structure.speciesList[s].cartCoord)
                    mlab.triangular_mesh(ch.points[:, 0], ch.points[:, 1], ch.points[:, 2], ch.simplices,
                                         color = defaultAtomObjects[self.structure.speciesList[s].element].color)
                except:
                    print 'no convex hull for ', self.structure.speciesList[s].__dict__
                                                                          
