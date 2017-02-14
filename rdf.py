''' This gives more functionality to the structure.radialDistributionFunction tool
    '''
import numpy as np
from coreClasses import Structure

class RDF():
    ''' Static methods at the moment
        Can make it into more useful object if need more processing/pickling '''
    def __init__():
        pass

    @staticmethod
    def calculateRDF(structure,
                     centralAtomIndices,
                     returnUnbinnedData = False,
                     nBins              = 100,
                     cutoffRadius       = 5.,
                     targetSpecies = None
                     ):

        # parallelize if needed with map async
        ''' calculate structure.radialDistributionFunction and bin the results 
            Always assume that r is the centre of the bin '''
        totalDisplacements = [structure.radialDistributionFunction(speciesListIndex = i,
                                                                   targetSpecies    = targetSpecies,
                                                                   cutoffRadius     = cutoffRadius)
                              for i in list(centralAtomIndices)]

        #flatten array - but flatten() requires same number elements in each list
        totalDisplacements = np.concatenate(totalDisplacements, axis = 0)
        
        if returnUnbinnedData:
            return totalDisplacements

        # bin data --- kernel density ?? may not be appropriate if -> 1
        histogram = np.histogram(totalDisplacements,
                                 bins = nBins)
        '''
        have option to return histogram (and if called rdfData, plot it thusly)
        hist, bins = rdfData
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        plt.show()
        '''

        # continue to turn histogram into g(r,dr)= (V/N)(/4(pai)r^2dr) 

#        dr = bins[1] - bins[0],
#        V  = 4. * np.pi * bins[-1] ** 3.

        
        return RDF.histogramToRDF(histogram)

    @staticmethod
    def histogramToRDF(hist):
        ''' hist[0] is set of counts 
            hist[1] is set of limits to bins (one more than counts) '''
        
        N  = np.sum(hist[0])
        dr = hist[1][1] - hist[1][0]
        r  = np.array([(hist[1][i+1] + hist[1][i]) / 2 for i in xrange(hist[0].shape[0])])
        V  = 4. * np.pi * hist[1][-1] ** 3.
#        print N, dr, r.shape, V

        #return g(r, dr) -- should I return histogram? or r values??
        return ((V/(N * 4. * np.pi * dr)) * hist[0] / r**2, hist[1])
