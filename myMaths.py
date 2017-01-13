import numpy as np
from scipy.special import sph_harm
# import this if you need to
# from sympy.physics.wigner import wigner_3j

## Create a sphere
#r = 0.3
#pi = np.pi
#cos = np.cos
#sin = np.sin
#phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]

#x = r * sin(phi) * cos(theta)
#y = r * sin(phi) * sin(theta)
#z = r * cos(phi)

def polarToCart(v):
    ''' r, theta, phi -> x,y,z '''
    return np.array([v[0] * np.sin(v[2]) * np.cos(v[1]),
                     v[0] * np.sin(v[2]) * np.sin(v[1]),
                     v[0] * np.cos(v[2])])

def cartToPolar(v):
    ''' x,y,z -> r, theta, phi '''
    r     = np.linalg.norm(v)
    phi   = np.arccos(np.dot(v, np.array([0., 0., 1.])) / r)
    theta = np.arcsin(v[1] / (r*np.sin(phi))) 
    return np.array([r, theta, phi])

def ql_norm(l, positions):
    ''' positions is n,3 matrix 
        calculate ql (l is ang mo) , then return something to do with norm (check norm is fine for complex) '''

    dummyTot = 0.
    for p in positions:
        polarP = cartToPolar(p)
        qVec = np.array([sph_harm(m, l, polarP[1], polarP[2]) for m in xrange(-l, l + 1)])
        dummyTot += np.dot(qVec, np.conjugate(qVec))
    return dummyTot.real / float(positions.shape[1])

assert(np.linalg.norm(cartToPolar(polarToCart(np.array([2.3, np.pi/3.2, np.pi/1.2]))) - np.array([2.3, np.pi/3.2, np.pi/1.2])) < 0.01)
#assert(np.dot

#print phi.shape, phi[3]

#sph_harm(n, m, theta, phi) -n, m is usually m, l -- usually returns complex value
print sph_harm(1, 1, np.pi/2., np.pi/4.)
tetrahedralPoints = np.array([[1., 1., 1.],
                              [-1., -1., 1.],
                              [1., -1., -1.],
                              [-1., 1., -1.]])
print ql_norm(4, tetrahedralPoints)
print ql_norm(6, tetrahedralPoints)
