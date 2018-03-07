import numpy as np
import math
from jacobi import jacobi

def HG2LG(n,m):
    """A function for Matlab which returns the coefficients and mode indices of
    the LG modes required to create a particular HG mode.
    Usage: coefficients,ps,ls = HG2LG(n,m)
    
    n,m:          Indces of the HG mode to re-create.
    coeffcients:  Complex coefficients for each order=n+m LG mode required to
                  re-create HG_n,m.
    ps,ls:        LG mode indices corresponding to coefficients.
    """
    # Mode order
    N = n+m;
    
    # Create empty vectors for LG coefficients/ indices
    coefficients = np.linspace(0,0,N+1,dtype=np.complex_)
    ps = np.linspace(0,0,N+1)
    ls = np.linspace(0,0,N+1)
    
    # Calculate coefficients
    for j in np.arange(0,N+1):
        
        # Indices for coefficients
        l = 2*j-N
        p = int((N-np.abs(l))/2)
        
        ps[j] = p
        ls[j] = l
        
        signl = np.sign(l)
        if (l==0):
            signl = 1.0

        # Coefficient
        c = (signl*1j)**m * np.sqrt(math.factorial(N-m)*math.factorial(m)/(2**N * math.factorial(np.abs(l)+p)*math.factorial(p)))
        c = c * (-1.0)**p * (-2.0)**m * jacobi(m,np.abs(l)+p-m,p-m,0.0)

        coefficients[j] = c
        
    return coefficients, ps, ls 

    
def LG2HG(p,l):
    """ Function to compute the amplitude coefficients
    of Hermite-Gauss modes whose sum yields a Laguerre Gauss mode
    of order n,m.
    Usage: coefficients, ns, ms = LG2HG(p,l)
    p:    Radial LG index
    l:    Azimuthal LG index
    The LG mode is written as u_pl with 0<=|l|<=p.
    The output is a series of coefficients for u_nm modes,
    given as complex numbers and respective n,m indices
    coefficients (complex array): field amplitude for mode u_nm
    ns (int array): n-index of u_nm
    ms (int array): m-index of u_nm

    
    The formula is adpated from M.W. Beijersbergen et al 'Astigmatic
    laser mode converters and transfer of orbital angular momentum',
    Opt. Comm. 96 123-132 (1993)
    We adapted coefficients to be compatible with our
    definition of an LG mode, it differs from
    Beijersbergen by a (-1)^p factor and has exp(il\phi) rather
    than exp(-il\phi).  Also adapted for allowing -l.
    Andreas Freise, Charlotte Bond    25.03.2007"""

    # Mode order
    N=2*p+np.abs(l)

    # Indices for coefficients
    n = np.abs(l)+p
    m = p

    # create empty vectors
    coefficients = np.linspace(0,0,N+1,dtype=np.complex_)
    ns = np.linspace(0,0,N+1)
    ms = np.linspace(0,0,N+1)

    # l positive or negative
    signl = np.sign(l)
    if (l==0):
        signl = 1.0
    
    # Beijersbergen coefficients
    for j in np.arange(0,N+1):
        ns[j]=N-j
        ms[j]=j

        c=(-signl*1j)**j * math.sqrt(math.factorial(N-j)*math.factorial(j)/(2**N * math.factorial(n)*math.factorial(m)))
        coefficients[j] = c * (-1.0)**p * (-2)**j * jacobi(j,n-j,m-j,0.0)
    
    return coefficients, ns, ms

x = LG2HG(2,0)
print(x)

y = HG2LG(2,0)
z = HG2LG(0,2)
print(y)
print(z)



