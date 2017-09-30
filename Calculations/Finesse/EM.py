import numpy as np

class Plane(object):
	"""
	Create Compled Amplitude  objects either plane wave or gaussian beam

	"""
	def __init__(self):
		return none	
	
	def ComplexAmplitude( A, z, phase=0):
		wavelength = 1064e-9
		k = 2*np.pi/wavelength
		E = A * np.exp( -1j*(k*z +phase) )
		return E


class Gaussian(object):
	"""
        Gaussian beam

	"""
	def __init_(self):
		return none
		
	def ComplexAmplitude(z,zR,r):
		wavelength = 1064e-9
		k = 2*np.pi/wavelength
		A = np.sqrt( (k*zR) / np.pi)
		q = z + 1j*zR
		E = (A/q) * np.exp(-1j*(k*z)) * np.exp(-1j*k*r**2 / (2*q)) * np.exp(-1j*np.arctan(z/zR))
		return E

