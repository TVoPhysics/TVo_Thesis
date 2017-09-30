import math
import numpy as np
#### Dynamic Response of a FP Cavity 


def Delta_L(freq, r_a, r_b, L):
	c = 3e8
	T = L/c
	tf = (1-r_a*r_b) / ( 1 -r_a*r_b*np.exp(-2*1j*freq*T) )
	return tf

def Delta_F(freq, r_a, r_b, L):
	c = 3e8
	T = L/c
	tf = (  (1-np.exp(-2*1j*freq*T)) * (1-r_a*r_b)  ) / (2*1j*freq*T * (1-r_a*r_b*np.exp(-2*1j*freq*T)) )
	return tf
