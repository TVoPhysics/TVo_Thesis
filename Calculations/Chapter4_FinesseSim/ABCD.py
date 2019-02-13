import numpy as np

def space(d):
	m = np.matrix( ((1,d),(0,1)) )
	return m
def mirror(r):
	m = np.matrix( ((1,0),(-2/r,1)) )
	return m
def lens(f):
	if str(f) == 'inf':
		m = np.matrix( ((1,0),(0,1)) )
	else:
        	m = np.matrix( ((1,0),(-1/f,1)) )
	return m

def find_f(space,mirror,lens):
	return none
