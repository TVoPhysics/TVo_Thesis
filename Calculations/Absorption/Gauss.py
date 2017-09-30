import numpy as np
### Gaussian Beam Properties

def Power(w,lamda,r1,r2,theta1,theta2):
	###Power of a gaussian beam with arbitrary angles and distances.
	###W is the beam size,not the waist size.

	A = 1/(2*np.pi)
	del_theta = theta2-theta1

	###Difference beween integrating to a value and integrating to inf
	if np.equal(str(r2),'inf'):
		power = (A*del_theta)*(np.exp(-2*(r1**2/w**2)))	
		print("Normalized Power is "+ str(power))
	else:
		power = (A*del_theta)*( np.exp(-2*(r1**2/w**2)) -np.exp(-2*(r2**2/w**2)) )
		print("Normalized Power is " + str(power))


	return power

def Overlap(q1,q2):
	overlap = abs(4*q1.imag * q2.imag)/abs(q1.conjugate()-q2)**2
	return overlap
