import numpy as np
import pykat
import pykat.optics.ABCD as abcd

#
#
#    q1   d1   f1  d2   f2  d3   q2
#    |   <-->  |  <-->  | <-->  |
#    | <---------- D ---------> |
#
#
#

def modematch(q_in, q_out, d1, d2, d3):
	d1_m = abcd.space(1.0,d1) 	
	d2_m = abcd.space(1.0,d2)
	d3_m = abcd.space(1.0,d3)
	d3_m_backwards = abcd.space(1.0,-d3)

	q1 =  abcd.apply(d1_m, q_in, 1, 1).q
	q4 = abcd.apply(d3_m_backwards, q_out, 1, 1).q

	alpha_1 = np.imag(1/q1)
	alpha_4 = np.imag(1/q4)
	  
	beta_2 = (np.sqrt(alpha_1/alpha_4 - (d2*alpha_1)**2) - 1.0)/d2	
	if isinstance(beta_2,complex):
		print("Warning! Real part of q_2 has imaginary part")


	alpha_2 = alpha_1	
	iq2 = complex(beta_2, alpha_2)
	q2 = 1/iq2
	q3 = abcd.apply(d2_m, q2, 1, 1).q

	f1 =float(np.real( ( (1/q1) - (1/q2) )**(-1)))	
	f2 =float(np.real( ( (1/q3) - (1/q4) )**(-1)))
	
	f1_m = abcd.lens(f1)
	f2_m = abcd.lens(f2)
	
	full_abcd = d3_m * f2_m * d2_m * f1_m * d1_m
	q_check = abcd.apply(full_abcd,q_in,1,1)


	overlap = pykat.BeamParam.overlap(q_check,q_out)
	
	if overlap>.9999:
		print("Successfully Mode-Matched! >99.9% overlap")
	else:
		print("Something went wrong!")
		print(q_check)
		print(q_out)


	print('f1 is ' +str(f1))
	print('f2 is ' +str(f2))	
	#print('Output_q is '+ str(q_out))
	
	return float(np.real(f1)), float(np.real(f2))
'''
#Testing the function
q_in  = complex(-3.53011974824173,1.11302760222413)
q_out = complex(-3.03860374824173,1.11302760222413)
D  = 0.491516 # total distance from SRC to OFI
d1 = 0.25
d2 = 0.0001
d3 = D-d1-d2

modematch(q_in,q_out,d1,d2,d3)

print(q_out)
'''
