import math

c = 3e8
L_FC = 16 

w0 = 1.78e15

IMFC_R = 1-(6.1e-6)
EMFC_R = 1.0


FSR = (math.pi * c) / L_FC

for n in range(0,int(1e8),int(1e6)):
	omega = n * FSR
	x = "{0:.4e}".format(omega)
	print(x,n)



