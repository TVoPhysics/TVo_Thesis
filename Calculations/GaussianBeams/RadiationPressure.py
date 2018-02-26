# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:13:35 2018
Calculate the radiation pressure noise in strain/rthz
@author: tvo88
"""

import numpy as np

h = 6.63E-34
hbar = h/(2*np.pi)
c = 3E8
lam = 1064E-9
k = 2*np.pi/lam
omega = 2*np.pi*c/lam

E = h*c/lam

Pin = 125

m = 40
l=4000

x = (1/(2*m*c*np.pi**2))*np.sqrt(E*Pin)

print(x)

dphi= np.sqrt( (hbar*omega) /Pin)

y_RP = dphi /(2*k*l)
print(y_RP)