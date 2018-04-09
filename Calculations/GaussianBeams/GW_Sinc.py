# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 09:56:02 2018

@author: tvo88
"""

import numpy as np
import math
import matplotlib.pyplot as pl

c = 3e8
lam = 5e-7
L = 4000

tau = 2*L/c   

def Sinc_Func(f, tau):
    out = f.copy()
    for i in range(len(out)):
        x = f[i]*tau
        out[i] =  tau/lam * 2 * np.pi *c * (math.sin(np.pi *x)/(np.pi*x)) * np.exp(1j*np.pi * f[i] * tau)
    return out

f = np.arange(1.,100000.,10.)

phis = Sinc_Func(f,tau)
print(phis)

pl.figure(figsize=(9, 5), dpi=80)

pl.loglog(f,np.absolute(phis), '-',linewidth=3)

pl.xlabel("Frequency (MHz)")
pl.ylabel("Normalized Power (watts)")

###Main plot formatting
pl.tick_params(labelsize=12)
pl.tick_params(which='both', width=2)
pl.tick_params(which='major', length=7)
pl.tick_params(which='minor', length=4)
pl.grid(True, zorder=-1)
pl.legend(fontsize=20)

pl.show()