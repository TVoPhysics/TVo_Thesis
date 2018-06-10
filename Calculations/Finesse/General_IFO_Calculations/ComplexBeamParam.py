# -*- coding: utf-8 -*-
"""
For a given q-parameter, plot the complex phase space.
"""
import numpy as np
import matplotlib.pyplot as pl

class PhaseSpace:
    
    def __init__(self,real,imag):
        self.z  = real
        self.zR = imag
    
    def Overlap(self,x,y):
        q1 = np.complex(self.z,self.zR)
        q2 = np.complex(self.z+x,self.zR+y)
        overlap = abs(4*q1.imag*q2.imag) / abs(q1.conjugate()-q2)**2
        return overlap
    
    def PlotContour(self):
        x = np.linspace(-100,100,20)
        y = np.linspace(-100,100,20)
        s = (np.size(x),np.size(y))
        z = np.zeros(s)
        for i in range(len(x)):
            for j in range(len(y)):
                z[i][j] = self.Overlap(x[i],y[j])
        
        pl.figure(figsize=(12,8))
        pl.contour(x,y,z)
        pl.colorbar()
        pl.show()

        return None
    
### Testing and Examples
x = PhaseSpace(1355.156,225.975)
x.PlotContour()