import numpy as np
import matplotlib.pyplot as plt
import cmath as cm

eps1=1
mu1=1
eps2=2
mu2=1

E2p=1
E2m=1

d=10

def Interfaz(eps1, mu1, eps2, mu2):
    Y1=cm.sqrt(eps1/mu1)
    Y2=cm.sqrt(eps2/mu2)
    A = np.matrix([ [Y1+Y2, Y1-Y2],[Y1-Y2, Y1+Y2] ]) 
    T = (1/(2*Y1)) * A
    return T 

E2 = np.matrix([ [E2p],[E2m] ])                  # Creates Electric Field.
T=Interfaz(eps1,mu1,eps2,mu2)                    # Creates a Transference matrix.

print T.T                                    # Transpose of A.
print T*E2                                   # Matrix multiplication of A and x.
print T.I                                    # Inverse of A.
print np.linalg.solve(T,E2)     # Solve the linear equation system.
