import numpy as np
import matplotlib.pyplot as plt
import cmath as cm

#eps1=1
#mu1=1
#eps2=2
#mu2=1

eps0=1. 
mu0=1.

E2p=1.
E2m=1.

d=10.

c=1.
wp=1.

def Interfaz(eps1, mu1, eps2, mu2):
    Y1=cm.sqrt(eps1/mu1)
    Y2=cm.sqrt(eps2/mu2)
    A = np.matrix([ [Y1+Y2, Y1-Y2],[Y1-Y2, Y1+Y2] ]) 
    T = (1/(2*Y1)) * A
    return T 

def Propagacion(d, w):
    k=(1/c)*cm.sqrt((w*w)-(wp*wp))
    D = np.matrix([ [cm.exp(-1j*k*d),0],[0, cm.exp(1j*k*d)] ]) 
    return D

def eps(w):
	eps= eps0 * (1.-(wp*wp)/(w*w))
	return eps
	

E2 = np.matrix([ [E2p],[E2m] ])                  # Creates Electric Field.

w=np.linspace(wp/10,2*wp, num = 1000)
x=np.linspace(-d,2*d,num=100)

eps_plasma= eps(w)
#plt.plot(w, eps_plasma, 'o')
#plt.show()

#I=Interfaz(eps0,mu0,eps_plasma,mu0)
#print I
E1=[]
E1pn=[] #Norma E1+
E1pp=[] #Fase E1+
E1mn=[] #Norma E1-
E1mp=[] #Fase E1-


for f in w:
	Ef = Interfaz(eps0,mu0,eps(f),mu0)*Propagacion(d,f)*Interfaz(eps(f),mu0,eps0,mu0)*E2
	Ef0 = np.array(Ef[0])[0].tolist()
	Ef1 = np.array(Ef[1])[0].tolist()
	E1pn.append(cm.polar(Ef0[0])[0]) #Norma E1+
	E1pp.append(cm.polar(Ef0[0])[1]) #Fase E1+
	E1mn.append(cm.polar(Ef1[0])[0]) #Norma E1-
	E1mp.append(cm.polar(Ef1[0])[1]) #Fase E1-
	#print cm.polar(E_f[0])
	E1.append(Ef)

#print E1

print E1pn
#plt.plot(w, E1pn)
#plt.plot(w, E1pp)
#plt.plot(w, E1mn)
#plt.plot(w, E1mp)
#plt.show()
#plt.plot(x1, y, 'o')
#[<matplotlib.lines.Line2D object at 0x...>]
#>>> plt.plot(x2, y + 0.5, 'o')
#[<matplotlib.lines.Line2D object at 0x...>]
#>>> plt.ylim([-0.5, 1])
#(-0.5, 1)
#>>> plt.show()

#print T.T                                    # Transpose of A.
print E2                                   # Matrix multiplication of A and x.
#print E1                                    # Inverse of A.
#print np.linalg.solve(T,E2)     # Solve the linear equation system.
