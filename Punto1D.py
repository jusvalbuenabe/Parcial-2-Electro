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
E2m=0. #Incide desde la izquierda

E2 = np.matrix([ [E2p],[E2m] ])                  # Creates Electric Field.

c=1.
wp=1.

lp= 2*np.pi*c/wp #\lambda_p
d=6*lp

def Interfaz(eps1, mu1, eps2, mu2):
    Y1=cm.sqrt(eps1/mu1)
    Y2=cm.sqrt(eps2/mu2)
    A = np.matrix([ [Y1+Y2, Y1-Y2],[Y1-Y2, Y1+Y2] ]) 
    T = (1/(2*Y1)) * A
    return T 

def Propagacion(d, w, eps):
    k=(1./c)*(w*w)*cm.sqrt(eps)
    D= np.matrix([ [cm.exp(-1j*k*d),0],[0, cm.exp(1j*k*d)] ]) 
    return D

def eps(w):
	eps= eps0 * ((w*w-wp*wp)/(w*w))
	return eps
	
#w=np.linspace(wp/10,2*wp, num = 1000)


#eps_plasma= eps(w)

E1p=[]
E1pNorm=[] #Norma E1+
E1pPhase=[] #Fase E1+
E1pReal=[] #Parte Real
E1m=[]
E1mNorm=[] #Norma E1-
E1mPhase=[] #Fase E1-

x1=np.linspace(-d,0,num=1000)
x2=np.linspace(0,d,num=1000)
x3=np.linspace(d,2*d,num=1000)

f=wp*0.7 #prueba
for z in x3:
	e=eps(f)
	A= Propagacion(z,f,eps0)
        E1 =A*E2
	E1p = np.array(E1[0])[0].tolist()
	E1m = np.array(E1[1])[0].tolist()
	#print E1p[0]
	#print E1m[0]
	#print E1m[0]/E1p[0]
	E1pNorm.append(cm.polar(E1p[0])[0]) #Norma E1+
	E1pPhase.append(cm.polar(E1p[0])[1]) #Fase E1+
	E1pReal.append(E1p[0].real) #Fase E1+
	E1mNorm.append(cm.polar(E1m[0])[0]) #Norma E1-
	E1mPhase.append(cm.polar(E1m[0])[1]) #Fase E1-
	#print cm.polar(E_f[0])
       
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#Grafica de Norma t	
plt.plot(x3, E1pNorm, label='$|E1^+|$')
plt.plot(x3, E1pPhase, label='Phase')
plt.plot(x3, E1pReal, label='Real')
plt.grid(True)
plt.legend()
plt.show()
for z in x2:
	e=eps(f)
	A= Propagacion(d,f,e)*Interfaz(eps0,mu0,e,mu0)*Propagacion(d,f,e)
        
for z in x1:
	e=eps(f)
	A=Propagacion(d,f,e)*Interfaz(eps0,mu0,e,mu0)*Propagacion(d,f,e)*Interfaz(e,mu0,eps0,mu0)
        
        #Campos#
	#E1 =A*E2
	#E1p = np.array(E1[0])[0].tolist()
	#E1m = np.array(E1[1])[0].tolist()
	#print E1p[0]
	#print E1m[0]
	#print E1m[0]/E1p[0]
	#E1pNorm.append(cm.polar(E1p[0])[0]) #Norma E1+
	#E1pPhase.append(cm.polar(E1p[0])[1]) #Fase E1+
	#E1mNorm.append(cm.polar(E1m[0])[0]) #Norma E1-
	#E1mPhase.append(cm.polar(E1m[0])[1]) #Fase E1-
	#print cm.polar(E_f[0])
