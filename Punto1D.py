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
d=2*lp

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
E1pReal=[] #Parte Real E+

E1m=[]
E1mNorm=[] #Norma E1-
E1mPhase=[] #Fase E1-
E1mReal=[] #Parte Real E-

E1TN=[] #Norma del campo Neto
E1TPhase=[]#Fase Neto
E1T=[] #Campo Neto REAL

x1=np.linspace(-d,0,num=1000)
x2=np.linspace(0,d,num=1000)
x3=np.linspace(d,2*d,num=1000)


x12=np.concatenate((x1,x2))
x123=np.concatenate((x12,x3))

f=wp*1.5 #prueba

for z in x1:
    e=eps(f)
    A=Propagacion(z,f,eps0)*Interfaz(eps0,mu0,e,mu0)*Propagacion(0-d,f,e)*Interfaz(e,mu0,eps0,mu0)*Propagacion(d-2*d,f,eps0)
    E1 =A*E2
    E1p = np.array(E1[0])[0].tolist()
    E1m = np.array(E1[1])[0].tolist()

    E1pNorm.append(cm.polar(E1p[0])[0]) #Norma E1+
    E1pPhase.append(cm.polar(E1p[0])[1]) #Fase E1+
    E1pReal.append(E1p[0].real) #Parte Real de  E1+
    
    E1mNorm.append(cm.polar(E1m[0])[0]) #Norma E1-
    E1mPhase.append(cm.polar(E1m[0])[1]) #Fase E1-
    E1mReal.append(E1m[0].real) #Parte Real de  E1-


    E1TN.append(cm.polar(E1p[0]+E1m[0])[0]) #Norma de E1+ + E1-
    E1TPhase.append(cm.polar(E1p[0]+E1m[0])[1]) #Fase de E1+ + E1-
    E1T.append((E1p[0]+ E1m[0]).real) #Parte Real de  E1+
    #print cm.polar(E_f[0])
       
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

	
plt.plot(x1, E1pNorm,'--', label='$|E_1^+|$')
plt.plot(x1, E1pReal, label='$E_1^+$')

plt.plot(x1, E1mNorm,'--',label='$|E_1^-|$')
plt.plot(x1, E1mReal, label='$E_1^-$')

plt.plot(x1, E1TN,'--', label='$|E_T$')
plt.plot(x1, E1T, label='$E_T$')

plt.grid(True)
plt.legend()
plt.show()

for z in x2:
    e=eps(f)
    A= Propagacion(z-d,f,e)*Interfaz(e,mu0,eps0,mu0)*Propagacion(d-2*d,f,eps0)
    E1 =A*E2
    E1p = np.array(E1[0])[0].tolist()
    E1m = np.array(E1[1])[0].tolist()
    #print (E1p[0])
    #print E1m[0]
    #print E1m[0]/E1p[0]
    E1pNorm.append(cm.polar(E1p[0])[0]) #Norma E1+
    E1pPhase.append(cm.polar(E1p[0])[1]) #Fase E1+
    E1pReal.append(E1p[0].real) #Parte Real de  E1+
    
    E1mNorm.append(cm.polar(E1m[0])[0]) #Norma E1-
    E1mPhase.append(cm.polar(E1m[0])[1]) #Fase E1-
    E1mReal.append(E1m[0].real) #Parte Real de  E1-

    ET=(E1p[0]+ E1m[0]).real
    E1T.append(ET) #Parte Real de  E1+
    E1TN.append(cm.polar(E1p[0]+E1m[0])[0])

	
plt.plot(x12, E1pNorm,'--', label='$|E_1^+|$')
plt.plot(x12, E1pReal, label='$E_1^+$')

plt.plot(x12, E1mNorm,'--',label='$|E_1^-|$')
plt.plot(x12, E1mReal, label='$E_1^-$')

plt.plot(x12, E1TN,'--', label='$|E_T$')
plt.plot(x12, E1T, label='$E_T$')

plt.grid(True)
plt.legend()
plt.show()


for z in x3:
    e=eps(f)
    A= Propagacion(2*d-z,f,eps0)
    E1 =A*E2
    E1p = np.array(E1[0])[0].tolist()
    E1m = np.array(E1[1])[0].tolist()
    #print E1p[0]
    #print E1m[0]
    #print E1m[0]/E1p[0]
    E1pNorm.append(cm.polar(E1p[0])[0]) #Norma E1+
    E1pPhase.append(cm.polar(E1p[0])[1]) #Fase E1+
    E1pReal.append(E1p[0].real) #Parte Real de  E1+
    
    E1mNorm.append(cm.polar(E1m[0])[0]) #Norma E1-
    E1mPhase.append(cm.polar(E1m[0])[1]) #Fase E1-
    E1mReal.append(E1m[0].real) #Parte Real de  E1-

    ET=(E1p[0]+ E1m[0]).real
    E1T.append(ET) #Parte Real de  E1+
    E1TN.append(cm.polar(E1p[0]+E1m[0])[0])

	
plt.plot(x123, E1pNorm,'--', label='$|E_1^+|$')
plt.plot(x123, E1pReal, label='$E_1^+$')

plt.plot(x123, E1mNorm,'--',label='$|E_1^-|$')
plt.plot(x123, E1mReal, label='$E_1^-$')

plt.plot(x123, E1TN,'--', label='$|E_T$')
plt.plot(x123, E1T, label='$E_T$')

plt.grid(True)
plt.legend()
plt.title('\omega=1.5')
plt.show()

    
