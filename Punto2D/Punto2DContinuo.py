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
    k=(1./c)*(w)*cm.sqrt(eps)
    D= np.matrix([ [cm.exp(-1j*k*d),0],[0, cm.exp(1j*k*d)] ]) 
    return D

def eps(w):
	eps= eps0 * ((w*w-wp*wp)/(w*w))
	return eps
	
#w=np.linspace(wp/10,2*wp, num = 1000)


#eps_plasma= eps(w)


n=np.arange(0.8,1.25,0.1)
for ii in n:
    if ii == 1.:
        continue
    f=wp*ii #prueba
    ns=str(ii)

    E1pNorm=[] #Norma E1+

    E1T=[] #Campo Neto REAL
    E1TN=[] #Norma del Campo
    E1TPhase=[]
    
    x1=np.linspace(-d,0,num=1000)
    x2=np.linspace(0,d,num=1000)
    x3=np.linspace(d,2*d,num=1000)

    x12=np.concatenate((x1,x2))
    x123=np.concatenate((x12,x3))
    for z in x1:
        e=eps(f)
        A=Propagacion(z,f,eps0)*Interfaz(eps0,mu0,e,mu0)*Propagacion(0-d,f,e)*Interfaz(e,mu0,2*eps0,mu0)*Propagacion(d-2*d,f,2*eps0)
        E1 =A*E2
        E1p = np.array(E1[0])[0].tolist()
        E1m = np.array(E1[1])[0].tolist()

        E1pNorm.append(cm.polar(E1p[0])[0]) #Norma E1+

        E1TN.append(cm.polar(E1p[0]+E1m[0])[0]) #Norma de E1+ + E1-
        E1TPhase.append(cm.polar(E1p[0]+E1m[0])[1]) #Fase de E1+ + E1-
        E1T.append((E1p[0]+ E1m[0]).real) #Parte Real de  E1+
    
    E1Ta=np.asarray(E1pNorm)
    E1MAX=max(np.abs(E1Ta))    
    
    for z in x2:
        e=eps(f)
        A= Propagacion(z-d,f,e)*Interfaz(e,mu0,2*eps0,mu0)*Propagacion(d-2*d,f,2*eps0)
        E1 =A*E2
        E1p = np.array(E1[0])[0].tolist()
        E1m = np.array(E1[1])[0].tolist()
       
        ET=(E1p[0]+ E1m[0]).real
        E1T.append(ET) #Parte Real de  E1+
        E1TN.append(cm.polar(E1p[0]+E1m[0])[0])
               
        
    for z in x3:
        e=eps(f)
        A= Propagacion(z-2*d,f,2*eps0)
        E1 =A*E2
        E1p = np.array(E1[0])[0].tolist()
        E1m = np.array(E1[1])[0].tolist()
       
        ET=(E1p[0]+ E1m[0]).real
        E1T.append(ET) #Parte Real de  E1+
        E1TN.append(cm.polar(E1p[0]+E1m[0])[0])

    E1Tb=np.asarray(E1T)
    
    ##GRAFICA DEL CAMPO
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
        
    plt.plot(x123, E1Tb/E1MAX, label=ns+'$\omega_p$')
    
plt.grid(True)
# plt.axis([-d-0.2,2*d+0.2,-1.5,1.5])

plt.legend(bbox_to_anchor=(1.05, 1.0), loc=2, borderaxespad=0.,fontsize=10)
   
    
#labelsx=[r'$-\lambda_p$','$0$','$\lambda_p$','$2\lambda_p$']
labelsx=[r'$-d$',r'$-\frac{d}{2}$',r'$0$',r'$\frac{d}{2}$',r'$d$',r'$\frac{3d}{2}$',r'$2d$']
plt.xticks(np.arange(-d,2*d+1,d/2), labelsx, fontsize=16)
    
plt.xlabel(r'$z$', fontsize=16)
plt.ylabel(r'$\frac{E}{|E_1^+|}$', fontsize=16)
ttitle='$\omega=N\omega_p$'
plt.title(ttitle)
#ssave2='E'+ns+'w.pdf'

plt.savefig('Continuo.pdf')
plt.show()
   
    
    


    
