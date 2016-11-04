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

N=21
a=1.5
w=a*wp
f=w
T=2.*np.pi/w
t=np.linspace(0.,1.0,num=N)

nf=str(a)

for ii in t:
    #print ii
    #Fase= np.matrix([ [cm.exp(-1j*w*T/N),0],[0, cm.exp(-1j*w*T/N)] ]) 
    E2=np.matrix([ [E2p*cm.exp(1j*w*ii*T)],[E2m*cm.exp(1j*w*ii*T)] ])
   
    #Fase= np.matrix([ [cm.exp(-1j*w*T/N),0],[0, cm.exp(-1j*w*T/N)] ]) 
    # E2=Fase*E2

    # print E2

    ns= "%.2f"%ii


    E1pNorm=[] #Norma E1+
  
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
        
	
        #plt.plot(x1, E1pNorm,'--', label='$|E_1^+|$')
        #plt.plot(x1, E1pReal, label='$E_1^+$')
        
        #plt.plot(x1, E1mNorm,'--',label='$|E_1^-|$')
        #plt.plot(x1, E1mReal, label='$E_1^-$')
        
        #plt.plot(x1, E1TN,'--', label='$|E_T$')
        #plt.plot(x1, E1T, label='$E_T$')
        
        #plt.grid(True)
        #plt.legend()
        #plt.show()
    
    E1Ta=np.asarray(E1pNorm)
    E1MAX=max(np.abs(E1Ta))   

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
        
	
        #plt.plot(x12, E1pNorm,'--', label='$|E_1^+|$')
        #plt.plot(x12, E1pReal, label='$E_1^+$')
        
        #plt.plot(x12, E1mNorm,'--',label='$|E_1^-|$')
        #plt.plot(x12, E1mReal, label='$E_1^-$')
        
        #plt.plot(x12, E1TN,'--', label='$|E_T$')
        #plt.plot(x12, E1T, label='$E_T$')
        
        #plt.grid(True)
        #plt.legend()
        #plt.show()
        
        
    for z in x3:
        e=eps(f)
        A= Propagacion(z-2*d,f,eps0)
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

   

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(x123, E1pNorm/E1MAX,'--', label='$|E_1^+|$')
    plt.plot(x123, E1pReal/E1MAX, label='$E_1^+$')
    
    plt.plot(x123, E1mNorm/E1MAX,'--',label='$|E_1^-|$')
    plt.plot(x123, E1mReal/E1MAX, label='$E_1^-$')
    
    plt.plot(x123, E1TN/E1MAX,'--', label='$|E_T|$')
    plt.plot(x123, E1T/E1MAX, label='$E_T$')
    
    plt.grid(True)

    plt.axis([-d-0.2,2*d+0.2,-2.5,2.5])
  # plt.axis([-d-0.2,2*d+0.2,-1.5,1.5])
    
    
    labelsx=[r'$-\lambda_p$','$0$','$\lambda_p$','$2\lambda_p$']
    labelsx=[r'$-d$',r'$-\frac{d}{2}$',r'$0$',r'$\frac{d}{2}$',r'$d$',r'$\frac{3d}{2}$',r'$2d$']
    plt.xticks(np.arange(-d,2*d+1,d/2), labelsx, fontsize=16)
    plt.xlabel(r'$z$', fontsize=16)
    plt.ylabel(r'$\frac{E}{|E_1^+|}$', fontsize=16)

    ttitle='$t=$'+ns+'$T$'
    plt.title(ttitle)
    plt.subplots_adjust(right=0.8)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc=2, borderaxespad=0.)
    #ssave=ns+'w.pdf'
    #plt.savefig(ssave)
    ssave='w'+nf+'t'+ns+'T.pdf'
    plt.savefig(ssave)
    #plt.show()
    plt.close()
      
    
    ##GRAFICA DEL CAMPO
    plt.plot(x123, E1T/E1MAX, label='$E_T$')
    
    plt.grid(True)
    plt.axis([-d-0.2,2*d+0.2,-2.5,2.5])

    
    labelsx=[r'$-\lambda_p$','$0$','$\lambda_p$','$2\lambda_p$']
    labelsx=[r'$-d$',r'$-\frac{d}{2}$',r'$0$',r'$\frac{d}{2}$',r'$d$',r'$\frac{3d}{2}$',r'$2d$']
    plt.xticks(np.arange(-d,2*d+1,d/2), labelsx, fontsize=16)
    
    plt.xlabel(r'$z$', fontsize=16)
    plt.ylabel(r'$\frac{E}{|E_1^+|}$', fontsize=16)
    ttitle='$t=$'+ns+'$T$'
    plt.title(ttitle)
    ssave2='w'+nf+'t'+ns+'u.pdf'
    plt.savefig(ssave2)
    #plt.show()
    plt.close()
    #print(ns)

#plt.show()
    
    
    


    
