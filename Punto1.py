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

def Propagacion(d, w):
    #if(w<wp):
    #    k=(1./c)*cm.sqrt((wp*wp)-(w*w))
    #    kd=k*d
    #    D = np.matrix([ [cm.exp(kd),0],[0, cm.exp(-kd)] ]) 
    #else:
    #    k=(1./c)*cm.sqrt((w*w)-(wp*wp))
    #    kd=(k*d)%(2*np.pi)
    #    D = np.matrix([ [cm.exp(-1.0*1j*kd),0],[0, cm.exp(1.0*1j*kd)] ]) 
   
    k=(1./c)*cm.sqrt((w*w)-(wp*wp))
    kd=(k*d)#%(2*np.pi)
    D = np.matrix([ [cm.exp(-1.0*1j*kd),0],[0, cm.exp(1.0*1j*kd)] ]) 
    return D 

def eps(w):
	eps= eps0 * (((w*w)-(wp*wp))/(w*w))
	return eps
	

E2 = np.matrix([ [E2p],[E2m] ])                  # Creates Electric Field.

w=np.linspace(wp/10.,2.*wp, num = 1000)
x=np.linspace(-d,2*d,num=100)

eps_plasma= eps(w)

rNorm=[] #Norma r
rPhase=[] #Fase r
tNorm=[] #Norma t
tPhase=[] #Fase t-

rpNorm=[] #Norma rp
rpPhase=[] #Fase rp
tpNorm=[] #Norma tp
tpPhase=[] #Fase tp

R=[]
T=[]
RT=[]
Rp=[]
Tp=[]
RTp=[]
for f in w:
	e=eps(f)
	A= Interfaz(eps0,mu0,e,mu0)*Propagacion(d,f)*Interfaz(e,mu0,eps0,mu0)
	
	#t
	t=1/A[0,0]
	tNorm.append(np.abs(t))
	tPhase.append(cm.phase(t))
	#T
	tt=np.abs(t)*np.abs(t) #*sqrt(Y2/Y1) 
	T.append(tt) 
	
	#r
	r=A[1,0]*t
	rNorm.append(np.abs(r))
	rPhase.append(cm.phase(r))
	#R
	rr = np.abs(r)*np.abs(r)
	R.append(rr)
        
        #R+T
        RT.append(rr+tt)
	
	#rp
	rp=-1.*A[0,1]*t
	rpNorm.append(np.abs(rp))
	rpPhase.append(cm.phase(rp))
	#Rp
	rrp = np.abs(rp)*np.abs(rp)
	Rp.append(rrp)
	
	
	#tp
	tp=(A[0,0]*A[1,1]-A[0,1]*A[1,0])/A[0,0] #t-r*rp/t
	tpNorm.append(np.abs(tp))
	tpPhase.append(cm.phase((A[0,0]*A[1,1]-A[0,1]*A[1,0])/(A[0,0])) )
	#Tp	
        ttp=np.abs(tp)*np.abs(tp) #*sqrt(Y2/Y1) 
	Tp.append(ttp) 
	#Rp+Tp
        RTp.append(rrp+ttp)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')



#Grafica de Norma t	
plt.plot(w, tNorm, label='$|t|$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel('$|t|$', fontsize=16)
plt.title('$|t|$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('t_N.pdf')
#plt.show()


#Grafica de Fase t	
plt.plot(w, tPhase, label='$\phi_t$')
plt.legend()
plt.xlabel('w')
plt.ylabel('Fase')
plt.title('Fase de t')
plt.show()

#Grafica de r
plt.plot(w, rNorm, label='r_N')
plt.legend()
plt.xlabel('w')
plt.ylabel('Norma')
plt.title('Norma de r')    
plt.show()
plt.plot(w, rPhase, label='r_phase')
plt.legend()
plt.xlabel('w')
plt.ylabel('Fase')
plt.title('Fase de r')
plt.show()


#Grafica de tp	
plt.plot(w, tpNorm, label='tp_N')
plt.legend()
plt.xlabel('w')
plt.ylabel('Norma')
plt.title('Norma de tp')    
plt.show()
plt.plot(w, tpPhase, label='tp_phase')
plt.legend()
plt.xlabel('w')
plt.ylabel('Fase')
plt.title('Fase de tp')
plt.show()

#Grafica de rp
plt.plot(w, rpNorm, label='rp_N')
plt.legend()
plt.xlabel('w')
plt.ylabel('Norma')
plt.title('Norma de rp')    
plt.show()
plt.plot(w, rpPhase, label='rp_phase')
plt.legend()
plt.xlabel('w')
plt.ylabel('Fase')
plt.title('Fase de rp')
plt.show()


#Grafica de T
plt.plot(w, T, label='T')
plt.legend()
plt.xlabel('w')
plt.ylabel('T')
plt.title('T')    

#Grafica de R
plt.plot(w, R, label='R')
plt.legend()
plt.xlabel('w')
plt.ylabel('R')
plt.title('R')    
#plt.show()


#Grafica de T+R
plt.plot(w, RT, label='R+T')
plt.legend()
plt.xlabel('w')
plt.ylabel('R+T')
plt.title('R+T')   
plt.show()

plt.plot(w,np.log(RT))
plt.show()


#Grafica de Tp
plt.plot(w, Tp, label='Tp')
plt.legend()
plt.xlabel('w')
plt.ylabel('Tp')
plt.title('Tp')    

#Grafica de Rp
plt.plot(w, Rp, label='Rp')
plt.legend()
plt.xlabel('w')
plt.ylabel('Rp')
plt.title('Rp')    
plt.plot(w,RTp)
plt.show()

