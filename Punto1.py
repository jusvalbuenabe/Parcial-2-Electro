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
E1pNorm=[] #Norma E1+
E1pPhase=[] #Fase E1+
E1mNorm=[] #Norma E1-
E1mPhase=[] #Fase E1-

#E1=[]
r=[]
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
Rp=[]
Tp=[]
RT=[]
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
	
	#RT
	#print rr
	#print tt
	#print rr+tt
	RT.append(rr+tt)
	
	#rp
	rp=-1.*A[0,1]*t
	rpNorm.append(np.abs(rp))
	rpPhase.append(cm.phase(rp))
	#Rp
	Rp.append(np.abs(rp)*np.abs(rp))
	
	#tp
	tp=A[1,1] + r*rp/t
	tpNorm.append(np.abs(tp))
	tpPhase.append(cm.phase(tp))
	#Tp
	Tp.append(np.abs(tp)*np.abs(tp))
	
	
	
		
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

#Grafica de t	
plt.plot(w, tNorm, label='t_N')
plt.legend()
plt.xlabel('w')
plt.ylabel('Norma')
plt.title('Norma de t')    
#plt.show()
plt.plot(w, tPhase, label='t_phase')
plt.legend()
plt.xlabel('w')
plt.ylabel('Fase')
plt.title('Fase de t')
#plt.show()

#Grafica de r
plt.plot(w, rNorm, label='r_N')
plt.legend()
plt.xlabel('w')
plt.ylabel('Norma')
plt.title('Norma de r')    
#plt.show()
plt.plot(w, rPhase, label='r_phase')
plt.legend()
plt.xlabel('w')
plt.ylabel('Fase')
plt.title('Fase de r')
#plt.show()


#Grafica de tp	
plt.plot(w, tpNorm, label='tp_N')
plt.legend()
plt.xlabel('w')
plt.ylabel('Norma')
plt.title('Norma de tp')    
#plt.show()
plt.plot(w, tpPhase, label='tp_phase')
plt.legend()
plt.xlabel('w')
plt.ylabel('Fase')
plt.title('Fase de tp')
#plt.show()

#Grafica de rp
plt.plot(w, rpNorm, label='rp_N')
plt.legend()
plt.xlabel('w')
plt.ylabel('Norma')
plt.title('Norma de rp')    
#plt.show()
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
#print RT
#print E1
#print E1pn
#plt.plot(w, E1pNorm)
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
#plt.show()
