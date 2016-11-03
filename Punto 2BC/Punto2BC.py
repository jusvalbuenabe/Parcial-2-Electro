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
#x=np.linspace(-d,2*d,num=100)

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
    A= Interfaz(eps0,mu0,e,mu0)*Propagacion(d,f)*Interfaz(e,mu0,2*eps0,mu0)
    
    #t
    t=1/A[0,0]
    tNorm.append(np.abs(t))
    tPhase.append(cm.phase(t))
    #T
    tt=np.abs(t)*np.abs(t)*np.sqrt(2.) #*sqrt(Y2/Y1) 
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
    ttp=np.abs(tp)*np.abs(tp) *(np.sqrt(2.)/2.)#*sqrt(Y2/Y1) 
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
plt.show()

#Grafica de Fase t	
plt.plot(w, tPhase, label='$\phi_t$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel(r'$\phi_t$', fontsize=16)
plt.title(r'$\phi_t$', fontsize=18)
labels=[r'$-2\pi$', r'$- \frac{3 \pi}{2}$', r'$- \pi$', r'$- \frac{ \pi}{2}$', r'$0$', r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2 \pi$']
plt.yticks(np.arange(-2*np.pi,2*np.pi, (np.pi)/2), labels)
plt.grid(True)
plt.axis([0.1*wp,2*wp,-np.pi-0.1,np.pi+0.1])
plt.savefig('t_f.pdf')
plt.show()

#Grafica de Norma de r
plt.plot(w, rNorm, label='|r|')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel('$|r|$', fontsize=16)
plt.title('$|r|$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('r_N.pdf')
plt.show()

#Grafica de Fase r
plt.plot(w, rPhase, label='$\phase_r$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel(r'$\phi_r$', fontsize=16)
plt.title(r'$\phi_r$', fontsize=18)
labels=[r'$-2\pi$', r'$- \frac{3 \pi}{2}$', r'$- \pi$', r'$- \frac{ \pi}{2}$', r'$0$', r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2 \pi$']
plt.yticks(np.arange(-2*np.pi,2*np.pi, (np.pi)/2), labels)
plt.grid(True)
plt.axis([0.1*wp,2*wp,-np.pi-0.1,np.pi+0.1])
plt.savefig('r_f.pdf')
plt.show()

#Grafica de Norma tp	
plt.plot(w, tpNorm, label=r'$t^\prime$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel(r'$|t^\prime|$', fontsize=16)
plt.title(r'$|t^\prime|$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.3])
plt.grid(True)
plt.savefig('tp_N.pdf')
plt.show()

#Grafica de Fase de tp
plt.plot(w, tpPhase, label=r'$\phi_{t^\prime}$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel(r'$\phi_{t^\prime}$', fontsize=16)
plt.title(r'$\phi_{t^\prime}$', fontsize=18)
labels=[r'$-2\pi$', r'$- \frac{3 \pi}{2}$', r'$- \pi$', r'$- \frac{ \pi}{2}$', r'$0$', r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2 \pi$']
plt.yticks(np.arange(-2*np.pi,2*np.pi, (np.pi)/2), labels)
plt.grid(True)
plt.axis([0.1*wp,2*wp,-np.pi-0.1,np.pi+0.1])
plt.savefig('tp_f.pdf')
plt.show()


#Grafica de Norma rp
plt.plot(w, rpNorm, label=r'$|r^\prime|$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel(r'$|r^\prime|$', fontsize=16)
plt.title(r'$|r^\prime|$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('rp_N.pdf')
plt.show()

#Grafica de Fase de rp
plt.plot(w, rpPhase, label=r'$\phi_{r^\prime}$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel(r'$\phi_{r^\prime}$', fontsize=16)
plt.title(r'$\phi_{r^\prime}$', fontsize=18)
labels=[r'$-2\pi$', r'$- \frac{3 \pi}{2}$', r'$- \pi$', r'$- \frac{ \pi}{2}$', r'$0$', r'$\frac{\pi}{2}$',r'$\pi$',r'$\frac{3\pi}{2}$',r'$2 \pi$']
plt.yticks(np.arange(-2*np.pi,2*np.pi, (np.pi)/2), labels)
plt.grid(True)
plt.axis([0.1*wp,2*wp,-np.pi-0.1,np.pi+0.1])
plt.savefig('rp_f.pdf')
plt.show()


#Grafica de T
plt.plot(w, T, label='$T$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel('$T$', fontsize=16)
plt.title('$T$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('T.pdf')
plt.show()

#Grafica de R
plt.plot(w, R, label='$R$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel('$R$', fontsize=16)
plt.title('$R$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('R.pdf')
plt.show()

#Grafica de T+R
plt.plot(w, T, label=r'$T$')
plt.plot(w, R, label=r'$R$')
plt.plot(w, RT, label=r'$T+R$')
plt.legend(loc=2)
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
#plt.ylabel('$T+R$', fontsize=16)
plt.title('$T+R$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('RT.pdf')
plt.show()


#Grafica de Tp
plt.plot(w, Tp, label='$T^\prime$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel('$T^\prime$', fontsize=16)
plt.title('$T^\prime$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('Tp.pdf')
plt.show()

#Grafica de Rp
plt.plot(w, Rp, label='$R^\prime$')
#plt.legend()
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
plt.ylabel('$R^\prime$', fontsize=16)
plt.title('$R^\prime$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('Rp.pdf')
plt.show()

#Grafica de Tp+Rp
plt.plot(w, Tp, label=r'$T^\prime$')
plt.plot(w, Rp, label=r'$R^\prime$')
plt.plot(w, RTp, label=r'$T^\prime+R^\prime$')
plt.legend(loc=2)
plt.xlabel(r'$\frac{\omega}{\omega_p}$', fontsize=16)
#plt.ylabel('$Tp+Rp$', fontsize=16)
plt.title(r'$T^\prime+R^\prime$', fontsize=18) 
plt.axis([0.1*wp,2*wp,-0.1,1.1])
plt.grid(True)
plt.savefig('RTp.pdf')
plt.show()
