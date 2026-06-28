import numpy as np
import matplotlib.pyplot as plt

def SSA_lc(t,p,F_tp,tp):

	m=0.88
	n=10.3
	a=2.*m+0.5
	b=(p+5.-6.*m)/2.
	v=1.
	vp=1.

	F_t=F_tp*1.582*((t/tp)**a)*((v/vp)**(5/2))*(1-np.exp(-((t/tp)**(-(a+b)))*((v/vp)**(-(p+4)/2))))

	return F_t

def SSA_sed(v,p,F_vp,vp):

	m=0.88
	n=10.3
	a=2.*m+0.5
	b=(p+5.-6.*m)/2.
	t=1.
	tp=1.

	F_t=F_vp*1.582*((t/tp)**a)*((v/vp)**(5/2))*(1-np.exp(-((t/tp)**(-(a+b)))*((v/vp)**(-(p+4)/2))))

	return F_t

def SSA_lc_fixp(t,F_tp,tp):

	m=0.88
	n=10.3
	p=2.8
	a=2.*m+0.5
	b=(p+5.-6.*m)/2.
	v=1.
	vp=1.

	F_t=F_tp*1.582*((t/tp)**a)*((v/vp)**(5/2))*(1-np.exp(-((t/tp)**(-(a+b)))*((v/vp)**(-(p+4)/2))))

	return F_t

#def SSA_FFA_lc(t,p,F_tp,tp,M_dot,Zavg,mu_e,R,v_w,T):
#
#	m=0.88
#	n=10.3
#	a=2.*m+0.5
#	b=(p+5.-6.*m)/2.
#	v=1.
#	vp=1.
#
#	F_t_SSA=F_tp*1.582*((t/tp)**a)*((v/vp)**(5/2))*(1-np.exp(-((t/tp)**(-(a+b)))*((v/vp)**(-(p+4)/2))))
#	tau_ff=(2.021e25*(M_dot**2)*Zavg)/((mu_e**2)*(v**2.1)*(R**3)*(v_w**2)*(T**1.35))
#	F_t_FFA=F_t_SSA*np.exp(-1.*tau_ff)
#
#	return F_t_FFA
#
#def SSA_FFA_sed(v,p,F_vp,vp,M_dot,Zavg,mu_e,R,v_w,T):
#
#	m=0.88
#	n=10.3
#	a=2.*m+0.5
#	b=(p+5.-6.*m)/2.
#	t=1.
#	tp=1.
#
#	F_t_SSA=F_vp*1.582*((t/tp)**a)*((v/vp)**(5/2))*(1-np.exp(-((t/tp)**(-(a+b)))*((v/vp)**(-(p+4)/2))))
#	tau_ff=(2.021e25*(M_dot**2)*Zavg)/((mu_e**2)*(v**2.1)*(R**3)*(v_w**2)*(T**1.35))
#	F_t_FFA=F_t_SSA*np.exp(-1.*tau_ff)
#	print(tau_ff)
#
#	return F_t_FFA

def SSA_FFA_lc(t,p,F_tp,tp,tcff):

	m=0.88
	n=10.3
	a=2.*m+0.5
	b=(p+5.-6.*m)/2.
	
	F_t_SSA=F_tp*1.582*((t/tp)**a)*(1.-np.exp(-(t/tp)**(-(a+b))))
	F_t_FFA=F_t_SSA*np.exp(-1.*((t/tcff)**(-3.*m)))

	return F_t_FFA

def SSA_FFA_sed(v,p,F_vp,vp):

	m=0.88
	n=10.3
	a=2.*m+0.5
	b=(p+5.-6.*m)/2.

	F_t_SSA=F_vp*1.582*((v/vp)**(5./2.))*(1-np.exp(-((v/vp)**(-(p+4.)/2.))))
	F_t_FFA=F_t_SSA*np.exp(-1.*((v/vp)**(-2.1)))

	return F_t_FFA


def new_L(L_orig,r_orig,r_new):
	
	L_new=(10**L_orig)*((r_new/r_orig)**2)
	log_Lnew=np.log10(L_new)
	
	return log_Lnew


def new_R(L_new,T):
	
	sig=5.67e-8
	L_sol=3.828e26
	R=np.sqrt(((10**L_new)*L_sol)/(4*np.pi*sig*(10**T)**4))
	
	return R


def R(L,T):
	
	lum=(10**L)*3.828e26
	temp=10**T
	sig=5.67e-8
	rad=np.sqrt(lum/(4*np.pi*sig*temp**4))
	
	return rad


def R_sol(L,T):
	
	Rsol=6.9634e8
	lum=(10**L)*3.828e26
	temp=10**T
	sig=5.67e-8
	rad=np.sqrt(lum/(4*np.pi*sig*temp**4))
	rad_sol=rad/Rsol
	
	return rad_sol


def v_esc(M_sol,R_m):
	
	mass=M_sol*1.989e30
	G=6.674e-11
	v_esc=np.sqrt((2*G*mass)/R_m)
	
	return v_esc


def mod2dist_err(err_mod,D):
	
	errD=0.461*D*err_mod
	
	return errD

