import numpy as np
import matplotlib.pyplot as plt
import json
from astropy import units as u
from astropy.modeling import models
import scipy.stats as sc
import sys
import os
import time
import pickle


'''
Inputs:

sys.argv[1] = LC data file (string)
sys.argv[2] = SN info file (string)

'''


def P15(Mc,Me,Re,kappa,Esn,t):
	'Apply the Piro15 model --> calculate luminosity from inputs and output temperature and radius'
	'                           as a function of time' 
	'Input: Core mass (Mc) and envelope mass (Me) in solar mass units'
	'       Envelope radius (Re) in solar radius units'
	'       Opacity (kappa) in cm^2/g units'
	'       Kinetic energy of the supernova (Esn) in ergs'
	'       Time array (t) in seconds'
	'Returns: Temperature and radius as a function of time'
	
	#constants:
	sigma_SB=5.67051e-5  #erg cm^-2 s^-1 K^-4
	Rsol2cm=6.9599e10  #cm/Rsol
	
	#calculated values from input parameters:
	E51=Esn/1e51   #unitless
	ve=2e9*(E51**.5)*(Mc**-.35)*((Me/.01)**-.15) #cm/s
	te_s=(Re*Rsol2cm)/ve    #s
	Ee=4e49*E51*(Mc**-.7)*((Me/.01)**.7)  #ergs
	tp_d=.9*((kappa/0.34)**.5)*(E51**(-1./4.))*(Mc**.17)*((Me/.01)**.57)  #days
	tp_s=tp_d*60.*60.*24.  #s
	
	L=((te_s*Ee)/(tp_s**2))*np.exp(-1.*((t*(t+2*te_s))/(2*(tp_s**2))))  # ergs/s
	R=(Re*Rsol2cm)+(ve*t)  #cm
	T=(L/(4.*np.pi*(R**2)*sigma_SB))**(1./4.)  #K
	
	return T,R,L


def synth_phot(filter_file,zero_pt,Re,D,T,R):
	'Synthetic photometry --> for each time step calculate blackbody function, scale it, integrate for flux'
	'                         within the filter, convert flux to absolute magnitude'
	'Inputs: File with optical filter info (filter_file)'
	'        Zero point value in Jy'
	'        Envelope radius (Re) in solar radius units'
	'        Temperature as a function of time in units Kelvin'
	'Returns: Absolute magnitude as function of time'
	
	#constants:
	pc2cm=10.*3.086e18  #cm
	
	filt=np.loadtxt(filter_file) #first column = wavelength (angstrom), second column = filter amplitude
	mag=np.empty(len(T))
	for i in range(len(T)):
		model=models.BlackBody(temperature=T[i]*u.K)
		bb=model(filt[:,0]*u.angstrom)
		scaled_bb=((4*(np.pi*u.steradian)*np.pi*(R[i]*u.cm)**2)/(4*np.pi*(pc2cm*u.cm)**2))*bb
		flux=np.trapz(scaled_bb*filt[:,1],filt[:,0]*u.angstrom,axis=-1)/np.trapz(filt[:,1],filt[:,0]*u.angstrom,axis=-1)
		mag[i]=-2.5*np.log10((flux)/(zero_pt*((u.erg)/(u.cm**2*u.s*u.Hz))*1e-23))
	
	return mag


def synth_phot_sing(filter_file,zero_pt,Re,D,T,R):
	'Synthetic photometry --> for each time step calculate blackbody function, scale it, integrate for flux'
	'                         within the filter, convert flux to absolute magnitude'
	'Inputs: File with optical filter info (filter_file)'
	'        Zero point value in Jy'
	'        Envelope radius (Re) in solar radius units'
	'        Temperature as a function of time in units Kelvin'
	'Returns: Absolute magnitude as function of time'
	
	#constants:
	pc2cm=10.*3.086e18  #cm
	
	filt=np.loadtxt(filter_file) #first column = wavelength (angstrom), second column = filter amplitude
	model=models.BlackBody(temperature=T*u.K)
	bb=model(filt[:,0]*u.angstrom)
	scaled_bb=((4*(np.pi*u.steradian)*np.pi*(R*u.cm)**2)/(4*np.pi*(pc2cm*u.cm)**2))*bb
	flux=np.trapz(scaled_bb*filt[:,1],filt[:,0]*u.angstrom,axis=-1)/np.trapz(filt[:,1],filt[:,0]*u.angstrom,axis=-1)
	mag=-2.5*np.log10((flux)/(zero_pt*((u.erg)/(u.cm**2*u.s*u.Hz))*1e-23))
	
	return mag


print('Reading in data and info files...')

#read in LC data file:
with open('LC_files/SN2016gkg.json','r') as read_LCs:
    data=json.load(read_LCs)


#read in SN info file:
with open('LC_files/SN2016gkg_info.json','r') as read_info:
    info=json.load(read_info)


SNe=os.path.split('LC_files/SN2016gkg.json')[1].split('.')[0]
band=info['Band']

if not os.path.exists('grid_search/{}_manual'.format(SNe)):
    os.makedirs('grid_search/{}_manual'.format(SNe))


#Extract LC data from json file
print('Collecting SN data...')
band_data=[]
errband=[]
mjd=[]
for i in range(len(data[SNe]['photometry'])):
	if data[SNe]['photometry'][i]['band'] == band:
		srcs=np.array([int(s) for s in data[SNe]['photometry'][i]['source'].split(',')])
		if 1 in srcs or 4 in srcs:
			band_data.append(np.float(data[SNe]['photometry'][i]['magnitude']))
			mjd.append(np.float(data[SNe]['photometry'][i]['time']))
			try:
				errband.append(np.float(data[SNe]['photometry'][i]['e_magnitude']))
			except:
				errband.append(0.0)
band_data=np.array(band_data)
errband=np.array(errband)
mjd=np.array(mjd)

#Extract SN info from json file:
D=np.float(info['Distance'])*1e6
Derr=np.float(info['Error in distance'])*1e6
exp_date=np.float(info['Lit explosion date'])
stop_pt=np.int(info['Data range'])

if band=='R':
	ext_corr=2.119*np.float(info['E(B-V)'])
	ext_corr_err=np.float(info['dE(B-V)'])
elif band=='V':
	ext_corr=2.682*np.float(info['E(B-V)'])
	ext_corr_err=np.float(info['dE(B-V)'])
else:
	print('Band is not R or V.. cannot correct for extinction.')
	


#convert SN data to absolute magnitude and MJD time to days post-explosion:
band_absmag_orig=band_data-(5*np.log10(D/10.))-ext_corr
errband_mag=np.sqrt((errband)**2+(-5.*(1/(D*np.log(10.)))*Derr)**2+(ext_corr_err)**2)
up_lim=17.360-(5*np.log10(D/10.))-ext_corr
lim_date=57651.165-exp_date
timeband=mjd-exp_date


# make 1000 samples of the data from random normal distribution
band_absmag=np.random.normal(loc=band_absmag_orig,scale=errband_mag,size=(1000,len(band_absmag_orig)))
print(band_absmag.shape)


# exp_check=-3.
# if (lim_date-3.)<=0:
# 	exp_check=-(lim_date-lim_date%5)
# print('Min time shift:',exp_check)

# exp_check=-3.
# if (lim_date)>0:
# 	exp_check=round(lim_date*2)/2
# print('Min time shift:',exp_check)




#grid search for optimal parameters:
Mc=np.float(info['Lit Mc'])
Me_range=np.arange(.008,.013+.0001,.0001)
Re_range=np.arange(10.,300.+5.,5.)
Esn=np.float(info['Lit Esn'])
# exp_range=np.arange(exp_check,3.+.5,.5)
exp_range=np.arange(round(lim_date*2)/2,round(lim_date*2)/2+1.,.5)
print('Limit date:',lim_date)
print('Min time shift:',exp_range[0])
print('Max time shift:',exp_range[-1])

t=timeband*(60.*60.*24.)
t=t[:stop_pt]
tday=t/(60.*60.*24.)

print('Beginning grid search...')
start=time.time()
chi_matrix=np.empty((len(band_absmag[:,0]),len(exp_range),len(Re_range),len(Me_range)))
chi=np.array([1e20])
dist=np.array([1e20])
pos=np.empty(3,dtype='int')
for e in range(len(exp_range)):
	print(e)
	for c in range(len(Re_range)):
# 		if c%10==0:
# 			print(c)
		for b in range(len(Me_range)):
			t_shift=t+(exp_range[e]*60*60.*24.)
			T,R,L=P15(Mc,Me_range[b],Re_range[c],0.34,Esn,t_shift)
			mag=synth_phot(info['Filter file'],info['Zero point'],Re_range[c],D,T,R)

			lim_shift=(lim_date+exp_range[e])*60*60.*24.
			Tlim,Rlim,Llim=P15(Mc,Me_range[b],Re_range[c],0.34,Esn,lim_shift)
			lim_mag=synth_phot_sing(info['Filter file'],info['Zero point'],Re_range[c],D,Tlim,Rlim)

			for i in range(len(band_absmag[:,0])):
				if lim_mag>=up_lim-0.3:
					chi_sq=np.sum((band_absmag[i][:stop_pt]-mag)**2/(errband_mag[:stop_pt]**2))
				else:
					chi_sq=1e20

				chi_matrix[i][e][c][b]=chi_sq

# 			if np.abs(1-chi_sq) < dist and chi_sq > 0.:
# 				dist[0]=np.abs(1-chi_sq)
# 				chi[0]=chi_sq
# 				pos[:]=np.array([e,c,b])


pickle.dump(chi_matrix,open('grid_search/{}_manual/{}_manual_MC_chi_square_matrix.pkl'.format(SNe,SNe),'wb'))
# np.savetxt('grid_search/{}_manual/{}_manual_optimal_params.txt'.format(SNe,SNe),np.array([exp_range[pos[0]],Re_range[pos[1]],Me_range[pos[2]]]))
end=time.time()
print('Finished grid search in',(end-start)/60.,'minutes.')
# print('Minimum chi squared value:',chi,'at',pos)
# print('Shift:',exp_range[pos[0]],', Re:',Re_range[pos[1]],', Me:',Me_range[pos[2]])


print('Extracting best fit values...')
min_chis=np.empty((len(band_absmag[:,0]),2))
for i in range(len(band_absmag)):
	min_pos=np.where(chi_matrix[i]==np.min(chi_matrix[i]))
	min_chis[i]=np.array([Re_range[min_pos[1][0]],Me_range[min_pos[2][0]]])


print('Calculating integrals...')
tday=np.linspace(0.,10.,100)
ts=tday*24.*60.*60.
integrals=np.empty(len(min_chis[:,0]))
for i in range(len(min_chis[:,0])):
	T,R,lum=P15(Mc,min_chis[i][1],min_chis[i][0],0.34,Esn,ts)
# 	print(lum.shape)
	integrals[i]=np.trapz(lum,ts,axis=-1)


print('Mean integral:',np.nanmean(integrals))
print('Std of integrals:',np.nanstd(integrals))


