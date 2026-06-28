import numpy as np
import matplotlib.pyplot as plt
import json
from astropy import units as u
from astropy.modeling import models
import sys
import os


'''
Inputs:

sys.argv[1] = LC data file (string)
sys.argv[2] = SN info file (string)
sys.argv[3] = Envelope mass in solar mass units (float)
sys.argv[4] = Envelope radius in solar radius units (float)
sys.argv[5] = Time shift (float)

'''


def P15(Mc,Me,Re,kappa,Esn,t):
	'Apply the Piro15 model --> calculate luminosity from inputs and output temperature and radius'
	'                           as a function of time' 
	'Input: Core mass (Mc) and envelope mass (Me) in solar mass units'
	'       Envelope radius (Re) in solar radius units'
	'       Opacity (kappa) in cm^2/g units'
	'       Kinetic energy of the supernova (Esn) in ergs'
	'	Time array (t) in seconds'
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
	
	return T,R


def synth_phot(filter_file,zero_pt,Re,D,T,R):
	'Synthetic photometry --> for each time step calculate blackbody function, scale it, integrate for flux'
	'                         within the filter, convert flux to absolute magnitude'
	'Inputs: File with optical filter info (filter_file)'
	'        Zero point value in Jy'
	'        Envelope radius (Re) in solar radius units'
	'        Temperature as a function of time in units Kelvin'
	'Returns: Absolute magnitude as function of time'
	
	#parsec to centimetre constant:
	pc2cm=10.*3.086e18  #cm
	
	filt=np.loadtxt(filter_file) #first column = wavelength (angstrom), second column = filter amplitude
	mag=np.empty(len(T))
	for i in range(len(T)):
		#make blackbody model and render with angstrom range from filter file
		model=models.BlackBody(temperature=T[i]*u.K)
		bb=model(filt[:,0]*u.angstrom)

		#scale the blackbody
		scaled_bb=((4*(np.pi*u.steradian)*np.pi*(R[i]*u.cm)**2)/(4*np.pi*(pc2cm*u.cm)**2))*bb

		#integrate for the flux
		flux=np.trapz(scaled_bb*filt[:,1],filt[:,0]*u.angstrom,axis=-1)/np.trapz(filt[:,1],filt[:,0]*u.angstrom,axis=-1)

		#calculate magnitude
		mag[i]=-2.5*np.log10((flux)/(zero_pt*((u.erg)/(u.cm**2*u.s*u.Hz))*1e-23))
	
	return mag


#read in LC data file:
with open(sys.argv[1],'r') as read_LCs:
    data=json.load(read_LCs)


#read in SN info file:
with open(sys.argv[2],'r') as read_info:
    info=json.load(read_info)
	

SNe=os.path.split(sys.argv[1])[1].split('.')[0]
band=info['Band']

band_data=[]
errband=[]
mjd=[]
for i in range(len(data[SNe]['photometry'])):
	if data[SNe]['photometry'][i]['band'] == band:
		band_data.append(np.float(data[SNe]['photometry'][i]['magnitude']))
		mjd.append(np.float(data[SNe]['photometry'][i]['time']))
		try:
			errband.append(np.float(data[SNe]['photometry'][i]['e_magnitude']))
		except:
			errband.append(0.0)
band_data=np.array(band_data)
errband=np.array(errband)
mjd=np.array(mjd)


#convert SN data to absolute magnitude and MJD time to days post-explosion:
D=np.float(info['Distance'])*1e6
Derr=np.float(info['Error in distance'])*1e6
exp_date=np.float(info['Lit explosion date'])

if band=='R':
	ext_corr=2.119*np.float(info['E(B-V)'])
	ext_corr_err=np.float(info['dE(B-V)'])
elif band=='V':
	ext_corr=2.682*np.float(info['E(B-V)'])
	ext_corr_err=np.float(info['dE(B-V)'])
elif band=='B':
	ext_corr=3.641*np.float(info['E(B-V)'])
	ext_corr_err=3.641*np.float(info['dE(B-V)'])
else:
	print('Band is not R, V, or B.. cannot correct for extinction.')

band_absmag=band_data-(5*np.log10(D/10.))-ext_corr
errband_mag=np.sqrt((errband)**2+(-5.*(1/(D*np.log(10.)))*Derr)**2+(ext_corr_err)**2)
time_shift=np.float(sys.argv[5])
timeband=mjd-exp_date+time_shift
print(timeband)

#find where SN data reaches 10 days and min and max magnitude for x and y lim when plotting:
stop_pt=np.int(info['Data range'])
#stop_pt=-1
minmag=np.min(band_absmag[:stop_pt])
maxmag=np.max(band_absmag[:stop_pt])


#make time array for range of early peak (~2 weeks):
#t=np.linspace(100.,12.*60.*60.*24.,100)
t=timeband*(60.*60.*24.)
t=t[:stop_pt]
t=np.linspace(0,t[-1],100)
tday=t/(60.*60.*24.)


#run Piro15 model and synthetic photometry:
Mcore=np.float(info['Lit Mc'])
Esn=np.float(info['Lit Esn'])
T,R=P15(Mcore,np.float(sys.argv[3]),np.float(sys.argv[4]),0.34,Esn,t)
mag=synth_phot(info['Filter file'],info['Zero point'],np.float(sys.argv[4]),D,T,R)

# print(17.360-(5*np.log10(D/10.))-ext_corr)


#Plot supernova data with model:
plt.figure(figsize=(10,10))
plt.errorbar(timeband[:stop_pt],band_absmag[:stop_pt],yerr=errband_mag[:stop_pt],fmt='o',color='green',label=SNe)
plt.plot(0+time_shift,17.360-(5*np.log10(D/10.))-ext_corr,'v',color='green')
plt.plot(tday,mag,'--',color='black',label='Piro15')
plt.xlabel('Time (days post-explosion)')
plt.ylabel('Absolute magnitude {} band'.format(band))
#plt.ylim(minmag-.5,maxmag+.5)
#plt.xlim(0,tday[-1]+2)
plt.gca().invert_yaxis()
plt.legend()
plt.show()







