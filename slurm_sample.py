# Hello! This python script will be more of what you normally see research wise for this project. It will work much like our sampling we did in the notebooks.

#import all neccesary modules
import numpy as np                
from rbvfit import model as m 
import dynesty

wave = np.loadtxt('wave_cloudy.csv',delimiter=',')
fnorm = np.loadtxt('flux_cloudy.csv',delimiter=',')  #loads the flux and wavelength files we made in QAS_Fitting_w_CLOUDY.ipynb 

enorm = enorm = 0.02 * np.ones((len(wave))) # fake error
redshift = np.array([0.0])
zabs = redshift

lambda_rest  = np.array([1845.5202, 1526.7066,2026.4768,2796.352])  #rest wavelengths for SiI,SiII,MgI,MgII respectively

elmnt = np.array(['Si','Si','Mg','Mg']) #These two arrays must line up 1 to 1. Together they read as: Si I, Si II, Mg I, Mg II
ion = np.array([1,2,1,2]) 




def frac_from_dens(redshift,density,temp,species,ion):  #useful function to produce ionization fraction
    
    cloudy = np.load('cloudy_table.npz')   #load table

    table = cloudy['table']
    #[redshift,density,temperature,species,ion]#
    
    rrange = np.arange(7,-1,-1) #redshift range from S Bird Cloudy tables
    trange = np.arange(3.,8.6,0.05) #temp range from S Bird Cloudy tables
    drange = np.arange(-7.,4.,0.2) #drange range from S Bird Cloudy tables
    
    spec = np.array(["H", "He", "C", "N", "O", "Ne", "Mg", "Si", "Fe"])
    
    trange = trange.round(3) #formatting necessary to read table

    r = np.where(rrange == redshift)
    t = np.where(trange == temp)
    el = np.where(spec == species)
    
    r = int(r[0])
    t = int(t[0])
    el = int(el[0])
    
    ion = ion-1 #account for formatting of np arrays
    
    frac = []
    for i in range (len(drange)):
        frac.append(table[r,i,t,el,ion])
    frac = np.array(frac)
    frac = 10**(frac)
    
    frac_dens = np.interp(density,drange,frac)
    
    return frac_dens



def gen_logN_arr(log_rho,log_los,log_temp,Met,elmnt,ion):
    """
    This function generates a logN array based off the set of samples dynesty will feed it. We will call it in our gen_model function 
    ____________________________________________________________________________________________________________________________________
    
        log(rho) = array of log 3D density values
        log_los = array of log of the length of line of sight values
        log_temp = array of log temperature values
        Met = array of metallicity values
        elmnt = array of element names. Must be given as strings in their atomic symbols. i.e. elmnt = ['Si','Mg'] for silicon and magnesium
        ion = array of corresponding ion numbers 
    
    """
    N_tot = 10**(log_rho)*10**(log_los)
    log_temp = np.round(log_temp,1)
    f_solar_H = 1
    f_solar_He = 8.33*10**(-2)
    f_solar_C =  10**(-3.45)
    f_solar_N = 10**(-3.95)
    f_solar_O = 10**(-3.28)
    f_solar_Ne = 10**(-3.95)
    f_solar_Mg = 10. ** (-9.22)
    f_solar_Si = 3.24 * 10.**(-4.48)
    f_solar_Fe = 10.**(-4.53)
    
    logN_arr = []
    for i in range (len(elmnt)):
        N = 0
        
        if (elmnt[i] == 'H'):
            f_ion = f_solar_H
        elif (elmnt[i] == 'He'):
            f_ion = f_solar_He
        elif (elmnt[i] == 'C'):
            f_ion = f_solar_C
        elif (elmnt[i] == 'N'):
            f_ion = f_solar_N
        elif (elmnt[i] == 'O'):
            f_ion = f_solar_O
        elif (elmnt[i] == 'Ne'):
            f_ion = f_solar_Ne
        elif (elmnt[i] == 'Mg'):
            f_ion = f_solar_Mg
        elif (elmnt[i] == 'Si'):
            f_ion = f_solar_Si
        elif (elmnt[i] == 'Fe'):
            f_ion = f_solar_Fe
        else:
            continue
        for j in range (n_clouds):
            f = frac_from_dens(0,log_rho[j],log_temp[j],elmnt[i],ion[i]) * f_ion
            N += N_tot * f * Met
        lgN = np.log10(N)
        logN_arr.append(lgN)
    
    logN_arr = np.array(logN_arr)
    logN_arr = logN_arr.ravel() #makes array 1-D
    logN_arr = logN_arr.reshape(len(lambda_rest),n_clouds)
    return logN_arr





def compile_model(log_rho,log_los,v,log_temp,b,Met,elmnt,ion):
    
    """ 
    This function takes the logN array as well as velo and doppler b to create our flux model at each rest wavelength
    """
    
    logN_arr = gen_logN_arr(log_rho,log_los,log_temp,Met,elmnt,ion)
    outflx = []
    for i in range (len(lambda_rest)):
        for j in range (n_clouds):
            theta = np.array([logN_arr[i,j],b[j],v[j]])
            nclump = n_clouds
            s= m.create_voigt(zabs,[lambda_rest[i]],nclump,1,FWHM = '6.5',verbose=True)
            flx = s.model_flux(theta,wave)
            outflx.append(flx)
    tau = -np.log(outflx)
    op_depth_arr = np.sum(tau,axis =0)
    model = np.exp(-op_depth_arr)
    
    return model

def prior_transform(utheta):  #defines our prior range and converts our paramters to (log_rho,log_lost,v,log_temp,b,Met)
    urho = utheta[0:n_clouds]
    ulos = utheta[n_clouds:2*n_clouds]
    uv = utheta[2*n_clouds:3*n_clouds]
    ut = utheta[3*n_clouds:4*n_clouds]
    ub = utheta[4*n_clouds:5*n_clouds]
    um = utheta[5*n_clouds:6*n_clouds]
    
    rho = 8.*urho -3.    #log rho range is from -3 to 5
    los = 7.*ulos + 18.  #log los range is from 18 to 25
    v = 400.* uv - 200.  #velo range is from -200 to 200
    t = 5.*ut + 3.       #log temp range is from 3 to 8
    b = 20.*ub + 5       # doppler b range is from 5 to 25
    m = 1. * um          # mettalicity range is from 0 to 1
    
    
 
    return np.concatenate([rho, los, v, t, b, m])



def lnlikelihood(theta):  #checks the samples with known values
    model = compile_model(theta[0:n_clouds],theta[n_clouds:2*n_clouds],theta[2*n_clouds:3*n_clouds],theta[3*n_clouds:4*n_clouds],theta[4*n_clouds:5*n_clouds],theta[5*n_clouds:6*n_clouds],elmnt,ion)
    inv_sigma2 = 1.0/(enorm**2 )
    all_c2_values = (fnorm - model) ** 2 * inv_sigma2  - np.log(inv_sigma2)
    return -0.5 * (np.sum( all_c2_values  )) 


# Now we are ready to sample. Lets start with just 1 cloud
n_clouds = 1 #1 cloud sample
nclump = n_clouds
ntransition=len(lambda_rest)
s= m.create_voigt(zabs,lambda_rest,nclump,ntransition,FWHM = '6.5',verbose=True)
dsampler = dynesty.DynamicNestedSampler(lnlikelihood, prior_transform, ndim=6*n_clouds, bound='multi', sample='auto')
dsampler.run_nested()
dres1 = dsampler.results
samp_1 = dres1.samples #best fit parameters
ev_1 = dres1.logz  #evidence

np.savetxt("slurm_sample_results_1cld_.csv",samp_1,delimiter =',') #saves samples
np.savetxt("slurm_sample_results_1cld_ev_.csv",ev_1,delimiter =',') #saves log evidence values


# This script should now be ready to be sent as a part of a slurm job - check out slurm_sample.sh now


