import numpy as np
import matplotlib.pyplot as plt
from importlib import reload
from rbvfit import model as m
reload(m)
from get_frac_from_dens import frac_from_dens

#log_rho = np.array([1.2])
#log_los = np.array([20.])
#logN = np.log10(10**(log_rho) * 10**(log_los))
#logT = np.array([4.6])
#Metallicity = np.array([1.0])
#
#lambda_rest = np.array([1206.500,1393.755,1526.7066,1562.002])
#elmnt = np.array(["Silicon",'Silicon','Silicon','Silicon'])
#ion = np.array([3,4,2,1]) in same order as lambda_rest wavelengths
#

#b_Si = np.array([12.])
#v_Si = np.array([50.])
#wave=np.arange(1200,1600,.01)
#
#zabs=np.array([0.0])

""" Above are example input paramters needed to generate spectra, rest of file should take care of the rest.
"""


def generate_params(n_clouds, log_rho,log_los,elmnt,ion,logT):  #just logN for now, will include b as well in future
    """
                This function will calculate LogN , b which will be necessary to produce our spectra, all arrays should be same length, length = n_clouds
                
        _____________________________________________________________________________________________
        
                n_clouds    :   how many clouds you define along this line of sight
                log_rho     :   array of log density of each cloud
                log_los     :   array of the log length of the cloud
                Metallicity :   array of the metallicity of each cloud
                elmnt       :   array of elements present
                ion         :   corresponding ionization states of each element in elmnt_arr
                                e.g.  element = ['Magnesium','Magnesium','Silicon','Silicon']
                                      ion     = [1,2,2,3]
                                        This means Mg-I , Mg-II, Si-II, and Si-III are all present
                zabs        :    redshift of each cloud
                logT        :    log temperature of each cloud
                
       ______________________________________________________________________________________________

    """
    f_solar_mg = 10. ** (-9.22)
    f_solar_si = 10.**(-4.46)
    N_Z = 10**(log_rho)*10**(log_los)*Metallicity
    model = []
    logN_arr = []
    for i in range (len(elmnt)):
            f = frac_from_dens(0,log_rho,logT,elmnt[i],ion[i]) * f_solar_si
            N = N_Z * f
            lgN = np.log10(N)
            logN_arr.append(lgN)
    logN_arr = np.array(logN_arr)
    logN_arr = logN_arr.ravel() #makes array 1-D
    return logN_arr


def compile_model(n_clouds,log_rho,log_los,v,log_temp,b):
    f_solar_si = 10.**(-4.46)
    log_temp = np.round(log_temp,1)
    logN_arr = generate_params(n_clouds,log_rho,log_los,elmnt,ion,log_temp)
    
    outflx = []
    for i in range (len(logN_arr)):
        theta = np.array([logN_arr[i],b_Si[0],v_Si[0]])
        n_clouds = 1
        nclump = n_clouds
        s= m.create_voigt(zabs,[lambda_rest[i]],nclump,1,FWHM = '6.5',verbose=True)
        flx = s.model_flux(theta,wave)
        outflx.append(flx)
    tau = -np.log(outflx)
    op_depth_arr = np.sum(tau,axis =0)
    model = np.exp(-op_depth_arr)
        
    
    return model

            
        
            
        
    
    



    
    
    
    
