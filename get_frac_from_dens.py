import numpy as np
import matplotlib.pyplot as plt




def frac_from_dens(redshift,denisty,temp,species,ion):
    
    cloudy = np.load('cloudy_table.npz')

    table = cloudy['table']
    #[redshift,density,temperature,species,ion]#
    
    rrange = np.arange(7,-1,-1) #redshift range from S Bird Cloudy tables
    trange = np.arange(3.,8.6,0.05) #temp range from S Bird Cloudy tables
    drange = np.arange(-7.,4.,0.2) #drange range from S Bird Cloudy tables
    
    spec = np.array(["Hydrogen", "Helium", "Carbon", "Nitrogen", "Oxygen", "Neon", "Magnesium", "Silicon", "Iron"])
    
    trange = trange.round(3)

    r = np.where(rrange == redshift)
    t = np.where(trange == temp)
    elmnt = np.where(spec == species)
    
    r = int(r[0])
    t = int(t[0])
    elmnt = int(elmnt[0])
    
    ion = ion-1 #account for formatting of np arrays
    
    frac = []
    for i in range (len(drange)):
        frac.append(table[r,i,t,elmnt,ion])
    frac = np.array(frac)
    frac = 10**(frac)
    
    frac_dens = np.interp(denisty,drange,frac)
    
    return frac_dens

def romanise_num(num):
    """Turn a number into a roman numeral (very badly)"""
    if 1 <= num <= 3:
        return "I"*num
    elif num == 4:
        return "IV"
    elif num == 5:
        return "V"
    elif num == 6:
        return "VI"
    elif num == 7:
        return "VII"
    elif num == 8:
        return "VIII"
    elif num == 9:
        return "IX"
    elif num== 10:
        return "X"
    
    else:
        raise RuntimeError("Not implemented")




def plot_frac_vs_dens(redshift,temp,species,ion):
    
    cloudy = np.load('cloudy_table.npz')

    table = cloudy['table']
    #[redshift,density,temperature,species,ion]#
    
    rrange = np.arange(7,-1,-1) #redshift range from S Bird Cloudy tables
    trange = np.arange(3.,8.6,0.05) #temp range from S Bird Cloudy tables
    drange = np.arange(-7.,4.,0.2) #drange range from S Bird Cloudy tables
    
    spec = np.array(["Hydrogen", "Helium", "Carbon", "Nitrogen", "Oxygen", "Neon", "Magnesium", "Silicon", "Iron"])
    
    trange = trange.round(3)

    r = np.where(rrange == redshift)
    t = np.where(trange == temp)
    elmnt = np.where(spec == species)
    
    r = int(r[0])
    t = int(t[0])
    elmnt = int(elmnt[0])
    
    ion = ion-1 #account for formatting of np arrays
    
    frac = []
    for i in range (len(drange)):
        frac.append(table[r,i,t,elmnt,ion])
    frac = np.array(frac)
    frac = 10**(frac)
    
    ion+= 1
    ion = romanise_num(ion)
    plt.plot(drange,frac)
    plt.xlabel(r"Log$\rho_\mathrm{Si}$ (cm$^{-3}$)")
    plt.ylabel("Ionization Fraction")
    plt.title("Ionization Fraction for "+species+"-"+ion+" at redshift "+str(redshift)+" and Temperature Log("+str(temp)+") K$")
    plt.show()
    
    
    
    
