import numpy as np
import matplotlib.pyplot as plt




def frac_from_temp(redshift,density,temp,species,ion):
    
    cloudy = np.load('cloudy_table.npz')

    table = cloudy['table']
    #[redshift,density,temperature,species,ion]#
    
    rrange = np.arange(7,-1,-1) #redshift range from S Bird Cloudy tables
    trange = np.arange(3.,8.6,0.05) #temp range from S Bird Cloudy tables
    drange = np.arange(-7.,4.,0.2) #drange range from S Bird Cloudy tables
    
    spec = np.array(["Hydrogen", "Helium", "Carbon", "Nitrogen", "Oxygen", "Neon", "Magnesium", "Silicon", "Iron"])
    
    drange = drange.round(3)
    
    r = np.where(rrange == redshift)
    d = np.where(drange == density)
    elmnt = np.where(spec == species)
    
    r = int(r[0])
    d = int(d[0])
    elmnt = int(elmnt[0])
    
    ion = ion-1 #account for formatting of np arrays
    
    frac = []
    for i in range (len(trange)):
        frac.append(table[r,d,i,elmnt,ion])
    frac = np.array(frac)
    frac = 10**(frac)
    
    frac_temp = np.interp(temp,trange,frac)
    
    return frac_temp

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


def plot_frac_vs_temp(redshift,density,species,ion):
    
    cloudy = np.load('cloudy_table.npz')

    table = cloudy['table']
    #[redshift,density,temperature,species,ion]#
    
    rrange = np.arange(7,-1,-1) #redshift range from S Bird Cloudy tables
    trange = np.arange(3.,8.6,0.05) #temp range from S Bird Cloudy tables
    drange = np.arange(-7.,4.,0.2) #drange range from S Bird Cloudy tables
    
    spec = np.array(["Hydrogen", "Helium", "Carbon", "Nitrogen", "Oxygen", "Neon", "Magnesium", "Silicon", "Iron"])
    
    drange = drange.round(3)
    
    r = np.where(rrange == redshift)
    d = np.where(drange == density)
    elmnt = np.where(spec == species)
    
    r = int(r[0])
    d = int(d[0])
    elmnt = int(elmnt[0])
    
    ion = ion-1 #account for formatting of np arrays
    
    frac = []
    for i in range (len(trange)):
        frac.append(table[r,d,i,elmnt,ion])
    frac = np.array(frac)
    frac = 10**(frac)
    
    ion+= 1
    ion = romanise_num(ion)
    plt.plot(trange,frac)
    plt.xlabel(r"Log(T) (K)")
    plt.ylabel("Ionization Fraction")
    plt.title("Ionization Fraction for "+species+"-"+ion+" at redshift "+str(redshift)+" and denisty "+str(density)+" cm$^{-3}$")
    plt.show()
    
    
    
