import numpy as np


def cd_cl(NACA, angle):
    naca = np.loadtxt('./NACA_values/NACA_'+str(NACA)+'.txt', skiprows=12, usecols=(0,1,2))
    x=np.where(naca==angle)
    #print(x[0][0])
    # Cd -> index 2
    # Cl -> index 1
    #print(naca[x[0][0]][2])
    Cd = naca[x[0][0]][2]
    Cl = naca[x[0][0]][1]
    return (Cd, Cl)

print(cd_cl('0009', 5)[0])