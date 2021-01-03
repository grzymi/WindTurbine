import numpy as np


def cd_cl(NACA, angle):
    naca = np.loadtxt('./NACA_values/NACA_'+str(NACA)+'.txt', skiprows=12, usecols=(0,1,2))
    maximum = naca[-1][0]
    minimum = naca[0][0]
    angles = naca[0:-1,0]
    indices = np.abs(np.subtract.outer(angles, angle)).argmin(0)
    print (angles[indices])
    #print (naca[0:-1,0])

    # if angle>minimum and angle<maximum:


    # try:
    #     x=np.where(naca==angle)
    #     Cd = naca[x[0][0]][2]
    #     Cl = naca[x[0][0]][1]
    #     return (Cd, Cl)
    # except:
    #     if angle < 0:
    #         x=np.where(naca==minimum)
    #         Cd = naca[x[0][0]][2]
    #         Cl = naca[x[0][0]][1]
    #         print ('Angle is out of range. The minimum value was taken:', minimum)
    #         return (Cd, Cl)
    #     else:
    #         x=np.where(naca==maximum)
    #         Cd = naca[x[0][0]][2]
    #         Cl = naca[x[0][0]][1]
    #         print ('Angle is out of range. The maximum value was taken:', maximum)
    #         return (Cd, Cl)            
    #print(x[0][0])
    # Cd -> index 2
    # Cl -> index 1
    #print(naca[x[0][0]][2])
    

print(cd_cl('0009', -18))