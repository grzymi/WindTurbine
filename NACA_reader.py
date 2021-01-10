import numpy as np
from bs4 import BeautifulSoup
import requests

r = requests.get('http://airfoiltools.com/polar/details?polar=xf-n0009sm-il-1000000')
soup = BeautifulSoup(r.text, "html.parser") #Parsowanie strony
table = soup.find('table', attrs={'class': 'tabdata'}) #znalezienie odpowiedniej tabeli,
                                                        #po uprzednim sprawdzeniu nazwy tabeli w źródle strony
no_rows = table.find_all('tr')
#print ('Row number:', len(no_rows)-1) #Wyświetla liczbę wierszy, umniejszoną o wiersz nagłówkowy

cd_cl_data = np.zeros((len(no_rows), 3)) #Tworzy tablicę o żądanym rozmiarze
#print (cd_cl_data)
for i, row in enumerate(no_rows):
#for row in no_rows:
    data = row.find_all('td') #znalezienie kolumn
    if len(data)==0:
        continue   #for i in data:
    #    print(i.getText())
    cd_cl_data[i][0] = data[0].getText()
    cd_cl_data[i][1] = data[1].getText()
    cd_cl_data[i][2] = data[2].getText()
    #print(data)
    #print(len(data))
    #print(data[0])

print(cd_cl_data)

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
    

#print(cd_cl('0009', -18))