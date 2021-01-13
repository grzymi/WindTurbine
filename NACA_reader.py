import numpy as np
from bs4 import BeautifulSoup
import requests

def cd_cl(NACA, angle):
    '''
    NACA - put the website with profile from airfoiltools.com
    Examples:
        http://airfoiltools.com/polar/details?polar=xf-naca0012h-sa-1000000
        http://airfoiltools.com/polar/details?polar=xf-naca0006-il-1000000
    '''
    r = requests.get(NACA)
    soup = BeautifulSoup(r.text, "html.parser") #Parsowanie strony
    table = soup.find('table', attrs={'class': 'tabdata'}) #znalezienie odpowiedniej tabeli,
                                                            #po uprzednim sprawdzeniu nazwy tabeli w źródle strony
    no_rows = table.find_all('tr')
    #print ('Row number:', len(no_rows)-1) #Wyświetla liczbę wierszy, umniejszoną o wiersz nagłówkowy

    cd_cl_data = np.zeros((len(no_rows), 3)) #Tworzy tablicę o żądanym rozmiarze

    for i, row in enumerate(no_rows):
        data = row.find_all('td') #znalezienie kolumn
        if len(data)==0:
            continue   #for i in data:
        cd_cl_data[i][0] = data[0].getText()
        cd_cl_data[i][1] = data[1].getText()
        cd_cl_data[i][2] = data[2].getText()
    
    cd_cl_data = np.delete(cd_cl_data, 0, 0) #Usunięcie pierwszego wiersza z tablicy
    angles = cd_cl_data[0:-1,0]
    indices = np.abs(np.subtract.outer(angles, angle)).argmin(0)
    x=np.where(cd_cl_data==cd_cl_data[indices][0])
    Cd = cd_cl_data[x[0][0]][2]
    Cl = cd_cl_data[x[0][0]][1]
    return [Cd, Cl]

print (cd_cl('http://airfoiltools.com/polar/details?polar=xf-naca0006-il-1000000', 8))

