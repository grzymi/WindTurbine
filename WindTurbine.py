import numpy as np

class WindTurbine:
    def __init__(self, radius, air_velocity,
                 pressure=101325, temperature=15, time=1, rpm=20,
                 no_of_elem=15, NACA='0009', height=[10,100], aoa = 0, 
                 roughness_class=0, eta_a=0.95, eta_m=0.95, n_blades=3):
        self.radius = radius
        self.pressure = pressure
        self.temperature = temperature
        self.time = time
        self.rpm = rpm
        self.air_velocity = air_velocity
        self.no_of_elem = no_of_elem
        self.NACA = NACA
        self.height = height
        self.aoa = aoa
        self.roughness_class = roughness_class
        self.eta_a = eta_a
        self.eta_m = eta_m
        self.n_blades = n_blades
   
    def velocities(self):
        #v0-before turbine, v1=in turbine, v2-after turbine
        #v0 = self.air_velocity
        v0 = self.velocity_profile()
        v2 = (v0/3)
        v1 = (v0+v2)/2
        return (v0, v1, v2)

    def velocity_profile(self):
        #Profil wiatru. Dodac wykladnik jakosci powierzchni do danych wejsciowych.
        #Jako lista lub slownik. Dane sa w notebooku, ktory robilem wczesniej.
        r_class = np.array([0, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4])
        exponent_value = np.array([0.0815, 0.1, 0.1377, 0.1499, 0.165, 0.1786, 0.2139, 0.2677, 0.305])
        indices = np.abs(np.subtract.outer(r_class, self.roughness_class)).argmin(0)
        #print (indices)
        #print (exponent_value[indices])
        v = self.air_velocity*pow((self.height[1]/self.height[0]),exponent_value[indices])
        if v>25:
            print ('Wind speed is to high. The turbine has stopped.')
            v = 0
            return v
        else:
            return v
    
    def area(self):
        A = np.pi*(pow(self.radius, 2))
        return A
    
    def density(self):
        rho = self.pressure/(287*(self.temperature+273.15))
        return rho
    
    def air_mass(self):
        m = self.density()*self.area()*self.velocities()[1]
        return m
    
    def kinetic_energy (self):
        E = (self.air_mass()*pow(self.velocities()[1],2)/2)
        return E
    
    def power_theoretic(self):
        P_t = self.kinetic_energy()/self.time
        return P_t
    
    def rot_speed(self):
        omega = self.rpm*(np.pi/30)
        return omega
    
    def hsd(self):
        #high-speed discriminant
        if self.velocities()[1] == 0:
            hisd = 0
        else:
            hisd = (self.rot_speed()*self.radius)/self.velocities()[1]
        return hisd
    
    def characteristic(self):
        #Characteristic of the wind turbine
        A0, A1, A2, A3, A4 = 0.093368, -0.1838, 0.118605, -0.01773, 0.000756
        cp = (A0*pow(self.hsd(),0)+A1*pow(self.hsd(),1)
              +A2*pow(self.hsd(),2)+A3*pow(self.hsd(),3)
              +A4*pow(self.hsd(),4))
        return cp
    
    def lenght(self):
        l = (self.radius-self.radius*0.1)/self.no_of_elem
        return l
    
    def radius_l(self):
        #radius for selected element on lenght
        #k - number of element
        r_k = np.zeros(self.no_of_elem)
        for i in range(1, self.no_of_elem+1):
            r_k[i-1] = (self.radius*0.1)+self.lenght()*i
        '''
        if k < 1:
            k=1
            print ("One element was take into account.")
        elif k > self.no_of_elem:
            k=self.no_of_elem
            print (self.no_of_elem, 'element was take into account.')
        r_k = (self.radius*0.1)+self.lenght()*k
        '''
        return r_k
    
    def circum_vel(self):
        #circumferential velocity
        u_k = np.zeros(self.no_of_elem)
        for i,n in enumerate(self.radius_l()):
            u_k[i] = (np.pi*n*self.rpm)/30
        return u_k
    
    def relative_vel(self):
        w_k = np.zeros(self.no_of_elem)
        for i,n in enumerate(self.circum_vel()):
            w_k[i] = np.sqrt(pow(self.velocities()[1],2)+pow(n,2))
        return w_k

    def cd_cl(self):
        naca = np.loadtxt('./NACA_values/NACA_'+str(self.NACA)+'.txt', skiprows=12, usecols=(0,1,2))
        angles = naca[0:-1,0]
        indices = np.abs(np.subtract.outer(angles, self.aoa)).argmin(0)
        #print (indices)
        x=np.where(naca==naca[indices][0])
        #print (x)
        Cd = naca[x[0][0]][2]
        Cl = naca[x[0][0]][1]
        #Przyjęto współczynniki dla kąta natarcia równego... - problem z tłumaczeniem
        #Zrobić tak, aby wyświetlało sie tylko raz przy wywołaniu!!!
        #print ('Cd and Cl was taken for AoA equal', naca[indices][0])
        return [Cd, Cl]
    
    def blade_width(self):
        #Cały moduł do sprawdzenia pod kątem poprawności matematycznej
        s_k = np.zeros(self.no_of_elem)
        Cd = self.cd_cl()[0]
        Cl = self.cd_cl()[1]
        radius = self.radius_l()
        rel_vel = self.relative_vel()
        for i in range(self.no_of_elem):
            #licznik = 4*np.pi*self.velocities()[1]*(self.velocities()[0]-self.velocities()[2])
            licznik = 4*np.pi*radius[i]*self.velocities()[1]*(self.velocities()[0]-self.velocities()[2])
            mianownik = self.n_blades*rel_vel[i]*(self.circum_vel()[i]*Cd+self.velocities()[0]*Cl)
            s_k[i] = licznik/mianownik
        return s_k

    def forces(self):
        #Podstawowe siły działające na łopatę
        Cd = self.cd_cl()[0]
        Cl = self.cd_cl()[1]        
        lift_force = np.zeros(self.no_of_elem)
        drag_force = np.zeros(self.no_of_elem)
        aerodynamic_force = np.zeros(self.no_of_elem)
        for i,n in enumerate (self.blade_width()):
            lift_force[i] = Cl*self.density()*n*self.lenght()*(self.relative_vel()[i]/2)
            drag_force[i] = Cd*self.density()*n*self.lenght()*(self.relative_vel()[i]/2)
            aerodynamic_force[i] = np.sqrt(pow(lift_force[i],2)+pow(drag_force[i],2))
        return [lift_force, drag_force, aerodynamic_force]

    def axial_force(self):
        #Wynikowe siły działające na łopatę
        #Siła osiowa
        axial_f = np.zeros(self.no_of_elem)
        #Siła hamująca
        breaking_f = np.zeros(self.no_of_elem)
        #Siła napędzająca
        driving_f = np.zeros(self.no_of_elem)
        #Siła wynikowa, powodująca obrót turbiny
        rotation_f = np.zeros(self.no_of_elem)
        #wywolanie funkcji
        rel_vel = self.relative_vel()
        forces = self.forces()
        for i,n in enumerate(self.circum_vel()):
            axial_f[i] = ((n/rel_vel[i])*forces[0][i])+((self.velocities()[0]/rel_vel[i])*forces[1][i])
            breaking_f[i] = ((n/rel_vel[i])*forces[1][i])
            driving_f[i] = ((self.velocities()[0]/rel_vel[i])*forces[0][i])
            rotation_f[i] = driving_f[i]-breaking_f[i]
        return [axial_f, breaking_f, driving_f, rotation_f]

    def power(self):
        suma = 0
        predkosc = self.circum_vel()
        sila = self.axial_force()[3]
        for i in range(self.no_of_elem):
            suma = suma + (sila[i]*predkosc[i])
        #power = self.characteristic()*self.eta_m*self.n_blades*suma
        power = self.eta_m*self.n_blades*suma
        return power




import time

podzial_lopaty = np.array([15,20,25,30,35,40,45,50,55,60,100])
katy = np.arange(-7,7.25,0.25)
predkosci = np.array([5,6,7,8,9,10,11,12,13])
moc_wiatru = np.zeros((len(predkosci),len(katy)))
moc_turbiny = np.zeros((len(predkosci),len(katy)))
czas_dzialania = np.zeros((len(predkosci),len(katy)))
for i,n in enumerate(predkosci):
    start_outer = time.time()
    for j,m in enumerate(katy):
        start_inner = time.time()
        turbina = WindTurbine(50, n, aoa=m, NACA='0012', no_of_elem=60, roughness_class=2, rpm=15)
        #print (turbina.power_theoretic()/1000000)
        #print (turbina.power()/1000000)
        moc_wiatru[i][j] = round((turbina.power_theoretic()/1000000),2)
        moc_turbiny[i][j] = round((turbina.power()/1000000),2)
        end_inner = time.time()
        czas_dzialania[i][j] = round(end_inner-start_inner,2)
        print ('\tInner iteration', j+1, 'is finished with time', round(end_inner-start_inner,2), 's', j,m)
    end_outer = time.time()
    print ('Outer iteration',i+1,'is finished within',round(end_outer-start_outer,2),'s',i,n)

sprawnosc = (moc_turbiny/moc_wiatru)
print (moc_turbiny)
print (moc_wiatru)
print (sprawnosc)
#print (czas_dzialania)


# tt = WindTurbine(50, 5, aoa=10, NACA='0009', no_of_elem=15, roughness_class=2, rpm=15)
# print (tt.power_theoretic()/1000000)
# print (tt.power()/1000000)

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')

X, Y = np.meshgrid(katy, predkosci)

surf = ax.plot_surface(X, Y, sprawnosc, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()