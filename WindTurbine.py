import numpy as np

class WindTurbine:
    def __init__(self, radius, air_velocity,
                 pressure=101325, temperature=15, time=1, rpm=20,
                 no_of_elem=15, Cd=0.2, Cl=0.8):
        self.radius = radius
        self.pressure = pressure
        self.temperature = temperature
        self.time = time
        self.rpm = rpm
        self.air_velocity = air_velocity
        self.no_of_elem = no_of_elem
        self.Cd = Cd
        self.Cl = Cl
   
    def velocities(self):
        #v0-before turbine, v1=in turbine, v2-after turbine
        v0 = self.air_velocity
        v2 = (v0/3)
        v1 = (v0+v2)/2
        return (v0, v1, v2)
    
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
    
    #def blade_width(self):
        #s_k = 