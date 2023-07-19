import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc 
class wakefield:

    def __init__(self):
        self.e0 = sc.epsilon_0
        self.u0 = sc.mu_0
        self.conduct = 250 *10**7 #conductivity (1/sigma*m)
        self.b = 0.044 #beam pipe radius (m)
        self.q = 1.6 * 10**11 * sc.elementary_charge #charge of beam: number of protons multiplied by charge (C)
        self.c = sc.speed_of_light #speed of light (m/s)
        self.x = np.linspace(0, 550, 10000)
        self.A = self.q/(2*np.pi*self.b)*np.sqrt(self.c/self.conduct)
        self.y = np.zeros_like(self.x)
        dx = self.x[1] - self.x[0]
        for i in range(10):
            mu = 7.49 * (i + 1) # defining centres of proton bunch
            sigma = 0.09 # std
            self.y = self.y + self.gaussian(self.x, mu, sigma) 

        alpha = 0.05
        self.yw_32 = np.zeros_like(self.y)
        self.yw_52 = np.zeros_like(self.y)
        for i in range(len(self.y)):
            q = self.y[i] * dx # q is charge over time need to be multiplied by time diff
            self.yw_32 = self.yw_32 + q * self.z_32(self.x - self.x[i], alpha)
            self.yw_52 = self.yw_52 + q * self.z_52(self.x - self.x[i], alpha)
        
    def E_r(self, x, y, z):
        return -3/4 * self.A * np.sqrt(x**2 + y**2) * self.yw_52 * 1/np.sqrt(4*np.pi*self.e0)
    
    def E_s(self, z):
        return self.A * self.yw_32

    def B_theta(self, x, y, z):
        return -3/4*self.A*np.sqrt(x**2 + y**2) * self.yw_52 * np.sqrt(self.u0 / (4*np.pi))

    def E_x(self, x, y, z):
        return self.E_r(x, y, z) * np.cos(np.arctan2(y, x)) 

    def E_y(self, x, y, z):
        return self.E_r(x, y, z) * np.sin(np.arctan2(y, x)) 

    def B_x(self, x, y, z):
        return self.B_theta(x, y, z) * -np.sin(np.arctan2(y, x))

    def B_y(self, x, y, z):
        return self.B_theta(x, y, z) * np.cos(np.arctan2(y, x))
    
    def z_52(self, z, alpha):
        pw = np.piecewise(z, [z <= alpha], [1./alpha**(5./2.), lambda x: 1./x**(5./2.)])
        pw[z <= 0] = 0
        return pw
    
    def z_32(self, z, alpha):
        pw = np.piecewise(z, [z <= alpha], [1./alpha**(3./2.), lambda x: 1./x**(3./2.)])
        pw[z <= 0] = 0
        return pw
    
    def gaussian(self, x, mu, sigma):
        pw = np.piecewise(x, [np.abs((x-mu)/sigma) > 5], [0, lambda x: 1/(np.sqrt(2*np.pi) * sigma) * np.exp(-(x-mu)**2/(2*sigma**2))])
        return pw