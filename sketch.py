import matplotlib.pyplot as plt
import numpy as np

class Wavefunction:
    def __init__(self, position_grid):
        self.values = np.ones_like(position_grid, dtype=np.complex128)
        self.normalize()
        self.mass = 1.

    def normalize(self):
        self.norm = np.linalg.norm(self.values, 2)
        if self.norm == 1.:
            return
        self.values /= self.norm
        self.norm = 1.

class Potential:
    def __init__(self, type='tophat', **kwargs):
        if type=='tophat':
            self.V = self.tophat
            self.start = kwargs['start']
            self.end = kwargs['end']

    def tophat(self, r):
        if r>=self.start and r<=self.end:
            return 1.

class System:
    def __init__(self):
        self.grid = np.linspace(0., 1., 100)
        self.potential = Potential(type='tophat', start=self.grid[0], end=self.grid[-1])
        self.psi = Wavefunction(self.grid)
        self.hbar = 1.
    
    def Hamiltonian(self):
        # Kinetic energy contribution
        kinetic_energy = np.zeros_like(self.grid, dtype=np.complex128)
        dx_inverse_sq = (self.grid[1]-self.grid[2])**(-2)
        coefficient = self.hbar**2/(2.*self.psi.mass)
        # Interior contribution
        for i in range(1, self.grid.shape[0]-1):
            psi_xx = (self.psi.values[i-1]-2.*self.psi.values[i]+self.psi.values[i+1])*dx_inverse_sq
            kinetic_energy[i] = coefficient*psi_xx
        # Boundary contribution
        psi_xx = (self.psi.values[0]-2.*self.psi.values[1]+self.psi.values[2])*dx_inverse_sq
        kinetic_energy[0] =  coefficient*psi_xx
        psi_xx = (self.psi.values[-3]-2.*self.psi.values[-2]+self.psi.values[-1])*dx_inverse_sq
        kinetic_energy[-1] = coefficient*psi_xx

        # Potential Energy contribution
        potential_energy = np.zeros_like(self.grid, dtype=np.complex128)
        for i in range(self.grid.shape[0]):
            potential_energy[i] = self.potential.V(self.grid[i])
        
        return kinetic_energy + potential_energy
        

    
if __name__ == '__main__':
    sim = System()
    print(sim.Hamiltonian())
