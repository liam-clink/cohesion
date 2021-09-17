import matplotlib.pyplot as plt
import numpy as np
from itertools import permutations
from itertools import combinations

class Wavefunction:
    def __init__(self, position_grid, particle_number):
        self.values = np.ones((particle_number, *position_grid.shape), \
            dtype=np.complex128)
        self.normalize()
        self.mass = 1.

    def normalize(self):
        self.norm = np.linalg.norm(self.values, ord=2)
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
    def __init__(self, part_num):
        self.grid = np.linspace(0., 1., 100)
        self.potential = Potential(type='tophat', start=self.grid[0], end=self.grid[-1])
        self.particle_number = part_num
        self.psi = Wavefunction(self.grid, self.particle_number)
        self.hbar = 1.
    
    def Hamiltonian(self):
        # Kinetic energy contribution
        kinetic_energy = np.zeros_like(self.grid, dtype=np.complex128)
        dx_inverse_sq = (self.grid[1]-self.grid[2])**(-2)
        coefficient = self.hbar**2/(2.*self.psi.mass)
        
        for n in range(self.particle_number):
            # Interior contribution
            for i in range(1, self.grid.shape[0]-1):
                psi_xx = (self.psi.values[n,i-1]-2.*self.psi.values[n,i]+self.psi.values[n,i+1])*dx_inverse_sq
                kinetic_energy[i] = coefficient*psi_xx
            # Boundary contribution
            psi_xx = (self.psi.values[n,0]-2.*self.psi.values[n,1]+self.psi.values[n,2])*dx_inverse_sq
            kinetic_energy[0] =  coefficient*psi_xx
            psi_xx = (self.psi.values[n,-3]-2.*self.psi.values[n,-2]+self.psi.values[n,-1])*dx_inverse_sq
            kinetic_energy[-1] = coefficient*psi_xx

            # Potential Energy contribution
            potential_energy = np.zeros_like(self.grid, dtype=np.complex128)
            for i in range(self.grid.shape[0]):
                potential_energy[i] = self.potential.V(self.grid[i]) # TODO: multiply wavefunction
        
        return kinetic_energy + potential_energy
        
def permutation_parity(permutation):
    '''
    Given a permutation of the digits 0..N in order as a list, 
    returns its parity (or sign): +1 for even parity; -1 for odd.
    '''
    permutation = list(permutation)
    parity = 1
    N = len(permutation)
    for i in range(0,N-1):
        if permutation[i] != i:
            parity *= -1
            mn = min(range(i,N), key=permutation.__getitem__)
            permutation[i], permutation[mn] = \
                permutation[mn], permutation[i]
    return parity

def antisymmetrize(sys):
    N = sys.particle_number
    # Calculate all the *unique* combinations of position coordinates
    # This is equivalent to getting all possible argument combinations
    indices = np.mgrid[0:sys.grid.shape[0]+.1:1].reshape(len(sys.grid.shape),-1).T
    combs = combinations(indices, N)

    # For each combination, sum over the permutations of arguments
    perms = permutations(range(N))
    parities = [permutation_parity(perm) for perm in perms]
    
    for combination in combs:
        temp = 0.
        for j in range(len(parities)):
            # Sum over parity times permuted wavefunction
            temp += parities[j]*sys.psi.values[combination[perms[j]]]
            
        # Now set the values to temp times the parity
    
if __name__ == '__main__':
    sim = System(2)
    print(antisymmetrize(sim))
