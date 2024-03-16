import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def quantum_fusion_reactor(n_particles, initial_positions, initial_momenta, hbar=1):
    """
    Simulates a Quantum Fusion Reactor (QFR) using the Schrödinger equation.
    
    Parameters:
    n_particles (int): Number of particles in the system.
    initial_positions (list or array): Initial positions of the particles.
    initial_momenta (list or array): Initial momenta of the particles.
    hbar (float): Reduced Planck constant (default 1).
    
    Returns:
    tuple: (final_positions, final_momenta, energy_history)
    """
    # Initialize positions and momenta
    positions = np.array(initial_positions)
    momenta = np.array(initial_momenta)
    
    # Initialize energy history
    energy_history = []
    
    # Time evolution parameters
    dt = 0.01
    n_steps = 1000
    
    # Time evolution loop
    for step in range(n_steps):
        # Calculate Hamiltonian matrix
        hamiltonian_matrix = np.zeros((n_particles, n_particles))
        for i in range(n_particles):
            for j in range(n_particles):
                if i != j:
                    r_ij = np.linalg.norm(positions[i] - positions[j])
                    hamiltonian_matrix[i, j] = 1 / r_ij
        
        # Calculate Hamiltonian and energy
        hamiltonian = np.eye(n_particles) - hbar**2 * np.dot(momenta, momenta.T) / 2 + hbar**2 * np.dot(positions, positions.T)
        energy = np.sum(np.diagonal(hamiltonian
