from amuse.lab import *
import numpy as np

M_SN = 30 | units.MSun  
supernova = Particles(1)
supernova.mass = M_SN
supernova.position = (50,0,0) | units.AU
supernova.velocity = (0,0,0) | units.kms
supernova.radius = 10 | units.RSun  

E_SN = 1*10**51 | units.erg 
v_ejecta = np.sqrt(2*E_SN/M_SN).in_(units.kms)  

N_shell = 1000
shell = Particles(N_shell)
shell_mass = M_SN / N_shell

for i, particle in enumerate(shell):
    theta = np.random.uniform(0, np.pi)
    phi = np.random.uniform(0, 2*np.pi)
    r = np.random.uniform(1, 2) | units.RSun
    particle.mass = shell_mass
    particle.position = supernova.position + r*(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi),np.cos(theta))
    particle.velocity = v_ejecta * (particle.position-supernova.position).normalized()

gravity = Hermite()
gravity.particles.add_particles(shell)

t_end = 10 | units.yr 
dt = 0.1 | units.yr

while gravity.model_time < t_end:
    gravity.evolve_model(gravity.model_time + dt)

gravity.stop()

