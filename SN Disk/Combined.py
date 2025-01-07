from amuse.lab import *
from amuse.couple import bridge
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from Disk import star_1, star_2, planets, disk
from SN import supernova, shell
def disk_pos(disk, inclination_angle):
    rotation_matrix = np.array([
        [1,0,0],
        [0,np.cos(inclination_angle),-np.sin(inclination_angle)],
        [0,np.sin(inclination_angle),np.cos(inclination_angle)]])
    for p in disk:
        p.position = rotation_matrix @ p.position.value_in(units.AU) | units.AU
        p.velocity = rotation_matrix @ p.velocity.value_in(units.kms) | units.kms


def supernova_pos(supernova, distance, angle):
    supernova.position = (distance*np.cos(angle), distance*np.sin(angle), 0) | units.AU

disk_pos(disk, np.radians(30))
supernova_pos(supernova, 50 | units.AU, np.radians(45))

system = star_1+star_2+planets+disk+shell

gravity = Hermite()
gravity.particles.add_particles(system)

disk_hydro = Fi()
disk_hydro.particles.add_particles(disk)

system_bridge = bridge.Bridge(use_threading=False)
system_bridge.add_system(gravity, (disk_hydro,))
system_bridge.add_system(disk_hydro, (gravity,))

def toomre_Q(disk):
    G = constants.G
    Q_values = []
    for particle in disk:
        c = np.sqrt(particle.u) 
        kappa = np.sqrt(G*central_star.mass/particle.position.length()**3)  
        sigma = particle.mass/(np.pi*(particle.position.length()**2))  
        Q = c*kappa/(np.pi*G*sigma)
        Q_values.append(Q)
    return np.mean(Q_values)

def disk_params(disk):
    surface_density = np.mean([p.mass/(np.pi*(p.position.length()**2))for p in disk])
    surface_temperature = np.mean([p.u for p in disk])
    t_Q = toomre_Q(disk)
    print(f"Surface Density: {surface_density.in_(units.g / units.cm**2)}")
    print(f"Surface Temperature: {surface_temperature.in_(units.K)}")
    print(f"Toomre Q: {t_Q}")

t_end = 10 | units.yr  
dt = 0.1 | units.yr

while system_bridge.model_time < t_end:
    system_bridge.evolve_model(system_bridge.model_time+dt)
    disk_params(disk)

gravity.stop()
disk_hydro.stop()
