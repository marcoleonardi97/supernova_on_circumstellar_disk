from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.couple import bridge
import numpy as np

M_1 = 1.0 | units.MSun
star_1 = Particles(1)
star_1.mass = M_1
star_1.position = (0,0,0) | units.AU
star_1.velocity = (0,0,0) | units.kms
star_1.radius = 1 | units.RSun

#The binary companion:
M_2 = 0.3 | units.MSun  
star_2 =Particles(1)
star_2.mass = M_perturber
star_2.position = (5,0,0) | units.AU  # Initial distance from star
star_2.velocity = (0,5,0) | units.kms  # Circular orbit
star_2.radius = 0.5 | units.RSun


#Planets in the disk:
planets = Particles(4)
planet_mass = [0.055,0.815,1.0,0.107] | units.MEarth  # Mercury, Venus, Earth, Mars
planet_position = [(0.39,0,0),(0.72,0,0),(1.0,0,0),(1.52,0,0)] | units.AU
planet_velocity = [(0,47.87,0),(0,35.02,0),(0,29.78,0),(0,24.07,0)] | units.kms

for i, (mass, pos, vel) in enumerate(zip(planet_mass, planet_position, planet_velocity)):
    planets[i].mass = mass
    planets[i].position = pos
    planets[i].velocity = vel
    planets[i].radius = 0.1 | units.REarth

Ndisk=1000 
R_in= 0.1 | units.AU,
R_out= 10 | units.AU
Mdisk=0.01 | units.MSun 

disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=None, densitypower=1.5, Rmin=R_in, Rmax=R_out, q_out=1.0, Mdisk=Mdisk).result
disk.position+=star_1.position

total= star_1+star_2+planets+disk

gravity = Hermite()
gravity.particles.add_particles(total)

disk_hydro = Fi()
disk_hydro.particles.add_particles(disk)
system_bridge = bridge.Bridge(use_threading=False)
system_bridge.add_system(gravity, (disk_hydro,))
system_bridge.add_system(disk_hydro, (gravity,))

t_end = 100 | units.yr
dt = 1 | units.yr

while system_bridge.model_time < t_end:
    system_bridge.evolve_model(system_bridge.model_time + dt)

gravity.stop()
disk_hydro.stop()
