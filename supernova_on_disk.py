# -*- coding: ascii -*-
"""
Creates a protoplanetary disk around a sun-like star and nearby supernova to simulate the collision
(project part 2)
"""



import os
import numpy as np
import glob

from amuse.units import units
from amuse.plot import scatter
from amuse.units import nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particles

from amuse.ext.star_to_sph import (pickle_stellar_model, convert_stellar_model_to_SPH)
from amuse.test.amusetest import get_path_to_results
from amuse.community.evtwin.interface import EVtwin
from amuse.community.fi.interface import Fi

from matplotlib import pyplot as plt
import imageio.v2 as imageio
from matplotlib.animation import FuncAnimation, PillowWriter


# Making the disk

N = 1000
tend = 2. | units.yr
Mstar = 1. | units.MSun

convert = nbody_system.nbody_to_si(Mstar, 1. | units.AU)
proto = ProtoPlanetaryDisk(
    N, convert_nbody=convert, densitypower=1.5, Rmin=4, Rmax=20, q_out=1.)
gas = proto.result
gas.h_smooth = 0.06 | units.AU

sun = Particles(1)
sun.mass = Mstar
sun.radius = 2. | units.AU
sun.x = 30. | units.AU
sun.y = 0. | units.AU
sun.z = 0. | units.AU
sun.vx = 0. | units.kms
sun.vy = 0. | units.kms
sun.vz = 0. | units.kms

gas.position += sun.position
gas.velocity += sun.velocity


# Hydro code for the disk
sph = Fi(convert, mode='openmp')

sph.parameters.use_hydro_flag = True
sph.parameters.radiation_flag = False
sph.parameters.self_gravity_flag = True
sph.parameters.gamma = 1.
sph.parameters.isothermal_flag = True
sph.parameters.integrate_entropy_flag = False
sph.parameters.timestep = 0.125 | units.yr

sph.particles.add_particles(gas)
sph.dm_particles.add_particles(sun)


# Making the supernova


def setup_stellar_evolution_model(time):
    out_pickle_file = os.path.join(get_path_to_results(), 
                                   "test.pkl")

    stellar_evolution = EVtwin(redirection="none")  # Using EVtwin instead of MESA
    stars = Particles(1)
    stars.mass = 10 | units.MSun  # Massive star setup for supernova
    stellar_evolution.particles.add_particles(stars)
    stellar_evolution.commit_particles()
    tt = 1 | units.Myr

    print(
        "Evolving a EVtwin star with mass:",
        stellar_evolution.particles[0].mass
    )

    while stellar_evolution.model_time < time:  
        stellar_evolution.evolve_model()
        if stellar_evolution.model_time >= tt:
            tt += 1 | units.yr
            print("Star:", stellar_evolution.particles[0].stellar_type, stellar_evolution.model_time.in_(units.Myr))
            #print(stellar_evolution.particles[0].luminosity)
    
    print("Evolved star to", stellar_evolution.particles[0].age)
    print("Radius:", stellar_evolution.particles[0].radius)
    
    # Save the stellar model structure for later use
    pickle_stellar_model(stellar_evolution.particles[0], out_pickle_file)
    stellar_evolution.stop()
    print("Done")
    return out_pickle_file

# Use the stellar evolution setup with EVtwin
t = 1 | units.yr

path_to_pickle = "tmpow7eucrj/test.pkl" # Change to your path. If you can't run EVTwin, download the test.pkl pickle file from this github.
if os.path.exists(path_to_pickle):
    pickle_file = path_to_pickle
    print(f"Using existing file {path_to_pickle}")
else:
    print("Setting up new star...")
    pickle_file = setup_stellar_evolution_model(t) # Assuming the star model hasnt already been created


# Modelling core and gas particles

number_of_sph_particles = 1000
print(pickle_file)
print("Creating initial conditions from a EVTwin stellar evolution model...")
model = convert_stellar_model_to_SPH(
        None,
        number_of_sph_particles,
        seed=12345,
        pickle_file=pickle_file,
        #        base_grid_options = dict(type = "glass", target_rms = 0.01),
        with_core_particle=True,
        target_core_mass = 1.4|units.MSun
    )
print("model=", model)
core, gas_without_core, core_radius = \
        model.core_particle, model.gas_particles, model.core_radius
print("Created", len(gas_without_core),
       "SPH particles and one 'core-particle':\n", core)
print("Setting gravitational smoothing to:", core_radius.in_(units.km))

def inject_supernova_energy(gas_particles, 
                            explosion_energy=1.0e+51|units.erg,
                            exploding_region=10|units.RSun):
    inner = gas_particles.select(
        lambda pos: pos.length_squared() < exploding_region**2,
        ["position"])
    print(len(inner), "innermost particles selected.")
    print("Adding", explosion_energy / inner.total_mass(), "of supernova " \
        "(specific internal) energy to each of the n=", len(inner), "SPH particles.")
    inner.u += explosion_energy / inner.total_mass()
    
inject_supernova_energy(gas_without_core, exploding_region=1|units.RSun)

converter = nbody_system.nbody_to_si(10|units.MSun, core_radius)

# Hydro code for the supernova
hydro_code = Fi(converter,mode='openmp',redirection='none')
hydro_code.parameters.epsilon_squared = core_radius**2
hydro_code.parameters.n_smooth_tol = 0.01
hydro_code.gas_particles.add_particles(gas_without_core)
hydro_code.dm_particles.add_particle(core)



# Simulation

tend = 10 | units.yr
timestep = 1 | units.day
time = 0 | units.yr

while time < tend:

    #hydro_code.evolve_model(sph.model_time + timestep)
    sph.evolve_model(time)
    time += sph.parameters.timestep
    print(f"Evolved to: {sph.model_time}")
    plt.figure()
    l = 55
    plt.xlim(-l,l)
    plt.ylim(-l,l)
    scatter(hydro_code.gas_particles.x.in_(units.au), hydro_code.gas_particles.y.in_(units.au), c='orange', label="SN")
    scatter(hydro_code.dm_particles.x.in_(units.au), hydro_code.dm_particles.y.in_(units.au), marker='*',c='yellow', label="Core")
    scatter(sph.particles.x.in_(units.au), sph.particles.y.in_(units.au), c='blue', label="Disk")
    scatter(sun.x.in_(units.au), sun.y.in_(units.au), marker='*', c='yellow')
    plt.legend(loc='upper right')
    plt.title(f"{time}")
    plt.savefig(f"nova_{time.number:.3f}.png")

sph.stop()



files = glob.glob("*.png")

numbers = []
for f in files:
    numbers.append(float(f[5:10]))
ff = [y for _,y in sorted(zip(numbers, files))]

with imageio.get_writer("nova_attempt.gif", mode='I', duration=0.1) as writer:
    for frame in ff:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory 
for frame in ff:
    os.remove(frame)



