import numpy as np
from amuse.units import units, constants
from amuse.plot import plot, scatter
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.lab import nbody_system
from matplotlib import pyplot as plt
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.community.ph4.interface import Ph4
from amuse.community.fi.interface import Fi
from amuse.couple import bridge
from amuse.ext.composition_methods import *



### Setting up the binary system and the disk ###

def orbital_period(Mtot, a):
    return (((4 * numpy.pi**2) * a**3)/(constants.G * Mtot)).sqrt()


# We can play around with these
M1 = 1.0 | units.MSun
M2 = 0.5 | units.MSun 
semi_a = 15 | units.au
ecc = 0.6
bodies = new_binary_from_orbital_elements(M1, M2, semi_a, ecc,
                                         G=constants.G)
bodies[0].name = "primary"
bodies[1].name = "secondary"
primary = bodies[bodies.name=="primary"]
secondary = bodies[bodies.name=="secondary"]
RH = semi_a * (1.0 - ecc) * (M2 / (3 * M1))**(1./3.)
R = 1 | units.au

converter = nbody_system.nbody_to_si(M1, R)
Ndisk = 1000
Rin = 1 | units.AU
Rout = 5 | units.AU
Pinner = orbital_period(M1, Rin)
Mdisk = 0.01 * M1

disk = ProtoPlanetaryDisk(Ndisk,
                          convert_nbody=converter,
                          Rmin=Rin/R,
                          Rmax=Rout/R,
                          q_out=10.0,
                          discfraction=Mdisk/M1).result
disk.name = "disk "
disk.move_to_center()
disk.position += primary.position
disk.velocity += primary.velocity

masses = Mdisk / float(Ndisk)
disk.mass = masses
rho = 3.0 | (units.g / units.cm**3)
disk.radius = (disk.mass / (4 * rho))**(1./3.)

bodies.add_particles(disk)
bodies.move_to_center()
star.move_to_center()


def plot_system(save=False, save_name=None, time=None):
    """
    Simple 2D plotting of the bodies object to visualize it over time.
    I'm currently using this function to save each plot as a frame and then stitch them together as a .gif
    This works but I'm guessing is way more time consuming than having an update func() in pyplot animation. - Marco
    """
    l = 10
    plt.xlim(-l,l)
    plt.ylim(-l,l)
    scatter(bodies.x.in_(units.AU), bodies.y.in_(units.AU), c='orange', alpha=0.5, s=10)
    scatter(primary.x, primary.y, marker="*",c='r', s=120,label="Primary Star")
    scatter(secondary.x, secondary.y, marker='*', c='y',s=120, label="Secondary Star")
    plt.legend()
    if time is not None:
        plt.title(f"Disk at {time:.2f} yr")

    if save:
        plt.savefig(save_name)

### Setting up the simulation workers (Ph4 and Fi) ###

gravity = Ph4(converter)
gravity.particles.add_particles(bodies)
channel = {"from stars": bodies.new_channel_to(gravity.particles),
            "to_stars": gravity.particles.new_channel_to(bodies)}

hydro = Fi(converter, mode="openmp")
hydro.parameters.use_hydro_flag = True
hydro.parameters.radiation_flag = False
hydro.parameters.gamma = 1
hydro.parameters.isothermal_flag = True
hydro.parameters.integrate_entropy_flag = False
hydro.parameters.timestep = 0.01*Pinner 
hydro.parameters.verbosity = 0
hydro.parameters.eps_is_h_flag = True  # h_smooth is NOT constant
eps = 0.1 | units.au
hydro.parameters.gas_epsilon = eps
hydro.parameters.sph_h_const = eps * 5  # This works with 1000 particles, but it can run into problems in long simulation. Keep an eye on this if we get '-1' errors in Fi. -Marco
hydro.particles.add_particles(disk)
hydro.dm_particles.add_particles(secondary.as_set())
channel.update({"from_disk": disk.new_channel_to(hydro.particles)})
channel.update({"to_disk": hydro.particles.new_channel_to(disk)})
channel.update({"from_moon": secondary.new_channel_to(hydro.dm_particles)})
channel.update({"to_moon": hydro.dm_particles.new_channel_to(secondary)})
channel.update({"pts": star.new_channel_to(planet)})
channel.update({"stp": planet.new_channel_to(star)})

### Setting up Bridge ###

gravhydro = bridge.Bridge(use_threading=False)  # , method=SPLIT_4TH_S_M4)
gravhydro.add_system(gravity, (hydro,))
gravhydro.add_system(hydro, (gravity,))
gravhydro.timestep = 0.001 * Pinner # This is important and will cause errors in long simulation if not set right.. with 0.001 and dt of 7*Pinner the simulation crashed after 9 years -Marco


# Actual simulation code

#Notes:
# written like this, the simulation runs fairly quickly for 1 years, but it will take around 1h ~ 10 years. 
# i'm assuming it's because i'm plotting the positions at each time step, maybe i'll just save them in a list and try to make it faster
# especially if we increase the number of particles in the disk. 

dtt = Pinner.number | units.s
t_end = 100.0 | units.yr
model_time = 0 | units.yr
dt = 1. * dtt  # Usually would be lower, but i'm increasing it to run the simulation for 100 years. Once i'm sure it works properly i'll try to do ~Myr
plot_counter = 0

while model_time < t_end:

    model_time += dt
    plot_counter += 1
    
    print("Time:", model_time.in_(units.yr))  

    gravhydro.evolve_model(model_time)
    channel["to_stars"].copy()
    channel["to_disk"].copy()
    channel["to_moon"].copy()
    channel["pts"].copy()
    channel["stp"].copy()

    fig = plt.figure()
    plot_system(save=True, save_name=f"time_{model_time.in_(units.day)}.png", time=model_time.number)
    

gravity.stop()
hydro.stop()


# Animation 

import glob
import os
import imageio.v2 as imageio
from matplotlib.animation import FuncAnimation, PillowWriter

# -------------------- This code is only needed if we don't name the plots with plot_counter
files = glob.glob("*.png")

numbers = []
for f in files:
    numbers.append(float(f[5:12]))
ff = [y for _,y in sorted(zip(numbers, files))]
# ---------------------




with imageio.get_writer("perturbed_disk_long.gif", mode='I', duration=0.1) as writer:
    for frame in ff:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory 
for frame in ff:
    os.remove(frame)
