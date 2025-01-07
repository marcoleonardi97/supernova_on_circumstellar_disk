"""
Example of bridging two different hydro codes (Fi) for two different protoplanetary disks.
"""
import numpy as np

from matplotlib import pyplot as plt

from amuse.community.fi.interface import Fi

from amuse.units import units
from amuse.plot import scatter
from amuse.units import nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particles
import glob
import os
import imageio.v2 as imageio
from matplotlib.animation import FuncAnimation, PillowWriter


# Make disk 1 ------

N = 1000
Mstar = 1. | units.MSun

convert = nbody_system.nbody_to_si(Mstar, 1. | units.AU)
proto = ProtoPlanetaryDisk(
    N, convert_nbody=convert, densitypower=1.5, Rmin=1, Rmax=5, q_out=1.)
gas = proto.result
gas.h_smooth = 0.06 | units.AU

sun = Particles(1)
sun.mass = Mstar
sun.radius = 2. | units.AU
sun.x = -7. | units.AU
sun.y = 0. | units.AU
sun.z = 0. | units.AU
sun.vx = 0. | units.kms
sun.vy = 0. | units.kms
sun.vz = 0. | units.kms

gas.position += sun.position
gas.velocity += sun.velocity

sph = Fi(convert)

sph.parameters.use_hydro_flag = True
sph.parameters.radiation_flag = False
sph.parameters.self_gravity_flag = True
sph.parameters.gamma = 1.
sph.parameters.isothermal_flag = True
sph.parameters.integrate_entropy_flag = False
sph.parameters.timestep = 0.125 | units.yr

sph.gas_particles.add_particles(gas)
sph.particles.add_particles(sun)

# Make disk 2 ------

N2 = 1000
Mstar2 = 1. | units.MSun

convert = nbody_system.nbody_to_si(Mstar, 1. | units.AU)
proto2 = ProtoPlanetaryDisk(
    N2, convert_nbody=convert, densitypower=1.5, Rmin=1.5, Rmax=4, q_out=1.)
gas2 = proto.result
gas2.h_smooth = 0.06 | units.AU

sun2 = Particles(1)
sun2.mass = Mstar
sun2.radius = 2. | units.AU
sun2.x = 7. | units.AU
sun2.y = 0. | units.AU
sun2.z = 0. | units.AU
sun2.vx = 0. | units.kms
sun2.vy = 0. | units.kms
sun2.vz = 0. | units.kms

gas2.position += sun2.position
gas2.velocity += sun2.velocity

sph2 = Fi(convert)

sph2.parameters.use_hydro_flag = True
sph2.parameters.radiation_flag = False
sph2.parameters.self_gravity_flag = True
sph2.parameters.gamma = 1.
sph2.parameters.isothermal_flag = True
sph2.parameters.integrate_entropy_flag = False
sph2.parameters.timestep = 0.125 | units.yr

sph2.gas_particles.add_particles(gas2)
sph2.particles.add_particles(sun2)


# bridge
from amuse.couple import bridge

all_particles = Particles()
all_particles.add_particles(gas)
all_particles.add_particles(gas2)

hydro = bridge.Bridge()
hydro.add_system(sph, (sph2,))
hydro.add_system(sph2, (sph,))


t = 0 | units.yr
tend = 20. | units.yr
while t < tend:

    hydro.evolve_model(t)
    print(f"t: {t.in_(tend.unit)}")
    t += sph2.parameters.timestep
    
    plt.figure()
    l = 25
    plt.xlim(-l,l)
    plt.ylim(-l,l)
    scatter(sph.gas_particles.x.in_(units.au), sph.gas_particles.y.in_(units.au), c='blue',alpha=0.4,label="Disk 1")
    scatter(sun.x.in_(units.au), sun.y.in_(units.au), marker='*', c='yellow')

    scatter(sph2.gas_particles.x.in_(units.au), sph2.gas_particles.y.in_(units.au), c='orange',alpha=0.4,label="Disk 2")
    scatter(sun2.x.in_(units.au), sun2.y.in_(units.au), marker='*', c='yellow')
    plt.xlabel('AU')
    plt.legend(loc='upper right')
    plt.title(f"{t}")
    plt.savefig(f'disk_{t.number:.3f}.png')
    plt.close()
print("Done")





files = glob.glob("*.png")

numbers = []
for f in files:
    numbers.append(float(f[5:10]))
ff = [y for _,y in sorted(zip(numbers, files))]

with imageio.get_writer("disk_test.gif", mode='I', duration=0.1) as writer:
    for frame in ff:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory 
for frame in ff:
    os.remove(frame)
