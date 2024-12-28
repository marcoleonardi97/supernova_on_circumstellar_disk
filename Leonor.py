import numpy
from amuse.units import (units, constants)
from amuse.lab import Particles
system= Particles(1)
star=system[0]
star.mass=1|units.Msun
star.position = (0, 0, 0) | units.au
#star.velocity = (0, 0, 0) | units.kms
from amuse.lab import Particles, units, nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.units.constants import G

def make_disk_around_star(star, Ndisk=1000, Rin=0.1|units.au, Rout=100|units.au, Mdisk_fraction=0.01):
    """
    Create a protoplanetary disk around a star.
    
    Parameters:
        star: Particle
            The central star particle.
        Ndisk: int
            Number of particles in the disk.
        Rin: Quantity
            Inner radius of the disk.
        Rout: Quantity
            Outer radius of the disk.
        Mdisk_fraction: float
            Disk mass as a fraction of the star's mass.
    
    Returns:
        disk: Particles
            The generated protoplanetary disk particles.
        converter: nbody_system.nbody_to_si
            Unit converter for the disk.
    """
    Mstar = star.mass
    converter = nbody_system.nbody_to_si(Mstar, Rout)

    Mdisk = Mdisk_fraction * Mstarc

    disk = ProtoPlanetaryDisk(Ndisk,
                              convert_nbody=converter,
                              Rmin=Rin / Rout,
                              Rmax=Rout / Rout,
                              q_out=10.0,
                              discfraction=Mdisk / Mstar).result
    disk.name = "disk"
    disk.move_to_center()

    # Assign masses to disk particles
    disk.mass = Mdisk / float(Ndisk)

    # Estimate particle radius based on density
    rho = 3.0 | (units.g / units.cm**3)
    disk.radius = (disk.mass / (4 * rho))**(1./3.)
    return disk, converter

# Generate a disk around the star
disk, converter = make_disk_around_star(star)

# Combine the star and disk into a single system
system.add_particles(disk)

print(f"Star mass: {star.mass.in_(units.MSun)}")
print(f"Disk has {len(disk)} particles with total mass {disk.mass.sum().in_(units.MSun)}")

from amuse.plot import scatter
def plot(star, disk):
    scatter(disk.x, disk.y, c='r')
    scatter(star.x, star.y,color='blue')
plot(star, disk)
from amuse.community.huayno.interface import Huayno
from amuse.community.fi.interface import Fi
gravity = Huayno(converter)
gravity.particles.add_particles(system)
channel = {"from star+disk": system.new_channel_to(gravity.particles),
            "to stars+disk": gravity.particles.new_channel_to(system)}
hydro = Fi(converter, mode="openmp")
hydro.parameters.use_hydro_flag = True
hydro.parameters.radiation_flag = False
hydro.parameters.gamma = 1
hydro.parameters.isothermal_flag = True
hydro.parameters.integrate_entropy_flag = False
hydro.parameters.timestep = 0.01|units.Myr
hydro.parameters.verbosity = 0
hydro.parameters.eps_is_h_flag = False  # h_smooth is constant
eps = 0.1 | units.au
hydro.parameters.gas_epsilon = eps
hydro.parameters.sph_h_const = eps

# Add disk (gas particles) to the hydrodynamic solver
hydro.particles.add_particles(disk)
channel.update({"from_disk": disk.new_channel_to(hydro.particles)})
channel.update({"to_disk": hydro.particles.new_channel_to(disk)})

# No need to add system to hydro.dm_particles, as Fi doesn't support dark matter
# Instead, ensure the system and disk interact via gravity:
system.add_particles(disk)

from amuse.couple import bridge
from amuse.ext.composition_methods import *
gravhydro = bridge.Bridge(use_threading=False)  # , method=SPLIT_4TH_S_M4)
gravhydro.add_system(gravity, (hydro,))
gravhydro.add_system(hydro, (gravity,))
gravhydro.timestep = 0.2

from amuse.ext.composition_methods import *
def gravity_hydro_bridge(gravity, hydro, gravhydro, system
                         t_end):

    gravity_initial_total_energy = gravity.get_total_energy() + hydro.get_total_energy()
    model_time = 0 | units.Myr
    dt = 0.5 | units.Myr  # 1.0*Pinner
    while model_time < t_end:

        model_time += dt

        dE_gravity = gravity_initial_total_energy / (
            gravity.get_total_energy() + hydro.get_total_energy()
        )
        print("Time:", model_time.in_(units.day), \
              "dE=", dE_gravity)  # , dE_hydro

        gravhydro.evolve_model(model_time)
        channel["to star+disk"].copy()
        print("g=", gravity.particles)
        print(gravity.particles.y.in_(units.au))

    gravity.stop()
    hydro.stop()

t_end = 10 | units.Myr
gravity_hydro_bridge(gravity, hydro, gravhydro,system,t_end)

