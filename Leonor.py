import numpy
from amuse.units import (units, constants)
from amuse.lab import Particles

from amuse.lab import Particles, units, nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.units.constants import G
def orbital_period(Mtot, a):
    return (((4 * numpy.pi**2) * a**3)/(constants.G * Mtot)).sqrt()

from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary

def get_binary():
    Mstar = 1.0 | units.MSun
    Mstar2 = 0.5| units.Msun
    a_binary = 5.2 | units.au
    e_binary = 0.6
    bodies = new_binary_from_orbital_elements(Mstar, Mstar2, a_binary, e_binary,
                                             G=constants.G)
    bodies[0].name = "primary_star"
    bodies[1].name = "secondary_star"
    primary = bodies[bodies.name=="primary_star"]
    secondary = bodies[bodies.name=="secondary_star"]
    RH_planet = a_binary * (1.0 - e_binary) * (Mstar2 / (3 * Mstar))**(1./3.)
   
    return bodies
bodies = get_binary()
print(bodies)

from amuse.lab import Particles, units, nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.units.constants import G

def make_disk_around_star(primary, Ndisk=1000, Rin=0.1|units.au, Rout=0.7|units.au, Mdisk_fraction=0.01):
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
    Mstar = primary.mass
    converter = nbody_system.nbody_to_si(Mstar, Rout)

    Mdisk = Mdisk_fraction * Mstar
    Ndisk = 1000
    Rin = 0.1|units.au
    Rout = 0.7|units.au
    Pinner = orbital_period(Mstar, Rin)
    Mdisk = 0.01 * Mstar


    disk = ProtoPlanetaryDisk(Ndisk,
                              convert_nbody=converter,
                              Rmin=Rin / Rout,
                              Rmax=Rout / Rout,
                              q_out=10.0,
                              discfraction=Mdisk / Mstar).result
    disk.name = "disk"
    disk.move_to_center()
    disk.position += primary.position
    disk.velocity += primary.velocity

    # Assign masses to disk particles
    disk.mass = Mdisk / float(Ndisk)

    # Estimate particle radius based on density
    rho = 3.0 | (units.g / units.cm**3)
    disk.radius = (disk.mass / (4 * rho))**(1./3.)
    return disk, converter,Pinner
    bodies.add_particles(disk)
    bodies.move_to_center()
    primary.move_to_center()

# Generate a disk around the star
primary = bodies[bodies.name == "primary_star"][0]
secondary = bodies[bodies.name == "secondary_star"]
disk, converter,Pinner= make_disk_around_star(primary)

# Combine the star and disk into a single system

print(f"Star mass: {primary.mass.in_(units.MSun)}")
print(f"Disk has {len(disk)} particles with total mass {disk.mass.sum().in_(units.MSun)}")

from amuse.plot import scatter
def plot(primary, disk):
    scatter(disk.x.in_(units.au), disk.y.in_(units.au), c='r')
    scatter(primary.x.in_(units.au), primary.y.in_(units.au),color='blue')
    #scatter(secondary.x.in_(units.au), secondary.y.in_(units.au),color='blue')
plot(primary, disk)
from amuse.community.huayno.interface import Huayno
from amuse.community.fi.interface import Fi
gravity = Huayno(converter)
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
hydro.parameters.eps_is_h_flag = False  # h_smooth is constant
eps = 0.1 | units.au
hydro.parameters.gas_epsilon = eps
hydro.parameters.sph_h_const = eps*5

hydro.particles.add_particles(disk)
hydro.dm_particles.add_particles(secondary.as_set())
channel.update({"from_disk": disk.new_channel_to(hydro.particles)})
channel.update({"to_disk": hydro.particles.new_channel_to(disk)})
channel.update({"from_star2": secondary.new_channel_to(hydro.dm_particles)})
channel.update({"to_star2": hydro.dm_particles.new_channel_to(secondary)})

bodies.add_particles(disk)

# No need to add system to hydro.dm_particles, as Fi doesn't support dark matter
# Instead, ensure the system and disk interact via gravity:


from amuse.couple import bridge
from amuse.ext.composition_methods import *
gravhydro = bridge.Bridge(use_threading=False)  # , method=SPLIT_4TH_S_M4)
gravhydro.add_system(gravity, (hydro,))
gravhydro.add_system(hydro, (gravity,))
gravhydro.timestep = 0.02*Pinner

from amuse.ext.composition_methods import *
def gravity_hydro_bridge(gravity, hydro, gravhydro, bodies,t_end):

    gravity_initial_total_energy = gravity.get_total_energy() + hydro.get_total_energy()
    model_time = 0 | units.yr
    dtt = Pinner.number | units.s
    dt = 1.0*Pinner
    time=[]
    while model_time < t_end:

        model_time += dt

        dE_gravity = gravity_initial_total_energy / (
            gravity.get_total_energy() + hydro.get_total_energy()
        )
        print("Time:", model_time.in_(units.yr), \
              "dE=", dE_gravity)  # , dE_hydro

        gravhydro.evolve_model(model_time)
        channel["to_stars"].copy()
        channel["to_disk"].copy()
        channel["to_star2"].copy()
        print("g=", gravity.particles)
        print(gravity.particles.y.in_(units.au))
        

    gravity.stop()
    hydro.stop()

t_end = 50| units.yr
gravity_hydro_bridge(gravity, hydro, gravhydro,bodies,t_end)

