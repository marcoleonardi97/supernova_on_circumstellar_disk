import numpy as np
from matplotlib import pyplot as plt

from amuse.community.fi.interface import Fi
from amuse.units import units, constants
from amuse.units import nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particles

from parameters_to_check import *
from plotting import *
from supernova import *



N = 20000
time = 0 |units.yr
tend = 50. | units.yr
Mstar = 15. | units.MSun
Q_history = [] # Q Parameter at each dt
Z_history = [] # Metallicities at each dt

convert = nbody_system.nbody_to_si(Mstar, 1. | units.AU)
proto = ProtoPlanetaryDisk(
    N, convert_nbody=convert, densitypower=1.5, Rmin=5, Rmax=20, q_out=1.)
gas = proto.result
gas.h_smooth = 0.06 | units.AU

sun = Particles(1)
sun.mass = Mstar
sun.radius = 2. | units.AU
sun.x = 0. | units.AU
sun.y = 0. | units.AU
sun.z = 0. | units.AU
sun.vx = 0. | units.kms
sun.vy = 0. | units.kms
sun.vz = 0. | units.kms

# External supernova setup
supernova_position = [0.5, 0.0, 0.0] | units.parsec  # Position relative to the disk
shock_speed = 1e4 | units.kms  # Typical supernova shock speed
explosion_energy = 1.0e+51 | units.erg
ejecta_mass = 10 | units.MSun
heavy_element_mass = 0.1 | units.MSun

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

supernova_time = 10. |units.yr
supernova_triggered = False

while time < tend:
    print(f"Triggering supernova at {time}")
    heavy_elements_mass = 0.01 | units.MSun # for example

    # This setup is for the host star going SN
        #ejecta_mass = 10 | units.MSun  # Total ejecta mass 
        #inject_supernova_energy(gas, heavy_elements_mass,exploding_region=10 | units.AU)
        #sun.mass -= ejecta_mass  # Mass lost from the star

    # This setup is for an external star going SN, which is what we want.
    # We can trigger this immediately, as the shockwave still needs time to reach the disk 
    # ps. not sure if this function works, probably not yet
    inject_external_supernova(
        gas_particles=gas,
        supernova_position=supernova_position,
        time=time,
        explosion_energy=explosion_energy,
        ejecta_mass=ejecta_mass,
        heavy_element_mass=heavy_element_mass,
        shock_speed=shock_speed
    )

    
    Q_values = compute_toomre_q(sph, gas, sun.mass)  # probably not very accurate
    Q_history.append(np.mean(Q_values))
    radii, metallicities = compute_metallicity_profile(gas) # this works but maybe it's better to just focus on 1 or 2 heavy elements, like 60F or some typical SN injected elements
    Z_history.append(metallicities)
    
    sph.evolve_model(time)
    print(f"Disk evolved to {time}")
    time += 1. | units.yr

    L = 50
    rhod, rho = make_density_map(sph, N=200, L=L)
    plt.figure(figsize=(8, 8))
    plt.imshow(np.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
                  extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15)
    plt.title(f'Time: {time}')
    plt.xlabel('AU')
    plt.savefig(f'{time}.png')

frames = [f'{time} yr' for time in range(int(round(tend.number)))]
animate_frames(frames, save_as='protodisk_sn')
animate_2d_plot(range(len(50)), Q_history, save_as='disk_stability')
