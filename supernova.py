import os
import numpy
from matplotlib import pyplot as plt
from amuse.units import units
from amuse.ext.star_to_sph import (pickle_stellar_model, convert_stellar_model_to_SPH)
from amuse.test.amusetest import get_path_to_results
from amuse.community.evtwin.interface import EVtwin
from amuse.datamodel import Particles
from amuse.plot import scatter
import imageio.v2 as imageio
from matplotlib.animation import FuncAnimation, PillowWriter

# Setting up the star

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

path_to_pickle = "tmpnvo376xy/test.pkl" # Change to your path 
if os.path.exsits(path_to_pickle):
  pickle_file = path_to_pickle
else:
  pickle_file = setup_stellar_evolution_model(t) # Assuming the star model hasnt already been created


# Modelling core and gas particles

number_of_sph_particles = 1000
print(pickle_file)
print("Creating initial conditions from a MESA stellar evolution model...")
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

from amuse.community.fi.interface import Fi
from amuse.lab import nbody_system
converter = nbody_system.nbody_to_si(10|units.MSun, core_radius)

hydro_code = Fi(converter,mode='openmp',redirection='none')
hydro_code.parameters.epsilon_squared = core_radius**2
hydro_code.parameters.n_smooth_tol = 0.01
hydro_code.gas_particles.add_particles(gas_without_core)
hydro_code.dm_particles.add_particle(core)



def plot_star_from_code(code, save = False, name=None, time=None):
    l = 100 
    plt.xlim(-l,l)
    plt.ylim(-l,l)
    gas = code.gas_particles
    core = code.dm_particles
    scatter(gas.x.value_in(units.au), gas.y.value_in(units.au), c='orange', label="Gas")
    scatter(core.x.value_in(units.au), core.y.value_in(units.au), marker='*',c='r', label="Core")
    if time is not None:
        plt.title(f"{time.number:.2f} days")
    if save:
        plt.savefig(f"sn_{name}.png")

tend = 100 | units.day
timestep = 1 | units.hour
frames = 0
while hydro_code.model_time < tend:
    frames += 1
    hydro_code.evolve_model(hydro_code.model_time + (timestep))
    print(f"Evolved to: {hydro_code.model_time.in_(units.day)}, Frame: {frames}")
    plot_star_from_code(hydro_code, save=True, name=f"50_{frames}", time=hydro_code.model_time.in_(units.day))

hydro_code.stop()


frames = [f"sn_50_{j}.png" for j in range(1,frames)]

with imageio.get_writer("sn_test.gif", mode='I', duration=0.1) as writer:
    for frame in frames:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory 
for frame in frames:
    os.remove(frame)
