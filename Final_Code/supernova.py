import os
import numpy
from amuse.community.fi.interface import Fi
from amuse.couple import bridge

from amuse.units import units
from amuse.plot import scatter
from amuse.units import nbody_system
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.support import options
from amuse.datamodel import Particles
from matplotlib import pyplot as plt
from amuse.ext.star_to_sph import (pickle_stellar_model, convert_stellar_model_to_SPH)
from amuse.test.amusetest import get_path_to_results
from amuse.community.evtwin.interface import EVtwin

import imageio.v2 as imageiox
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D

os.environ["OMPI_MCA_rmaps_base_oversubscribe"]="true"
options.GlobalOptions.instance().override_value_for_option("polling_interval_in_milliseconds", 10)
class Supernova(object):
    def __init__(self, pickle=None, mass=10|units.MSun, nparticles=1000, external_object=None):
        self.pickle = pickle
        self.mass = mass
        self.nparticles = nparticles
        self.external_object = external_object
        self.particles = Particles()
        self.gas_particles = Particles()

        if pickle is None:
            print("Creating initial conditions from a EVTwin stellar evolution model...")
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
        
            while stellar_evolution.model_time < 1 | units.yr:  
                stellar_evolution.evolve_model()
            
            print("Evolved star to", stellar_evolution.particles[0].age)
            print("Radius:", stellar_evolution.particles[0].radius)
            
            # Save the stellar model structure for later use
            pickle_stellar_model(stellar_evolution.particles[0], out_pickle_file)
            stellar_evolution.stop()
            print("Done")
            self.pickle = out_pickle_file

        model = convert_stellar_model_to_SPH(
        None,
        self.nparticles,
        seed=12345,
        pickle_file=self.pickle,
        #        base_grid_options = dict(type = "glass", target_rms = 0.01),
        with_core_particle=True,
        target_core_mass = 1.4|units.MSun)
        self.core, self.gas_without_core, self.core_radius = \
                model.core_particle, model.gas_particles, model.core_radius
        print("Created", len(self.gas_without_core),
               "SPH particles and one 'core-particle':\n", self.core)
        print("Setting gravitational smoothing to:", self.core_radius.in_(units.km))

        self.core.position = (0,0,0) | units.au
        self.gas_without_core.h_smooth = 2 | units.AU
        self.gas_without_core.name = "nova"
        self.explode()

        self.converter = nbody_system.nbody_to_si(self.mass, self.core_radius)
        self.hydro = Fi(self.converter,mode='openmp',redirection='none')
        self.hydro.parameters.epsilon_squared = self.core_radius**2
        self.hydro.parameters.timestep = 1. | units.hour
        self.hydro.parameters.n_smooth_tol = 0.1
        self.hydro.parameters.isothermal_flag = False
        self.hydro.parameters.integrate_entropy_flag = True
        self.hydro.gas_particles.add_particles(self.gas_without_core)
        self.hydro.dm_particles.add_particle(self.core)

        self.particles.add_particle(self.core)
        self.particles.add_particles(self.gas_without_core)
        self.gas_particles.add_particles(self.gas_without_core)

        if self.external_object is not None:
            self.external_code = Fi(self.converter, mode='openmp', redirection='none')
            self.external_code.parameters.use_hydro_flag = True
            self.external_code.parameters.radiation_flag = False
            self.external_code.parameters.self_gravity_flag = True
            self.external_code.parameters.gamma = 1.
            self.external_code.parameters.isothermal_flag = True
            self.external_code.parameters.integrate_entropy_flag = False
            self.external_code.parameters.timestep = 0.1 | units.yr
            self.external_code.particles.add_particles(self.external_object)

    def explode(self, explosion_energy=1.0e+51|units.erg, exploding_region=1|units.RSun):
        inner = self.gas_without_core.select(lambda pos: pos.length_squared() < exploding_region**2, ["position"])
        print(len(inner), "innermost particles selected.")
        print("Adding", explosion_energy / inner.total_mass(), "of supernova " \
            "(specific internal) energy to each of the n=", len(inner), "SPH particles.")
        inner.u += explosion_energy / inner.total_mass()

    def bridge(self):
        bridged = bridge.Bridge(use_threading=False, method=SPLIT_4TH_S_M4)
        bridged.add_system(self.hydro, (self.external_code,))
        bridged.add_system(self.external_code, (self.hydro,))
        bridged.timestep = self.hydro.parameters.timestep
        return bridged

    def evolve(self, tend, plot=False, plot3d = False, verbose=False):
        time = 0 | units.day
        
        if self.external_object is None:
            code = self.hydro
        else:
            code = self.bridge()
            
        while time < tend:
            
            code.evolve_model(time)
            print(f"Evolution: {time.in_(tend.unit)}")
            time += self.hydro.parameters.timestep

            # i give up Bridge won't transfer energy so i will (still particles don't move accordingly though)
            if self.external_object is not None:
                f = 0.5
                for particle in self.external_object:
                    closest = sn.gas_without_core.find_closest_particle_to(*particle.position)
                    if (closest.position - particle.position).length().in_(units.au) < 1 | units.au:
                        particle.u += f * closest.u
                        particle.vx += f  * closest.vx
                        particle.vy += f * closest.vy
                        particle.vy += f * closest.vz
                        closest.u -= f * particle.u
                        closest.vx -= f * particle.vx
                        closest.vy -= f * particle.vy
                        closest.vz -= f * particle.vz

                
                
                self.external_code.particles.new_channel_to(self.external_object).copy()
                self.hydro.particles.new_channel_to(self.gas_without_core).copy()
                self.gas_without_core.new_channel_to(self.external_object).copy_attributes(["u", "u", "vx", "vy",
                                                                                           "vz", "mass"])
                
                
            if plot:
                self.plot_system(show=False, save=True, time=time)
    
            if plot3d:
                self.plot3d(show=False, save=True, time=time)
    
            if verbose:
                print("Hydro particles: ", np.mean(self.hydro.particles.x.in_(units.au)))
                print("Local self.particles (core): ", np.mean(self.particles.x.in_(units.au)))
                print("Local self.gas_particles (gas): ", np.mean(self.gas_particles.x.in_(units.au)))
                 
        print("Done.")

    def evolve_external_only(self, tend, plot=False, plot3d = False, verbose=False):
        time = 0 | units.yr
        
        if self.external_object is None:
            print("No external object assigned")
            return
            
        while time < tend:
            #self.hydro.evolve_model(time)
            self.external_code.evolve_model(time)
            print(f"Evolution: {time.in_(tend.unit)}")
            time += self.external_code.parameters.timestep

            if plot:
                plt.figure()
                scatter(self.external_code.particles.x.in_(units.au), self.external_code.particles.y.in_(units.au))
                plt.xlim(-10,10)
                plt.ylim(-10,10)
                plt.title(f"{time.number:.3f}{time.unit}")
                plt.savefig(f"disk_{time.number:.3f}.png")
                
            if plot3d:
                self.plot3d(show=False, save=True, time=time)
    
            if verbose:
                print(f"External object mean energy: {np.mean(self.external_object.u)}")
                 
        print("Done.")
    
    def plot_system(self, show=False, save=False, time=None):
        plt.figure()
        if self.external_object is None:
            l = 10
        else: 
            l = max(self.external_object.x.in_(units.au)).number * 2
        plt.xlim(-l,l)
        plt.ylim(-l,l)
    
        scatter(self.hydro.gas_particles.x.in_(units.au), self.hydro.gas_particles.y.in_(units.au), c='blue', alpha=0.1, label="SN")
        scatter(self.hydro.dm_particles.x.in_(units.au), self.hydro.dm_particles.y.in_(units.au), marker='*', c='yellow', label="Core")
        if self.external_object is not None:
            try:
                scatter(self.external_code.particles.x.in_(units.au), self.external_code.particles.y.in_(units.au), c='orange', alpha=0.1, label=self.external_object.name[0])
            except:
                scatter(self.external_code.particles.x.in_(units.au), self.external_code.particles.y.in_(units.au), c='orange', alpha=0.1, label="External Object")
         
        plt.legend(loc='upper right')
        
        if time is not None:
            plt.title(f"{time.number:.3f} {time.unit}")
        if show:
            plt.show()
        if save:
            plt.savefig(f"nova_{time.number:.3f}.png")
            plt.close() 


            

    def plot3d(self, show=False, save=False, time=None):
        x = self.hydro.gas_particles.x.in_(units.au)
        y = self.hydro.gas_particles.y.in_(units.au)
        z = self.hydro.gas_particles.z.in_(units.au)
        
        # Create a 3D plot
        fig = plt.figure()
        l = 10
        ax = fig.add_subplot(111, projection='3d')
    
        ax.set_xlim(-l,l)
        ax.set_ylim(-l,l)
        ax.set_zlim(-l,l)
    
        # Scatter plot of the particles
        ax.scatter(x.number, y.number, z.number, s=1, c='blue', marker='o', label="SN")  # Adjust size and color as needed
        if self.external_object is not None:
            ex = self.external_object.x.in_(units.au)
            ey = self.external_object.y.in_(units.au)
            ez = self.external_object.z.in_(units.au)
            try:
                ax.scatter(ex.number, ey.number, ez.number,c='orange', alpha=0.1, label=self.external_object.name[0])
            except:
                 ax.scatter(ex.number, ey.number, ez.number, c='orange', alpha=0.1, label="External Object")
     
        plt.legend(loc='upper right')
        # Set labels
        ax.set_xlabel('X au')
        ax.set_ylabel('Y au')
        ax.set_zlabel('Z au')
        if time is not None:
            plt.title(f"{time.number:.3f}{time.unit}")
        if show:
            plt.show()
        if save:
            plt.savefig(f"nv3d_{time.number:.3f}.png")
            plt.close()

    def __str__(self):
        print("Supernova class WIP")
        return ""
    

        if self.external_object is not None:
          print(f"The supernova is hitting a body at position:")
          print(self.external.object.position)
        return ""
