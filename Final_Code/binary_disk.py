import glob
import os
import numpy as np

from matplotlib import pyplot as plt
from amuse.plot import plot, scatter
from amuse.units import units, constants
import imageio.v2 as imageio
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D


from amuse.lab import Particles
from amuse.units import nbody_system

from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

from amuse.community.ph4.interface import Ph4
from amuse.community.fi.interface import Fi
from amuse.couple import bridge
from amuse.ext.composition_methods import *
from amuse.support import options

from amuse.io import write_set_to_file
from amuse.io import read_set_from_file

from numba import jit # Currently not used



# In order to debug the individual codes, i added the option to create the single binary pair or the single disk,
# and to use only the gravity code or only the hydro code if needed. 
# When you initiate the class, just call: components = "all" or "disk" or "stars".


# the following fixes are highly recommended
#allow oversubscription for openMPI
os.environ["OMPI_MCA_rmaps_base_oversubscribe"]="true"
# use lower cpu resources for idle codes
options.GlobalOptions.instance().override_value_for_option("polling_interval_in_milliseconds", 10)



class BinaryDisk(object):
    def __init__(self, from_set=None, rin= 1 | units.au, rout= 5 | units.au, semimaj = 15 | units.au, ndisk=1000, components = "all"):
        self.components = components
        self.m1 = 1 | units.MSun
        self.m2 = 0.5 | units.MSun
        self.semimaj = semimaj
        self.ecc = 0.6
        self.rin = rin
        self.rout = rout
        self.ndisk = ndisk
        self.pinner = (((4 * np.pi**2) * self.rout**3)/(constants.G * (self.m1+self.m2))).sqrt()
        self.particles = Particles()
        self.gas_particles = Particles()
        self.all_particles = Particles()
        self.system_time = 0 | units.yr
        self.converter = nbody_system.nbody_to_si(self.m1, 1 | units.au)

        if from_set is None:

            if components == "all": 
                self.setup()
    
            elif components == "stars":
                self.make_stars()
    
            elif components == "disk":
                self.make_disk()
                
            elif components not in {"all", "stars", "disk"}:
                raise ValueError(f"Invalid value for parameter: {components}. Must be one of 'all', 'stars', 'disk'.")
        else:
            self.all_particles.add_particles(from_set)
            self.gas_particles.add_particles(from_set[from_set.name=="disk"])
            self.particles.add_particles(from_set[from_set.name!="disk"])

        
        # Grav
        self.gravity = Ph4(self.converter)
        self.gravity.particles.add_particles(self.all_particles)
        self.gravity.parameters.timestep_parameter = 0.07
        self.gravity.parameters.epsilon_squared = 1./self.ndisk**(2./3) | units.au**2  # 0.01 au **2 with 1000 particles, sometimes it works sometimes it doesn't

        # Hydro
        
        self.hydro = Fi(self.converter, mode="openmp", redirection='none')

        self.hydro.parameters.use_hydro_flag = True
        self.hydro.parameters.radiation_flag = False
        self.hydro.parameters.self_gravity_flag = True
        self.hydro.parameters.gamma = 1
        self.hydro.parameters.isothermal_flag = True
        self.hydro.parameters.integrate_entropy_flag = False
        self.hydro.parameters.timestep = self.gravity.parameters.timestep_parameter | units.yr #0.125 | units.yr
        self.hydro.parameters.verbosity = 1
        self.hydro.parameters.eps_is_h_flag = True  # True = h_smooth is not constant, False=constant
        eps = 0.1 | units.au

        self.hydro.parameters.gas_epsilon = eps 
        #self.hydro.parameters.n_smooth_tol = 0.1
        #self.hydro.parameters.sph_h_const = eps * 5


        # Disk
        self.hydro.particles.add_particles(self.gas_particles)
        if self.components == "disk":
            self.hydro.dm_particles.add_particles(self.particles)



    
    def make_stars(self):
        stars = new_binary_from_orbital_elements(self.m1, self.m2, self.semimaj, self.ecc, G=constants.G)
        stars[0].name = "Primary"
        stars[1].name = "Secondary"
        self.particles.add_particles(stars)
        self.all_particles.add_particles(stars)

    def make_disk(self):
        if self.components == "disk":
            s = Particles(1)
            s.mass = self.m1
            s.position = (0, 0, 0) | units.au
            s.velocity = (0, 0, 0) | units.kms
            self.particles.add_particles(s)
            
        R = 1 | units.au
        Mdisk = self.m1 * 0.01
        disk = ProtoPlanetaryDisk(self.ndisk,
                          convert_nbody=self.converter,
                          Rmin=self.rin/R,
                          Rmax=self.rout/R,
                          q_out=1.0,
                          discfraction=Mdisk/self.m1).result
        disk.name = "disk"
        disk.move_to_center()
        disk.position += self.particles[0].position
        disk.velocity += self.particles[0].velocity
        
        masses = Mdisk / float(self.ndisk)
        disk.mass = masses
        rho = 3.0 | (units.g / units.cm**3)
        disk.radius = (disk.mass / (4 * rho))**(1./3.)
        #disk.h_smooth = self.rin ? 
        self.gas_particles.add_particles(disk)
        self.all_particles.add_particles(disk)

        
    def setup(self):
        self.make_stars()
        self.make_disk()



    def plot_system(self, show=False, part="all", save=False, save_name=None, time=None):
        """
        show: bool; show the plot, used when you're just looking at the system without evolving.
        part: list or str; select which source of particles to plot. Accepted inputs
        are: combinations of ["hydro", "gravity", "gas", "stars"], "all".
        save: bool; flag to save plots
        save_name: str; file names
        time: str; plot titles. 
        """
        plt.figure()
        star = self.all_particles[self.all_particles.name=="Primary"]
        planet = self.all_particles[self.all_particles.name=="Secondary"]
        l = 2 * self.semimaj.number
        plt.xlim(-l,l)
        plt.ylim(-l,l)
        if "hydro" in part:
            scatter(self.hydro.particles.x.in_(units.AU), self.hydro.particles.y.in_(units.AU), c='blue', alpha=0.5, s=10)
        if "gravity" in part:
            scatter(self.gravity.particles.x.in_(units.AU), self.gravity.particles.y.in_(units.AU), c='yellow')
        if "gas" in part:
            scatter(self.gas_particles.x.in_(units.AU), self.gas_particles.y.in_(units.AU), c='blue', alpha=0.5, s=10)
        if "stars" in part:
            scatter(star.x, star.y, marker="*",c='r', s=120,label="Primary Star")
            scatter(planet.x, planet.y, marker='*', c='y',s=120, label="Secondary Star")
        if part == "all":
            scatter(self.all_particles.x.in_(units.AU), self.all_particles.y.in_(units.AU), c='orange', alpha=0.5, s=10)
            scatter(star.x, star.y, marker="*",c='r', s=120,label="Primary Star")
            scatter(planet.x, planet.y, marker='*', c='y',s=120, label="Secondary Star")
            
        plt.legend(loc='upper right')
        plt.close()
        if time is not None:
            plt.title(f"Disk at {time.number:.2f} {time.unit}")

        if show:
            plt.show()
    
        if save:
            plt.savefig(save_name)
            plt.close()

       def plot3d(self, show=False, save=False, savename=None, time=None):
    
        x = self.all_particles.x.in_(units.au)
        y = self.all_particles.y.in_(units.au)
        z = self.all_particles.z.in_(units.au)
        star = self.all_particles[self.all_particles.name=="Primary"]
        planet = self.all_particles[self.all_particles.name=="Secondary"]

        # Create a 3D plot
        fig = plt.figure()
        l = 2 * self.semimaj.number
        ax = fig.add_subplot(111, projection='3d')
    
        ax.set_xlim(-l,l)
        ax.set_ylim(-l,l)
        ax.set_zlim(-l,l)
    
    
        ax.scatter(x.number, y.number, z.number, s=1, c='blue', marker='o', alpha=0.5)  # Adjust size and color as needed
        ax.scatter(star.x.number, star.y.number, star.z.number, marker='*', c='r', s=120, label='Primary Star')
        ax.scatter(planet.x.number, planet.y.number, planet.z.number, marker='*', c='y', s=120, label='Secondary Star')
        plt.legend(loc='upper right')
        # Set labels
        ax.set_xlabel('X au')
        ax.set_ylabel('Y au')
        ax.set_zlabel('Z au')
        
        if time is not None:
            plt.title(f"Disk at: {time.number:.3f} {time.unit}")
        if show:
            plt.show()
        if save:
            plt.savefig(savename)
            plt.close()


    def channel(self):
        channel = {"from_stars": self.all_particles.new_channel_to(self.gravity.particles),
                  "to_stars": self.gravity.particles.new_channel_to(self.all_particles),
                  "from_disk": self.all_particles.new_channel_to(self.hydro.particles),
                  "to_disk": self.hydro.particles.new_channel_to(self.gas_particles),
                  "h_to_all": self.hydro.particles.new_channel_to(self.all_particles)}

        return channel

    #@jit(nopython=True)
    def bridge(self):
        gravhydro = bridge.Bridge(use_threading=False, method=SPLIT_4TH_S_M4)
        gravhydro.add_system(self.gravity, (self.hydro,))
        gravhydro.add_system(self.hydro, (self.gravity,))
        gravhydro.timestep = self.hydro.parameters.timestep
        return gravhydro
        
    #@jit(nopython=True)
    def evolve(self, tend, display="all", plot=False, backup = False, 
               backup_file = f"simulation_backup.hdf5", backup_dt = 10, verbose=False):

        if tend < self.system_time:
            print(f"System is already evolved to {self.system_time}")
            print(f"(Plots will looks still until the code evolves past that.)")
            
        if self.components == "stars":
            self.evolve_gravity_only(tend, display, plot, verbose)
        elif self.components == "disk":
            self.evolve_hydro_only(tend, display, plot, verbose)
        else:
            channel = self.channel()
            code = self.bridge()
            time = self.system_time
            dt =  code.timestep
            loops = 0
            
            print("Starting simulation...")
            while time < tend:
                code.evolve_model(time)
                time += dt
                loops += 1
                
                print(f"System evolved to: {time.in_(tend.unit)}")
                channel["to_stars"].copy()
                channel["to_disk"].copy()

                if backup == True and loops % backup_dt == 0:
                    write_set_to_file(self.all_particles.savepoint(time), backup_file, 'amuse', overwrite_file=True)
                
                if verbose:
                    print("------------------", "\n")
                    print("Gravity particles (x pos): ", np.mean(self.gravity.particles.x.in_(units.au)))
                    print("Hydro particles: ", np.mean(self.hydro.particles.x.in_(units.au)))
                    print("Local self.particles (stars): ", np.mean(self.particles.x.in_(units.au)))
                    print("Local self.gas_particles (disk): ", np.mean(self.gas_particles.x.in_(units.au)))
                    print("Local self.all_particles: ", np.mean(self.all_particles.x.in_(units.au)))
                    print("\n")
                    
                if plot:
                    self.plot_system(part=display, save=True, save_name=f"disk_{time.number:.3f} {time.unit}.png", time=time)
            print("Done.")
            #self.gravity.stop() better to stop them manually when you want
            #self.hydro.stop()

    def evolve_gravity_only(self, tend, part=["stars", "gravity"], plot=False, verbose=False):
        channel = self.channel()
        code = self.gravity
        time = self.system_time
        dt = self.gravity.parameters.timestep_parameter | units.yr
    
        print("Starting simulation...")
        while time < tend:
            time += 1 | units.yr
            code.evolve_model(time)
            print(f"System evolved to: {time}")
            channel["to_stars"].copy()
            channel["from_stars"].copy()
            if verbose:
                print("------------------", "\n")
                print("Gravity particles (x pos): ", np.mean(self.gravity.particles.x.in_(units.au)))
                print("Local self.particles (stars): ", np.mean(self.particles.x.in_(units.au)))
                print("Local self.all_particles: ", np.mean(self.all_particles.x.in_(units.au)))
                print("\n")
                
            if plot:
                self.plot_system(part=part, save=True, save_name=f"disk_{time.number:.3f} {time.unit}.png", time=time)
        print("Done.")
        

    def evolve_hydro_only(self, tend, part=["gas","hydro"], plot=False, verbose=False):
        channel = self.channel()
        code = self.hydro
        time = self.system_time
        dt = code.parameters.timestep
        
        print("Starting simulation...")
        while time < tend:
            time += dt
            code.evolve_model(time)
            print(f"System evolved to: {time}")
            channel["h_to_all"].copy()

            if verbose:
                print("------------------", "\n")
                print("Hydro particles: ", np.mean(self.hydro.particles.x.in_(units.au)))
                print("Local self.gas_particles (disk): ", np.mean(self.gas_particles.x.in_(units.au)))
                print("Local self.all_particles: ", np.mean(self.all_particles.x.in_(units.au)))
                print("\n")
                
            if plot:
                self.plot_system(part=part, save=True, save_name=f"disk_{time.number:.3f} {time.unit}.png", time=time)
        print("Done.")


    def evolve_without_bridge(self, dt, tend, part="all",plot=False, verbose=False):
        channel = self.channel()
        time = self.system_time
        print("Starting simulation without bridge...")
        
        while time < tend:
            time += dt
            self.gravity.evolve_model(time)
            print(f"Ph4 evolved to: {time}")
            self.hydro.evolve_model(time)
            print(f"Fi evolved to: {time}")
            channel["to_stars"].copy()
            channel["h_to_all"].copy()

            if plot:
                self.plot_system(part=part, save=True, save_name=f"disk_{time.number:.3f} {time.unit}.png", time=time)

    
    
    def save_particles(self, filename, memory="all", ow=False):
        try:
            if memory == "all":
                write_set_to_file(self.all_particles, f'{filename}.hdf5', 'amuse', overwrite_file=ow)
                print(f"Saved the system at time {self.system_time} in {filename}.hdf5")
            elif memory == "gas":
                write_set_to_file(self.gas_particles, f'{filename}.hdf5', 'amuse', overwrite_file=ow)
                print(f"Saving the system at time {self.system_time} in {filename}.hdf5")
            elif memory == "particles":
                write_set_to_file(self.particles, f'{filename}.hdf5', 'amuse', overwrite_file=ow)
                print(f"Saving the system at time {self.system_time} in {filename}.hdf5")
                
        except:
            print("This file already exists, but overwrite is set to False.","\n",
                 "Please set a new file name, or set ow=True to overwrite.")

    def __str__(self):
        d = {"all": "protoplanetary disk in a binary star system",
            "disk": "protoplanetary disk",
            "stars": "binary star system"}
        print (f"This class represents a {d[self.components]}. Currently at time {self.system_time}", "\n"
        "The local particles storages are:", "\n", 
        f"1. Stars (self.particles): {len(self.particles)}", "\n",
        f"2. Gas particles (self.gas_particles): {len(self.gas_particles)}", "\n",
        f"3. All particles (self.all_particles): {len(self.all_particles)}", "\n",
        f"Particles in the gravity code: {len(self.gravity.particles)}", "\n"
        f" Particles in the hydro code: {len(self.hydro.particles)}")
        return ""

