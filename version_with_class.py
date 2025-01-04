import glob
import os
import numpy as np

from matplotlib import pyplot as plt
from amuse.plot import plot, scatter
from amuse.units import units, constants
import imageio.v2 as imageio
from matplotlib.animation import FuncAnimation, PillowWriter

from amuse.lab import Particles
from amuse.units import nbody_system

from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

from amuse.community.ph4.interface import Ph4
from amuse.community.fi.interface import Fi
from amuse.couple import bridge
from amuse.ext.composition_methods import *

from numba import jit



# This version uses a class to create and evolve the system.
# In order to debug the individual codes, i added the option to create the single binary pair or the single disk,
# and to use only the gravity code or only the hydro code if needed. 
# When you initiate the class, just call: components = "all" or "disk" or "stars".

class BinaryDisk(object):
    def __init__(self, ndisk=1000, components = "all"):
        self.components = components
        self.m1 = 1 | units.MSun
        self.m2 = 0.5 | units.MSun
        self.semimaj = 15 | units.au
        self.ecc = 0.6
        self.ndisk = ndisk
        self.converter = nbody_system.nbody_to_si(self.m1, 1 | units.au)
        self.rin = 1 | units.au   # Disk inner radius
        self.rout = 5 | units.au  # Disk outer radius
        self.pinner = (((4 * np.pi**2) * self.semimaj**3)/(constants.G * (self.m1+self.m2))).sqrt() # Rev. Period of a disk particle in the innner disk
        self.particles = Particles()
        self.gas_particles = Particles()
        self.all_particles = Particles()

        if components == "all": 
            self.setup()

        elif components == "stars":
            self.make_stars()

        elif components == "disk":
            self.make_disk()
            
        elif components not in {"all", "stars", "disk"}:
            raise ValueError(f"Invalid value for parameter: {components}. Must be one of 'all', 'stars', 'disk'.")


       
        # This is a little messy, initialising all of these outside of a function, but i didn't want to call them by accident and restart them
      
        # Grav
        self.gravity = Ph4(self.converter)
        self.gravity.particles.add_particles(self.all_particles)

        # Hydro
        self.hydro = Fi(self.converter, mode="openmp")
        self.hydro.parameters.use_hydro_flag = True
        self.hydro.parameters.radiation_flag = False
        self.hydro.parameters.self_gravity_flag = True
        self.hydro.parameters.gamma = 1
        self.hydro.parameters.isothermal_flag = True
        self.hydro.parameters.integrate_entropy_flag = False
        self.hydro.parameters.timestep = 0.01 * self.pinner / 8
        self.hydro.parameters.verbosity = 1
        self.hydro.parameters.eps_is_h_flag = True  # True = h_smooth is not constant, False=constant
        eps = 0.06 | units.au
        self.hydro.parameters.gas_epsilon = eps
        self.hydro.parameters.sph_h_const = eps 
        self.hydro.particles.add_particles(self.gas_particles)
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
                          q_out=10.0,
                          discfraction=Mdisk/self.m1).result
        disk.name = "disk"
        disk.move_to_center()
        disk.position += self.particles[0].position
        disk.velocity += self.particles[0].velocity
        
        masses = Mdisk / float(self.ndisk)
        disk.mass = masses
        rho = 3.0 | (units.g / units.cm**3)
        disk.radius = (disk.mass / (4 * rho))**(1./3.)
        self.gas_particles.add_particles(disk)
        self.all_particles.add_particles(disk)

    def setup(self):
        self.make_stars()
        self.make_disk()

    def plot_system(self, save=False, save_name=None, time=None):
      """
      save: bool 
      save_name: str
      time: str
      """
        plt.figure()
        star = self.all_particles[self.all_particles.name=="Primary"]
        planet = self.all_particles[self.all_particles.name=="Secondary"]
        l = 2 * self.semimaj.number
        plt.xlim(-l,l)
        plt.ylim(-l,l)
        if not self.components == "stars":
            scatter(self.hydro.particles.x.in_(units.AU), self.hydro.particles.y.in_(units.AU), c='blue', alpha=0.5, s=10)
        if not self.components == "disk":
            scatter(self.all_particles.x.in_(units.AU), self.all_particles.y.in_(units.AU), c='orange', alpha=0.5, s=10)
        scatter(star.x, star.y, marker="*",c='r', s=120,label="Primary Star")
        scatter(planet.x, planet.y, marker='*', c='y',s=120, label="Secondary Star")
        plt.legend(loc='upper right')
        if time is not None:
            plt.title(f"Disk at {time.number:.2f} {time.unit}")
    
        if save:
            plt.savefig(save_name)


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
        gravhydro.timestep = 2 * self.hydro.parameters.timestep
        return gravhydro
        
    #@jit(nopython=True)
    def evolve(self, tend, plot=False, verbose=False):
        if self.components == "stars":
            self.evolve_gravity_only(tend, plot, verbose)
        elif self.components == "disk":
            self.evolve_hydro_only(tend, plot, verbose)
        else:
            channel = self.channel()
            code = self.bridge()
            time = 0 | units.yr
            dt = code.timestep
            print("Starting simulation...")
            while time < tend:
                time += dt
                code.evolve_model(time)
                print(f"System evolved to: {time}")
                channel["from_stars"].copy()
                channel["from_disk"].copy()
                if verbose:
                    print("------------------", "\n")
                    print("Gravity particles (x pos): ", np.mean(self.gravity.particles.x.in_(units.au)))
                    print("Hydro particles: ", np.mean(self.hydro.particles.x.in_(units.au)))
                    print("Local self.particles (stars): ", np.mean(self.particles.x.in_(units.au)))
                    print("Local self.gas_particles (disk): ", np.mean(self.gas_particles.x.in_(units.au)))
                    print("Local self.all_particles: ", np.mean(self.all_particles.x.in_(units.au)))
                    print("\n")
                    
                if plot:
                    self.plot_system(save=True, save_name=f"disk_{time.number:.3f} {time.unit}.png", time=time)
            print("Done.")
            self.gravity.stop()
            self.hydro.stop()

    def evolve_gravity_only(self, tend, plot=False, verbose=False):
        channel = self.channel()
        code = self.gravity
        time = 0 | units.yr
        dt = 1 | units.yr
        print("Starting simulation...")
        while time < tend:
            time += dt
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
                self.plot_system(save=True, save_name=f"disk_{time.number:.3f} {time.unit}.png", time=time)
        print("Done.")
        

    def evolve_hydro_only(self, tend, plot=False, verbose=False):
        channel = self.channel()
        code = self.hydro
        time = 0 | units.yr
        dt = code.parameters.timestep
        print("Starting simulation...")
        while time < tend:
            time += dt
            code.evolve_model(time)
            print(f"System evolved to: {time}")
            #scatter(self.hydro.particles.x, self.hydro.particles.y)
            #plt.savefig(f"{time}.png")
            if verbose:
                print("------------------", "\n")
                print("Hydro particles: ", np.mean(self.hydro.particles.x.in_(units.au)))
                print("Local self.gas_particles (disk): ", np.mean(self.gas_particles.x.in_(units.au)))
                print("Local self.all_particles: ", np.mean(self.all_particles.x.in_(units.au)))
                print("\n")
                
            if plot:
                self.plot_system(save=True, save_name=f"disk_{time.number:.3f} {time.unit}.png", time=time)
        print("Done.")


    def __str__(self):
        d = {"all": "protoplanetary disk in a binary star system",
            "disk": "protoplanetary disk",
            "stars": "binary star system"}
        return f"This class represents a {d[self.components]}. "



# To evolve this system you can simply call evolve():

system = BinaryDisk(1000, components="stars")  # Running this code with the system as 'stars' will just animate the stars orbit.
t_end = 100 | units.yr
system.evolve(t_end, plot=True, verbose=True)



# -------------------- This code is only needed if we don't name the plots with plot_counter
files = glob.glob("*.png")

numbers = []
for f in files:
    numbers.append(float(f[5:10]))  # change according to whatever name you use for the frames, this work with f"disk_{time:.3f}.png"
ff = [y for _,y in sorted(zip(numbers, files))]

with imageio.get_writer("disk.gif", mode='I', duration=0.1) as writer:
    for frame in ff:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory 
for frame in ff:
    os.remove(frame)
