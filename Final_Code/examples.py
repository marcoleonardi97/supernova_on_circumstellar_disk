from binary_disk import BinaryDisk # %run on jupyter
from supernova import Supernova

# PART 1 examples ---------------------------------------------------------------------------
# Make a system using the class and choosing the components (="all" by default).
# you can evolve using system.evolve(tend, display="all", plot=True, verbose=False) or system.evolve_without_bridge(tend, plot=True).
# system.plot_system(show=True) will show a snapshot of the system
# print(system) will show the model_time of the system and where the particles are currently stored.

tend = 10 | units.yr

system = BinaryDisk(rin = 0.5 | units.au, rout = 1 | units.au, components = "all") # Simple small disk
system.plot_system(show=True)
system.evolve(tend, plot=True, verbose=False)

system = BinaryDisk(rin = 0.5 | units.au, rout = 1 | units.au, components = "all") # 3D Plot! Just swap plot=True for plot3d=True. Don't activate both at once.
system.plot3d(show=True)
system.evolve(tend, plot3d=True, verbose=False)

system = BinaryDisk(rin = 1 | units.au, rout = 10 | units.au, components = "disk") # Disk only without the binary pair
system.plot_system(show=True)
system.evolve(tend, plot=True, verbose=False)

system = BinaryDisk(rin = 4 | units.au, rout = 5 | units.au, components = "stars") # Binary orbit, no disk
system.plot_system(show=True)
system.evolve(tend, plot=True, verbose=False)


#------- How to export and import existing systems ---------#

# You can save a system using the class function save_particles(filename, memory=default, ow(overwrite)=False)
system = BinaryDisk(rin = 1 | units.au, rout = 5 | units.au, semimaj = 15 | units.au, components="all")
system.evolve(tend, plot=False, verbose=False)
system.save_particles("test") # This will output the system.all_particles into an hdf5 file in your directory.

# To reopen a saved system as a BinaryDisk class, you can call the parameter from_set at initialisation:
saved_system = read_set_from_file('test.hdf5', 'amuse', close_file=True)
new_system = BinaryDisk(from_set = saved_system)
new_system.plot_system(show=True)

# If you want to save during a simulation because it may crash, you can use the backup parameter in evolve():
# backup_dt: the system will be backed up to the file every backup_dt timesteps. 
system.evolve(tend, plot=True, verbose=False, backup=True, backup_file="backup.hdf5", backup_dt = 10)


# PART 2 examples --------------------------------------------------------------------------------

# Make a system and move it at the desired distance from the supernova.
# The supernova will always explode at position (0,0,0).
system = BinaryDisk(components="disk")
p = (10, 0, 0) | units.au
system.move_system(p)

# Insert the system's particles in a supernova object using the external_object parameter
pickle = "tmpow7eucrj/test.pkl"
sn = Supernova(pickle=pickle, nparticles=10000, external_object=system.all_particles) # maybe we can also just use system.gas_particles here
sn.evolve(3| units.day, plot=True)

# Evolve the object to look at the effects
system.evolve(5 | units.yr, plot=True)

# Notes: the class object remembers if you evolved a system before, so it will keep simulating from the last timestep unless you re-initialise.
# There's no problem with this if you call components='all' and evolve with bridge, but there will be some strange stuff happening if you evolve a disk only (part2)
# nothing breaking but you may get some still frames in your gifs.




# All of these simulations will fill your directory with images called "disk{time:.3f}.png", run this code to glue them together in a gif:

import glob
import os
import imageio.v2 as imageio
from matplotlib.animation import FuncAnimation, PillowWriter

# Make sure you don't have other .png files in your directory before running this code.
files = glob.glob("*.png") # This will select all of your pngs

numbers = []
for f in files:
    numbers.append(float(f[5:10]))
ff = [y for _,y in sorted(zip(numbers, files))]

with imageio.get_writer("name.gif", mode='I', duration=0.1) as writer:
    for frame in ff:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory 
for frame in ff:
    os.remove(frame)
