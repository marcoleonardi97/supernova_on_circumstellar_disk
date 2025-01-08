from binary_disk import BinaryDisk

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

with imageio.get_writer("15au_10years_5disk.gif", mode='I', duration=0.1) as writer:
    for frame in ff:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory 
for frame in ff:
    os.remove(frame)
