from binary_disk import BinaryDisk

# Make a system using the class and choosing the components (="all" by default).
# you can evolve using system.evolve(tend, display="all", plot=True, verbose=False) or system.evolve_without_bridge(tend, plot=True).
# system.plot_system(show=True) will show a snapshot of the system
# print(system) will show the model_time of the system and where the particles are currently stored.

tend = 10 | units.yr

system = BinaryDisk(rin = 4 | units.au, rout = 5 | units.au, components = "all") # Cool ring figure
system.plot_system(show=True)
system.evolve(tend, plot=True, verbose=False)

system = BinaryDisk(rin = 0.5 | units.au, rout = 1 | units.au, components = "all") # Small disk
system.plot_system(show=True)
system.evolve(tend, plot=True, verbose=False)

system = BinaryDisk(rin = 1 | units.au, rout = 10 | units.au, components = "disk") # Disk only without the binary pair
system.plot_system(show=True)
system.evolve(tend, plot=True, verbose=False)

system = BinaryDisk(rin = 4 | units.au, rout = 5 | units.au, components = "stars") # Binary orbit, no disk
system.plot_system(show=True)
system.evolve(tend, plot=True, verbose=False)

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
