from binary_disk import BinaryDisk  # %run on juputer lab
from supernova import Supernova
import glob
import os
import imageio.v2 as imageio
from matplotlib.animation import FuncAnimation, PillowWriter

def make_animation(name, files):
    numbers = []
    for f in files:
        numbers.append(float(f[5:10]))
    
    ff = [y for _,y in sorted(zip(numbers, files))]
    
    with imageio.get_writer(f"{name}.gif", mode='I', duration=0.1) as writer:
        for frame in ff:
            image = imageio.imread(frame)
            writer.append_data(image)
    
    # Clean up the directory 
    for frame in ff:
        os.remove(frame)

# import each system that we evolved in part 1
imported = read_set_from_file('sim1_50_yr.hdf5', 'amuse')
system = BinaryDisk(from_set=imported)

# Let's have the supernova hit it
pickle = "tmpow7eucrj/test.pkl"
sn = Supernova(pickle=pickle, nparticles=10000, external_object=system.all_particles)
sn.evolve(4.8 | units.day, plot3d=True) # sometimes it crashes if you evolve for more than 5 days...


files = glob.glob("*.png") # This will select all of your pngs
make_animation('sim1_nova3d', files)


# Now let's look at the effects on the disk...
system.evolve(20 | units.yr, display=["gas", "stars"], plot=True) # you have to display gas and stars instead of "all" for energy transfer for some reason # i think this is fixed now and you don't have to do it

files = glob.glob("*.png") # This will select all of your pngs
make_animation('sim1_after_nova', files)
