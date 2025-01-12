from binary_disk import BinaryDisk  # %run on juputer lab
from supernova import Supernova

# Problem with this code: you still have to change the position of the disk manually from the class code (don't move the supernova or it will fuck up)
# i'll implement a way to call the disk position from the initialisation

# Let's evolve a binary disk before 
system = BinaryDisk(components="all")
system.evolve(50 | units.yr, plot=True)

# or import the saved one (50.hdf5 is a system evolved for 55 or 56 yrs)
imported = read_set_from_file('50.hdf5', 'amuse')
system = BinaryDisk(from_set=imported)

# Let's have the supernova hit it
pickle = "tmpow7eucrj/test.pkl"
sn = Supernova(pickle=pickle, nparticles=10000, external_object=system.all_particles)
sn.evolve(2 | units.day, plot=True)

# Now let's look at the effects on the disk...
system.evolve(10 | units.yr, plot=True)


# Make the animation (in jupyter lab you can run this between each evolution, i could also just make this a function)
files = glob.glob("*.png") # This will select all of your pngs

numbers = []
for f in files:
    numbers.append(float(f[5:10]))
ff = [y for _,y in sorted(zip(numbers, files))]

with imageio.get_writer("disk_before.gif", mode='I', duration=0.1) as writer:
    for frame in ff:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory 
for frame in ff:
    os.remove(frame)
