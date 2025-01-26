from binary_disk import BinaryDisk  # %run on juputer lab
from supernova import Supernova


# Let's evolve a normal disk before (no binary)
system = BinaryDisk(components="disk")
system.move_system((10, 0, 0) | units.au)
system.evolve(40 | units.yr, plot=True)

# Let's have the supernova hit it
pickle = "tmpow7eucrj/test.pkl"
sn = Supernova(pickle=pickle, nparticles=10000, external_object=system.all_particles)
sn.evolve(4.8 | units.day)

# Now let's look at the effects on the disk
system.evolve(40 | units.yr, plot=True)


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
