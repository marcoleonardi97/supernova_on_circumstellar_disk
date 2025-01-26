
from binary_disk import BinaryDisk  # %run binary_disk.ipynb on jupyter lab


# To evolve this system with bridge you can simply call system.evolve():
# To check the particles and model_time of the system you can print the class: print(system)

# This is an example of the real part 1 code, run for however long you want.
# The 4 simulations to present are saved in the backups folder (evolved for 50 years)
system = BinaryDisk(rin = 1 | units.au, rout = 5 | units.au, semimaj = 25 | units.au, density=10e-10 | (units.g / units.cm**3) components="all")
system.plot_system(show=True) # This will show a snapshot of the system
t_100 | units.yr
system.evolve(t_end, plot=True, backup=True)


# Animation - this method only works if you don't have other pngs in your directory.
files = glob.glob("*.png")

numbers = []
for f in files:
    numbers.append(float(f[5:10]))  # change according to whatever name you use for the frames, this work with f"disk_{time:.3f}.png"
ff = [y for _,y in sorted(zip(numbers, files))]

with imageio.get_writer("disk_your_parameters.gif", mode='I', duration=0.1) as writer:
    for frame in ff:
        image = imageio.imread(frame)
        writer.append_data(image)

# Clean up the directory - this should not remove all of your pngs, but be careful anyways
for frame in ff:
    os.remove(frame)
