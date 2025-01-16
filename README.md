The final code uses two different custom classes to simulate all the systems we need: BinaryDisk() and Supernova().
These allow for good customisation and debugging (i hope), simple evolution and for different methods of plotting.
They are kind of documented in the examples folder.  
Presentation WIP: https://docs.google.com/presentation/d/1C42A9Nemq37XhAS2Ti2RfmwunVOVZkHnxbtuOxaauZY/edit?usp=sharing

-------------------------------

Simulating a supernova explosion in the vicinity of a proto-planetary disk, with and without an additional perturbing object. 

1. Make a protoplanetary disk around a 1 MSun star, orbited by a smaller (< 0.5 MSun) object that inclines and flares the disk. 
2. Explode a ($\approx 30 M_\odot$) star in the vicinity (30 AU) without the perturbing object
3. Explode the same star in the vicinity of the perturbed disk, monitoring the effect on the flares and structures of the disk. 

<ins>Minimum passing grade plan</ins>

Evolve each step as its own system and comment on the results:

1. Simulation (1) should clearly show what kind of perturbations and structures arise from having a perturbing object orbiting the disk.
  - This will include bridging gravity with hydrodynamics, understanding the orbital parameters of the binary, plotting and animating the dynamics of the system over an extended period of time (~ 100 years).
2. Simulation (2) will show the effects of a nearby supernova on a vanilla disk, for reference.
  - This also needs a bridge between gravity and hydrodynamics and good control over the initial conditions (disk face-on, inclination, distance, disk outer density...)
3. Simulation (3) will merge the two above results showing the effects of the SN blast wave on the flares and protuberances of the perturbed disk.
  - Depending on the cost of the calculation we shall run the perturbed disk system again or use the results from simulation (1). This is the most important part of the project, so we will focus on showing the results via plots and animations wherever possible, but also on monitoring important parameters like surface temperature, surface density, Toomre Q parameter ...)




Results:

Disk in a binary system made with the code in "final_code/main_part_1.py" over 56 years (crashed at 56). It looks a little slow so i might change the timesteps a bit
Parameters: disk radius 1-5 au, semi major axis: 25 (i think). 

![15au_55years_5disk](https://github.com/user-attachments/assets/bc597912-7fc9-4ee1-adf1-152972bb3cd5)




Fast simulation of the disk without bridging hydro!

![no_bridge](https://github.com/user-attachments/assets/d8c20d91-d780-4060-b3c8-ef337908a3a5)



Supernova explosion using the code "supernova.py", looks like with the current parameters it reaches max distances after about 10 days:

![sn_test](https://github.com/user-attachments/assets/e4461af8-fdc9-4e71-9692-e943b66d6b77)

3d plot of a supernova hitting a vanilla disk. Still problems with the interaction of the two systems (this is using the supernova

![nova3d (1)](https://github.com/user-attachments/assets/e74787a7-98d6-4b4c-98c5-f2a799b438d2)
_on_disk file)


The disk does work but it runs on its own timestep... not sure how to fix this yet

![disk3d](https://github.com/user-attachments/assets/4dbe4a79-f261-4b9f-8d63-c1ecce0a7181)


