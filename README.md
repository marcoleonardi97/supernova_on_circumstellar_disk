Outline (old):
https://docs.google.com/document/d/1dHtwYFsnluikh4RLamNygjsehVWr4a8ct8pp-gWv71c/edit?tab=t.0

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


