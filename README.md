Outline (old):
https://docs.google.com/document/d/1dHtwYFsnluikh4RLamNygjsehVWr4a8ct8pp-gWv71c/edit?tab=t.0

Simulating a supernova explosion in the vicinity of a proto-planetary disk, with and without an additional perturbing object. 

1. Make a protoplanetary disk around a 1 MSun star, orbited by a smaller (< 0.5 MSun) object that inclines and flares the disk. 
2. Explode a ($\approx 30 M_\odot$) star in the vicinity (50 AU) without the perturbing object
3. Explode the same star in the vicinity of the perturbed disk, monitoring the effect on the flares and structures of the disk. 

<ins>Minimum passing grade plan</ins>

Evolve each step as its own system and comment on the results:

1. Simulation (1) should clearly show what kind of perturbations and structures arise from having a perturbing object orbiting the disk.
  - This will include bridging gravity with hydrodynamics, understanding the orbital parameters of the binary, plotting and animating the dynamics of the system over an extended period of time (~ 100 years).
2. Simulation (2) will show the effects of a nearby supernova on a vanilla disk, for reference.
  - This also needs a bridge between gravity and hydrodynamics and good control over the initial conditions (disk face-on, inclination, distance, disk outer density...)
3. Simulation (3) will merge the two above results showing the effects of the SN blast wave on the flares and protuberances of the perturbed disk.
  - Depending on the cost of the calculation we shall run the perturbed disk system again or use the results from simulation (1). This is the most important part of the project, so we will focus on showing the results via plots and animations wherever possible, but also on monitoring important parameters like surface temperature, surface density, Toomre Q parameter ...)




Example animations:
50yr evolution around a $15M_\odot$ star |  Toomre Q Parameter
:-------------------------:|:-------------------------:
![protodisk_with_big_SN](https://github.com/user-attachments/assets/7783ab68-66ac-4596-9d19-d4a923509a1b)  |  ![diskstability](https://github.com/user-attachments/assets/f80d806e-2a47-4e18-8688-a009c60bab45)

Disk in a binary system made with the code in "perturbed_disk.py" over 1 year.
![perturbed_disk_wide](https://github.com/user-attachments/assets/091c945b-b11d-4d12-b5ed-2e863b4afe3e)




