N = 20000
time = 0 |units.yr
tend = 50. | units.yr
Mstar = 15. | units.MSun
Q_history = []

convert = nbody_system.nbody_to_si(Mstar, 1. | units.AU)
proto = ProtoPlanetaryDisk(
    N, convert_nbody=convert, densitypower=1.5, Rmin=5, Rmax=20, q_out=1.)
gas = proto.result
gas.h_smooth = 0.06 | units.AU

sun = Particles(1)
sun.mass = Mstar
sun.radius = 2. | units.AU
sun.x = 0. | units.AU
sun.y = 0. | units.AU
sun.z = 0. | units.AU
sun.vx = 0. | units.kms
sun.vy = 0. | units.kms
sun.vz = 0. | units.kms

sph = Fi(convert)
sph.parameters.use_hydro_flag = True
sph.parameters.radiation_flag = False
sph.parameters.self_gravity_flag = True
sph.parameters.gamma = 1.
sph.parameters.isothermal_flag = True
sph.parameters.integrate_entropy_flag = False
sph.parameters.timestep = 0.125 | units.yr

sph.gas_particles.add_particles(gas)
sph.particles.add_particles(sun)

supernova_time = 10. |units.yr
supernova_triggered = False

while time < tend:
    if time >= supernova_time and not supernova_triggered:
        print(f"Triggering supernova at {time}")
        
        ejecta_mass = 10 | units.MSun  # Total ejecta mass
        inject_supernova_energy(gas, exploding_region=10 | units.AU)
        sun.mass -= ejecta_mass  # Mass lost from the star
        supernova_triggered = True

    Q_values = compute_toomre_q(sph, gas, sun.mass)  
    Q_history.append(np.mean(Q_values))
    
    sph.evolve_model(time)
    print(f"Disk evolved to {time}")
    time += 1. | units.yr

    L = 50
    rhod, rho = make_map(sph, N=200, L=L)
    plt.figure(figsize=(8, 8))
    plt.imshow(numpy.log10(1.e-5 + rho.value_in(units.amu / units.cm**3)),
                  extent=[-L / 2, L / 2, -L / 2, L / 2], vmin=10, vmax=15)
    plt.title(f'Time: {time}')
    plt.xlabel('AU')
    plt.savefig(f'{time}.png')
