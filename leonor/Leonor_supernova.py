def inject_supernova_energy(gas_particles, 
                            explosion_energy=1.0e+51 | units.erg,
                            exploding_region=10 | units.RSun):
    inner = gas_particles.select(
        lambda pos: pos.length_squared() < exploding_region**2,
        ["position"])
    print(len(inner), "innermost particles selected.")
    print("Adding", explosion_energy / inner.total_mass(), "of supernova " \
        "(specific internal) energy to each of the n=", len(inner), "SPH particles.")
    inner.u += explosion_energy / inner.total_mass()

# Make disk 1 ------

N = 1000
Mstar = 1. | units.MSun

convert = nbody_system.nbody_to_si(Mstar, 1. | units.AU)
proto = ProtoPlanetaryDisk(
    N, convert_nbody=convert, densitypower=1.5, Rmin=1, Rmax=5, q_out=1.)
gas = proto.result
gas.h_smooth = 0.06 | units.AU

sun = Particles(1)
sun.mass = Mstar
sun.radius = 2. | units.AU
sun.x = -7. | units.AU
sun.y = 0. | units.AU
sun.z = 0. | units.AU
sun.vx = 0. | units.kms
sun.vy = 0. | units.kms
sun.vz = 0. | units.kms

gas.position += sun.position
gas.velocity += sun.velocity

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


# Replace Disk 2 with a Supernova Explosion ------

# Generate a spherical distribution of particles for the supernova
from amuse.ic.plummer import new_plummer_model

N_sn = 1000  # Number of particles in the supernova
M_sn = 10. | units.MSun  # Supernova mass

supernova = new_plummer_model(N_sn, convert_nbody=nbody_system.nbody_to_si(M_sn, 1. | units.AU))

# Position the supernova near the disk
supernova.x += 10. | units.AU
supernova.y += 0. | units.AU
supernova.z += 0. | units.AU
gas_sn = proto.result
gas_sn.h_smooth = 0.06 | units.AU
gas.position += supernova.position

# Add the supernova particles to an SPH solver
sph_sn = Fi(convert)

sph_sn.parameters.use_hydro_flag = True
sph_sn.parameters.radiation_flag = False
sph_sn.parameters.self_gravity_flag = True
sph_sn.parameters.gamma = 1.
sph_sn.parameters.isothermal_flag = True
sph_sn.parameters.integrate_entropy_flag = False
sph_sn.parameters.timestep = 0.125 | units.yr

sph_sn.gas_particles.add_particles(gas_sn)
sph_sn.particles.add_particles(supernova)

# Inject supernova energy into the central region
inject_supernova_energy(sph_sn.particles, 
                        explosion_energy=1.0e+51 | units.erg, 
                        exploding_region=1 | units.RSun)

# Bridge the disk and supernova
hydro = bridge.Bridge()
hydro.add_system(sph, (sph_sn,))
hydro.add_system(sph_sn, (sph,))

# Update the evolution loop
t = 0 | units.yr
tend = 20. | units.yr
while t < tend:

    hydro.evolve_model(t)
    print(f"t: {t.in_(tend.unit)}")
    t += sph_sn.parameters.timestep
    
    plt.figure()
    l = 50
    plt.xlim(-l, l)
    plt.ylim(-l, l)
    scatter(sph.gas_particles.x.in_(units.au), sph.gas_particles.y.in_(units.au), c='blue', alpha=0.4, label="Disk")
    scatter(sun.x.in_(units.au), sun.y.in_(units.au), marker='*', c='yellow')

    scatter(sph_sn.gas_particles.x.in_(units.au), sph_sn.gas_particles.y.in_(units.au), c='red', alpha=0.4, label="Supernova")
    plt.xlabel('AU')
    plt.legend(loc='upper right')
    plt.title(f"{t}")
    plt.savefig(f'supernova_disk_{t.number:.3f}.png')
    plt.close()
print("Done")
