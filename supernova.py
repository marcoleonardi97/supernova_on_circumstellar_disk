from amuse.units import units

def inject_supernova_energy(gas_particles, heavy_elements_mass,
                            explosion_energy=1.0e+51|units.erg,
                            exploding_region=10|units.RSun):
    inner = gas_particles.select(
        lambda pos: pos.length_squared() < exploding_region**2,
        ["position"])
    print(len(inner), "innermost particles selected.")
    print("Adding", explosion_energy / inner.total_mass(), "of supernova " \
        "(specific internal) energy to each of the n=", len(inner), "SPH particles.")
    inner.u += explosion_energy / inner.total_mass()
    inner.metallicity += heavy_elements_mass / inner.total_mass()
