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


def inject_external_supernova(gas_particles, supernova_position, time, explosion_energy, ejecta_mass, heavy_element_mass, shock_speed):
    """
    Injects energy, momentum, and heavy elements into the gas particles
    as a shockwave from an external supernova impacts them.

    Parameters:
    - gas_particles: Particles object, the gas particles of the disk.
    - supernova_position: VectorQuantity, position of the exploding star.
    - time: ScalarQuantity, current simulation time.
    - explosion_energy: ScalarQuantity, total energy released by the supernova.
    - ejecta_mass: ScalarQuantity, total mass of the ejecta.
    - heavy_element_mass: ScalarQuantity, mass of heavy elements in the ejecta.
    - shock_speed: ScalarQuantity, speed of the shockwave.
    """
    # Calculate distance from the supernova to each particle
    distances = (gas_particles.position - supernova_position).lengths()
    shock_radius = shock_speed * time

    # Identify particles within the shockwave radius
    affected_particles = gas_particles.select(
        lambda r: r < shock_radius, ["position"]
    )
    if len(affected_particles) == 0:
      print("The SN shockwave didn't interact with any gas particle")
      return 
    
    # Inject energy and momentum
    specific_energy = explosion_energy / affected_particles.total_mass()
    specific_momentum = (ejecta_mass * shock_speed / affected_particles.total_mass()).value_in(units.kms)

    direction_to_particles = (affected_particles.position - supernova_position).normalized()
    affected_particles.u += specific_energy
    affected_particles.vx += specific_momentum * direction_to_particles[:, 0]
    affected_particles.vy += specific_momentum * direction_to_particles[:, 1]
    affected_particles.vz += specific_momentum * direction_to_particles[:, 2]

    # Inject heavy elements
    specific_heavy_elements = heavy_element_mass / affected_particles.total_mass()
    affected_particles.metallicity += specific_heavy_elements

    print(f"Injected supernova effects into {len(affected_particles)} particles at time {time}.")
