def compute_toomre_q(sph, gas_particles, star_mass, gamma=1.0):
    """
    Compute the Toomre Q parameter for the disk. 
    Defines stability of the disk (stable for Q > 1)
    """
    G = constants.G  # Gravitational constant
    Q_values = []

    # Get positions, velocities and densities using hydro state from SPH
    for gas in gas_particles:
        # Get position and density of the gas particle
        r = gas.position.length()  # Radial distance from the center
        
        # Extract velocity components (use these for the calculation)
        vx, vy, vz = gas.velocity.x, gas.velocity.y, gas.velocity.z
        
        # Use sph to get the density, sound speed, etc., at the particle position
        x = gas.position.x
        y = gas.position.y
        z = gas.position.z
        vx = vx
        vy = vy
        vz = vz

        # Retrieve the hydro state at the particle's position
        rho, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(
            x, y, z, vx, vy, vz)
        
        # Sound speed (assuming an isothermal equation of state)
        cs = np.sqrt(gamma * rhoe / rho)  # Assuming an ideal gas equation of state
        
        # Epicyclic frequency (approximation for circular orbit)
        kappa = np.sqrt(G * star_mass / r**3)
        
        # Surface density (assuming disk is thin)
        sigma = rho * gas.h_smooth  # rho * smoothing length to estimate surface density
        
        # Compute Toomre Q
        Q = (cs * kappa) / (np.pi * G * sigma)
        Q_values.append(Q)

    return np.array(Q_values)
