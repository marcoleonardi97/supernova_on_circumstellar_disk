def make_density_map(sph, N=100, L=1):
    """
    Makes density map for animation
    """
    x, y = numpy.indices((N + 1, N + 1))
    x = L * (x.flatten() - N / 2.) / N
    y = L * (y.flatten() - N / 2.) / N
    z = x * 0.
    vx = 0. * x
    vy = 0. * x
    vz = 0. * x

    x = units.AU(x)
    y = units.AU(y)
    z = units.AU(z)
    vx = units.kms(vx)
    vy = units.kms(vy)
    vz = units.kms(vz)

    rhod, rhovx, rhovy, rhovz, rhoe = sph.get_hydro_state_at_point(
        x, y, z, vx, vy, vz)
    rho = rhod.reshape((N + 1, N + 1))

    return rhod, numpy.transpose(rho)
