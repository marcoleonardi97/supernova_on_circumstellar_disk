%matplotlib inline
import os
import numpy as np
from matplotlib import pyplot
from amuse.units import units
from amuse.units import units
from amuse.units.constants import G
from amuse.couple import bridge
from amuse.ext.composition_methods import *
from amuse.lab import nbody_system
from amuse.ic.kingmodel import new_king_model
from amuse.community.ph4.interface import Ph4
from amuse.community.bhtree.interface import BHTree
class GalacticCenterGravityCode(object):
    def __init__(self,R, M, alpha):
        self.radius=R
        self.mass=M
        self.alpha=alpha
    def get_gravity_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        m=self.mass*(r/self.radius)**self.alpha
        fr=G*m/r2
        ax=-fr*x/r
        ay=-fr*y/r
        az=-fr*z/r
        return ax,ay,az
    def circular_velocity(self,r):
        m=self.mass*(r/self.radius)**self.alpha
        vc=(G*m/r)
        return vc**0.5
Rgal=1*10**3|units.pc
alpha=1.2
Mgal=1.6*10**10|units.Msun
galactic_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)
galactic_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)
N=100
W0 = 3.0
Mcluster=50000|units.Msun
Rcluster=0.8|units.pc
converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.parsec)
gravity = Ph4(converter)
def make_king_model_cluster(nbodycode, N, W0, Mcluster,
    Rcluster,parameters=[]):
    converter = nbody_system.nbody_to_si(Mcluster, Rcluster)
    bodies = new_king_model(N, W0, convert_nbody=converter)
    code = nbodycode(converter)
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    return code

#cluster_code= make_king_mo
Rgal = 1. | units.kpc
Mgal = 1.6e10 | units.MSun
alpha = 1.2
R_initial=50|units.pc
galaxy_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)
cluster_code = make_king_model_cluster(Ph4, N, W0, Mcluster, Rcluster,
parameters=[("epsilon_squared",
(0.01 | units.parsec)**2)])
stars = cluster_code.particles.copy()
stars.x += R_initial
stars.vy = 0.8*galaxy_code.circular_velocity(R_initial)
channel = stars.new_channel_to(cluster_code.particles)
channel.copy_attributes(["x","y","z","vx","vy","vz"])
system = bridge.Bridge(verbose=False)
system.add_system(cluster_code, (galaxy_code,))
tend=50
timestep=0.25
times = np.arange(0,tend,timestep)
for i,t in enumerate(times):
        system.evolve_model(t|units.Myr,timestep=timestep|units.Myr)
x = system.particles.x.value_in(units.parsec)
y = system.particles.y.value_in(units.parsec)
cluster_code.stop()

from matplotlib import pyplot as plt
plt.scatter(x,y)
plt.xlabel('X (pc)')
plt.ylabel('Y (pc)')
plt.show()
