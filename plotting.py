import imageio.v2 as imageio
import os 
import numpy 
from amuse.units import units
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib import pyplot as plt


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


def animate_frames(name, max, save_as):
    """
    Makes a .gif animation from a list of .png images.
    @ name (str): generator object that lists all of the frames
    @ max (int): number of frames
    @ save_as (str): name of the final .gif
    """
    frames = [f"{name}.png" for i in range(1,max+1)]
    with imageio.get_writer(f"{save_as}.gif", mode='I', duration=0.1) as writer:
        for frame in frames:
            image = imageio.imread(frame)
            writer.append_data(image)
    
    # Optional: Remove individual frame files after the GIF is created
    for frame in frames:
        os.remove(frame)

def animate_2d_plot(x, y, save_as, xlab=None, ylab=None, title=None):
    """
    Animates a standard 2d plot
    @ x, y (arrays) to plot
    @ save_as (str): name of the final .gif
    @ xlab, ylab, title (str)
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_ylabel(xlab)
    ax.set_xlabel(ylab)
    ax.set_title(title)
    plt.legend()
    line, = ax.plot([], [], color="orange")
    

    def update(frame):
        real_line.set_data(x[:frame], y[:frame])
        return real_line,
    
    # Create animation
    anim = FuncAnimation(fig, update, frames=len(x), interval=50, blit=True)
    gif_writer = PillowWriter(fps=30)
    anim.save(f"{save_as}.gif", writer=gif_writer)
