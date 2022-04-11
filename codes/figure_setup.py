import numpy as np
import matplotlib.pyplot as plt

def dual_half_sphere(ax=None):
    """
    Plot the dual half sphere.
    """


    theta, phi = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
    r = 1

    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)

    if ax is None:
        fig = plt.figure(figsize=(6, 6))
        ax = plt.axes(projection='3d')
    # Plot the half-sphere
    ax.plot_surface(x, y, z, color='w', alpha=0.5, rstride=1, cstride=1,
                    linewidth=0, antialiased=False)
    # Plot the half-sphere's boundary
    #ax.plot_wireframe(x, y, z, color='k', rstride=1, cstride=1)
    return ax

def setup_fig(xlim=(-2, 2), ylim=(-2, 2), zlim=(-2, 2), xlabel=r'X [GSM, $R_\oplus$]',
              ylabel=r'Y [GSM, $R_\oplus$]', zlabel=r'Z [GSM, $R_\oplus$]',):
    """
    Set up the figure for plotting the field lines with real time data.
    """
    fig = plt.figure(figsize=(6, 6))
    ax = plt.axes(projection='3d')
    ax.axvline(0, ls=':', color='k')
    ax.axhline(0, ls=':', color='k')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    ax.set_aspect('auto')
    dual_half_sphere(ax=ax)
    return ax