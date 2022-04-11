import datetime
import itertools
import multiprocessing as mp
import sys
import time
import importlib

import geopack.geopack as gp
import matplotlib.pyplot as plt
import numpy as np

import get_sw_params as get_sw_params
import figure_setup as figure_setup

importlib.reload(get_sw_params)
importlib.reload(figure_setup)

start = time.time()
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{color}')

today_date = datetime.datetime.today().strftime('%Y-%m-%d')


def trace_lines(*args):
    """
    Trace lines from the field line to the surface of the Earth.
    The code needs two other functions: "setup_fig" and "dual_half_sphere" to function properly.
    """

    param_time = 1640115443.165065
    # Get the parameters for the field line tracing
    # Convert the param time to a datetime object in gmt time
    #param_time_dto_strt = datetime.datetime.utcfromtimestamp(param_time)
    #param_time_dto_end = datetime.datetime.utcfromtimestamp(param_time - 300)
    #param_time_str_strt = param_time_dto_strt.strftime('%Y-%m-%d %H:%M:%S')
    #param_time_str_end = param_time_dto_end.strftime('%Y-%m-%d %H:%M:%S')
    #sw_params = get_sw_params(trange = [param_time_str_strt, param_time_str_end])

    sw_params = {}
    sw_params["p_dyn"] = 1e-5
    sw_params["by_gsm"] = -3
    sw_params["bz_gsm"] = -2

    # Define the parameter for computing the total magnetic field from Tsyganenko model (T04)
    # NOTE: For the Tsyganenko model, the elemnts of 'param' are solar wind dynamic pressure, DST,
    # y-component of IMF, z-component of IMF, and 6 zeros.
    param = [sw_params["p_dyn"], -30, sw_params["by_gsm"], sw_params["bz_gsm"], 0, 0, 0, 0, 0, 0]

    # Compute the unix time and the dipole tilt angle
    #t_unix = datetime.datetime(1970, 1, 1)
    #time_dipole = (df_dsco.index[indx_min] - t_unix).total_seconds()
    ps = gp.recalc(param_time)
    theta = args[0][0]
    phi = args[0][1]
    x_gsm = np.sin(theta) * np.cos(phi)
    y_gsm = np.sin(theta) * np.sin(phi)
    z_gsm = np.cos(theta)
    _, _, _, xx1, yy1, zz1 = gp.trace(x_gsm, y_gsm, z_gsm, dir=-1, rlim=30, r0=.99999, parmod=param,
                                      exname='t96', inname='igrf', maxloop=10000)
    _, _, _, xx2, yy2, zz2 = gp.trace(x_gsm, y_gsm, z_gsm, dir=1, rlim=30, r0=.99999, parmod=param,
                                      exname='t96', inname='igrf', maxloop=10000)
    return xx1, yy1, zz1, xx2, yy2, zz2, ps


def line_trace(
    obs_time=None,
    theta_arr=None,
    phi_arr=None,
    ax=None):
    """
    Create a line trace of the field lines for a given date.
    """

    try:

        if obs_time is None:
            obs_time = 1640115443.165065

        if theta_arr is None:
            theta_arr = np.linspace(0, np.pi, 5)

        if phi_arr is None:
            phi_arr = np.array([0])

        print(f"Code execution started at (UTC):" +
              f"{datetime.datetime.utcfromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')}")

        try:
            plt.close("all")
        except:
            pass

        theta_arr = np.linspace(0, 2 * np.pi, 5)
        phi_arr = np.array([0])  # np.linspace(-11 * np.pi/180, 11*np.pi/180, 2)

        p = mp.Pool()
        input = ((i, j) for i, j in itertools.product(theta_arr, phi_arr))
        res = p.map(trace_lines, input)
        p.close()
        p.join()

        ax=figure_setup.setup_fig()
        for r in res:
            xx1 = r[0]
            yy1 = r[1]
            zz1 = r[2]
            xx2 = r[3]
            yy2 = r[4]
            zz2 = r[5]
            ps = r[6]
            ax.plot(xx1, yy1, zz1)
            ax.plot(xx2, yy2, zz2)

        figure_time = f"{datetime.datetime.utcfromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')}"

        ax.text(0.01, 0.01, f'Figure plotted on {figure_time[0:10]} at {figure_time[11:]} UTC',
                ha='left', va='bottom', transform=ax.transAxes, fontsize=12)
        ax.text(0.99, 0.99, f'Real-time T-96 model', ha='right', va='top', transform=ax.transAxes, 
                fontsize=12)
        ax.text(0.01, 0.99, f'Dipole Tilt: {np.round(np.rad2deg(ps), 2)}$^\circ$', 
                ha='left', va='top', transform=ax.transAxes,
                fontsize=12)
        fig_name = f"../figures/Earths_magnetic_field.png"
        fig_name_gdr = f"../figures/Earths_magnetic_field.png"

        t = int(datetime.datetime.today().replace(tzinfo=datetime.timezone.utc).timestamp())
        fig_name_hist = f"../figures/Earthsmagnetic_field_{t}.png"

        plt.savefig(fig_name, bbox_inches='tight', pad_inches=0.05, format='png', dpi=300)
        #plt.savefig(fig_name_hist, bbox_inches='tight', pad_inches=0.05, format='png', dpi=300)
        #plt.savefig(fig_name_gdr, bbox_inches='tight', pad_inches=0.05, format='png', dpi=300)
        #plt.close("all")

        print(f"Code execution finished at (UTC):" +
              f"{datetime.datetime.utcfromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')}\n")
        print(f"Figure saved at (UTC):{figure_time}\n")
        print(f"Waiting for about 15 minutes before running the code again.\n")
    except:
        print("Error:", sys.exc_info()[0])
        print(f"Code execution finished at (UTC):" +
              f"{datetime.datetime.utcfromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')}\n")
        print(f"Waiting for about 15 minutes before running the code again.\n")
        #s.enter(1000, 1, line_trace, (sd,))

    return ax


axs = line_trace()