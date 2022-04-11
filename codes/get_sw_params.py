import datetime

import warnings

import geopack.geopack as gp
import numpy as np
import pyspedas as spd
import pytplot as ptt
from tabulate import tabulate

def get_sw_params(
    probe=None,
    omni_level="hro",
    time_clip=True,
    trange=None,
    mms_probe_num=None,
    verbose=False
    ):
    r"""
    Get the solar wind parameters from the OMNI database.

    Parameters
    ----------
    probe : str
        The probe to use. Default is 'None'.
    omni_level : str
        The omni data level to use. Options are 'hro' and 'hro2'. Default is 'hro'.
    time_clip : bool
        If True, the data will be clipped to the time range specified by trange. Default is True.
    trange : list or an array of length 2
        The time range to use. Should in the format [start, end], where start and end times should
        be a string in the format 'YYYY-MM-DD HH:MM:SS'.
    mms_probe_num : str
        The MMS probe to use. Options are '1', '2', '3' and '4'. Default is None.
    verbose : bool
        If True, print out a few messages and the solar wind parameters. Default is False.

    Raises
    ------
    ValueError: If the probe is not one of the options.
    ValueError: If the trange is not in the correct format.

    Returns
    -------
    sw_params : dict
        The solar wind parameters.
    """

    if trange is None:
        raise ValueError("trange must be specified as a list of start and end times in the format 'YYYY-MM-DD HH:MM:SS'.")

    # Check if trange is either a list or an array of length 2
    if not isinstance(trange, (list, np.ndarray)) or len(trange) != 2:
        raise ValueError(
            "trange must be specified as a list or array of length 2 in the format 'YYYY-MM-DD HH:MM:SS.")

    # Download the OMNI data (default level of 'hro_1min') for the specified timerange.
    omni_varnames = ['BX_GSE', 'BY_GSM', 'BZ_GSM', 'proton_density', 'Vx', 'Vy', 'Vz', 'SYM_H']
    omni_vars = spd.omni.data(trange=trange, varnames=omni_varnames, level=omni_level,
                             time_clip=time_clip)

    omni_time = ptt.get_data(omni_vars[0])[0]
    print(f'omni_time: {len(omni_time)}')
    omni_bx_gse = ptt.get_data(omni_vars[0])[1]
    omni_by_gsm = ptt.get_data(omni_vars[1])[1]
    omni_bz_gsm = ptt.get_data(omni_vars[2])[1]
    omni_np = ptt.get_data(omni_vars[3])[1]
    omni_vx = ptt.get_data(omni_vars[4])[1]
    omni_vy = ptt.get_data(omni_vars[5])[1]
    omni_vz = ptt.get_data(omni_vars[6])[1]
    omni_sym_h = ptt.get_data(omni_vars[7])[1]

    # Get mms postion in GSM coordinates for the specified time range

    if (mms_probe_num is not None):
        mms_varnames = [f'mms{mms_probe_num}_mec_r_gsm']
        mms_vars = spd.mms.mec(trange=trange, varnames=mms_varnames, probe=mms_probe_num,
                               data_rate='srvy', level='l2', time_clip=time_clip,
                               latest_version=True)
        mms_time = ptt.get_data(mms_vars[0])[0]
        # Position of MMS in GSM coordinates in earth radii (r_e) units
        r_e = 6378.137  # Earth radius in km
        mms_sc_pos = ptt.get_data(mms_vars[0])[1:3][0]/r_e

        # TODO: Find out why adding 'mms_fgm_varnames' as a variable causes the code to give out no data.
        mms_fgm_varnames = [f'mms{mms_probe_num}_fgm_b_gsm_srvy_l2_bvec']
        mms_fgm_vars = spd.mms.fgm(trange=trange, probe=mms_probe_num, time_clip=time_clip,
                                   latest_version=True)
        mms_fgm_time = ptt.get_data(mms_fgm_varnames[0])[0]
        mms_fgm_b_gsm = ptt.get_data(mms_fgm_varnames[0])[1:4][0]
    else:
        mms_time = None
        mms_sc_pos = None
        mms_fgm_time = None
        mms_fgm_b_gsm = None
        pass

    time_imf = np.nanmedian(omni_time)
    #print(time_imf, type(time_imf))
    b_imf_x = np.nanmedian(omni_bx_gse)
    b_imf_y = np.nanmedian(omni_by_gsm)
    b_imf_z = np.nanmedian(omni_bz_gsm)

    if (b_imf_z > 15 or b_imf_z < -18):
        warnings.warn(
        f"The given parameters produced the z-component of IMF field (b_imf_z) {b_imf_z} nT,"
        f"which is out of range in which model is valid (-18 nT < b_imf_z < 15 nT)"
        )

    time_imf_hrf = datetime.datetime.utcfromtimestamp(time_imf)
    np_imf = np.nanmedian(omni_np)
    vx_imf = np.nanmedian(omni_vx)
    vy_imf = np.nanmedian(omni_vy)
    vz_imf = np.nanmedian(omni_vz)
    sym_h_imf = np.nanmedian(omni_sym_h)
    v_imf = [vx_imf, vy_imf, vz_imf]
    b_imf = [b_imf_x, b_imf_y, b_imf_z]
    imf_clock_angle = np.arctan2(b_imf[1], b_imf[2]) * 180 / np.pi
    if imf_clock_angle < 0:
        imf_clock_angle += 360
    if mms_probe_num is not None:
        mean_mms_sc_pos = np.round(np.nanmean(mms_sc_pos, axis=0), decimals=2)
        mean_mms_fgm_b_gsm = np.round(np.nanmedian(mms_fgm_b_gsm, axis=0), decimals=2)
    else:
        mean_mms_sc_pos = None
        mean_mms_fgm_b_gsm = None

    print("IMF parameters found:")
    if (verbose):
        print(tabulate(
            [["Time of observation (UTC)", time_imf_hrf],
             ["IMF Magnetic field [GSM] (nT)", b_imf],
             ["IMF Proton density (1/cm^-3)", np_imf],
             ["IMF Plasma velocity (km/sec)", v_imf],
             ["IMF clock angle (degrees)", imf_clock_angle],
             ["IMF Sym H", sym_h_imf],
             ["MMS position (GSM) (R_E)", mean_mms_sc_pos]],
            headers=["Parameter", "Value"], tablefmt="fancy_grid", floatfmt=".2f",
            numalign="center"))

    # Check if the values are finite, if not then assign a default value to each of them
    if ~(np.isfinite(np_imf)):
        np_imf = 5
    if ~(np.isfinite(vx_imf)):
        vx_imf = -500
    if ~(np.isfinite(vy_imf)):
        vy_imf = 0
    if ~(np.isfinite(vz_imf)):
        vz_imf = 0
    if ~(np.isfinite(sym_h_imf)):
        sym_h_imf = -1

    m_proton = 1.672e-27  # Mass of proton in SI unit

    rho = np_imf * m_proton * 1.15

    #  Solar wind ram pressure in nPa, including roughly 4% Helium++ contribution
    p_dyn = 1.6726e-6 * 1.15 * np_imf * (vx_imf**2 + vy_imf**2 + vz_imf**2)

    if (p_dyn > 8.5 or p_dyn < 0.5):
        warnings.warn(
            f"The given parameters produced a dynamic pressure of {p_dyn} nPa which is out of"
            f" range in which model is valid (0.5 nPa < p_dyn < 8.5 nPa)",
        )
    param = [p_dyn, sym_h_imf, b_imf_y, b_imf_z, 0, 0, 0, 0, 0, 0]

    # Compute the dipole tilt angle
    ps = gp.recalc(time_imf)

    # Make a dictionary of all the solar wind parameters
    sw_dict = {}
    sw_dict['time'] = time_imf
    sw_dict['b_imf'] = b_imf
    sw_dict['rho'] = rho
    sw_dict['ps'] = ps
    sw_dict['p_dyn'] = p_dyn
    sw_dict['sym_h'] = sym_h_imf
    sw_dict['imf_clock_angle'] = imf_clock_angle
    sw_dict['param'] = param
    sw_dict['mms_time'] = mms_time
    sw_dict['mms_sc_pos'] = mms_sc_pos
    sw_dict['mms_b_gsm'] = mean_mms_fgm_b_gsm

    return sw_dict
