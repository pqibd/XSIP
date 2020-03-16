import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re
from pathlib import Path
import inspect
from scipy.interpolate import interp1d
import scipy.constants as C
import os


def focus_p( r, theta_b, chi,nu=0.20):
    """
    Calculator for Polychromatic Focus of asymmetric bent Laue crystal
    :return:
    """
    theta_b = np.radians(theta_b)
    chi = np.radians(chi)
    denom = 2*np.sin(chi+theta_b)+(1+nu)*np.sin(2*chi)*np.cos(chi+theta_b)
    f_p = r*np.sin(2*theta_b)/denom
    return f_p


def focus_g(f_s,r,theta_b,chi):
    """
    Calculator for Geometric Focus of asymmetric bent Laue crystal
    :param fs: Source to crystal distance (meter)
    :param r: Crystal bending radius (meter)
    :param theta_b: Bragg angle (degree)
    :param chi: Asymmetric angle (degree)
    :return:  Geometric Focus of asymmetric bent Laue crystal (meter)
    """
    theta_b = np.radians(theta_b)
    chi = np.radians(chi)
    denom = 2/r + np.cos(chi+theta_b)/f_s
    f_g = np.cos(chi-theta_b)/denom
    return f_g


def trans_rate(materials, energy):
    """
    Calculator for the total transmission rate of X-ray at certain energy through some material(s).
    Parameters
    ----------
    materials : 2-element tuple or list, or list of 2-element tuple or list. 
                      eg. ('air', 300) for 300cm air
                      eg. [('air', 300),('al',0.01)] for 300cm air and 0.01cm Aluminum.
    energy : float
        Energy in `keV`. 
    Returns
    -------
    The transmission rate.

    """
