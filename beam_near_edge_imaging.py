from near_edge_imaging import *
import time


def beam_near_edge_imaging(tomo_data, beam_parameters, materials, lowpass=False,
                           ct=False, side_width=0,flat_gamma=1.0, fix_vertical_motion=False,
                           fix_cross_over=False,  flip=False, width_factor=1,
                           use_sm_data=False, use_weights=False, use_measured_standard=False, *xover):
    '''
    This is the main calculation function.
    - We get get mu_rho values for every material at every y,x position on the detector (in the image).
    - We calculate the $\mu t$ for every y position (representing energy) at every x position
        (representing horizontal position in the sample, in every tomo image.
    - We calculate the $\rho t$ at every horizontal position for every material.
    In theory, if there is only one material, we can solve the $\rho t$ with the information at one energy
    position, by $(\mu t)/(\mu/\rho)$. When we have 3 materials, we can solve it with 3 energy points.
    In reality, we have sometimes about 900 energy points, so we use linear regression (or other algorithm)
    to solve the coefficient of every material.
    :param tomo_data: 3d array, measured values in [n_tomo,ny,nx]
    :param beam_parameters:
    :param materials: dictionary, with material names, and where we get murhos
    :return: mu_rhos, rts
    '''
    flat = beam_parameters.beam_files.flat
    dark = beam_parameters.beam_files.dark
    flat_dark = flat-dark
    nx = flat.shape[1];    ny = flat.shape[0]
    gaussian_energy_width = beam_parameters.e_width * width_factor  # gaussian edge width in terms of energy
    exy = beam_parameters.exy
    beam = beam_parameters.beam

    ##########  get murho values for every material at [y,x] position  ############
    print('(beam_near_edge_imaging) Calling "nei_determine_murhos"')
    mu_rhos = nei_determine_murhos(materials, exy, gaussian_energy_width=gaussian_energy_width,
                                  use_sm_data=use_sm_data, use_measured_standard=use_measured_standard)
    #############   Things in IDL, but have not been done in python   #############
    # Todo: nei_mean_square_murhos  (might not be needed)
    # Todo: use_weights
    # Todo: put_dark_back
    # Todo: fix_crossover
    # Todo: fix_vertical_motion
    # Todo: fix_skes_data_edges
    # Todo: general_nei_vector (might not be needed)

    ####################  calculate -ln(r)=  mu/rho * rho * t   #################
    mu_t = calculate_mut(tomo_data, beam_parameters,lowpass=lowpass,
                         ct=ct,side_width=side_width)

    ####################  Todo: something to reduce artifact   ############


    ####################          calculate rho*t               #################
    rts = calculate_rhot(mu_rhos,mu_t,beam)

    return mu_rhos, rts
    # mu_rhos: {name: [ny,nx],..}, mu_rho values for every material at every y,x position
    # rts:  {name:[n_tomo,nx],...}, rho_t values for every material at every x position in
    #       every tomo file
