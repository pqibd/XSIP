from near_edge_imaging import *
from nei_beam_parameters import nei_beam_parameters


def nei(materials='', path='', n_proj=900, algorithm='sKES_equation',
        slice=0, multislice=False, ct=False, side_width=0,
        display=True, e_range=0,lowpass=False,use_torch=True,
        fix_vertical_motion=False,  # maybe change default to True
        clip=False, flip=False, fix_cross_over=False,width_factor=1.0,
        use_sm_data=False, use_measured_standard=False,
        use_weights=False, energy_weights=0, flat_gamma=1.0,
        put_dark_back=False, fix_detector_factor=0, snr=False,
        Verbose=False):
    """
    Get beam_parameters.
    Get $\mu/\rho$ {n_materials:2d-array([n_energies,n_horizontal_positions]),...}
    Get $\mu t$ [n_projections,n_energies,n_horizontal_positions]
    Get $\rho t$ {material: 2d-array([n_projection,n_horizontal_positions]),...}
    :param materials:
    :param path:
    :param ct:
    :param side_width:
    :param n_proj:
    :param slice:
    :param multislice:
    :param e_range:
    :param fix_vertical_motion:
    :param fix_cross_over:
    :param flat_gamma:
    :param Verbose:
    :return: $\rho t$ {material: 2d-array([n_projection,n_horizontal_positions]),...}.
             Values in the dictionary are sinograms in 2d-arrays.
    """

    ###############   define materials       ######################
    if materials == '':
        names = ['K2SeO4', 'K2SeO3', 'Se-Meth', 'Water']
        sources = ['FILE', 'FILE', 'FILE', 'SYSTEM']
        materials = {}
        for i in range(len(names)):
            materials[names[i]] = sources[i]

    ##############   get path for experiment data file  ################
    if path == '':
        path = choose_path()
    print("Data directory: ",path)
    # start counting time
    start = time.clock()
    #############  get system setup info from arrangement.dat file ##########
    setup = nei_get_arrangement(path)
    detector = setup.detector
    # overwrite energy_range from arrangement file if needed
    if e_range != 0: setup.energy_range = e_range

    ########  get beam files: averaged flat, dark, and edge  ############
    print('\n(nei) Running "get_beam_files"')
    beam_files = get_beam_files(path=path, clip=clip, flip=flip, Verbose=Verbose)

    #####################  get tomo data  ########################
    print('\n(nei) Running "get_tomo_files"')
    tomo_data = get_tomo_files(path,multislice=multislice,slice=slice,n_proj=n_proj)

    #################### Get beam_parameters #####################
    print('\n(nei) Running "nei_beam_parameters"')

    beam_parameters = nei_beam_parameters(display=display, beam_files=beam_files,
                                              setup=setup, detector=detector,
                                              fix_vertical_motion=fix_vertical_motion,
                                              clip=clip,Verbose=Verbose)

    ####################  Main calculation  #################################
    '''
    The following is the main calculation for Energy Dispersive Xray Absorption Spectroscopy.
    
    - We get get mu_rho values for every material at every y,x position on the detector (in the image).
    - We calculate the $\mu t$ for every y position (representing energy) at every x position
        (representing horizontal position in the sample) in every tomo image.
    - We calculate the $\rho t$ at every horizontal position for every material.
    In theory, if there is only one material, we can solve the $\rho t$ with the information at one energy
    position, by $(\mu t)/(\mu/\rho)$. When we have 3 materials, we can solve it with 3 energy points.
    In reality, we have sometimes about 900 energy points, so we use linear regression (or other algorithm)
    to solve the coefficient of every material.
    '''
    ##########  get murho values for every material at [y,x] position  ############
    print('\n(nei) Running "nei_determine_murhos"')
    gaussian_energy_width = beam_parameters.e_width * width_factor  # gaussian edge width in terms of energy
    exy = beam_parameters.exy
    mu_rhos = nei_determine_murhos(materials, exy, gaussian_energy_width=gaussian_energy_width,
                                    use_measured_standard=use_measured_standard)

    ####################  calculate -ln(r)=  mu/rho * rho * t   #################
    print('\n(nei) Running "calculate_mut"')
    mu_t = calculate_mut(tomo_data, beam_parameters, lowpass=lowpass,
                         ct=ct, side_width=side_width)

    ####################  Todo: something to reduce artifact   ############

    ####################          calculate rho*t               #################
    beam = beam_parameters.beam
    print('\n(nei) Running "calculate_rhot"')
    rts = calculate_rhot(mu_rhos, mu_t, beam,algorithm=algorithm,use_torch=use_torch)

    ####################   get signal to noise ratio if needed  #################
    snrs=''
    if snr:
        print('(nei) Running "signal_noise_ratio"')
        snrs = signal_noise_ratio(mu_rhos,mu_t,rts,beam_parameters,tomo_data,use_torch)

    ####################   Wrap up results and return  ###########################
    beam_parameters.setup=setup
    class Result:
        def __init__(self, beam_parameters, mu_rhos, mu_t,rts):
            self.beam_parameters = beam_parameters
            self.mu_rhos = mu_rhos
            self.mu_t = mu_t
            self.rts = rts
    print('\n(nei) Total running time for "nei":'
          '\n     ',round(time.clock()-start,2),'seconds')
    return Result(beam_parameters, mu_rhos, mu_t,rts)
