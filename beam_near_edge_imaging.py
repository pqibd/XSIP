from near_edge_imaging import *
import time


def beam_near_edge_imaging(tomo_data, beam_parameters, materials,
                           projection=False, flat_gamma=1.0, fix_vertical_motion=False,
                           fix_cross_over=False, left=0, right=0, flip=False, width_factor=1,
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
    :param materials: dictionary, with material names, and where do we get murhos
    :return: mu_rhos, rts
    '''
    flat = beam_parameters.beam_files.flat
    dark = beam_parameters.beam_files.dark
    flat_dark = flat-dark
    nx = flat.shape[1];    ny = flat.shape[0]
    gaussian_energy_width = beam_parameters.e_width * width_factor  # gaussian edge width in terms of energy
    exy = beam_parameters.exy
    beam = beam_parameters.beam

    ########################  get murho values for every material at [y,x] position  ############
    print('(beam_near_edge_imaging) Calling "nei_determine_murhos"')
    mu_rhos = nei_determine_murhos(materials, exy, gaussian_energy_width=gaussian_energy_width,
                                  use_sm_data=use_sm_data, use_measured_standard=use_measured_standard)
    ########################   Things in IDL, but have not been done in python   ################
    # nei_mean_square_murhos  (might not be needed)
    # use_weights
    # adjustment and entrance for CT recon
    # put_dark_back
    # fix_crossover
    # fix_vetical_motion
    # fix_skes_data_edges
    # general_nei_vector (might not be needed)

    ####################  calculate -ln(r)=  mu/rho * rho * t   #################
    # tomo_data.shape is [n_tomo,ny,nx]
    n_tomo = tomo_data.shape[0]
    mu_t = tomo_data*0.0
    for i in range(mu_t.shape[0]):
        mu_t[i]= -np.log((tomo_data[i]-dark)/flat_dark)

    ####################  something to make mu_t less noisy   ########




    #################### calculate rho*t  ##############################
    # import the regression tool we want to use
    import scipy.optimize.nnls as nnls

    # turn dictionary into array
    names = list(mu_rhos.keys())
    murhos_array = np.array(list(mu_rhos.values()))

    rho_t = np.zeros(shape=(n_tomo,nx,len(names)))
    counter=0
    start_time = time.clock()
    print('(beam_near_edge_imaging) Started calculating RHO_T with linear regression')
    print('                         ',end='')
    for t in range(n_tomo):
        for x in range(nx):
            df = pd.DataFrame(murhos_array[:,:,x]).T
            df.columns=names
            df['Mu_t']= mu_t[t,:,x] #mu_t [n_tomo,ny,nx]
            # Only use the part in the 'beam' range
            beam_range = beam[:,x]>0
            # do linear regression
            coef = nnls(df[beam_range][names],df[beam_range]['Mu_t'])[0]
            rho_t[t,x]=coef
            # print progress
            print('>'*(counter%int(n_tomo*nx/48)==0),end='')
            counter += 1
    print('\n                         Finished calculation for'
          '\n                        ',n_tomo,' tomo files in',
          round(time.clock()-start_time,2),'seconds')

    rts={}
    for i in range(len(names)):
        rts[names[i]]=rho_t[:,:,i]

    return mu_rhos, rts
    # mu_rhos: {name: [ny,nx],..}, mu_rho values for every material at every y,x position
    # rts:  {name:[n_tomo,nx],...}, rho_t values for every material at every x position in
    #       every tomo file
