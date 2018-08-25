from near_edge_imaging import *
import math
import scipy.constants as C


def nei_beam_parameters(display, beam_files, setup, detector, fix_vertical_motion,
                        clip, no_fit=False, Verbose=False, poly_degree=5):
    '''
    Unit of exy in the returned parameters: keV
    :param display:
    :param sub_dir:
    :param setup:
    :param detector:
    :param fix_vertical_motion:
    :param clip:
    :param no_fit:
    :param Verbose:
    :param poly_degree:
    :return:
    '''
    hkl = setup.hkl
    chi_degrees = setup.chi_degrees
    energy = setup.energy  # edge energy
    dist_fd = setup.dist_fd
    diffraction_plane = setup.diffaction_plane
    e_range = setup.energy_range

    flip = detector.flip
    pixel = detector.pixel
    det_type = detector.type
    phperapu = detector.phperapu
    dispxd = detector.disp_x_demag
    pct_max = detector.pct_max

    ############### physics in the crystal #################
    chi = math.radians(chi_degrees)
    a0 = 5.4305  # Angstroms ##silicon crystal unit cell length at 300K.
    # This is usually used as the internal standard for silicon
    e_edge = energy
    d_hkl = a0 / math.sqrt((np.array(hkl) ** 2).sum())
    # C.c: speed of light, C.h: planck constant, C.eV: eV to Joer
    lamb = (C.c * C.h / C.eV) / (e_edge * 1000) * (10 ** 10)  # WaveLength, Unit: Angstroms
    theta_b = math.asin(lamb / (2 * d_hkl))  # lambda = 2dsin(theta)
    print('(nei_beam_parameters) Bragg angle in degree:\n'
          '                     ', theta_b * 180 / math.pi)

    # What is this?
    bfact = math.cos(theta_b + chi) / math.cos(theta_b - chi)

    #######           Get beam_files       ################################
    flat = beam_files.flat
    dark = beam_files.dark
    edge = beam_files.edge
    beam_shape = flat.shape
    ny = beam_shape[0]
    nx = beam_shape[1]
    x_range = np.arange(nx)
    y_range = np.arange(ny)

    median_width = 5

    edge_dark = edge - dark  # dark corrected edge
    flat_dark = flat - dark  # dark corrected flat

    ##########   Determine the beam size we want to keep  ################################
    thresh = pct_max / 100.0
    print('(nei_beam_parameters) Running "beam_edges"')
    beam_position = beam_edges(flat_dark, thresh, no_fit=no_fit, poly_deg=poly_degree)
    beam = beam_position.beam  # 2D array (image), trimmed top and bottom to remove the darker part vertically
    beam_top = beam_position.top
    beam_bot = beam_position.bot
    beam_peak = beam_position.peak

    # flat_dark_filtered = median_filter(flat_dark, size=median_width)
    # flat_max = flat_dark_filtered.max()

    if diffraction_plane.lower() == 'horizontal':
        raise Exception('"Diffraction plane: horizontal" has '
                        'not been set up for nei_beam_parameters. Come back in future.')

    # # a filter to amplify the center of every spectrum
    # median_flat_dark = np.median(median_filter(flat_dark, 5), axis=1)
    # filter = median_flat_dark / median_flat_dark.max()

    ################ Find absorption edge y_positions ############################
    r = edge_dark / flat_dark
    mu_t = -np.log(r)
    mu_t_median = np.median(mu_t, axis=1)
    deriv_med = abs(np.gradient(mu_t_median))
    mu_t_smooth = median_filter(mu_t, [10, 5])
    deriv_all = abs(np.gradient(mu_t_smooth, axis=0))
    deriv_fwhm = mphy.fwhm(np.arange(ny), deriv_med)
    conv_filter = deriv_med[deriv_fwhm[1]:deriv_fwhm[2]]  # fwhm[1]: left side of fwhm; fwhm[2]: right side of fwhm
    deriv_conv = []
    for x in range(nx):
        deriv_x = deriv_all[:, x]
        conv = np.convolve(deriv_x, conv_filter, 'same')  # Convolution filter used here. See Notes for details.
        deriv_conv.append(conv)
    deriv_conv = np.array(deriv_conv).T
    edge_positions_origin = deriv_conv.argmax(axis=0)
    ###### polynomial fit for edge_positions
    edge_positions = mphy.polyfit(np.arange(nx), median_filter(edge_positions_origin, 3), degree=5)
    edge_positions = np.round(edge_positions).astype(int)

    ##############  flat beam as filter to amplify center#############
    # mu_t = -np.log(r*beam+(1.0-beam))
    # deriv_mut = []
    # for i in range(nx):
    #     top_value = mu_t[beam_top[i], i]
    #     bot_value = mu_t[beam_bot[i], i]
    #     mu_t[0:beam_bot[i], i] = bot_value
    #     mu_t[beam_top[i]:, i] = top_value
    #     deriv = median_filter(np.abs(np.gradient(mu_t[:, i])), 5)
    #     deriv = np.array(deriv) * filter # beam center amplified here
    #     deriv_mut.append(deriv)
    # deriv_mut = np.array(deriv_mut).T
    # deriv_mut = median_filter(deriv_mut, [5, 20])
    # edge_positions = deriv_mut.argmax(axis=0)
    ##################################################################
    if Verbose == True:
        plt.plot(deriv_all[:, 0], label='One spectrum derivative')
        plt.plot(deriv_med, label='Median derivative')
        plt.title('nei_beam_parameters: $\mu t$ & derivatives')
        plt.twinx()
        plt.plot(mu_t[:, 0], color='y', alpha=0.4, label='$\mu t$')
        plt.show()
        plt.figure()
        plt.imshow(mu_t)
        plt.plot(edge_positions, label='edge_positions')
        plt.show()

    ################ Find fwhm and then gaussian with as the edge width##
    ################ To be fixed not sure what is the right way to get the correct real values
    # get edge_widths by calculating fwhm for derivative of each spectrum
    # It could have false values on left and right end when noise is too much
    # fwhms = []
    # for i in range(nx):
    #     fwhms.append(mphy.fwhm(np.arange(ny),deriv_mut[:,i])[0])
    # fwhms = np.array(fwhms)
    # edge_widths = fwhms/(2*np.sqrt(2*np.log(2)))

    # get gaussian_edge_width by calculating the median fwhm for the whole mu_t image
    fw = mphy.fwhm(np.arange(ny), deriv_med)[0]
    edge_width = fw / (2 * np.sqrt(2 * np.log(2)))

    # fit the edge and peak with a 1st order polynomial to get slope of beam
    # for chi and bragg angle corrections.
    # If it is bad, then we cannot change it in software, we have to change the hardware
    # setting on the beamline. So that these variable are only used when we are doing experiment
    edge_slope = mphy.polyfit(x_range, edge_positions, degree=1)  # K-Edge y values
    peak_slope = mphy.polyfit(x_range, beam_peak, degree=1)  # Peak of beam y values

    # get mean values for edge and peak
    edge_mean = edge_positions.mean()
    peak_mean = beam_peak.mean()

    #########   align pixel with energy values      ####################
    '''
        edge_positions-y_index: get the relative position to edge. Use matrix to do it
        exy is the energy(eV) at every [y,x] location
        10**10 is used to line up the unit to the magnitude of Angstrom
    '''
    y_relative = edge_positions - y_range.reshape((ny, 1))
    exy = (C.h * C.c / C.eV) * 10 ** 10 / (2 * d_hkl * np.sin(theta_b + 0.5 * np.arctan(y_relative * pixel / dist_fd)))
    exy = exy / 1000  # change the unit to keV

    ########### If e_range is set, use the set values  #################
    if (type(e_range) == list) & (e_range[0] >= 0) & (e_range[1] > e_range[0]):
        # if e_range is set with reasonable number, then use the e_range to define beam top and bottom
        print('(nei_beam_parameters) Energy range is manually set to: \n'
              '                      ', e_range)
        beam[:, :] = 0.0
        range_index = np.where((exy >= (e_range[0])) & (exy <= (e_range[1])))
        if len(range_index) > 0:
            beam[range_index] = 1.0
            top = []
            bot = []
            for x in x_range:
                beam_inrange = np.where(beam[:, x] > 0)[0]
                top.append(beam_inrange.max())
                bot.append(beam_inrange.min())
            beam_top = np.array(top)
            beam_bot = np.array(bot)
        else:
            raise ValueError('(nei_beam_parameters) The wanted energy range is not available. '
                             'Please change to a reasonable range, or reset range to [0,0] '
                             'to use the whole available energy range')

    ############   calculate gaussian width in terms of Energy  #############
    # energy at the absorption edge
    edge_energies = np.array([exy[edge_positions[i], i] for i in x_range])
    # energy of one pixel away
    edge1_energies = np.array([exy[edge_positions[i] + 1, i] for i in x_range])
    e_per_pixel = abs(edge_energies - edge1_energies).mean()
    e_width = edge_width * e_per_pixel  # gaussian edge width in terms of ENERGY
    # In IDL, we also had the std from gaussian width
    print('(nei_beam_parameters) Gaussian Width measured from Se metal film: ')
    print('                      Energy Width(eV) = ', e_width * 1000)
    print('                      Pixel Width      = ', edge_width)

    #######################  Fix Vertical Motion  ############################
    '''
    ;since beam seems to move vertically, a non-vertical motion affected flat can be created if keyword FIX_VERTICAL_MOTION is set
    ;  this flat is used to I/Io correct the data
    if keyword_set( FIX_VERTICAL_MOTION ) then flt = find_best_average_flat(flat_path, dark)
;'''

    ##################### wrap up things to return  #########################
    class Edges:
        def __init__(self, beam_top, beam_bot, beam_peak, edge_positions):
            self.top = beam_top
            self.bot = beam_bot
            self.peak = beam_peak
            self.edge = edge_positions

    edges = Edges(beam_top, beam_bot, beam_peak, edge_positions)

    class Parameters:
        def __init__(self, beam_files, beam, edges, mu_t, edge_width, edge_slope, peak_slope, exy,
                     e_per_pixel, e_width):
            self.beam_files = beam_files
            self.beam = beam
            self.edges = edges
            self.mu_t = mu_t
            self.pixel_edge_width = edge_width
            self.e_per_pixel = e_per_pixel
            self.e_width = e_width
            self.edge_slope = edge_slope
            self.peak_slope = peak_slope
            self.exy = exy

    beam_parameters = Parameters(beam_files, beam, edges, mu_t, edge_width, edge_slope, peak_slope, exy,
                                 e_per_pixel, e_width)
    print('(nei_beam_parameters) Finished "nei_beam_parameters"')

    return (beam_parameters)


def get_beam_parameters(path='', e_range=0, Verbose=False):
    ##############    get  Path for experiment data file  #######################
    if path == '':
        path = choose_path()
    print("Data directory: ", path)

    #############  get system setup info from arrangement.dat   ##########
    setup = nei_get_arrangement(setup_type='File',path=path)
    detector = setup.detector
    # redefine energy_range if needed
    if e_range != 0: setup.energy_range = e_range

    ########  get beam files: averaged flat, dark, and edge  ############
    beam_files = get_beam_files(path=path, Verbose=Verbose)

    ########  get beam parameters    ####################################
    beam_parameters = nei_beam_parameters(display=False, beam_files=beam_files,
                                          setup=setup, detector=detector,
                                          fix_vertical_motion=False,
                                          clip=False, Verbose=Verbose)
    return beam_parameters
