import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
import math_physics as mphy

# __all__ = ['NeiSubDir','file_search',
#            'nei_get_arrangement','read_average_tifs','get_tomo_files','nei_determine_murhos',
#            'get_beam_files','nei','beam_near_edge_imaging', 'nei_beam_parameters',]

class NeiSubDir:
    def __init__(self, path, After=False, EdgeB=False):
        self.DarkBefore = path + r'/DarkBefore'
        self.FlatBefore = path + r'/FlatBefore'
        self.EdgeABefore = path + r'/EdgeABefore'
        self.Tomo = path + r'/Tomo'
        if After == True:
            self.DarkAfter = path + r'/DarkAfter'
            self.FlatAfter = path + r'/FlatAfter'
            self.EdgeAAfter = path + r'/EdgeAAfter'
        else:  # use 'Before' for all 'After'
            self.DarkAfter = path + r'/DarkBefore'
            self.FlatAfter = path + r'/FlatBefore'
            self.EdgeAAfter = path + r'/EdgeABefore'
        if EdgeB == True:
            self.EdgeBBefore = path + r'/EdgeBBefore'
            if After == True:
                self.EdgeBAfter = path + r'/EdgeBAfter'
            else:
                self.EdgeBAfter = path + r'/EdgeBBefore'
        else:  # use EdgeA for all EdgeB
            self.EdgeBBefore = path + r'/EdgeABefore'
            if After == True:
                self.EdgeBAfter = path + r'/EdgeAAfter'
            else:
                self.EdgeBAfter = path + r'/EdgeABefore'


def file_search(path, file_filter='*'):
    """
    :param path: search files in this path
    :param file_filter: change the pattern of filter to match the wanted files
    :return:
    """
    import os
    import fnmatch
    files = fnmatch.filter(os.listdir(path), file_filter)
    if len(files) == 0:
        print('------------------------'
              '\nWarning: No files found\n'
              '------------------------')
        return files
    files = [path + r'/' + file for file in files]
    return files


def nei_get_arrangement(setup_type, path):
    class get_arrangement:
        def __init__(self, path):
            filename = path + r'\arrangement.dat'
            try:
                data = pd.read_csv(filename, index_col=0,
                                   header=None, sep=r',\s+', engine='python').T
                for i in range(data.shape[1]):  # remove the ' sign in some strings
                    data.iloc[0, i] = data.iloc[0, i].replace("'", "")

                k = int(data.loc[1, 'k'])
                l = int(data.loc[1, 'l'])
                h = int(data.loc[1, 'h'])
                energy_range_low = float(data.loc[1, 'energy_range_low'])
                energy_range_high = float(data.loc[1, 'energy_range_high'])

                self.diffaction_plane = data.loc[1, 'diffraction_plane']
                self.type = data.loc[1, 'type']
                self.chi_degrees = float(data.loc[1, 'chi_degrees'])
                self.hkl = [h, k, l]
                self.energy = float(data.loc[1, 'energy'])
                self.energy_range = [energy_range_low, energy_range_high]
                self.dist_fd = float(data.loc[1, 'dist_fd'])
                self.detector = self.detector(data)
            except:
                print('Error: Something is wrong when reading in the arrangement file')

        class detector:
            def __init__(self, data):
                self.type = data.loc[1, 'det_type']
                self.pixel = float(data.loc[1, 'det_pixel'])
                self.flip = float(data.loc[1, 'det_flip'])
                self.phperapu = float(data.loc[1, 'det_phperapu'])
                self.disp_x_demag = float(data.loc[1, 'det_disp_x_demag'])
                self.pct_max = float(data.loc[1, 'det_pct_max'])

    arrangement = get_arrangement(path)
    arrangement_parameters = {'diffaction_plane': ' DIFFRACTION PLANE:',
                              'type': ' TYPE:',
                              'chi_degrees': ' ASYMETTRY ANGLE (CHI):',
                              'hkl': ' HKL:',
                              'energy': ' ENERGY (keV):',
                              'energy_range': ' ENERGY RANGE (keV):',
                              'dist_fd': ' DISTANCE FOCUS-DETECTOR (mm):',
                              'detector': '\n DETECTOR PARAMETERS:'}
    detector_parameters = {'type': ' DETECTOR TYPE:',
                           'pixel': ' DETECTOR PIXEL (mm):',
                           'flip': ' DETECTOR FLIP:',
                           'disp_x_demag': ' DETECTOR DISPLAY X DEMAGNIFICATION:',
                           'pct_max': ' DETECTOR PERCENT THRESHOLD:',
                           'phperapu': ' phperapu'}

    for name, value in vars(arrangement).items():
        if name == 'detector':  # print detector information in a seperate section
            print(arrangement_parameters[name])
            for name_2, value_2 in vars(arrangement.detector).items():
                print(detector_parameters[name_2], value_2)
        else:
            print(arrangement_parameters[name], value)

    return arrangement


def read_average_tifs(files,flip=False,xlow=0,xhigh=0,
                      rotate_90=False,twelve_bit=0):
    # read all the image files into a 3D array[n_images,rows,columns],
    # take the average along the images, so that we get an average image
    # Returned value is 2d array
    from PIL import Image

    # twelve_bit=
    n_files = len(files)
    image_array = []
    for i in range(n_files):
        image_array.append(np.array(Image.open(files[i])))
    image_array = np.array(image_array)
    average = image_array.mean(axis=0)
    return(average)

def get_beam_files(path,Verbose=False,clip=False, flip=False):
    '''
    return raw_images, only averaged, nothing else.
    :param path:
    :param Verbose:
    :param clip:
    :param flip:
    :return:
    '''

    sub_dir = NeiSubDir(path, After=False, EdgeB=False)

    flat_path = sub_dir.FlatBefore
    dark_path = sub_dir.DarkBefore
    edge_path = sub_dir.EdgeABefore

    flat_files = file_search(flat_path,'*.tif')
    dark_files = file_search(dark_path,'*.tif')
    edge_files = file_search(edge_path,'*.tif')

    # read and average the raw data files - flats, darks, edges
    # 2d numpy.array(image) is returned by read_average_tifs
    flat = read_average_tifs(flat_files)
    dark = read_average_tifs(dark_files)
    edge = read_average_tifs(edge_files)

    image_size = flat.shape
    n_horizontal = image_size[1] # Number of horizontal positions(pixels)
    n_vertical   = image_size[0] # Number of energy points (the vertical direction on detector)

    class OriginBeamFiles:
        def __init__(self,flat,dark,edge,n_horizontal,n_vertical):
            self.flat        = flat
            self.dark        = dark
            self.edge        = edge
            self.n_horizontal= n_horizontal
            self.n_vertical  = n_vertical
    origin_beam_files = OriginBeamFiles(flat,dark,edge,n_horizontal,n_vertical)

    horizontal_low = 0; horizontal_high= n_horizontal  # Or n_horizontal-1?????

    if clip: # keyword clip is used only to trim off the left and right black area out of the beam,
             # This is usually not necessary, because the beam is supposed to fill the image all the
             # way horizontally.
        t = flat-dark
        t_smooth = median_filter(t,5) # median smooth
        t_smooth = t_smooth/(t_smooth.max())
        t_average= t_smooth.mean(axis=0) # each element is the average of the vertical direction
        bright_ratio  = t_average/(t_average.max())
        # find where(bright_deriv >= pct_max/100.0)
        bright_deriv  = np.gradient(bright_ratio) # the derivative for every element.
                                                  # Considering the brightness center in the beam,
                                                  # it should be positive on one side, and negative
                                                  # on the other side
        fn = bright_deriv/bright_ratio
        x = np.arange(0.0,float(n_horizontal))
        # plt.ion()
        if Verbose==True:
            plt.plot(bright_ratio)
            plt.plot(bright_deriv)
            plt.plot(fn)
            plt.legend(['Flat-Dark','Derivative(Flat-Dark)','Deriv/(Flat-Dark)'])
            plt.title('Brightness along horizontal axis, and its derivative')
            plt.show()

        idx_left = (fn>0); idx_right = (fn<0)
        fwhm_left,left_low,left_high = mphy.fwhm(x[idx_left],fn[idx_left],Verbose=Verbose)
        fwhm_right,right_low,right_high=mphy.fwhm(x[idx_right],-fn[idx_right],Verbose=Verbose)
        horizontal_low=left_high
        horizontal_high=right_low
        flat = flat[:,horizontal_low:horizontal_high]
        dark = dark[:,horizontal_low:horizontal_high]
        edge = edge[:,horizontal_low:horizontal_high]

    class BeamFiles:
        def __init__(self,flat,dark,edge,horizontal_low,
                     horizontal_high,origin_beam_files):
            self.flat = flat
            self.dark = dark
            self.edge = edge
            self.horizontal_low = horizontal_low
            self.horizontal_high=horizontal_high
            self.origin_beam_files= origin_beam_files
    return BeamFiles(flat,dark,edge,horizontal_low,horizontal_high,origin_beam_files)


def get_tomo_files(path, Verbose=False, After=False, EdgeB=False):
    sub_dir = NeiSubDir(path, After=After, EdgeB=EdgeB)
    tomo_files = file_search(sub_dir.Tomo, '*tif')
    return (tomo_files)


def nei_determine_murhos(material_datasource, exy, gaussian_energy_width, interpol_kind='linear',
                         use_sm_data=False, use_measured_standard=False):
    '''
    For every compound, every horizontal position, get the murho value for that
    compound at every energy point (y position on the detector)
    ; mats structure
    ; tags:
    ; names - what you call each material, i.e. 'K2SeO4'
    ; datasource - how we find the mu/rhos
    ;       - system   : get it from calculation
    ;       - FILE : get it from a file
            - standard: get murho from experiment standard. Stardard data are collected with
                        current experiment setting, and standard selenium compound solution ,etc.
    ; i.e. mats = { names:names, datasource:datasource }
    ; where: names = ['Water', 'Bone', 'Selenite', 'U'   ]
    ;        types = [  'system',  'system',     'FILE', 'system' ]
    :param material_datasource:
    :param exy:
    :param gaussian_energy_width:
    :param interpol_kind: default value is 'linear'
    :param use_sm_data:
    :param use_measured_standard:
    :return:
    '''
    nx = exy.shape[1]
    ny = exy.shape[0]  # number of energies
    emin = np.median(exy, axis=1).min()
    emax = np.median(exy, axis=1).max()
    e_range = emax - emin
    energies = np.linspace(emin - 2 * (e_range), emax + 2 * (e_range), 5 * ny)

    murhos = {}
    for name, datasource in material_datasource.items():
        print('(nei_determine_murhos) Geting murho data for ' + name)
        if datasource.lower() == 'file':
            # get murho from saved file
            mu_rho = mphy.murho_selenium_compounds(name, energies)
        elif datasource.lower() == 'system':
            # get murho by calculating it for every element
            mu_rho = mphy.murho(name, energies)
        elif datasource.lower() == 'standard':
            # get murho from experiment standard. Stardard data are collected with current experiment setting,
            # and standard selenium compound solution ,etc.
            pass
        else:
            raise Exception('Material murho data source code is invalid.\n '
                            'Please choose from ["SYSTEM", "FILE", "STANDARD"],\n'
                            'and redefine the "source" variable\n')
        murhos[name] = mu_rho

    ####################  Blur the edge if needed  ##################
    if gaussian_energy_width != 0:
        energies_stepsize = energies[1] - energies[0]
        width = gaussian_energy_width / energies_stepsize
        ###### auto gaussian blur
        for name, value in murhos.items():
            # truncate 3.0 means 3 sigma will be the width used for gaussian filter
            murhos[name] = gaussian_filter(value, sigma=width, mode='nearest', truncate=3.0)
    #     ###### manual gaussian blur
    #     width_int = round(width)
    #     sigma_3X = 3*width_int
    #     sigma_3X = max(sigma_3X,5)
    #     region_3sigma = np.linspace(-sigma_3X,sigma_3X,int(2*sigma_3X)+1)
    #     gauss_filter = np.exp((-region_3sigma**2)/(2*width**2))
    #     total_gf = gauss_filter.sum()
    #     murho1 = murhos.copy()
    #     for name,value in murho1.items():
    #         murho1[name] = np.convolve(value,gauss_filter,'same')/total_gf

    ############## Use interpol to get the murho value at every energy point in exy for every element  ####
    murhos_all = {}
    print('(nei_determine_murhos) Started doing interpolation to bundle MU_RHOS and EXY')
    for name in murhos.keys():  # literate over every compound, save each 2d array in dictionary
        # print('                       Started interpolation for ' + name)
        # print('                       ',end='')
        interpol = interp1d(energies, murhos[name], kind=interpol_kind)  # Get the interpol object for this compound
        murhos_exy = np.empty((ny, nx)) # to save the murhos for one compound
        for x in range(nx):
            murhos_exy[:, x] = interpol(exy[:, x])  # do the interpol prediction
            # print('#'*(x%int(nx/50)==0),end='')
            # if x==nx-1: print()
        murhos_all[name] = murhos_exy  # save to dict
        print('                       Finished interpolation for ' + name)
    print('(nei_determine_murhos) Finished "nei_determine_murhos"\n')

    return murhos_all

def beam_edges(flat_dark,threshold,no_fit=False,Verbose=False,poly_deg=5):
    '''
    :param flat_dark: flat-dark, 2D array
    :param threshold: determines where we cut the beam on the energy axis
    :param no_fit:
    :param Verbose:
    :param poly_deg:
    :return:
    '''
    shape = flat_dark.shape
    nx = shape[1]; ny = shape[0]
    top_positions=[]; bot_positions=[];peak_positions=[]
    for ind_x in range(nx): # loop through x axis
        spectrum = flat_dark[:,ind_x]
        spectrum = median_filter(spectrum,5)
        y_index  = np.arange(len(spectrum))
        peak = spectrum.max()
        peak_ind = spectrum.argmax()
        # Threshold 1: [-sigma,sigma] of normal distribution
        t1 = np.exp(-0.5)
        sigma_ind = y_index[spectrum>(t1*peak)] #The index for all values > t1
        sigma_width=0.5*(sigma_ind[-1]-sigma_ind[0])

        estimate = [peak,peak_ind,sigma_width] #estimate for the gaussian function parameters
                                               #This is also the plan B

        # Call Gaussfit to calculate center and sigma, then the width
        spectrum_gauss,gauss_popt = mphy.gaussfit(y_index,spectrum,estimate)
        gauss_center = gauss_popt[1]
        gauss_sigma  = gauss_popt[2]
        y_peak = int(round(gauss_center))
        if y_peak not in y_index:
            raise ValueError('The peak index is not in the reasonable range')
        half_width = int(np.rint(gauss_sigma*np.sqrt(-2*np.log(threshold))))
        y_top = round(y_peak+half_width)
        if y_top>y_index.max():
            raise ValueError('Impossible! The position of beam top is higher than the beam.')
        y_bot = round(y_peak-half_width)
        if y_bot<y_index.min():
            raise ValueError('Impossible! The position of beam bottom is lower than the beam.')
        top_positions.append(y_top)
        bot_positions.append(y_bot)
        peak_positions.append(y_peak)

    peak_positions= median_filter(peak_positions,5)
    peak_positions= np.rint(np.array(peak_positions)).astype(int)
    top_positions = (peak_positions+half_width).astype(int)
    bot_positions = (peak_positions-half_width).astype(int)

    beam  = flat_dark*0 # Create zero array with the same shape of flat_dark

    class BeamEdgesClass:
        def __init__(self,top,bot,peak,beam):
            self.top=top
            self.bot=bot
            self.peak=peak
            self.beam=beam

    if no_fit: # return with no polynomial fit
        for i in range(nx):
            beam[bot_positions[i]:top_positions[i], i] = 1.0
        s = BeamEdgesClass(top_positions, bot_positions, peak_positions, beam)
        return(s)

    ###############  Doing polynomial fit   ##############################
    #Approach 1 : polynomial fit
    # Only for the peak
    peak_poly_coef = np.polyfit(np.arange(nx), peak_positions,deg=poly_deg)
    p = np.poly1d(peak_poly_coef)
    peak_positions_poly = p(np.arange(nx))

    peak_positions_poly=np.rint(np.array(peak_positions_poly)).astype(int) #round up the positions to integer
    top_positions = (peak_positions_poly+half_width).astype(int)
    bot_positions = (peak_positions_poly-half_width).astype(int)

    for i in range(nx): # Light it up between top and bottom
        beam[bot_positions[i]:top_positions[i]+1,i]=1.0 # Note the +1 to include the top position line

    s = BeamEdgesClass(top_positions, bot_positions, peak_positions_poly, beam)

    if Verbose:
        plt.scatter(np.arange(len(peak_positions)),peak_positions,s=2,alpha=0.2,label='Before polynomial')
        plt.plot(peak_positions_poly,color='r',label='Poly degree=5')
        plt.title('peak_position from beam_edge.pro')
        plt.legend()
        plt.figure()
        plt.imshow(beam)
        plt.title('Selected Beam Area')
        plt.show()

    ########### poly degress study ####################################
    # peak_poly_coef = np.polyfit(np.arange(nx), peak_positions,deg=4)
    # p = np.poly1d(peak_poly_coef)
    # peak_positions_poly4 = p(np.arange(nx))
    # if Verbose:
    #     plt.plot(peak_positions_poly4,color='g',label='Poly degree=4')
    #
    # peak_poly_coef = np.polyfit(np.arange(nx), peak_positions,deg=6)
    # p = np.poly1d(peak_poly_coef)
    # peak_positions_poly6 = p(np.arange(nx))
    # if Verbose:
    #     plt.plot(peak_positions_poly6,color='y',label='Poly degree=6')
    #     plt.legend()
    #     plt.show()
    ###################################################################
    #Approach 2:
    ################################################################
    # Linear regression only gives a straight line
    # machine learning linear regression
    # from sklearn.linear_model import LinearRegression
    # linear_regression = LinearRegression().fit(pd.DataFrame(np.arange(nx)),peak_positions)
    # peak_positions_ml = linear_regression.predict(pd.DataFrame(np.arange(nx)))
    # if Verbose:
    #     plt.plot(peak_positions_ml,color='g',label='After linear regression',alpha=0.3)
    #     plt.legend()
    #     plt.show()
    ###############################################################

    #################################################################
    # # Tried RandomForest, easily overfitted.
    # from sklearn.ensemble import RandomForestRegressor
    # import pandas as pd
    # rf = RandomForestRegressor(n_estimators=100,random_state=42)
    # rf = rf.fit(pd.DataFrame(np.arange(nx)),peak_positions)
    # peak_rf_positions = rf.predict(pd.DataFrame(np.arange(nx)))
    # if Verbose:
    #     plt.plot(peak_rf_positions,color='g',label='After RandomForest',alpha=0.3)
    #     plt.legend()
    #     plt.show()
    #################################################################

    return s






