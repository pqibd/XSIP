import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import torch
from scipy.ndimage import median_filter
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
import math_physics as mphy
from toolkit import *

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


def nei_get_arrangement(path,setup_type='FILE'):
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
                raise Exception('Error: Something is wrong when reading in the arrangement file')

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
    # take the average along the images, so that we get an average image.
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
    return averaged flat, dark, edge images in form of 2d-arrays.
    :param path:
    :param Verbose:
    :param clip: Used only to trim off the left and right black area out of the beam.
                 This is usually not necessary, because the beam is supposed to fill
                 the image all the way horizontally
    :param flip:
    :return: beam_files: averaged flat, dark, edge images in form of 2d-arrays.
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
        """
        This is used to keep a copy of beam files with their original shape when clip is true.
        """
        def __init__(self,flat,dark,edge):
            self.flat        = flat
            self.dark        = dark
            self.edge        = edge

    origin_beam_files = OriginBeamFiles(flat,dark,edge)

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


def get_tomo_files(path, multislice=False, slice=0, n_proj=900, Verbose=False, After=False, EdgeB=False):
    """
    Read all the tomo files in tomo folder, or n_proj images for one slice of CT.
    :param path:
    :param multislice: If True, n_proj images will be read in for the wanted single slice.
    :param slice: Which slice do you want? Start from 0.
    :param n_proj: How many projection images are there in one slice
    :param Verbose:
    :param After:
    :param EdgeB:
    :return: 3d-array [n_tomo(n_proj), n_energies, n_horizontal_positions]
    """

    sub_dir = NeiSubDir(path, After=After, EdgeB=EdgeB)
    tomo_files = file_search(sub_dir.Tomo, '*tif')

    # get the tomo images for ONE slice when there are multislices of projections in tomo image folder.
    if multislice == True:
        # print('-----------------------------------------------------'
        #       '\nReminder: "Slice" starts from 1. There is NO 0th slice'
        #       '-----------------------------------------------------')
        # if slice < 1:
        #     print('-----------------------------------------------------'
        #           '\nWarning: "Slice" starts from 1. "slice=' + str(slice) + ' is entered\n'
        #           '-----------------------------------------------------')
        i_begin = n_proj * slice
        i_end = i_begin + n_proj
        tomo_files = tomo_files[i_begin:i_end]
        print('(get_tomo_files) Tomo files in 1 slice loaded')
    n_tomo = len(tomo_files)
    print('(get_tomo_files) Number of Tomo files: ', n_tomo)  # equal to n_projections

    # tomo files to data array
    from PIL import Image
    counter=0
    tomo_data = []
    if n_tomo >= 200:
        print('|--------------------------------------------------|\n|', end='')
    for i in range(n_tomo):
        tomo_data.append(np.array(Image.open(tomo_files[i])))
        if n_tomo >= 200:
            print('>' * (counter % int(n_tomo / 50) == 0), end='')
        counter+=1
    print()
    tomo_data = np.array(tomo_data)
    return (tomo_data)


def nei_determine_murhos(materials, exy, gaussian_energy_width, interpol_kind='linear',
                         use_file=True,use_sm_data=False, use_measured_standard=False):
    '''
    For every compound, every horizontal position, get the murho value for that
    compound at every energy point (y position on the detector). Ways to get murho values
    are defined by "source".
    materials structure:
        - names: what you call each material, i.e. 'K2SeO4'
        - sources: how we find the mu/rhos for that element or compound
            - system   : get it from calculation with known chemical information
            - FILE : get it from a file
            - standard: get murho from experiment standard. Stardard data are collected
                         with current experiment setting, and standard selenium compound solution ,etc.
    i.e. materials = { names:names, sources:sources }
    where: names = ['Water', 'Bone', 'Selenite', 'U'   ]
           sources = [  'system',  'system',     'FILE', 'system' ]
    :param materials:
    :param exy:
    :param gaussian_energy_width:
    :param interpol_kind: default value is 'linear'. See 'scipy.interpolate.interp1d' for other available
                          "interpol_kind"s.
    :param use_sm_data:
    :param use_measured_standard:
    :return:
    '''
    nx = exy.shape[1]
    ny = exy.shape[0]  # number of energies
    emin = np.median(exy, axis=1).min()
    emax = np.median(exy, axis=1).max()
    e_range = emax - emin
    # create an energies array with the size 5 * e_range
    energies = np.linspace(emin - 2 * (e_range), emax + 2 * (e_range), 5 * ny)

    murhos = {}
    # for name, source in materials.items():
    #     print('(nei_determine_murhos) Getting murho data for ' + name)
        # if source.lower() == 'file':
        #     # get murho from saved file
        #     mu_rho = mphy.murho_selenium_compounds(name, energies)
        # elif source.lower() == 'system':
        #     # get murho by calculating it for every element
        #     mu_rho = mphy.murho(name, energies)
        # elif source.lower() == 'standard':
        #     # Todo:
        #     # get murho from experiment standard. Stardard data are collected with current experiment setting,
        #     # and standard selenium compound solution ,etc.
        #     pass
        # else:
        #     raise Exception('Material murho data source code is invalid.\n '
        #                     'Please choose from ["SYSTEM", "FILE", "STANDARD"],\n'
        #                     'and redefine the "source" variable\n')
    # for name, source in materials.items():
    for name in materials:
        print('(nei_determine_murhos) Getting murho data for ' + name)
        mu_rho = mphy.murho(name, energies, use_file=use_file)
        murhos[name] = mu_rho

    ####################  Blur the edge if needed  ##################
    if gaussian_energy_width != 0:
        energies_stepsize = energies[1] - energies[0]
        width = gaussian_energy_width / energies_stepsize
        ###### auto gaussian blur
        for name, value in murhos.items():
            # truncate 3.0 means 3*sigma defines the width (number of pixels) to apply gaussian filter
            murhos[name] = gaussian_filter(value, sigma=width, mode='nearest', truncate=3.0)
    #     ###### manual gaussian blur
    #     width_int = round(width)
    #     sigma_3X = 3*width_int
    #     sigma_3X = max(sigma_3X,5) # we don't want the sigma_3X too small
    #     region_3sigma = np.linspace(-sigma_3X,sigma_3X,int(2*sigma_3X)+1)
    #     gauss_filter = np.exp((-region_3sigma**2)/(2*width**2))
    #     total_gf = gauss_filter.sum()
    #     murho1 = murhos.copy()
    #     for name,value in murho1.items():
    #         murho1[name] = np.convolve(value,gauss_filter,'same')/total_gf

    ############## Use interpol to get the murho value at every energy point in exy for every element  ####
    murhos_all = {}
    print('(nei_determine_murhos) Started doing interpolation to bundle MU_RHOS and EXY')
    for name in murhos.keys():  # iterate over every compound, save each 2d array in dictionary
        interpol = interp1d(energies, murhos[name], kind=interpol_kind)  # Get the interpol object for this compound
        murhos_exy = np.empty((ny, nx)) # to save the murhos for one compound
        for x in range(nx):
            murhos_exy[:, x] = interpol(exy[:, x])  # do the interpol prediction
        murhos_all[name] = murhos_exy  # save to dict
        print('                       Finished interpolation for ' + name)
    print('(nei_determine_murhos) Finished "nei_determine_murhos"')

    return murhos_all


def beam_edges(flat_dark,threshold,no_fit=False,Verbose=False,poly_deg=5):
    '''
    Locate the peak positions in the beam. Decide the useful region in the beam we are going to keep.
    :param flat_dark: flat-dark, 2D array
    :param threshold: percentage of the max value. Beam values higher than threshold will be kept.
    :param no_fit: No polynomial fit for the peak positions in the beam.
    :param Verbose:
    :param poly_deg: The degree parameter for doing polynomial fit.
    :return: Peak position in flat_dark image.
             Beam region we will use for analysis that follows.
             Top positions of the useful beam region.
             Bottom postions of the useful beam region.
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
            raise ValueError('The peak index is not in reasonable range')
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

    class BeamEdges:
        def __init__(self,top,bot,peak,beam):
            self.top=top
            self.bot=bot
            self.peak=peak
            self.beam=beam

    if no_fit:  # return with no polynomial fit
        for i in range(nx):
            beam[bot_positions[i]:top_positions[i], i] = 1.0
        s = BeamEdges(top_positions, bot_positions, peak_positions, beam)
        return s

    ###############  Doing polynomial fit   ##############################
    #Approach 1 : polynomial fit
    # Only for the peak
    print('(beam_edges) Doing polynomial fit for beam peak positions')
    peak_poly_coef = np.polyfit(np.arange(nx), peak_positions,deg=poly_deg)
    p = np.poly1d(peak_poly_coef)
    peak_positions_poly = p(np.arange(nx))

    peak_positions_poly=np.rint(np.array(peak_positions_poly)).astype(int) #round up the positions to integer
    top_positions = (peak_positions_poly+half_width).astype(int)
    bot_positions = (peak_positions_poly-half_width).astype(int)

    for i in range(nx): # Light it up between top and bottom
        beam[bot_positions[i]:top_positions[i]+1,i]=1.0 # Note the +1 to include the top position line

    s = BeamEdges(top_positions, bot_positions, peak_positions_poly, beam)

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


def idl_recon(sinogram,pixel,center=0):
    """
    CT reconstruction with the "normalized_fbp" function from IDL. Note: Licensed IDL software is required.
    :param sinogram: 3d or 2d-array [n_projections,n_horizontal_positions]
    :param pixel: pixel size in centimeter. For example: 0.0009 for 9um.
    :param center: Default 0. The sample rotation center during CT imaging.
    :return: reconstructed square 2d-array [n_horizontal_positions,n_horizontal_positions]
    """
    from idlpy import IDL
    dimensions = sinogram.ndim
    if dimensions ==3:
        recon=[]
        for i in range(sinogram.shape[0]):
            recon.append(IDL.normalized_fbp(sinogram[i],dx=center,pixel=pixel))
        recon = np.array(recon)
    elif dimensions==2:
        recon = IDL.normalized_fbp(sinogram,dx=center,pixel=pixel)
    else:
        raise Exception('The dimensions of input "sinogram" should be either 2 or 3.'
                        ,dimensions,'dimensions were given.')
    return recon


def skimage_recon(sinogram,n_proj,pixel_size=1.0,output_size=None,filter='ramp',center=0,circle=None):
    """
    CT reconstruction using Inverse Radon Transform, with Filtered Back Projection algorithm.
    See "radon_transform" in skimage.transform for more detail.


    :param sinogram:array_like, dtype=float
                    Image containing radon transform (sinogram).
                    If n_dimension==2, each row of the image corresponds to
                    a projection along a different angle. The tomography
                    rotation axis should lie at the pixel index.
                    If n_dimension>=3, the last two dimensions should be the
                    sinogram.
    :param n_proj: Number of projections taken for one slice of CT imaging
    :param pixel_size: Pixel size of the detector. Unit: cm.
    :param output_size: The width of the output reconstruction image.
                        If Output_size not specified, use the image horizontal width.
    :param filter: str, optional (default ramp)
                    Filter used in frequency domain filtering. Ramp filter used by default.
                    Filters available: ramp, shepp-logan, cosine, hamming, hann.
                    Assign None to use no filter.
    :param center: The rotation center of CT imaging.(default 0)
    :param circle: boolean, optional
                    Assume the reconstructed image is zero outside the inscribed circle.
                    The default behavior (None) is equivalent to False.
    :return: Reconstruction image arrays, with the same number of dimensions of the input
             sinogram array.
    """
    from skimage.transform import iradon
    # if sinogram.ndim>=3, which means there is more than one sinogram. Iterate them all.
    # The last two dimension should be the sinogram array.
    if sinogram.ndim>=3:
        recon = []
        for i in range(sinogram.shape[0]):
            recon.append(skimage_recon(sinogram[i],n_proj=n_proj,pixel_size=pixel_size,
                                       output_size=output_size,filter=filter,center=center,
                                       circle=circle))
        recon = np.array(recon)
        return recon
    # when sinogram.ndim==2, do the reconstruction.
    print('(skimage_recon) Started one CT reconstruction...',end='')
    if not output_size: # If Output_size not specified, use the image horizontal width
        output_size=sinogram.shape[1]
    sinogram = sinogram.transpose(1,0) # transpose row and column to meet the order in skimage.transform.
    theta = np.linspace(0,180,n_proj) # make the angles of projections from 0 to 180 degree
    recon = iradon(sinogram,theta=theta,output_size=output_size,filter=filter,
                   center_drift=center,circle=circle)
    # correct result with pixel size (cm).
    recon = recon/pixel_size
    print('...Finished')

    return recon



def calculate_mut(tomo_data, beam_parameters,lowpass=False,ct=False,side_width=0):
    """
    mu*t =  (mu/rho) * (rho * t) = -ln[(tomo-dark)/(flat-dark)]
    :param tomo_data:
    :param beam_parameters:
    :param lowpass: Use gaussian filter to smooth the mu_t spectrum along the energy axis
    :param ct: If ct, use the left and right of projection to remove air absorption
    :param side_width: the width used for the "ct" parameter
    :return: mu_t: 3d-array [n_tomo, ny, nx]
    """
    ####################  calculate -ln[(tomo-dark)/(flat-dark)]   #################
    # tomo_data.shape is [n_tomo,ny,nx]
    flat = beam_parameters.beam_files.flat
    dark = beam_parameters.beam_files.dark
    flat_dark = flat-dark
    beam = beam_parameters.beam
    nx = beam.shape[1]; ny = beam.shape[0]
    n_tomo = tomo_data.shape[0]

    mu_t = tomo_data*0.0 # make a 3_d array with the same shape of tomo_data, fill with zeros
    print('(calculate_mut) Started calculating MU_T in every tomo image at every [energy,horizontal] position')
    if n_tomo >= 200:
        print('|--------------------------------------------------|\n|', end='')
    counter=0
    for i in range(n_tomo):
        mu_t[i]= -np.log((tomo_data[i]-dark)/flat_dark)
        mu_t[i]= np.nan_to_num(mu_t[i])  #replace nan values with 0.

        if ct: # remove air absorption. Calculated from the left and right of the projection image,
               # where there should be only air, no sample.
            if side_width<=0: # making sure side_width is valid.
                raise Exception('When "CT" is True, "side_width" is used to remove air absorption, and '
                                '"side_width" has to be an integer >0')
            mut_left_total = (mu_t[i]*beam)[:,0:side_width].sum()
            mut_right_total= (mu_t[i]*beam)[:,-side_width:].sum()
            mut_left_count = beam[:,0:side_width].sum()
            mut_right_count= beam[:,-side_width:].sum()
            mut_left_avg   = mut_left_total/mut_left_count
            mut_right_avg  = mut_right_total/mut_right_count
            x_position = np.arange(nx)
            filter_1d = mut_left_avg+(mut_right_avg-mut_left_avg)*((x_position-side_width/2)/(nx-side_width))
            filter_2d = filter_1d*(beam*0.0+1.0)
            # remove air absorption from total mu_t
            mu_t[i] = mu_t[i]-filter_2d
        if n_tomo>=200:
            print('>' * (counter % int(n_tomo / 50) == 0), end='')
        counter += 1
    print('')
        # print("\r                  %d%%" % (round((i/n_tomo)*100)),end='')
    # use a lowpass filter (gaussian filter) to remove some high frequency noise
    # Todo: decide the default value for "lowpass"
    pixel_gaussian_width = beam_parameters.pixel_edge_width
    if lowpass:
        print('(calculate_mut) Applying lowpass filter along energy axis')
        mu_t = gaussian_filter(mu_t,[0,pixel_gaussian_width,0],truncate=3.0,mode='nearest')

    return mu_t


def calculate_rhot(mu_rhos,mu_t,beam,names,algorithm='',use_torch=True):
    """
    calculate the $\rho t$ for every material at every horizontal position in every projection
    :param mu_rhos: mu_rhos is obtained from "nei_determine_murhos". {material: [n_energies,,nx]...}
    :param mu_t: mu_t is obtained from "calculate_mut" .[n_projections,n_energies,nx]
    :param beam: beam is obtained from "beam_parameters.beam". [n_energies,nx]
    :param algorithm: The core algorithm to calculate $\rho t$.
                      Availabe options are ["nnls", "sKES_equation"] for now (Aug 27, 2018).
                      If "nnls": `scipy.optimize.nnls` will be used to perform linear regression
                                  for the spectrum at every horizontal position in every projection image.
                      If "sKES_equation": A equation derived with least-square approach is used.
                                           Because matrix operation is used here for calculation, it is much
                                           faster than doing all the iterations with "nnls".
                                           [Ref: Ying Zhu,2012 Dissertation]
    :return: 3d-array with shape [n_materials, n_projection,nx]}. For CT data, the last two dimensions
             form the sinogram.
    """

    if algorithm=='':
        algorithm=input('Choose algorithm from:  "nnls", "sKES_equation"\n'
                        '(type and enter): ')


    nm = mu_rhos.shape[0] # number of materials

    if algorithm == 'nnls':
        print('(calculate_rhot) Algorithm: "scipy.optimize.nnls"')
        import scipy.optimize.nnls as nnls
        n_tomo = mu_t.shape[0]
        nx = mu_t.shape[2]
        rho_t = np.zeros(shape=(nm,n_tomo, nx))
        counter = 0
        start_time = time.clock()
        print('(calculate_rhot) Started calculating RHO_T with linear regression')
        print('                 ', end='')
        for t in range(n_tomo):
            for x in range(nx):
                df = pd.DataFrame(mu_rhos[:, :, x]).T
                df.columns = names
                df['Mu_t'] = mu_t[t, :, x]  # mu_t [n_tomo,ny,nx]
                # Only use the part in the 'beam' range
                beam_range = beam[:, x] > 0
                # do linear regression. Todo: Find the best method to do non-negative linear regression
                coef = nnls(df[beam_range][names], df[beam_range]['Mu_t'])[0]
                rho_t[:, t, x] = coef
                # print progress
                print('>' * (counter % int(n_tomo * nx / 48) == 0), end='')
                counter += 1
        print('\n                 Finished calculation for'
              '\n                 ', n_tomo, ' tomo files in',
              round(time.clock() - start_time, 2), 'seconds')


    elif algorithm == 'sKES_equation':
        print('(calculate_rhot) Algorithm: "sKES_equation"')
        n_materials = mu_rhos.shape[0]
        nx = mu_rhos.shape[2]
        ny = mu_rhos.shape[1]
        n_proj = mu_t.shape[0]
        # mu_rhos.shape: [n_materials,ny,nx]  mu_t.shape:[n_proj,ny,nx]   beam.shape:[ny,nx]
        if use_torch:
            print('(calculate_rhot) Converting numpy.array to torch.tensor')
            mu_rhos = torch.from_numpy(mu_rhos)
            mu_t = torch.from_numpy(mu_t)
            beam = torch.from_numpy(beam)
            print('(calculate_rhot) Preparing matrix 1 of 2')
            inverted_mean_square_murhos = torch.zeros(size=(n_materials, n_materials, nx))
            for i in range(n_materials):
                for j in range(n_materials):
                    inverted_mean_square_murhos[i, j] = ((mu_rhos[i] * mu_rhos[j]) * beam).sum(dim=0) / (
                        beam.sum(dim=0))
            for x in range(nx):
                inverted_mean_square_murhos[:, :, x] = torch.inverse(inverted_mean_square_murhos[:, :, x])

            print('(calculate_rhot) Preparing matrix 2 of 2')
            if n_proj >= 200:
                print('|--------------------------------------------------|\n|', end='')
            sum_vector = torch.zeros(size=(n_materials, n_proj, nx))
            counter=0
            for i in range(n_materials):
                for j in range(n_proj):
                    sum_vector[i, j] = (mu_rhos[i] * mu_t[j] * beam).sum(dim=0) / beam.sum(dim=0)
                    if n_proj>=200:
                        print('>' * (counter % int(n_materials * n_proj / 50) == 0), end='')
                    counter += 1

            print('\n(calculate_rhot) Multiplying matrix 1 and 2')
            rho_t = sum_vector * 0.0
            for i in range(n_materials):
                rho_t[i] = (inverted_mean_square_murhos[i] * (sum_vector.transpose(1, 0))).sum(dim=1)
            rho_t=rho_t.numpy()

        else:
            print('(calculate_rhot) Using "numpy"')
            print('(calculate_rhot) Preparing matrix 1 of 2')
            inverted_mean_square_murhos = np.zeros(shape=(n_materials, n_materials, nx))
            for i in range(n_materials):
                for j in range(n_materials):
                    inverted_mean_square_murhos[i, j] = ((mu_rhos[i] * mu_rhos[j]) * beam).sum(axis=0) / (
                        beam.sum(axis=0))
            for x in range(nx):
                inverted_mean_square_murhos[:, :, x] = np.linalg.inv(inverted_mean_square_murhos[:, :, x])

            print('(calculate_rhot) Preparing matrix 2 of 2')
            sum_vector = np.zeros(shape=(n_materials, n_proj, nx))
            for i in range(n_materials):
                for j in range(n_proj):
                    sum_vector[i,j] = (mu_rhos[i] * mu_t[j] * beam).sum(axis=0) / (beam.sum(axis=0))

            print('(calculate_rhot) Multiplying matrix 1 and 2')
            rho_t = sum_vector * 0.0
            for i in range(n_materials):
                rho_t[i] = (inverted_mean_square_murhos[i] * (sum_vector.transpose(1, 0, 2))).sum(axis=1)

    else: raise Exception('"Algorithm" is not properly defined. Please choose from ["nnls", "sKES_equation"]')
    print('(calculate_rhot) Finished "calculate_rhot"')
    return rho_t


def signal_noise_ratio(mu_rhos,mu_t,rho_t,beam_parameters,tomo_data,use_torch=True):
    """
    Spectral KES for m-components SNR equation is used for SNR calculation.
    [Ref: Ying Zhu,2012 Dissertation]
    :param mu_rhos:
    :param mu_t:
    :param rho_t:
    :param beam_parameters:
    :param tomo_data:
    :param use_torch: Default True. Use torch.tensor for matrix operations. If False,
                      numpy.array is used instead, which is about half the speed of
                      tensor for the computation here.
    :return: snrs. Numpy array, in shape of [n_materials,n_projections,n_horizontal_positions]
    """
    if use_torch:
        print('(signal_noise_ratio) Preparing things for calculation')
        mu_rhos = torch.from_numpy(mu_rhos).float()
        mu_t = torch.from_numpy(mu_t).float()
        rho_t = torch.from_numpy(rho_t).float()
        beam = torch.from_numpy(beam_parameters.beam).float()
        flat = torch.from_numpy(beam_parameters.beam_files.flat).float()
        dark = torch.from_numpy(beam_parameters.beam_files.dark).float()
        beam_width = beam.sum(dim=0)
        n_materials = mu_rhos.size(0)
        n_energies = mu_rhos.size(1)
        nx = mu_rhos.size(2)
        n_proj = mu_t.size(0)
        # Todo: Use raw tomo_data? or filtered tomo_data from mu_t?
        matrix1 = 1/(tomo_data-dark)+1/(flat-dark)

        inverted_mean_square_murhos = torch.zeros(size=(n_materials,n_materials,nx),dtype=torch.float)
        for i in range(n_materials):
            for j in range(n_materials):
                inverted_mean_square_murhos[i,j] = (mu_rhos[i]*mu_rhos[j]*beam).sum(dim=0)/beam_width
        for x in range(nx):
            inverted_mean_square_murhos[:,:,x] = torch.inverse(inverted_mean_square_murhos[:,:,x])

        loops = n_materials*n_proj
        print('(signal_noise_ratio) Calculating SIGNAL to NOISE RATIO')
        if loops >= 200:
            print('|--------------------------------------------------|\n|', end='')
        counter = 0
        snrs = torch.zeros(size=(n_materials,n_proj,nx),dtype=torch.float)
        for m in range(n_materials):
            matrix_temp=torch.zeros(size=(n_materials,n_energies,nx),dtype=torch.float)
            for e in range(n_energies):
                matrix_temp[:,e,:] = inverted_mean_square_murhos[m]*mu_rhos[:,e,:]
            matrix2 = matrix_temp.sum(dim=0)**2  #[n_energies,nx]
            denom = torch.zeros(size=(n_proj,nx),dtype=torch.float)
            for p in range(n_proj):
                # Todo: make sure where the "beam_width" should be placed
                denom[p] = torch.sqrt((matrix2 * (matrix1[p] * beam)).sum(dim=0) / beam_width)
                if loops >= 200:
                    print('>' * (counter % int(loops / 50) == 0), end='')
                counter += 1
            snrs[m] = rho_t[m] / denom  # [n_proj,nx]
        print()
        snrs = snrs.numpy()

    else: # use numpy
        # mu_rhos = np.array(list(mu_rhos.values())).astype(float) #[n_materials,n_energies,nx]
        # rho_t = np.array(list(rho_t.values())).astype(float) #[n_materials,n_projections,nx]
        beam = beam_parameters.beam.astype(float)  #[n_energies,nx]
        flat = beam_parameters.beam_files.flat.astype(float)
        dark = beam_parameters.beam_files.dark.astype(float)
        beam_width = beam.sum(axis=0)
        n_materials = mu_rhos.shape[0]
        n_energies = mu_rhos.shape[1]
        nx = mu_rhos.shape[2]
        n_proj = mu_t.shape[0]

        matrix1 = 1/(tomo_data-dark)+1/(flat-dark) #[n_proj,n_energies,nx]

        # mu_rhos.shape: [n_materials,ny,nx]  mu_t.shape:[n_proj,ny,nx]   beam.shape:[ny,nx]
        inverted_mean_square_murhos = np.zeros(shape=(n_materials, n_materials, nx))
        for i in range(n_materials):
            for j in range(n_materials):
                inverted_mean_square_murhos[i, j] = ((mu_rhos[i] * mu_rhos[j]) * beam).sum(axis=0) / beam_width
        for x in range(nx):
            inverted_mean_square_murhos[:, :, x] = np.linalg.inv(inverted_mean_square_murhos[:, :, x])

        loops = n_materials*n_proj
        print('(signal_noise_ratio) Calculating SIGNAL to NOISE RATIO')
        if loops >= 200:
            print('|--------------------------------------------------|\n|', end='')
        counter = 0
        snrs = np.zeros((n_materials, n_proj, nx))
        for m in range(n_materials):
            matrix_temp = np.zeros((n_materials,n_energies,nx),float)
            for e in range(n_energies):
                matrix_temp[:, e, :] = inverted_mean_square_murhos[m] * mu_rhos[:, e, :]
            matrix2 = matrix_temp.sum(axis=0) ** 2  # [n_energies,nx]
            denom = np.zeros((n_proj, nx),float)
            for p in range(n_proj):
                # Todo: make sure where the "beam_width" should be placed
                denom[p] = np.sqrt((matrix2 * (matrix1[p] * beam)).sum(axis=0) / beam_width)
                if loops >= 200:
                    print('>' * (counter % int(loops / 50) == 0), end='')
                counter += 1
            snrs[m] = rho_t[m] / denom  # [n_proj,nx]
            print()

    return snrs # [n_materials,n_proj,nx]


def rho_in_ct(recon,names=None,center=[],width=0.0):
    """
    Calculate the average $\rho$ value in the target area in input recon image. If *center* and *width*
    are provided, target area is the wanted area. If not provided, this function will find the brightest
    area and calculate the average $\rho$ in it (actually only a square in the bright area)
    :param recon:Numpy array. Number of dimension >=2. If ndim==2,
    :param names: As reference for multi recon images, for the convenience of making plot labels.
    :param center: The center [x,y] pixel location of the square area to be calculated.
                   It can be either one center location [x,y], or a sequence of [x,y]s.
    :param width: The width of the square area to be calculated.
    :return:
    """

    # If more than 1 recon, use RECURSION to go through all of them.
    # The last two dimensions should be the recon image array
    if recon.ndim >= 3:
        if recon.ndim==3 and names and len(names)>1:
            # when there are more than one material, and there are no more than 3 dimensions,
            # the 1st dimension represents the materials. Use the names as the reference for
            # making plots.
            mean_rho=[]
            for i in range(recon.shape[0]):
                mean_rho.append(rho_in_ct(recon[i],names[i],center=center,width=width))
        else:
            mean_rho = []
            for i in range(recon.shape[0]):
                mean_rho.append(rho_in_ct(recon[i], names, center=center, width=width))
        return np.array(mean_rho)

    center = np.array(center).round().astype(int)
    if len(center)==0: # No center position is provided
        # find the bright area
        ct_filtered = median_filter(recon, 20)  # A big filter here, just to locate the center of target
        target = np.where(ct_filtered > ct_filtered.max() * 0.7) # Threshold is 0.7 to locate target
        xmin = target[1].min()
        xmax = target[1].max()
        ymin = target[0].min()
        ymax = target[0].max()
        x0 = round((xmin + xmax) / 2)
        y0 = round((ymin + ymax) / 2)
        center = np.array([x0, y0]).reshape((1,2))
        width = min([ymax - ymin, xmax - xmin]) / 2 # Use the 0.5 of maxs to be the width for square
    elif center.ndim==1:  # One center position is provided
        # y0=center[1]; x0=center[0]
        center = center.reshape((1,2))
    width = round(width)
    plt.figure()
    plt.imshow(recon*1000, cmap='gray_r')
    plt.colorbar().set_label('$mg/cm^3$',rotation=0,position=(1,1),labelpad=-5)
    mean_rho=[]
    for i in range(center.shape[0]):
        y0 = center[i][1]; x0 = center[i][0]
        width = round(width)
        x1 = int(x0 - 0.5 * width);
        x2 = int(x0 + 0.5 * width);
        y1 = int(y0 - 0.5 * width);
        y2 = int(y0 + 0.5 * width)
        mean_rho.append((recon[y1:y2, x1:x2].mean()*1000).round(2)) # change the unit to mg/cm^3
        draw_square([y0,x0], width, color='b')
    if names: # if we know the name of the material for the CT recon
        mean_molar_concentration = 1000 * np.array(mean_rho)/mphy.molar_mass(names) #mM
        plt.title('Concentration of ' + str(names))
        print('(rho_in_ct) Average density of '+str(names)+' in the square(s) is:\n'
          '          ', mean_rho,'mg/cm^3, Or',mean_molar_concentration.round(2),'mM.')
    else:
        print('(rho_in_ct) Average density in the square(s) is:\n'
          '          ', mean_rho,'mg/cm^3')

    return np.array(mean_rho)


def calculate_distance(path, x_location, proj, smooth_width=20):
    """
    Distance between beam focus and detector. This function is not required for spectral KES calculation.
    This function is only useful with Selenate absorption spectrum, in which there are two significant
    peaks near selenium k-edge energy.
    :return:
    """
    import math

    tomo = get_tomo_files(path)
    arrangement = nei_get_arrangement(path)
    beam_files = get_beam_files(path)
    flat = beam_files.flat
    dark = beam_files.dark
    mut = -np.log((tomo[proj]-dark)/(flat-dark))

    d_theta = math.radians(0.01999) # angle distance between the 2 peaks of selenate spectrum
    a0 = 5.4305  # Angstroms ##silicon crystal unit cell length at 300K.
    hkl = arrangement.hkl
    d_hkl = a0 / math.sqrt((np.array(hkl) ** 2).sum())
    pixel_size = arrangement.detector.pixel #mm

    mut = np.median(mut[:,x_location-round(smooth_width/2):x_location+round(smooth_width/2)],axis=1)











