import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import time
from scipy.ndimage import median_filter
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
import math_physics as mp
from toolkit import *


# __all__ = ['NeiSubDir','file_search',
#            'nei_get_arrangement','read_average_tifs','get_tomo_files','nei_determine_murhos',
#            'get_beam_files','nei','beam_near_edge_imaging', 'nei_beam_parameters',]

class NeiSubDir:
    """

    Parameters
    ----------
    path : str
        Path to the folder which contains one set of imaging data. Usually in this folder, there are subfolders named as 'DarkBefore', 'FlatBefore', 'EdgeABefore', 'Tomo', 'arrangement.dat'.
    After : bool, optional (default False). 
        Whether or not there are 'flat', 'dark', 'edge' collected after 'tomo' data collection.
    EdgeB : bool, optional (default False). 
        Whether or not there is a second 'edge' image folder.
    
    Returns
    -------
    A class object. 
        Contains paths of all the subfolders.
    """
    def __init__(self, path, After=False, EdgeB=False):
        self.DarkBefore = path /'DarkBefore'
        self.FlatBefore = path /'FlatBefore'
        self.EdgeABefore = path /'EdgeABefore'
        self.Tomo = path /'Tomo'
        if After == True:
            self.DarkAfter = path /'DarkAfter'
            self.FlatAfter = path /'FlatAfter'
            self.EdgeAAfter = path /'EdgeAAfter'
        else:  # use 'Before' for all 'After'
            self.DarkAfter = path /'DarkBefore'
            self.FlatAfter = path /'FlatBefore'
            self.EdgeAAfter = path /'EdgeABefore'
        if EdgeB == True:
            self.EdgeBBefore = path /'EdgeBBefore'
            if After == True:
                self.EdgeBAfter = path /'EdgeBAfter'
            else:
                self.EdgeBAfter = path /'EdgeBBefore'
        else:  # use EdgeA for all EdgeB
            self.EdgeBBefore = path /'EdgeABefore'
            if After == True:
                self.EdgeBAfter = path /'EdgeAAfter'
            else:
                self.EdgeBAfter = path /'EdgeABefore'


def nei_get_arrangement(path,setup_type='FILE'):
    """Reads in the arrangement for the experiment from the 'arrangement.dat' file. 

    Parameters
    ----------
    path : str
        Path to the folder which contains one set of imaging data. Usually in this folder, there are subfolders named as 'DarkBefore', 'FlatBefore', 'EdgeABefore', 'Tomo', 'arrangement.dat'.
    setup_type : str, optional,{'File','Manual'}. 
    Defines the way to input system arrangements.
        'File'. Find 'arrangement.dat' in the given path for arrangement. If 'arrangement.dat' is not found, it will automatically switch to 'Manual' mode.
        'Manual'. Popup a window to manually input arrangement. Todo.
    Returns
    -------
    arrangement : Class
        The system arrangement information
    """
    class get_arrangement:
        def __init__(self, path):
            filename = path/'arrangement.dat'
            data = pd.read_csv(filename, index_col=0, header=None, sep=r',\s+', engine='python')
            data=data[1] # series
            for i in range(len(data)):  # remove the ' sign in some strings
                data[i] = data[i].replace("'", "")
            data = data.to_dict()

            k = int(data['k'])
            l = int(data['l'])
            h = int(data['h'])
            energy_range_low = float(data['energy_range_low'])
            energy_range_high = float(data['energy_range_high'])

            self.diffaction_plane = data['diffraction_plane']
            self.type = data['type']
            self.chi_degrees = float(data['chi_degrees'])
            self.hkl = [h, k, l]
            self.energy = float(data['energy'])
            self.energy_range = [energy_range_low, energy_range_high]
            self.dist_fd = float(data['dist_fd'])
            self.detector = self.detector(data)

        class detector:
            def __init__(self, data):
                self.type = data['det_type']
                self.pixel = float(data['det_pixel'])
                self.flip = float(data['det_flip'])
                self.phperapu = float(data['det_phperapu'])
                self.disp_x_demag = float(data['det_disp_x_demag'])
                self.pct_max = float(data['det_pct_max'])
    try:
        arrangement = get_arrangement(path)
    except:
        print('(nei_get_arrangement)The "arrangement.dat" file either does not exist in the specified directory, or has error in it')
        #pop up a tkinter window for the 'arrangement' settings
        
    arrangement_parameters = {'diffaction_plane': ' DIFFRACTION PLANE:',
                              'type': ' TYPE:',
                              'chi_degrees': ' ASYMMETRY ANGLE (CHI):',
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
    """

    Parameters
    ----------
    files : `list` of `str`. 
        List of path+filenames of every image file.
    flip : 
        This is not used in the program for now (Nov. 15, 2018)
    xlow : 
        This is not used in the program for now (Nov. 15, 2018)
    xhigh : 
        This is not used in the program for now (Nov. 15, 2018)
    rotate_90 : 
        This is not used in the program for now (Nov. 15, 2018)
    twelve_bit : 
        This is not used in the program for now (Nov. 15, 2018)
    Returns
    -------
    average : array_like with shape [ny, nx].
        The average of input images.
    """

    n_files = len(files)
    image_array = []
    for i in range(n_files):
        image_array.append(np.array(Image.open(files[i])))
    image_array = np.array(image_array)
    average = image_array.mean(axis=0)
    return(average)


def get_beam_files(path,After=False,Verbose=False,clip=False, flip=False):
    '''
    return averaged flat, dark, edge images in form of 2d-arrays.
    Parameters
    ----------
    path : str 
        Path to the folder which contains one set of imaging data. Usually in this folder, there are subfolders named as 'DarkBefore', 'FlatBefore', 'EdgeABefore', 'Tomo', 'arrangement.dat'.
    After : 
        This is not used in the program for now (Nov. 15, 2018)
    clip : bool, optional (default False). 
        Whether to clip the side of the image view. This is needed only when the beam does not fill the whole view of image in the horizontal direction. Usually this is not needed.
    Verbose: bool, optional (default False). 
        Whether to print some information or show some images/plots for inspection during the running of this function.
    flip: 
        This is not used in the program for now (Nov. 15, 2018)
    Returns
    -------
    BeamFiles: Class object
        .flat: [ny, nx] array. Averaged flat image.
        .dark: [ny, nx] array. Averaged dark image.
        .edge: [ny, nx] array. Averaged edge image.
        .horizontal_low: int. Left border pixel index of useful beam region.
        .horizontal_high: int. Right border pixel index of useful beam region.
        .origin_beam_files: Class object. Contains beam files without clipping.

    '''
    path = Path(path)
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
        fwhm_left,left_low,left_high = mp.fwhm(x[idx_left],fn[idx_left],Verbose=Verbose)
        fwhm_right,right_low,right_high=mp.fwhm(x[idx_right],-fn[idx_right],Verbose=Verbose)
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
            self.horizontal_high= horizontal_high
            self.origin_beam_files= origin_beam_files
    return BeamFiles(flat,dark,edge,horizontal_low,horizontal_high,origin_beam_files)


def get_tomo_files(path, multislice=False, slice=0, n_proj=900, Verbose=False):
    """
    Read all the tomo files in tomo folder, or n_proj images for one slice of CT.
    Parameters
    ----------
    path : str
        Path to the folder which contains one set of imaging data. Usually in this folder, there are subfolders named as 'DarkBefore', 'FlatBefore', 'EdgeABefore', 'Tomo', 'arrangement.dat'.
    multislice : bool 
        Default False. Whether the 'tomo' data folder contains projection images for multislices.
    slice : int, non-negative (default 0). 
        If multislice = True, slice specifies which slice of all the projection images will be used for data analysis and CT reconstruction.
    n_proj : int, positive (default 900) 
        If multislice = True, n_proj tells the program how many projection images that one slice contains. It is used to determine which slice the projection images belongs to.
    Verbose : bool, optional (default False). 
        Whether to print some information or show some images/plots for inspection during the running of this function.    
    Returns
    -------
    tomo_data: float, [n_tomo_images, ny, nx] ndarray. 
        The tomo data in 3d array (of one slice, if multislice = True).
    """
    path = Path(path)
    sub_dir = NeiSubDir(path)# todo
    tomo_files_all = file_search(sub_dir.Tomo, '*tif')

    # get the tomo images for ONE slice when there are multislices of projections in tomo image folder.
    if multislice == True:
        i_begin = n_proj * slice
        i_end = i_begin + n_proj
        tomo_files = tomo_files_all[i_begin:i_end]
        print('(get_tomo_files) Tomo files in 1 slice loaded')
    else:
        tomo_files = tomo_files_all
    n_tomo = len(tomo_files)
    print('(get_tomo_files) Number of Tomo files: ', n_tomo)  # equal to n_projections

    # tomo files to data array
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
    return tomo_data


def nei_determine_murhos(materials, exy, gaussian_energy_width, interpol_kind='linear',
                         use_file=True,use_sm_data=False, use_measured_standard=False):
    '''
    Parameters
    ----------
    materials : `list` of `str`. 
        Example: materials = ['Water', 'Bone', 'Selenite', 'U' ]. Names of materials that we are investigating. They can be a standard compound name that gives all the information for the composition ('Na2SeO4'), or it can be an element name ('Se'), or it can be a nick name of a compound or material that we have stored the information in the 'COMPOSIT.DAT' file. The names are used to find the according μ/ρ of that material, so that they are very important.
    exy : float, [ny(n_energies), nx(n_horizontal_pixel)]  ndarray. 
        The energy 'map' on the detector. The value @ [iy, ix] location is the x-ray energy @ that pixel. The energy will be used to find according μ/ρ.
    gaussian_energy_width : float 
        The 'sigma' parameter of Gaussian function in terms of energy. It has a close relation to the system ENERGY RESOLUTION.
    interpol_kind : str, optional
        Ddefault 'linear'. This is used for the scipy.interpolate.interp1d function. See scipy documentation for detail (https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.interpolate.interp1d.html).
    use_file : bool
        Default True. Whether use the information in a stored file for the μ/ρ. If use_file = True, but actually there is no file for it, then it will automatically use computed value for \mu/\rhoμ/ρ.
    use_sm_data : bool
        Default False. Whether use a measured reference or a stored reference (in a file or by computation) from past for selenomethionine. This is not ready in the program yet. (Nov. 16, 2018)
    use_measured_standard: bool. 
        Default False. Whether use a measured reference or a stored reference (in a file or by computation)from past. This is not ready in the program yet. (Nov. 16, 2018)
    Returns
    -------
    murhos_all: float [n_materials, ny, nx] array. 
        For every material, the μ/ρ at every pixel on the detector.
    '''
    nx = exy.shape[1]
    ny = exy.shape[0]  # number of energies
    emin = np.median(exy, axis=1).min()
    emax = np.median(exy, axis=1).max()
    e_range = emax - emin
    # create an energies array with the size 5 * e_range
    energies = np.linspace(emin - 2 * (e_range), emax + 2 * (e_range), 5 * ny)

    murhos = {}
    # the following comment-out code is used in the past when we define materials together
    # with the source of murhos
    # for name, source in materials.items():
    #     print('(nei_determine_murhos) Getting murho data for ' + name)
        # if source.lower() == 'file':
        #     # get murho from saved file
        #     mu_rho = mp.murho_selenium_compounds(name, energies)
        # elif source.lower() == 'system':
        #     # get murho by calculating it for every element
        #     mu_rho = mp.murho(name, energies)
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
        mu_rho = mp.murho(name, energies, use_file=use_file)
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
    Locate the peak positions in the flat beam. Decide the useful region in the beam we are going to keep.
    Parameters
    ----------
    flat_dark : float, [ny, nx] ndarray. 
        average_flat - average_dark. Dark corrected average flat image.
    threshold : float
        In range of [0,100] as percentage of the max value of flat_dark. The threshold is specified in the 'arrangement.dat'. It is used to define the region of the beam (in the vertical direction) that we are going to use. When 'threshold' is 50, it means we use the fwhm as the useful beam region.
    no_fit : bool. Optional.
        Default False. Whether to do the polynomial fitting for peak positions over horizontal direction.
    poly_deg : int, positive
        Default 5. The degree of polynomial fitting. 5 is suggested.
    Verbose: bool, optional
        Default False. Whether to print some information or show some images/plots for inspection during the running of this function.    
    Returns
    -------
    BeamEdges: Class object
        Contains .top, .bot, .peak,.beam.
        .top: float, [nx] array. The top positions of the usable beam.
        .bot: float, [nx] array. The bottom positions of the usable beam.
        .peak: float, [nx] array. The peak positions of the flat_dark image. Usually polynomial fitted.
        .beam: float, [ny, nx] array. The usable beam region. This 2d array has the same shape flat and other image data array. And the usable region has is labeled with value 1, while the rest is labeled with value 0.
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
        spectrum_gauss,gauss_popt = mp.gaussfit(y_index,spectrum,estimate)
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

    ###################################################################
    #Approach 2:
    ################################################################
    # Linear regression of course only gives a straight line
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


def calculate_mut(tomo_data, beam_parameters,lowpass=False,ct=False,side_width=0):
    """
    mu*t =  (mu/rho) * (rho * t) = -ln[(tomo-dark)/(flat-dark)]
    Parameters
    ----------
    tomo_data : float, [n_tomo_images (n_angles), ny, nx] ndarray
        The tomo data in 3d array (of one slice, if multislice = True). It is returned by the function `get_tomo_file`.
    beam_parameters :  Class object. 
        Contains plenty of information, returned by the function nei_beam_parameters.
    lowpass: bool, optional
        Default False. Whether use a low_pass filter (Gaussian Filter)to smooth the signal.
    ct: bool, optional
        Default False. Whether use some pixels on the side of the image to calculate the absorption by air for two things.
            1. To remove the air absorption from the signal.
            2. To fix the ring artifect by normalizing the value of air absorption through the CT scan.
    side_width : int, required if `ct` is `True`
        Default 0. The number of pixels on the side of the image for calculating air absorption. Make sure it is only air in all the tomo images in these pixels.
    Returns
    -------
    mu_t: float, [n_tomo_images (n_angles), ny, nx] ndarray
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
                raise Exception('When `ct` is `True`, `side_width` needs to be defined for removing air absorption, and '
                                'it has to be an integer >0')
            # mut_left_total = (mu_t[i]*beam)[:,0:side_width].sum()
            # mut_right_total= (mu_t[i]*beam)[:,-side_width:].sum()
            mut_left_total = mu_t[i][:, 0:side_width].sum()
            mut_right_total = mu_t[i][:, -side_width:].sum()
            mut_left_count = beam[:,0:side_width].sum()
            mut_right_count= beam[:,-side_width:].sum()
            mut_left_avg   = mut_left_total/mut_left_count
            mut_right_avg  = mut_right_total/mut_right_count
            x_position = np.arange(nx)
            filter_1d = mut_left_avg+(mut_right_avg-mut_left_avg)*((x_position-side_width/2)/(nx-side_width))
            filter_2d = filter_1d*(beam*0.0+1.0)
            # remove air absorption from total mu_t
            # Todo: possible correction  mu_t[i] = mu_t[i]-filter_2d
            mu_t[i] = mu_t[i]-filter_2d
        if n_tomo>=200:
            print('>' * (counter % int(n_tomo / 50) == 0), end='')
        counter += 1
    print('')
    # use a lowpass filter (gaussian filter) to remove high frequency noise
    # Todo: decide the default sigma value for "lowpass" filter
    pixel_gaussian_width = beam_parameters.pixel_edge_width
    if lowpass:
        print('(calculate_mut) Applying lowpass filter along energy axis')
        mu_t = gaussian_filter(mu_t,[0,pixel_gaussian_width,0],truncate=3.0,mode='nearest')

    return mu_t


def calculate_rhot(mu_rhos,mu_t,beam,names,algorithm='',use_torch=True):
    """
    calculate the **ρ⋅t** for every material at every horizontal position in every projection
    Parameters
    ----------
    mu_rhos : dict, `{'material_name': [ny(n_energies),nx(n_horizontal_positions)] array, ...}` 
        `mu_rhos` is obtained from nei_determine_murhos.
    mu_t : float, [n_projections(n_angles),ny(n_energies),nx(n_horizontal_positions)] ndarray. 
        `mu_t` is obtained from calculate_mut.
    beam : float, [n_energies,nx] ndarray. 
        `beam` is obtained from nei_beam_parameters.
    algorithm : {"sKES_equation","nnls"}
        The algorithm used for calculating ρ⋅t. Available options are ["nnls", "sKES_equation"] for now (Aug 27, 2018).
            If "nnls": ``scipy.optimize.nnls`` will be used to perform linear regression for the spectrum at every horizontal position in every projection image.
            If "sKES_equation" (default) : An equation derived with least-square approach is used. Because the calculation here is done in matrix operation, it is much faster than iterating with "nnls". [Ref: Ying Zhu,2012 Dissertation]
    use_torch : bool, optional
        Default True. Whether use `torch.tensor` for matrix operation.
    Returns
    -------
    rho_t: float, [n_materials,n_tomo,nx] ndarray 
        Sinograms (the last two dimensions) of all the materials.
    """
    try:
        import torch
    except ModuleNotFoundError:
        print('(signal_noise_ratio) Module Pytorch is not found. Numpy will be used instead.')
        use_torch=False

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
    Parameters
    ----------
    mu_rhos : dict, {material_name: [ny(n_energies),nx(n_horizontal_positions)] array, ...} 
        `mu_rhos` (μ/ρ) is returned from nei_determine_murhos.
    mu_t : float, [n_projections(n_angles),ny(n_energies),nx(n_horizontal_positions)] array. 
        `mu_t` (μ⋅t) is returned from calculate_mut.
    rho_t: float, [n_materials,n_tomo,nx] array. 
        `rho_t` (ρ⋅t) is returned from calculate_rhot.
    Returns
    -------
    snrs: float, [n_materials,n_tomo,nx] array
    """
    try:
        import torch
    except:
        print('(signal_noise_ratio) Module pytorch is not available. Numpy will be used instead.')
        use_torch=False

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


def idl_recon(sinogram,pixel_size,center=0):
    """
    CT reconstruction with the "normalized_fbp" function from IDL. Note: Licensed IDL software is required.
    Parameters
    ----------
    sinogram: 3d or 2d-array, shape [(n_something),n_projections (n_angles),n_horizontal_positions]
    pixel_size: float. 
        Pixel size (resolution) of the detector in centimeter. For example: 0.0009 for 9um.
    center: int, default 0. 
        The rotation center in pixel during CT imaging. If 0, the horizontal center pixel in the image is used as the rotation center. Positive integer means the number of pixels to the right of the horizontal center of the image. Negative integer means the number of pixels to the left of the horizontal center of the image
    Returns
    -------
    recon : [(n_something),n_horizontal_positions,n_horizontal_positions] 3d or 2d array. 
        Reconstructed image(s) with a square shape.
    """
    from idlpy import IDL
    dimensions = sinogram.ndim
    if dimensions ==3:
        recon=[]
        for i in range(sinogram.shape[0]):
            recon.append(IDL.normalized_fbp(sinogram[i],dx=center,pixel=pixel_size))
        recon = np.array(recon)
    elif dimensions==2:
        recon = IDL.normalized_fbp(sinogram,dx=center,pixel=pixel_size)
    else:
        raise Exception('The dimensions of input "sinogram" should be either 2 or 3.'
                        ,dimensions,'dimensions were given.')
    return recon


def skimage_recon(sinogram,pixel_size=1.0,output_size=None,filter='ramp',center=0,degrees = 180,circle=True):
    """
    CT reconstruction using Inverse Radon Transform, with Filtered Back Projection algorithm.
    See "radon_transform" in skimage.transform for more detail.
    Parameters
    ----------
    sinogram : 3d or 2d-array [(n_something),n_projections (n_angles),n_horizontal_positions]
    pixel_size : float. 
        Pixel size (resolution) of the detector in centimeter. For example: 0.0009 for 9um.
    output_size : int, optional. 
        The width of the output reconstruction image. If output_size not specified, use the image horizontal width.
    filter : str, default 'ramp'. 
        The filter used for reconstruction. Please see "radon_transform" in skimage.transform for more detail. Todo: reference to radon_transform
    center : int, default 0. 
        The rotation center in pixel during CT imaging. If 0, the horizontal center pixel in the image is used as the rotation center. Positive integer means the number of pixels to the right of the horizontal center of the image. Negative integer means the number of pixels to the left of the horizontal center of the image
    degrees : int, {180,360}. 
        Either a 180 degree CT or a 360 degree CT.
    circle : bool, optional. 
        Assume the reconstructed image is zero outside the inscribed circle. The default behavior (None) is equivalent to False.
    Returns
    -------
    recon: 3d or 2d-array [(n_something),n_horizontal_positions,n_horizontal_positions].

    """
    from skimage.transform import iradon # The imported `iradon` here is a modified version
    # if sinogram.ndim>=3, which means there is more than one sinogram. Iterate them all.
    # The last two dimension should be the sinogram array.
    if sinogram.ndim>=3:
        recon = []
        for i in range(sinogram.shape[0]):
            recon.append(skimage_recon(sinogram[i],pixel_size=pixel_size,
                                       output_size=output_size,filter=filter,center=center,
                                       circle=circle))
        recon = np.array(recon)
        return recon
    # when sinogram.ndim==2, do the reconstruction.
    print('(skimage_recon) Started one CT reconstruction...',end='')
    n_proj = sinogram.shape[0]
    if not output_size: # If Output_size not specified, use the image horizontal width
        output_size=sinogram.shape[1]
    sinogram = sinogram.transpose(1,0) # transpose row and column to meet the order in skimage.transform.
    theta = np.linspace(0,degrees,n_proj) # make the angles of projections from 0 to 180 degree
    recon = iradon(sinogram,theta=theta,output_size=output_size,filter=filter,
                   center_drift=center,circle=circle)
    # correct result with pixel size (cm).
    recon = recon/pixel_size
    print('...Finished')

    return recon


def auto_center(data):
    """
    Use this function to auto locate the rotation center of CT reconstruction imaging.
    Parameters
    ----------
    data : 2-d array
        Dimension[0] refers to projections, and dimension[1] refers to positions in one projection.
    Returns
    -------
    drift : int 
        The relative pixel distance between the detected rotation center and the middle of the x-axis in the image. Negative value means rotation center is on the left side of the center of the image, and positive value means the opposite.
    """
    areas=[]
    center = round(data.shape[1]/2)
    drifts = np.arange(-center+1,center) # center is also the radius
    for d in drifts:
        if d <= 0:
            a1 = data[0]
            a2 = np.append(data[-1,2*d-1:0:-1],data[-1,2*d-1:])
        else:
            a1 = data[0]
            a2 = np.append(data[-1,0:2*d],data[-1,:2*d-1:-1])
    #     print(a1.shape,a2.shape,d)
        a = np.vstack((a1,a2)).min(axis=0)
    #     print(d,a.shape)
        area = (a-a.min()).sum()
        areas.append(area)
    drift = drifts[np.array(areas).argmax()]
    return(drift)


def rho_in_ct(recon,names=None,center=[],width=0.0,save_path=''):
    """
    Calculate the average $\rho$ value in the target area in input recon image. If *center* and *width*
    are provided, target area is the wanted area. If not provided, this function will find the brightest
    area and calculate the average $\rho$ in it (actually only a square in the bright area)
    Parameters
    ----------
    recon : n-dimension (n>=2) array [..,nx,nx]
        The last two dimensions form one image.
    names : list of strings, optional, default `None`. 
        If you want to label the calculated area with a name, pass in the list of names.
    center : 2-element list of `int`, or list of 2-element lists, Optional. 
        The center [x,y] pixel location of the square area to be calculated. It can be either one center location [x,y], or a sequence of [x,y]s. If this is specified, \rhoρ will be calculated for these locations.
    width: float, optional. 
        Use width to define edge length of the square area around the center locations for \rhoρ calculation.
    save_path: string, optional, default ''. 
        If a save_path is specified, recon-images marked with concentrations will be saved in the save_path.
    Returns
    -------
    np.array(mean_rhos): float array
        ρ values (mg/cm^3) in auto located bright areas or specified areas.
    """

    # If more than 1 recon, use RECURSION to go through all of them.
    # The last two dimensions should be the recon image array
    if recon.ndim >= 3:
        if recon.ndim==3 and names and len(names)>1:
            # when there are more than one material, and there are 3 dimensions,
            # the 1st dimension represents the materials. Use the names as the reference for
            # making plots.
            mean_rho=[]
            for i in range(recon.shape[0]):
                mean_rho.append(rho_in_ct(recon[i],names[i],center=center,width=width,save_path=save_path))
        else:
            mean_rho = []
            for i in range(recon.shape[0]):
                mean_rho.append(rho_in_ct(recon[i], names, center=center, width=width,save_path=save_path))
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
    plt.figure(figsize=(9,9))
    plt.imshow(recon*1000, cmap='gray_r')
    plt.colorbar().set_label('$mg/cm^3$',rotation=0,position=(1,1),labelpad=-5)
    mean_rhos=[]
    for i in range(center.shape[0]):
        y0 = center[i][1]; x0 = center[i][0]
        width = round(width)
        x1 = int(x0 - 0.5 * width);
        x2 = int(x0 + 0.5 * width);
        y1 = int(y0 - 0.5 * width);
        y2 = int(y0 + 0.5 * width)
        mean_rhos.append((recon[y1:y2, x1:x2].mean()*1000).round(2)) # change the unit to mg/cm^3
        draw_square([y0,x0], width, color='b')

    if save_path=='':
        save_path = choose_path('Please select directory to save reconstruction images:')
        if save_path =='': # No path is selected, just show the result and return.
            print('(rho_in_ct) Average density in the square(s) is:\n'
                  '          ', mean_rhos, 'mg/cm^3')
            return np.array(mean_rhos)

    save_path = str(save_path)
    figures = fnmatch.filter(os.listdir(save_path),'*.png')
    n_fig = len(figures)
    if names: # if we know the name of the material for the CT recon
        mean_molar_concentration = 1000 * np.array(mean_rhos)/mp.molar_mass(names) #mM
        plt.title('Concentration of ' + str(names))
        plt.savefig(save_path + str(names) + '.png')

        print('(rho_in_ct) Average density of '+str(names)+' in the square(s) is:\n'
          '          ', mean_rhos,'mg/cm^3, Or',mean_molar_concentration.round(2),'mM.')
    else:
        plt.savefig(save_path  + str(n_fig) + '.png')
        print('(rho_in_ct) Average density in the square(s) is:\n'
          '          ', mean_rhos,'mg/cm^3')

    return np.array(mean_rhos)


def calculate_distance(path, x_location, proj, smooth_width=20):
    """
    UNDER CONSTRUCTION

    Distance between beam focus and detector. This function is not required for spectral KES calculation.
    This function is only useful with Selenate absorption spectrum data, in which there are two significant
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


