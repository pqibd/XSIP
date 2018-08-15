import pandas as pd
from toolkit import *
import numpy as np
import math_physics as mphy
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter

class NEISubDir:
    def __init__(self,path,After=False,EdgeB=False):
        self.DarkBefore = path + r'/DarkBefore'
        self.FlatBefore = path + r'/FlatBefore'
        self.EdgeABefore= path + r'/EdgeABefore'
        self.Tomo       = path + r'/Tomo'
        if After == True:
            self.DarkAfter=path+ r'/DarkAfter'
            self.FlatAfter=path+ r'/FlatAfter'
            self.EdgeAAfter=path+r'/EdgeAAfter'
        else: # use 'Before' for all 'After'
            self.DarkAfter = path + r'/DarkBefore'
            self.FlatAfter = path + r'/FlatBefore'
            self.EdgeAAfter = path + r'/EdgeABefore'
        if EdgeB == True:
            self.EdgeBBefore=path + r'/EdgeBBefore'
            if After ==True:
                self.EdgeBAfter=path + r'/EdgeBAfter'
            else:
                self.EdgeBAfter=path + r'/EdgeBBefore'
        else: # use EdgeA for all EdgeB
            self.EdgeBBefore = path + r'/EdgeABefore'
            if After ==True:
                self.EdgeBAfter=path + r'/EdgeAAfter'
            else:
                self.EdgeBAfter=path + r'/EdgeABefore'

class NEIGetArrangement:
    #import pandas as pd
    def __init__(self, path):
        filename = path+r'\arrangement.dat'
        try:
            data = pd.read_csv(filename,index_col=0,
                           header=None, sep=r',\s+', engine='python').T
            print('arrangement.dat successfully loaded')
            # print(data.head())
            for i in range(data.shape[1]):# remove the ' sign in some strings
                data.iloc[0,i] = data.iloc[0,i].replace("'","")

            k                = int(data.loc[1, 'k'])
            l                = int(data.loc[1, 'l'])
            h                = int(data.loc[1, 'h'])
            energy_range_low = float(data.loc[1, 'energy_range_low'])
            energy_range_high= float(data.loc[1, 'energy_range_high'])

            self.diffaction_plane = data.loc[1, 'diffraction_plane']
            self.type             = data.loc[1, 'type']
            self.chi_degrees      = float(data.loc[1, 'chi_degrees'])
            self.hkl              = [h,k,l]
            self.energy           = float(data.loc[1, 'energy'])
            self.energy_range     = [energy_range_low,energy_range_high]
            self.dist_fd          = float(data.loc[1, 'dist_fd'])
            self.detector         = self.detector(data)
        except:
            print('Error: Something is wrong when reading in the arrangement file')

    class detector:
        def __init__(self,data):
            self.type = data.loc[1, 'det_type']
            self.pixel = float(data.loc[1, 'det_pixel'])
            self.flip = float(data.loc[1, 'det_flip'])
            self.phperapu = float(data.loc[1, 'det_phperapu'])
            self.disp_x_demag = float(data.loc[1, 'det_disp_x_demag'])
            self.pct_max = float(data.loc[1, 'det_pct_max'])

def nei_print_arrangement(arrangement):
    arrangement_parameters = {'diffaction_plane': ' DIFFRACTION PLANE:',
                              'type': ' TYPE:',
                              'chi_degrees': ' ASYMETTRY ANGLE (CHI):',
                              'hkl':' HKL:',
                              'energy': ' ENERGY (keV):',
                              'energy_range': ' ENERGY RANGE (keV):',
                              'dist_fd': ' DISTANCE FOCUS-DETECTOR (mm):',
                              'detector': '\n DETECTOR PARAMETERS:'}
    detector_parameters    = {'type': ' DETECTOR TYPE:',
                              'pixel': ' DETECTOR PIXEL (mm):',
                              'flip': ' DETECTOR FLIP:',
                              'disp_x_demag': ' DETECTOR DISPLAY X DEMAGNIFICATION:',
                              'pct_max': ' DETECTOR PERCENT THRESHOLD:',
                              'phperapu': ' phperapu'}
    for name, value in vars(arrangement).items():
        if name == 'detector': # print detector information in a seperate section
            print(arrangement_parameters[name])
            for name_2, value_2 in vars(arrangement.detector).items():
                print(detector_parameters[name_2], value_2)
        else:
            print(arrangement_parameters[name], value)

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

class GetBeamFiles:
    # return raw_images, only averaged, nothing else.

    def __init__(self,sub_dir,Verbose=False,clip=False,*,flip=False):
        flat_path = sub_dir.FlatBefore
        dark_path = sub_dir.DarkBefore
        edge_path = sub_dir.EdgeABefore

        flat_files = file_search(flat_path,'*.tif')
        dark_files = file_search(dark_path,'*.tif')
        edge_files = file_search(edge_path,'*.tif')

        # read and average the raw data files - flats, darks, edges
        # 2d numpy.array returned by read_average_tifs
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
            n_horizontal = flat.shape[1]

        self.flat = flat
        self.dark = dark
        self.edge = edge
        self.n_horizontal = n_horizontal
        self.n_vertical   = n_vertical
        self.horizontal_low = horizontal_low
        self.horizontal_high=horizontal_high
        self.origin_beam_files= origin_beam_files

def beam_edges(flat_dark,threshold,no_fit=False,Verbose=False,poly_deg=5):
    # flat_dark: flat-dark, 2D array
    # threshold: determines where we cut the beam on the energy axis
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

    # top_positions = median_filter(top_positions,5)
    # bot_positions = median_filter(bot_positions,5)
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

    return s
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

def nei_determine_murhos(materials, exy, gaussian_energy_width,use_sm_data=False,
                         use_measured_standard=False):
    '''
    ; mats structure
    ; tags:
    ; names - what you call each material, i.e. 'selenite'
    ; types - how we find the mu/rhos
    ;       - IDL   : get it from IDL
    ;       - FILE : get it from a file
    ; i.e. mats = { names:names, types:types }
    ; where: names = ['Water', 'Bone', 'Selenite', 'U'   ]
    ;        types = [  'IDL',  'IDL',     'FILE', 'IDL' ]
    :param materials:
    :param exy:
    :param gaussian_energy_width:
    :param use_sm_data:
    :param use_measured_standard:
    :return:
    '''
    n_energies = exy.shape[0]
    emin = exy.min()
    emax = exy.max()
    e_range=emax-emin
    e_range=np.linspace(emin-2*(e_range),emin+2*(e_range),5*n_energies)






