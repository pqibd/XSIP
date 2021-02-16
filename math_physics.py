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


class Constants:
    elements = ['H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE',
                'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA',
                'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
                'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR',
                'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN',
                'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND',
                'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',
                'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
                'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH',
                'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM']
    atom_weights = [ 1.0079000,      4.0026002,      6.9410000,      9.0121803,
                    10.8100004,     12.0109997,     14.0066996,     15.9994001,
                    18.9984035,     20.1790009,     22.9897709,     24.3050003,
                    26.9815407,     28.0855007,     30.9737606,     32.0600014,
                    35.4529991,     39.9480019,     39.0983009,     40.0800018,
                    44.9558983,     47.9000015,     50.9415016,     51.9959984,
                    54.9379997,     55.8470001,     58.9332008,     58.7000008,
                    63.5460014,     65.3799973,     69.7200012,     72.5899963,
                    74.9216003,     78.9599991,     79.9039993,     83.8000031,
                    85.4677963,     87.6200027,     88.9058990,     91.2200012,
                    92.9064026,     95.9400024,     98.0000000,    101.0699997,
                   102.9055023,    106.4000015,    107.8679962,    112.4100037,
                   114.8199997,    118.6900024,    121.7500000,    127.5999985,
                   126.9045029,    131.3000031,    132.9053955,    137.3300018,
                   138.9055023,    140.1199951,    140.9076996,    144.2400055,
                   145.0000000,    150.3999939,    151.9600067,    157.2500000,
                   158.9253998,    162.5000000,    164.9304047,    167.2599945,
                   168.9342041,    173.0399933,    174.9669952,    178.4900055,
                   180.9479065,    183.8500061,    186.2070007,    190.1999969,
                   192.2200012,    195.0899963,    196.9665070,    200.5899963,
                   204.3699951,    207.1999969,    208.9803925,    209.0000000,
                   210.0000000,    222.0000000,    223.0000000,    226.0254059,
                   227.0278015,    232.0381012,    231.0359039,    238.0290070,
                   237.0482025,    244.0000000,    243.0000000,    247.0000000,
                   247.0000000,    251.0000000,    252.0000000,    257.0000000]
    element_density = [0.000090, 0.000179, 0.530000, 1.850000,
                    2.340000, 2.620000, 0.001251, 0.001429,
                    0.001696, 0.000901, 0.970000, 1.740000,
                    2.700000, 2.330000, 1.820000, 2.070000,
                    0.003170, 0.001784, 0.860000, 1.550000,
                    3.000000, 4.500000, 5.800000, 7.190000,
                    7.430000, 7.860000, 8.900000, 8.900000,
                    8.960000, 7.140000, 5.910000, 5.320000,
                    5.720000, 4.800000, 3.120000, 0.003740,
                    1.530000, 2.600000, 4.500000, 6.490000,
                    8.550000, 10.200000, 11.500000, 12.200000,
                    12.400000, 12.000000, 10.500000, 8.650000,
                    7.310000, 7.300000, 6.680000, 6.240000,
                    4.920000, 0.005890, 1.870000, 3.500000,
                    6.700000, 6.780000, 6.770000, 7.000000,
                    6.475000, 7.540000, 5.260000, 7.890000,
                    8.270000, 8.540000, 8.800000, 9.050000,
                    9.330000, 6.980000, 9.840000, 13.100000,
                    16.600000, 19.299999, 21.000000, 22.400000,
                    22.500000, 21.400000, 19.299999, 13.530000,
                    11.850000, 11.400000, 9.800000, 9.400000,
                    0.000000, 0.009910, 0.000000, 5.000000,
                    10.070000, 11.700000, 15.400000, 18.900000,
                    20.400000, 19.799999, 13.600000, 13.511000,
                    0.000000, 0.000000, 0.000000, 0.000000]
    atom_mass_unit = 1.6606 * 10**(-24) # grams. atom mass unit
    electron_radius = 2.8179 * 10**(-13) # cm. Classical electron radius
    h = C.h # Planck constant. Joer*s
    c = C.c # speed of light. meter/s
    eV = C.eV # eV to Joer
    a0 = 5.4305  # Unit: Angstroms ##silicon crystal unit cell length at 300K.
                                                            ## This is usually used as the internal standard for silicon


def fwhm(x,y,Verbose=False):
    '''Calculate the Full Width Half Maximum (FWHM) of the input `x` for index and `y` for value.
    (0,0,0) will be returned if 'max_index' is too close to the side of the array, or if the left 
    or right region out of 'FWHM' has length of 0.
    Parameters
    ----------
    x : array
        Index array. It could be the array index of `y`, or energy values as index for `y`.

    y : array
        Value array.

    Verbose : bool, optional, default False
        Whether to print some information or show some images/plots for inspection during the running of this function

    Returns
    -------

    (fwhm, left_fwhm, right_fwhm) : tuple of float or tuple of int
        fwhm : fwhm using `x` as the index
        left_fwhm : the `x` value of the very **left** element in the `fwhm` region.
        right_fwhm : the `x` value of the very **right** element in the `fwhm` region.
    '''
    ## If max is too close to either end, raise error
    ind_max = y.argmax()

    if ind_max<=2 or ind_max>=len(y)-2:
        print('(fwhm)ind_max<=2 or ind_max>=len(y)-2')
        return(0,0,0)
    half_max  = y.max()/2
    ind_max   = y.argmax()
    low_index = np.where(y<half_max)[0] # index of all points lower than 'half_max'
    left      = low_index[low_index<ind_max] # index of points on the left side out of 'fwhm'
    right     = low_index[low_index>ind_max] # index of points on the right side out of 'fwhm'
    if len(left)*len(right)==0:
        print('(fwhm)len(left)*len(right)')
        return(0,0,0)
    left_fwhm = x[left[-1]] # The index of the point at the very right of the 'left', so that it is also the left limit of 'fwhm'
    right_fwhm= x[right[0]] # The index of the point at the very left of the 'right', so that it is also the right limit of 'fwhm'
    fwhm      = right_fwhm-left_fwhm # Width
    print(fwhm)

    # show the plot if Verbose==True
    if Verbose:
        plt.figure()
        plt.scatter(np.arange(len(y)),y,s=3)
        plt.plot([left_fwhm,left_fwhm],[0,y.max()],color='y')
        plt.plot([right_fwhm,right_fwhm],[0,y.max()],color='y')
        plt.title('FWHM')
        plt.show()
    return(fwhm,left_fwhm,right_fwhm)


def gaussfit(x, y, *estimate):
    '''Get the best fitted Gaussian curve for the input `x`(index) and `y`(value).

    Only 3 terms are used here for fitting the Gaussian function.
    >In mathematics, a Gaussian function, often simply referred to as a Gaussian, is a function of the form:
    $$f(x)=a\cdot e^{-{\frac {(x-b)^{2}}{2c^{2}}}}$$
    >for arbitrary real constants a, b and non zero c. It is named after the mathematician Carl Friedrich Gauss. The graph of a Gaussian is a characteristic symmetric "bell curve" shape. The parameter a is the height of the curve's peak, b is the position of the center of the peak and c (the standard deviation, sometimes called the Gaussian RMS width) controls the width of the "bell".

    Parameters
    ----------
    x : array
        Index array.

    y : array
        Value array. 

    Returns
    -------
    y_gauss : array
        The best fitted Gaussian curve.

    popt : list of float
        The values of the 3 terms to define the fitted curve.
    '''

    def gauss_func(x, a0, a1, a2):  # ,a3=0,a4=0,a5=0):
        # define the gaussian function with 3 terms
        return a0 * np.exp(-(x - a1) ** 2 / (2 * a2 ** 2))  # +a3+a4*x+a5*x**2

    # Call the curve_fit, imported from scipy.optimize
    if estimate:
        popt, pcov = curve_fit(gauss_func, x, y, estimate)
    else:
        popt, pcov = curve_fit(gauss_func, x, y)
    y_gauss = gauss_func(x, popt[0], popt[1], popt[2])
    return y_gauss, popt


def polyfit(x, y, degree):
    '''A wrapper of `numpy.polyfit` for least squares polynomial fit.
    
    `numpy.polyfit` returns a vector of coefficients p that minimises the squared error. This wrapper uses these coefficients and the input x-coordinates to produce the fitting polynomial curve. 
    > Fit a polynomial p(x) = p[0] * x**deg + ... + p[deg] of degree deg to points (x, y).

    Parameters
    ---
    x : array_like, shape(M,)
        x-coordinates of the M sample points (x[i], y[i])

    y : array_like, shape (M,) or (M, K)
        y-coordinates of the sample points. Several data sets of sample points sharing the same x-coordinates can be fitted at once by passing in a 2D-array that contains one dataset per column.

    degree : int
        Degree of the fitting polynomial.

    Returns
    ---
    y_poly : array, shape(M,) or (M, K)
        y-coordinates of the fitting polynomial curve, with the same shape of parameter `y`.
    '''
    coef = np.polyfit(x, y, degree)
    p = np.poly1d(coef)
    y_poly = p(x)
    return y_poly


def element_info(element_name,no_whine=False):
    '''
    search through all the elements list, to see if we can find our interested element
    'element_name' takes in the name of one element.
    function return one row in form of dataframe, containing some physical information
    of the element
    :param element_name:
    :param no_whine:
    :return:
    '''

    elements = ['H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F', 'NE',
                'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA',
                'SC', 'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN',
                'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR',
                'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN',
                'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND',
                'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',
                'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG',
                'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH',
                'PA', 'U', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM']

    lineno = [1, 31, 61, 90, 119, 148, 178, 208, 237, 266,
              297, 326, 355, 385, 414, 443, 472, 501, 530, 559,
              588, 616, 644, 673, 701, 729, 758, 787, 817, 846,
              877, 908, 939, 970, 1001, 1033, 1065, 1097, 1129, 1161,
              1192, 1224, 1256, 1288, 1320, 1352, 1384, 1416, 1449, 1482,
              1515, 1548, 1581, 1616, 1651, 1684, 1717, 1752, 1788, 1824,
              1859, 1895, 1928, 1964, 1999, 2034, 2069, 2104, 2139, 2171,
              2203, 2235, 2267, 2300, 2333, 2365, 2397, 2429, 2462, 2494,
              2526, 2559, 2592, 2628, 2661, 2694, 2727, 2760, 2793, 2826,
              2858, 2890, 2923, 2955, 2988, 3021, 3054, 3087, 3120, 3155]

    fullnames = [ 'Hydrogen'    , 'Helium'      , 'Lithium'     , 'Beryllium'   ,
                  'Boron'       , 'Carbon'      , 'Nitrogen'    , 'Oxygen'      ,
                  'Fluorine'    , 'Neon'        , 'Sodium'      , 'Magnesium'   ,
                  'Aluminum'    , 'Silicon'     , 'Phosphorus'  , 'Sulfur'      ,
                  'Chlorine'    , 'Argon'       , 'Potassium'   , 'Calcium'     ,
                  'Scandium'    , 'Titanium'    , 'Vanadium'    , 'Chromium'    ,
                  'Manganese'   , 'Iron'        , 'Cobalt'      , 'Nickel'      ,
                  'Copper'      , 'Zinc'        , 'Gallium'     , 'Germanium'   ,
                  'Arsenic'     , 'Selenium'    , 'Bromine'     , 'Krypton'     ,
                  'Rubidium'    , 'Strontium'   , 'Yttrium'     , 'Zirconium'   ,
                  'Niobium'     , 'Molybdenum'  , 'Technetium'  , 'Ruthenium'   ,
                  'Rhodium'     , 'Palladium'   , 'Silver'      , 'Cadmium'     ,
                  'Indium'      , 'Tin'         , 'Antimony'    , 'Tellurium'   ,
                  'Iodine'      , 'Xenon'       , 'Cesium'      , 'Barium'      ,
                  'Lanthanum'   , 'Cerium'      , 'Praseodymium', 'Neodymium'   ,
                  'Promethium'  , 'Samarium'    , 'Europium'    , 'Gadolinium'  ,
                  'Terbium'     , 'Dysprosium'  , 'Holmium'     , 'Erbium'      ,
                  'Thulium'     , 'Ytterbium'   , 'Lutetium'    , 'Hafnium'     ,
                  'Tantalum'    , 'Tungsten'    , 'Rhenium'     , 'Osmium'      ,
                  'Iridium'     , 'Platinum'    , 'Gold'        , 'Mercury'     ,
                  'Thallium'    , 'Lead'        , 'Bismuth'     , 'Polonium'    ,
                  'Astatine'    , 'Radon'       , 'Francium'    , 'Radium'      ,
                  'Actinium'    , 'Thorium'     , 'Protactinium', 'Uranium'     ,
                  'Neptunium'   , 'Plutonium'   , 'Americium'   , 'Curium'      ,
                  'Berkelium'   , 'Californium' , 'Einsteinium' , 'Fermium'     ]

    #The atomic weights of the elements
    atom_weights = [ 1.0079000,      4.0026002,      6.9410000,      9.0121803,
                    10.8100004,     12.0109997,     14.0066996,     15.9994001,
                    18.9984035,     20.1790009,     22.9897709,     24.3050003,
                    26.9815407,     28.0855007,     30.9737606,     32.0600014,
                    35.4529991,     39.9480019,     39.0983009,     40.0800018,
                    44.9558983,     47.9000015,     50.9415016,     51.9959984,
                    54.9379997,     55.8470001,     58.9332008,     58.7000008,
                    63.5460014,     65.3799973,     69.7200012,     72.5899963,
                    74.9216003,     78.9599991,     79.9039993,     83.8000031,
                    85.4677963,     87.6200027,     88.9058990,     91.2200012,
                    92.9064026,     95.9400024,     98.0000000,    101.0699997,
                   102.9055023,    106.4000015,    107.8679962,    112.4100037,
                   114.8199997,    118.6900024,    121.7500000,    127.5999985,
                   126.9045029,    131.3000031,    132.9053955,    137.3300018,
                   138.9055023,    140.1199951,    140.9076996,    144.2400055,
                   145.0000000,    150.3999939,    151.9600067,    157.2500000,
                   158.9253998,    162.5000000,    164.9304047,    167.2599945,
                   168.9342041,    173.0399933,    174.9669952,    178.4900055,
                   180.9479065,    183.8500061,    186.2070007,    190.1999969,
                   192.2200012,    195.0899963,    196.9665070,    200.5899963,
                   204.3699951,    207.1999969,    208.9803925,    209.0000000,
                   210.0000000,    222.0000000,    223.0000000,    226.0254059,
                   227.0278015,    232.0381012,    231.0359039,    238.0290070,
                   237.0482025,    244.0000000,    243.0000000,    247.0000000,
                   247.0000000,    251.0000000,    252.0000000,    257.0000000]
    #The density of the elements if needed
    element_dens = [0.000090,       0.000179,       0.530000,       1.850000,
                    2.340000,       2.620000,       0.001251,       0.001429,
                    0.001696,       0.000901,       0.970000,       1.740000,
                    2.700000,       2.330000,       1.820000,       2.070000,
                    0.003170,       0.001784,       0.860000,       1.550000,
                    3.000000,       4.500000,       5.800000,       7.190000,
                    7.430000,       7.860000,       8.900000,       8.900000,
                    8.960000,       7.140000,       5.910000,       5.320000,
                    5.720000,       4.800000,       3.120000,       0.003740,
                    1.530000,       2.600000,       4.500000,       6.490000,
                    8.550000,      10.200000,      11.500000,      12.200000,
                   12.400000,      12.000000,      10.500000,       8.650000,
                    7.310000,       7.300000,       6.680000,       6.240000,
                    4.920000,       0.005890,       1.870000,       3.500000,
                    6.700000,       6.780000,       6.770000,       7.000000,
                    6.475000,       7.540000,       5.260000,       7.890000,
                    8.270000,       8.540000,       8.800000,       9.050000,
                    9.330000,       6.980000,       9.840000,      13.100000,
                   16.600000,      19.299999,      21.000000,      22.400000,
                   22.500000,      21.400000,      19.299999,      13.530000,
                   11.850000,      11.400000,       9.800000,       9.400000,
                    0.000000,       0.009910,       0.000000,       5.000000,
                   10.070000,      11.700000,      15.400000,      18.900000,
                   20.400000,      19.799999,      13.600000,      13.511000,
                    0.000000,       0.000000,       0.000000,       0.000000]

    df = pd.DataFrame.from_dict({'Elements':elements,'FullNames':fullnames,'Lineno':lineno,
                                 'AtomWeights':atom_weights,'ElementDens':element_dens})

    name = element_name.replace(' ', '').upper()

    information = df[df.Elements==name]
    if len(information)==0:
        raise Exception('Element "'+name+'" is not found in the elements list.')
    class InfoClass:
        def __init__(self,information):
            self.element = information.Elements.values[0]
            self.full_name= information.FullNames.values[0]
            self.lineno  = information.Lineno.values[0]
            self.atom_weight= information.AtomWeights.values[0]
            self.density = information.ElementDens.values[0]
            self.atom_z  = information.index.values[0]+1
    info = InfoClass(information)

    return(info)


def molar_mass(name,Verbose=False):
    """
    Unit: g/mol
    :param name: Element name or compound name. Example: "Na2SeO4"
    :return: Molar mass of input element or compound. Unit: g/mol.
    """
    atom_weights = dict(zip(Constants.elements, Constants.atom_weights))
    try: # Find all the elements and the according atom counts by the input name of compound.
         # If fails, go to the "except".
        total_mass = 0
        # atom_weights = dict(zip(Constants.elements,Constants.atom_weights))
        all = re.findall(r'[A-Z][a-z]*|\d+|\(|\)',name)
        if Verbose:
            print(all)
        if len(all)==0:
            raise Exception('Name is not in standard form.')
        for i in range(len(all)):
            if not all[i].isdigit():
                try:
                    total_mass+= atom_weights[all[i].upper()]
                except:
                    raise Exception('No such element:'+str(all[i]))
            else:
                number= float(all[i])-1
                total_mass+=(atom_weights[all[i-1].upper()]*number)
    except: # Find the name of the compound in "composit.dat", and use the molecular information
            # provided there to calculate molar mass.
        total_mass = 0
        with open(Path('MU/COMPOSIT.DAT'), 'r') as file:
            content = file.read().upper()
        # Use pattern1 to find the composite we want, and get the number of types of element it has
        pattern1 = ' ' + name.upper() + '.+\n'
        thereitis = re.search(pattern1, content)
        if thereitis:
            n_elements = thereitis.group().split()[2]
        else:
            raise Exception(name.upper() + ' is not a legitimate composite name in "COMPOSIT.DAT"!')
        # Use pattern2 to pack up all the information for the composite we want
        pattern2 = ' ' + name.upper() + '.+(?:\n.+){' + n_elements + '}'
        # And put things into a list like this
        # [['SEO3', '1.00000', '3'],
        # ['NA', '2.0000000'],
        # ['SE', '1.0000000'],
        # ['O', '3.0000000']]
        composite_info = [i.split() for i in re.search(pattern2, content).group().split('\n')]
        if Verbose:
            print('Composite_info:\n',composite_info)
        for element in composite_info[1:]:  # composite_info[0] will not be used
            element_name = element[0]
            atom_count = float(element[1])
            atom_weight = atom_weights[element_name]
            total_mass+=(atom_count * atom_weight)
    return total_mass #g/mol


def density(name,Verbose=False):
    """
    Now[2018-Oct-09] this density calculator only works for element density from `Constants`.
     Unit: g/cm^3
    :param name: Element name.
    :param Verbose:
    :return:
    """
    element_dens = dict(zip(Constants.elements, Constants.element_density))
    return element_dens[name.upper()]


def read_absorber(element_name,Verbose=False):
    """

    :param element_name:
    :param Verbose:
    :return:
    """
    element_name = element_name.upper()
    if Verbose:
        print('(read_absorber) Looking for "'+element_name+'" in absorber.dat')
    # open 'absorber.dat' file, get to the interested element by regex,
    # get the ek, eta, ef..
    abs_path = Path(os.path.dirname(os.path.abspath(__file__)))
    with open(abs_path/Path('MU/ABSORBER.DAT'), 'r') as file:
        content = file.read()
    # create a pattern to find the element and the information we want
    pattern = '(?:\A|\n)(' + element_name + ' (?:.*\n\W)+.*)\n\w{1,2} ' # all the info for one element
    ms = re.search(pattern, content)
    if Verbose:
        print(ms.group())
    if ms: # if ms is not none (we found what we wanted), save the info into list of lists
        ms = ms.group(1)
        ms = ms.split('\n')# split different lines to a list
        lines = [] # save the values in each line in to a list
        for line in ms:
            lines.append(line.split())
        ms = lines
    else:
        raise Exception('element "'+ element_name+'" cannot be found in "absorber.dat" file.\n'
                        'Check the pattern of element name. For example, "SE" is correct pattern for Selenium.\n'
                        'If nothing is wrong with the name, you should go the "get_absorber.py" to '
                        'check \nthe regular expression searching pattern')
    # prepare de,c_a, c_b, e0, xj
    if Verbose:
        print(ms)
    n_edges = int(ms[0][1])
    de= np.zeros(3)
    c_a=np.zeros([3,8])
    c_b=np.zeros([3,16])
    e0 =np.zeros(3)
    xj =np.zeros([3,6])
    for i in range(3):
        de[i]= float(ms[n_edges+(n_edges-1)//6+2+6*i][0])
        c_a[i,:] = ms[n_edges+(n_edges-1)//6+2+6*i][1:9]
        c_b[i:,] = ms[n_edges+(n_edges-1)//6+3+6*i]+ms[n_edges+(n_edges-1)//6+4+6*i]+ms[n_edges+(n_edges-1)//6+5+6*i]\
                   +ms[n_edges+(n_edges-1)//6+6+6*i]
        e0[i]    = float(ms[n_edges+(n_edges-1)//6+7+6*i][0])
        xj[i,:]  = ms[n_edges+(n_edges-1)//6+7+6*i][1:7]
    c_a = np.array(c_a).astype(float)
    c_b = np.array(c_b).astype(float)
    xj  = np.array(xj).astype(float)
    edges=[]
    for i in range(1,(2+(n_edges-1)//6)):
        edges=edges+ms[i]
    edges = np.array(edges).astype(float)
    class Results:
        def __init__(self,ms,de,c_a,c_b,e0,xj):
            self.n_edges = int(ms[0][1])
            self.edges   = edges
            self.eta = float(ms[2+(n_edges-1)//6][0])
            self.ef  = float(ms[2+(n_edges-1)//6][1])
            self.ek  = float(ms[2+(n_edges-1)//6][2])
            self.za  = float(ms[2+(n_edges-1)//6][3])
            self.a   = np.array(ms[2+(n_edges-1)//6+1:n_edges+2+(n_edges-1)//6]).astype(float)
            self.de  = de
            self.c_a = c_a
            self.c_b = c_b
            self.e0  = e0
            self.xj  = xj
    return Results(ms,de,c_a,c_b,e0,xj)


def mu_calculator(element_name,energies):
    # what is amu?
    # what is emu?
    ####### get element absorption information first  ###############
    ms = read_absorber(element_name)

    ####### do the calculation for mu  ##############################
    try:
        n_energies = len(energies)
    except TypeError:
        energies = np.array([energies])
        n_energies = len(energies)
    amu = np.zeros(n_energies)
    emu = np.zeros(n_energies)

    for i in range(n_energies):
        vector = [energies[i]**(-(j+1)) for j in range(4)] #[^-1,^-2,^-3,^-4]
        # decide which row of a to use. If all the edges are lower than this energy, use the first row of a.
        # if not, use the lowest one that is higher than this energy(but not the last one)
        edges = ms.edges
        e     = energies[i]
        ic = len(edges[edges>=e])-1

        if ic==-1:
            tau = (ms.a[0]*vector).sum()
        elif ic >= len(edges)-1:
            tau = (ms.a[-2]*vector).sum()
        else:
            tau = (ms.a[ic]*vector).sum()

        amu[i] = tau
        if ms.ek <= energies[i]:
            emu[i] = tau * (1-ms.eta*ms.ef/energies[i])
        else:
            emu[i] = tau
    return(amu,emu)


def sigma_calculator(element_name,energies):
    # THIS FUNCTION CALCULATES SCATTERING CROSS SECTIONS
    # SIG(MA)COH(ERENT), SIG(MA)INC(OHERENT), AND SIG(MA)EN(ENERGY)
    # IN UNITS OF CM**2/GM
    # COEFFICIENTS OF EXPANSIONS GIVEN BY C"S
    # ZA=ATOMIC NUMBER/ATOMIC MASS
    # CALCULATION OF SIG(MA)"S AT ENERGIES E, UNIT IS KEV
    try:
        n_energies = len(energies)
    except TypeError:
        energies = np.array([energies])
        n_energies = len(energies)
    ms = read_absorber(element_name)
    sig_coh = np.zeros(n_energies)
    sig_incoh = np.zeros(n_energies)
    sig_energy = np.zeros(n_energies)

    for i in range(n_energies):
        e = energies[i]
        x = e / 511.006  # what is this?
        sig = ms.de * (ms.c_a[:, 0] + ms.c_a[:, 1] * x + ms.c_a[:, 2] * x ** 2 + ms.c_a[:, 3] * x ** 3) / (
                    1.0 + ms.c_a[:, 4] * x + ms.c_a[:, 5] * x ** 2 + ms.c_a[:, 6] * x ** 3 + ms.c_a[:, 7] * x ** 4)
        spln = np.zeros(3)
        for j in range(3):  # for coh, incoh, energy
            if e >= ms.e0[j]:
                spln[j] = ms.c_b[j, 10] + ms.c_b[j, 11] / e ** 2 + ms.c_b[j, 12] / e ** 3 \
                          + ms.c_b[j, 13] / (ms.c_b[j, 14] * e + ms.c_b[j, 15] * e ** 4)
            else:
                sum = 0
                for k in range(4, 10):
                    if e > ms.xj[j, k - 4]:
                        sum = sum + ms.c_b[j, k] * (e - ms.xj[j, k - 4]) ** 3
                spln[j] = ms.c_b[j, 0] + ms.c_b[j, 1] * e + ms.c_b[j, 2] * e ** 2 + ms.c_b[j, 3] * e ** 3 + sum
        sig_coh[i] = ms.za * sig[0] * spln[0]
        sig_incoh[i] = ms.za * sig[1] * spln[1]
        sig_energy[i] = ms.za * sig[2] * spln[2]

    return(sig_coh,sig_incoh,sig_energy)


def element_murho(name,energies):
    # return the murho value of the wanted element @ energy(s)
    # 	;input:		name:  		    element symbol
    # 	;			energy:		    photon energy [KeV]
    # 	;return:	mu_total:	    mass absorption coefficient - total [cm^2/g]
    # 	;output:	mu_energy:	    mass absorption coefficient - energy absorption [cm^2/g]
    # 	;			density:	    density [g/cm^3]
    # 	;			atom_weight:	standard atomic weight [g/mol]
    # Some constants
    atom_mass_unit = 1.6606 * 10**(-24) # grams. atom mass unit
    radius_e = 2.8179 * 10**(-13) # cm. Classical electron radius

    energies = np.array(energies).astype(float)
    info = element_info(name)

    density = info.density
    atom_weight=info.atom_weight
    full_name = info.full_name
    atom_z    = info.atom_z
    z_per_weight = atom_z/atom_weight

    # Note: rename these value to something readable in future
    amu,emu = mu_calculator(name,energies)
    sig_coh,sig_incoh,sig_energy = sigma_calculator(name,energies)
    photo_mu = amu
    mu_total = amu + sig_coh +sig_incoh
    mu_energy= emu + sig_energy

    if len(mu_total)<=1:
        mu_total=mu_total[0]
    return(mu_total)


def composite_murho(name,energies,use_file=True):
    # 	;input:		name:  		    composite symbol
    # 	;			energies:	    photon energy [KeV]
    # 	;return:	mu_total:	    mass absorption coefficient - total [cm^2/g]
    # 	;output:	mu_energy:	    mass absorption coefficient - energy absorption [cm^2/g]
    # 	;			density:	    density [g/cm^3]
    # 	;			atom_weight:	standard atomic weight [g/mol]
    abs_path = Path(os.path.dirname(os.path.abspath(__file__)))
    with open(abs_path/Path('MU/COMPOSIT.DAT'), 'r') as file:
        content = file.read().upper()
    # Use pattern1 to find the composite we want, and get the number of types of element it has
    pattern1 = ' '+name.upper()+' .+\n'
    thereitis = re.search(pattern1,content)
    if thereitis:
        composite_outline = thereitis.group().split()
        # eg. ['K2SEO4', '1.00000', '3', 'SEO4-PH7-1.CRS']
        n_elements = composite_outline[2] #todo: Bug fix. When the name is 'SeMet', 2 is out of index.
        if len(composite_outline)==4:
            # eg. ['K2SEO4', '1.00000', '3', 'SEO4-PH7-1.CRS']
            file_present=True
            file_name = composite_outline[3]
        else: # eg. ['WATER', '1.00000', '2']
            file_present=False
    else:
        raise Exception(name.upper()+' is not a legitimate composite name in "COMPOSIT.DAT"!')
    # Use pattern2 to pack up all the information for the composite we want
    pattern2 = ' '+name.upper()+'.+(?:\n.+){'+n_elements+'}'
    # And put things into a list like this
    #[['SEO3', '1.00000', '3'],
    # ['NA', '2.0000000'],
    # ['SE', '1.0000000'],
    # ['O', '3.0000000']]
    composite_info = [i.split() for i in re.search(pattern2,content).group().split('\n')]

    # calculate mu/rho for every element
    composite_mu=[]
    mol_weight = []
    for element in composite_info[1:]:  #composite_info[0] will not be used
        element_name = element[0]
        molecular_weight = float(element[1])
        atom_weight = element_info(element_name).atom_weight
        ele_murho = element_murho(element_name,energies)
        composite_mu.append(molecular_weight*atom_weight*ele_murho)
        mol_weight.append(molecular_weight*atom_weight)
    composite_mu = np.array(composite_mu).sum(axis=0)
    mol_weight   = np.array(mol_weight).sum(axis=0)
    composite_mu = composite_mu/mol_weight

    if use_file and file_present:
        composite_mu = murho_from_file(name,file_name,energies)

    return(composite_mu)


def murho(name,energies,use_file=True,Verbose=False):
    elements = Constants.elements
    if name.upper() in elements:
        mu_total = element_murho(name,energies)
    else:
        mu_total = composite_murho(name,energies,use_file=use_file)
    if Verbose:
        print(r'murho of "'+name+'" is :', mu_total,r'cm^2/g')
    return(mu_total)


def murho_selenium_compounds(name, energies, interpol_kind='linear'):
    '''
    read in xas data for selenium compounds, generate murho values from the values in files,
    interpol to make value for interested energies
    :param name:
    :param energies:
    :param interpol_kind:
    :return:
    '''
    # path = Path(r'C:\Users\qcyus\Dropbox (X-ray Imaging Group)\IDL procedures\Spectra Selenium compounds')
    path = Path('MU/LIB')
    name = name.upper()
    if name == 'SE-METH':
        file = path /'semet-solid.CRS'
        estop = -1
        estart = -1
    elif name == 'K2SEO3':
        file = path /'seo3-ph7.CRS'
        estop = -1
        estart = -1
    elif name == 'K2SEO4':
        file = path / r'seo4-ph7-1.CRS'
        estop = -1
        estart = -1
    else:
        raise Exception("(murho_selenium_compounds) Selenium compound '" + name +
                        r"' is not found in 'Spectra Selenium compounds'.")
    data = pd.read_csv(file, delimiter=r"\s+", skiprows=1,
                       names=['energies', 'cross_over', 'normalized'])
    e1 = energies
    e2 = data.energies / 1000
    a = data.cross_over / e2  # absorbance? absorption? attenuation?
    normalized_atten = data.normalized
    n_energies = len(e2)

    murho_e2 = murho(name, e2)
    # line up the first and last value of murho and a, then we use the a values to fake the murho values
    a = murho_e2[0] + (murho_e2[-1] - murho_e2[0]) * (a - a[0]) / (a.iloc[-1] - a[0])  # type(a): pandas...series
    a = a - 0.5 * (e2 - e2.min())  # what is this used for?
    dataframe_a = pd.DataFrame.from_dict({'energy': e2, 'murho': a}, )

    murho_e1 = murho(name, e1)
    try:
        dataframe_murho_e1 = pd.DataFrame.from_dict({'energy': e1, 'murho': murho_e1})
    except:
        dataframe_murho_e1 = pd.DataFrame.from_dict({'energy': [e1], 'murho': [murho_e1]})

    df1 = dataframe_murho_e1[dataframe_murho_e1.energy < e2[0]]  # where energy < e2[first]
    df3 = dataframe_murho_e1[dataframe_murho_e1.energy > e2.iloc[-1]]  # where energy > e2[last]
    df2 = dataframe_a
    murho_all = pd.concat([df1, df2, df3], ignore_index=True)
    # in IDL, an interpol(mall,eall,es) was used. 'LSQUADRATIC or QUADRATIC or SPLINE is not set,
    # the default is to use linear interpolation'
    # In python, 2 options here.
    ## numpy.interp(x, xp, fp, left=None, right=None, period=None): only do linear interpolation
    ## scipy.interpolate.interp1d(x, y, kind='linear', axis=-1, copy=True, bounds_error=None, fill_value=nan, assume_sorted=False)
    #### this one we can choose to use 'linear' or 'spline' or 'cubic' or something else.
    interpol = interp1d(murho_all.energy.values, murho_all.murho.values, kind=interpol_kind)
    murho_interpol = interpol(e1)

    if type(e1) is int or type(e1) is float:
        murho_interpol = murho_interpol.item()
    return murho_interpol


def murho_from_file(name,file_name, energies, interpol_kind='linear'):
    """
    read in xas data for selenium compounds, generate murho values from the values in files,
    interpol to make value for interested energies.
    Equivalent to "murho_file_compound"
    :param file_name:
    :param energies:
    :param interpol_kind:
    :return:
    """

    path = Path('MU/LIB')
    file = path/file_name
    data = pd.read_csv(file, delimiter=r"\s+", skiprows=1,
                       names=['energies', 'cross_over', 'normalized'])
    e1 = energies
    e2 = data.energies / 1000
    a = data.cross_over / e2
    normalized_atten = data.normalized
    n_energies = len(e2)

    murho_e2 = murho(name, e2, use_file=False)
    # line up the first and last value of murho and a, then we use the a values to simulate the murho values
    a = murho_e2[0] + (murho_e2[-1] - murho_e2[0]) * (a - a[0]) / (a.iloc[-1] - a[0])  # type(a): pandas...series
    a = a - 0.5 * (e2 - e2.min())  # Todo: what is this used for? It makes a slope for the a.
    dataframe_a = pd.DataFrame.from_dict({'energy': e2, 'murho': a}, )

    murho_e1 = murho(name, e1, use_file=False)
    try:
        dataframe_murho_e1 = pd.DataFrame.from_dict({'energy': e1, 'murho': murho_e1})
    except:
        dataframe_murho_e1 = pd.DataFrame.from_dict({'energy': [e1], 'murho': [murho_e1]})

    df1 = dataframe_murho_e1[dataframe_murho_e1.energy < e2[0]]  # where energy < e2[first]
    df3 = dataframe_murho_e1[dataframe_murho_e1.energy > e2.iloc[-1]]  # where energy > e2[last]
    df2 = dataframe_a
    murho_all = pd.concat([df1, df2, df3], ignore_index=True)
    # in IDL, an interpol(mall,eall,es) was used. 'LSQUADRATIC or QUADRATIC or SPLINE is not set,
    # the default is to use linear interpolation'
    # In python, 2 options here.
    ## numpy.interp(x, xp, fp, left=None, right=None, period=None): only do linear interpolation
    ## scipy.interpolate.interp1d(x, y, kind='linear', axis=-1, copy=True, bounds_error=None, fill_value=nan, assume_sorted=False)
    #### this one we can choose to use 'linear' or 'spline' or 'cubic' or something else.
    interpol = interp1d(murho_all.energy.values, murho_all.murho.values, kind=interpol_kind)
    murho_interpol = interpol(e1)

    if type(e1) is int or type(e1) is float:
        murho_interpol = murho_interpol.item()
    return murho_interpol


def bragg(hkl=[1,1,1],energy=None,theta=None):
    """
    Calculator for bragg angle and bragg energy. Requires the input of either energy or theta.
    If the input is ENERGY, output will be THETA_B.
    If the input is THETA_B, output will be ENERGY.
    :param hkl: reflection lattice. Default is [1,1,1]
    :param energy: Unit [keV]
    :param theta: Unit [degree]
    :return: Energy[keV] or Theta_b[degree]
    """
    a0 = Constants.a0  # Unit: Angstroms ##silicon crystal unit cell length at 300K.
                 # This is usually used as the internal standard for silicon
    c = Constants.c
    h = Constants.h
    eV = Constants.eV
    d_hkl = a0 / np.sqrt((np.array(hkl) ** 2).sum())


    if energy is not None:
        theta=[]
        for n in [1,2,3]:
            lamb = (c * h / eV) / (energy * 1000) * (10 ** 10)  # WaveLength, Unit: Angstroms. E = h*c/lambda
            theta.append(np.degrees(np.arcsin(n*lamb / (2 * d_hkl))))  # lambda = 2d*sin(theta)
        return theta # degree

    if theta is not None:
        energy = []
        for n in [1,2,3]:
            lamb = 2*d_hkl *np.sin(np.radians(theta))/n
            energy.append( (c*h/eV)*(10**10)/(1000*lamb))
        return energy # keV


# def magic(target='',theta=None, chi = None, R=None, nu=None,f_s=None):
#     def condition(theta, chi, R, nu, f_s):
#         ''' Return the absolute difference between `f_p` and `f_g`. When `condition == 0`, the magic condtion is met.
#         '''
#         # polychromatic focus function
#         f_p = (R*np.sin(2.0*theta))/(2.0*np.sin(chi+theta)+(1+nu)*np.sin(2.0*chi)*np.cos(chi+theta))
#         # geometric focus function
#         f_g = np.cos(chi-theta)/(np.cos(chi+theta)/f_s+2.0/R)
#         # magic condition
#         condition = abs(f_p-f_g)
#         return condition
#
#     if target =='theta':
#         theta_range_a = np.arange(0.0,360.0)
#     elif target =='chi':#todo
#         pass
#     scores_range_a = condition(np.radians(theta_range_a),chi,R,nu,f_s)
#     theta_a = theta_range_a[np.argmin(scores_range_a)]
#     # from PIL import Image
#     plt.plot(theta_range_a,scores_range_a)
#     plt.show()
#     # print('a', )
#     theta_range_b = np.linspace(theta_a-2,theta_a+2,401)
#     scores_range_b = condition(np.radians(theta_range_b), chi, R, nu, f_s)
#     theta_b = theta_range_b[np.argmin(scores_range_b)]
#     # print('b')
#     theta_range_c = np.linspace(theta_b-0.02,theta_b+0.02,401)
#     scores_range_c=condition(np.radians(theta_range_c), chi, R, nu, f_s)
#     theta_c = theta_range_c[np.argmin(scores_range_c)]
#     # print('c')
#     return(theta_c)


def magic_condition(target='', theta=None, chi=None, R=None, nu=None, f_s=None, verbose=False):
    """
    Martinson, M., Samadi, N., Bassey, B., Gomez, A., & Chapman, D. (2015). Phase-preserving beam expander for
    biomedical X-ray imaging. Journal of Synchrotron Radiation, 22(3), 801–806.
    https://doi.org/10.1107/S1600577515004695

    Example: magic_condition(target='theta', chi=3.33, R=-0.5, nu=0.22, f_s=22)

    :param target: The requested magic condition variable. Choose from ['theta', 'chi', 'r']
    :param theta: Bragg angle in degree for the interested energy and the used crystal diffraction lattice.
    :param chi: Asymmetry angle (between diffraction lattice plane and the crystal surface normal) in degree.
    :param R: Crystal bending radius in meter.
    :param nu: Poisson's ratio. Assumed to be an isotropic value in the crystal bending plane.
    :param f_s: Source distance in meter.
    :return: The requested magic condition variable in a list. The useful one is usually the one with smallest absolute
             value.
    """
    def condition(theta, chi, R, nu, f_s):

        theta = np.radians(theta)
        chi = np.radians(chi)
        # polychromatic focus function
        f_p = (R * np.sin(2.0 * theta)) / (
                    2.0 * np.sin(chi + theta) + (1 + nu) * np.sin(2.0 * chi) * np.cos(chi + theta))
        # geometric focus function
        f_g = np.cos(chi - theta) / (np.cos(chi + theta) / f_s + 2.0 / R)
        # magic condition
        condition = f_p - f_g
        return condition, f_p, f_g

    def theta_study():
        theta_range_a = np.arange(-95.0, 95.2, 0.2)
        scores_range_a, f_p, f_g = condition(theta_range_a, chi, R, nu, f_s)
        magic_thetas = []
        for i in range(len(scores_range_a) - 1):
            # find solutions for `condition=0` by looking for the point crosses x-axis
            if scores_range_a[i] == 0:  # todo
                pass
            if scores_range_a[i] * scores_range_a[i + 1] <= 0 and (
                    (scores_range_a[i] - scores_range_a[i - 1]) * (scores_range_a[i + 1] - scores_range_a[i])) > 0:
                # [1]. Two points are on different sides of the x-axis. [2]. The derivative is not changing sign (
                # stays positive/negative)
                # get condition value of 10000 points in this small range, then get the minimun absolute value
                thetas = np.linspace(theta_range_a[i], theta_range_a[i + 1], 10001)
                scores, _, _ = condition(thetas, chi, R, nu, f_s)
                magic_theta = thetas[abs(scores).argmin()]
                magic_thetas.append(magic_theta)
        note = '$\\chi(degree) = $' + str(chi) + '\nBent Radius(m) = ' + str(R) + '\n$\\nu = $' + str(
            nu) + '\nSource Distance(m) = ' + str(f_s)
        if verbose:
            make_plot(theta_range_a, scores_range_a, f_p, f_g, title='Hunting Magic $\\theta_B$', note=note)
        return magic_thetas

    def chi_study():
        chi_range_a = np.arange(-95.0, 96.0)
        scores_range_a, f_p, f_g = condition(theta, chi_range_a, R, nu, f_s)
        magic_chis = []
        for i in range(len(scores_range_a) - 1):
            # find solutions for `condition=0` by looking for the point crosses x-axis
            if scores_range_a[i] == 0:  # todo
                pass
            if scores_range_a[i] * scores_range_a[i + 1] <= 0 and (
                    (scores_range_a[i] - scores_range_a[i - 1]) * (scores_range_a[i + 1] - scores_range_a[i])) > 0:
                # [1]. Two points are in different sides on the x-axis. [2]. The derivative is not changing sign (
                # stay positive/negative)
                # get condition value of 10000 points in this small range, then get the minimun absolute value
                chis = np.linspace(chi_range_a[i], chi_range_a[i + 1], 10001)
                scores, _, _ = condition(theta, chis, R, nu, f_s)
                magic_chi = chis[abs(scores).argmin()]
                magic_chis.append(magic_chi)
        note = '$\\theta_B(degree) = $' + str(theta) + '\nBent Radius(m) = ' + str(R) + '\n$\\nu = $' + str(
            nu) + '\nSource Distance(m) = ' + str(f_s)
        if verbose:
            make_plot(chi_range_a, scores_range_a, f_p, f_g, title='Hunting Magic $\\chi$', note=note)
        return magic_chis

    def r_study():  # crystal bent radius
        r_range_a = np.append(np.arange(-10.0, -0.1, 0.1), np.arange(0.1, 10.1, 0.1))
        scores_range_a, f_p, f_g = condition(theta, chi, r_range_a, nu, f_s)
        magic_rs = []
        for i in range(len(scores_range_a) - 1):
            # find solutions for `condition=0` by looking for the point crosses x-axis
            if scores_range_a[i] == 0:  # todo
                pass
            if scores_range_a[i] * scores_range_a[i + 1] <= 0 and (
                    (scores_range_a[i] - scores_range_a[i - 1]) * (scores_range_a[i + 1] - scores_range_a[i])) > 0:
                # [1]. Two points are in different sides on the x-axis. [2]. The derivative is not changing sign (
                # stay positive/negative)
                # get condition value of 10000 points in this small range, then get the minimun absolute value
                rs = np.linspace(r_range_a[i], r_range_a[i + 1], 10001)
                scores, _, _ = condition(theta, chi, rs, nu, f_s)
                magic_r = rs[abs(scores).argmin()]
                magic_rs.append(magic_r)
        if len(
                magic_rs) > 1:  # Get rid of the match at 0. When the bent radius is 0, f_p and f_g are both 0. And
            # we dont need this.
            magic_rs = np.array(magic_rs)  # transition to numpy array
            to_delete = abs(magic_rs).argmin()  # get the index of the '0' value
            magic_rs = np.delete(magic_rs, to_delete).tolist()  # delete the '0' value and transition back to list

        note = '$\\theta_B(degree) = $' + str(theta) + '\n$\\chi$(degree) = ' + str(chi) + '\n$\\nu = $' + str(
            nu) + '\nSource Distance(m) = ' + str(f_s)  # the annotation content
        #         make_plot(r_range_a,scores_range_a,f_p,f_g,title='Hunting Best Bent Radius',note=note,ylim=10)
        if verbose:
            plt.style.use('seaborn')
            fig = plt.figure(figsize=(18, 9))
            plt.plot(r_range_a, scores_range_a, linestyle='--', color='r', label='Distance between P & G')
            #         plt.plot(r_range_a,np.zeros(len(xdata)),color='k',alpha=0.5)# straight line on x-axis
            plt.xticks(np.linspace(r_range_a[0], r_range_a[-1], 21), fontsize=12)
            plt.yticks(fontsize=12)
            plt.xlabel('Crystal Bent Radius', fontsize=14, position=(1.0, 1.0))
            plt.ylabel('Meter', fontsize=14, rotation=0, position=(1.0, 1.01))
            plt.plot(r_range_a, f_p, label='Polychromatic Focus', color='b')
            plt.plot(r_range_a, f_g, label='Geometric Focus', color='g')
            ylim = 1
            plt.ylim(-ylim, ylim)
            #         plt.xlim(-1,1)
            plt.legend(prop={'size': 16}, loc=1)
            plt.title('Hunting Magic Bent Radius', fontsize=18)
            xytext = [r_range_a[-1] * 2 / 3, -ylim * 1.6 / 2.0]
            plt.annotate(note, xy=(0, 0), xytext=xytext, size=16)

            # A small plot window
            subplot = fig.add_axes([0.17, 0.2, 0.2, 0.2])
            subplot.set_facecolor('white')
            subplot.plot(r_range_a, scores_range_a, color='r', linestyle='--', alpha=0.5, label='P to G')
            subplot.plot(r_range_a, f_p, color='b', alpha=0.5, label='Polychromatic')
            subplot.plot(r_range_a, f_g, color='g', alpha=0.5, label='Geometric')
            subplot.grid(color='gray', alpha=0.1)
            subplot.set_xticks(np.linspace(r_range_a[0], r_range_a[-1], 7))
            subplot.legend()
            subplot.set_title('Full View')
            plt.show()

        return magic_rs

    def make_plot(xdata, score, f_p, f_g, title='', note='', ylim=2):
        plt.style.use('seaborn')
        fig = plt.figure(figsize=(18, 9))
        plt.plot(xdata, score, linestyle='--', color='r', label='Distance between P & G')
        plt.plot(xdata, np.zeros(len(xdata)), color='k', alpha=0.5)  # straight line on x-axis
        plt.xticks(np.linspace(xdata[0], xdata[-1], 39), fontsize=12)
        plt.yticks(fontsize=12)
        plt.xlabel('Degree', fontsize=14, position=(1.0, 1.0))
        plt.ylabel('Meter', fontsize=14, rotation=0, position=(1.0, 1.01))
        plt.plot(xdata, f_p, label='Polychromatic Focus', color='b')
        plt.plot(xdata, f_g, label='Geometric Focus', color='g')
        plt.ylim(-ylim, ylim)
        plt.legend(prop={'size': 16})
        plt.title(title, fontsize=18)
        plt.annotate(note, xy=(0, 0), xytext=[60, -1.6], size=16)

        # A small plot window
        subplot = fig.add_axes([0.17, 0.2, 0.2, 0.2])
        subplot.set_facecolor('white')
        subplot.plot(xdata, score, color='r', linestyle='--', alpha=0.5, label='P to G')
        subplot.plot(xdata, f_p, color='b', alpha=0.5, label='Polychromatic')
        subplot.plot(xdata, f_g, color='g', alpha=0.5, label='Geometric')
        subplot.grid(color='gray', alpha=0.1)
        subplot.set_xticks(np.linspace(xdata[0], xdata[-1], 7))
        subplot.legend()
        subplot.set_title('Full View')
        plt.show()

    if target.upper() == 'THETA' or theta == None:
        target = 'theta'
        unit = 'deg'
        mc = theta_study()
    elif target.upper() == 'CHI' or chi == None:
        target = 'chi'
        unit = 'deg'
        mc = chi_study()
    elif target.upper() == 'R' or R == None:
        target = 'R'
        unit = 'm'
        mc = r_study()
    else:
        raise Exception('Choose your target from ["theta","chi","R"]')
    print('The magic condition is met when\n{} is {} {}.'.format(target,mc,unit))
    return mc


def relative_energy_resolution(theta,chi, R, nu, f_s, t=0.001):
    """
    This function is used to calculate the energy resolution of the magic condition monochromator. Note that this
    is not suitable for calculating general bent crystal energy resolution.
    :param theta: Bragg angle in degree for the interested energy and the used crystal diffraction lattice.
    :param chi: Asymmetry angle (between diffraction lattice plane and the crystal surface normal) in degree.
    :param R: Crystal bending radius in meter.
    :param nu: Poisson's ratio. Assumed to be an isotropic value in the crystal bending plane.
    :param f_s: Source distance in meter.
    :param t: Crystal thickness in meter.
    :return: dE/E, the relative energy resolution.
    """
    import blxo
    mono = blxo.mc.BentLaueMono(chi=np.radians(chi),
                                theta=np.radians(theta),
                                nu=nu,
                                t = t*1000,
                                r = R*1000,
                                p = f_s*1000)
    return mono.energy_resolution['de_all']


def quasi_mono_beam_width(theta,chi, R, nu, f_s, t):
    """
    This function is used to calculate the quasi-mono-beam-width of a bent Laue monochromator.
    Reference: Qi, P., Samadi, N. & Chapman, D. (2020). New insights into the focusing and energy dispersion properties
               of bent laue crystals. Manuscript in preparation.
    :param theta: Bragg angle in degree for the interested energy and the used crystal diffraction lattice.
    :param chi: Asymmetry angle (between diffraction lattice plane and the crystal surface normal) in degree.
    :param R: Crystal bending radius in meter.
    :param nu: Poisson's ratio. Assumed to be an isotropic value in the crystal bending plane.
    :param f_s: Source distance in meter.
    :param t: Crystal thickness in meter.
    :return: qmbw in mm.
    """
    import blxo
    mono = blxo.mc.BentLaueMono(chi=np.radians(chi),
                                theta=np.radians(theta),
                                nu=nu,
                                t = t*1000,
                                r = R*1000,
                                p = f_s*1000)
    return mono.qmb['width']
