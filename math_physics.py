import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re
from scipy.interpolate import interp1d


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


def fwhm(x,y,Verbose=False):
    '''
    :param x:
    :param y:
    :param Verbose:
    :return:
    '''
    ## If max is too close to either end, raise error
    ind_max = y.argmax()
    if ind_max<=2 or ind_max>=len(y)-2:
        # raise Exception('(In fwhm.py,) Peak occurs at extreme limits of x range ... returning')
        return(0,0,0)
    half_max  = y.max()/2
    ind_max   = y.argmax()
    low_index = np.where(y<half_max)[0]
    left      = low_index[low_index<ind_max]
    right     = low_index[low_index>ind_max]
    if len(left)*len(right)==0:
        return(0,0,0)
    left_fwhm = left[-1]
    right_fwhm= right[0]
    fwhm      = right_fwhm-left_fwhm

    # show the plot if Verbose==True
    if Verbose==True:
        plt.figure()
        plt.scatter(x,y,s=3)
        plt.plot([left_fwhm,left_fwhm],[0,y.max()],color='y')
        plt.plot([right_fwhm,right_fwhm],[0,y.max()],color='y')
        plt.title('FWHM')
        plt.show()
    return (fwhm,left_fwhm,right_fwhm)


def gaussfit(x, y, *estimate):
    '''

    :param x:
    :param y:
    :param estimate:
    :return:
    '''

    def gauss_func(x, a0, a1, a2):  # ,a3=0,a4=0,a5=0):
        # define the gaussian function
        # z  = (x-a1)/a2
        # ez = np.exp(-z**2/2.0) # Gaussian part
        return a0 * np.exp(-(x - a1) ** 2 / (2 * a2 ** 2))  # +a3+a4*x+a5*x**2

    # Call the curve_fit, imported from scipy.optimize
    if estimate:
        popt, pcov = curve_fit(gauss_func, x, y, estimate)
    else:
        popt, pcov = curve_fit(gauss_func, x, y)
    y_gauss = gauss_func(x, popt[0], popt[1], popt[2])
    return y_gauss, popt


def polyfit(x, y, degree):
    '''
    use the numpy.polyfit method, but directly returns the y_values for the fitted polynomial function
    :param x:
    :param y:
    :param degree:
    :return:
    '''
    coef = np.polyfit(x, y, degree)
    p = np.poly1d(coef)
    y_poly = p(x)
    return (y_poly)


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


def read_absorber(element_name,Verbose=False):
    element_name = element_name.upper()
    if Verbose:
        print('(read_absorber) Looking for "'+element_name+'" in absorber.dat')
    # open 'absorber.dat' file, get to the interested element by regex,
    # get the ek, eta, ef..
    with open(r'C:\Users\qcyus\Dropbox (X-ray Imaging Group)'
              r'\IDL procedures\Dean Procedures\local\MU\ABSORBER.DAT','r') as file:
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
        de[i]= float(ms[n_edges+n_edges//6+2+6*i][0])
        c_a[i,:] = ms[n_edges+n_edges//6+2+6*i][1:9]
        c_b[i:,] = ms[n_edges+n_edges//6+3+6*i]+ms[n_edges+n_edges//6+4+6*i]+ms[n_edges+n_edges//6+5+6*i]\
                   +ms[n_edges+n_edges//6+6+6*i]
        e0[i]    = float(ms[n_edges+n_edges//6+7+6*i][0])
        xj[i,:]  = ms[n_edges+n_edges//6+7+6*i][1:7]
    c_a = np.array(c_a).astype(float)
    c_b = np.array(c_b).astype(float)
    xj  = np.array(xj).astype(float)

    edges=[]
    for i in range(1,(2+n_edges//6)):
        edges=edges+ms[i]
    edges = np.array(edges).astype(float)
    class GotAbsorber:
        def __init__(self,ms,de,c_a,c_b,e0,xj):
            self.n_edges = int(ms[0][1])
            self.edges   = edges
            self.eta = float(ms[2+n_edges//6][0])
            self.ef  = float(ms[2+n_edges//6][1])
            self.ek  = float(ms[2+n_edges//6][2])
            self.za  = float(ms[2+n_edges//6][3])
            self.a   = np.array(ms[2+n_edges//6+1:n_edges+2+n_edges//6]).astype(float)
            self.de  = de
            self.c_a = c_a
            self.c_b = c_b
            self.e0  = e0
            self.xj  = xj

    return(GotAbsorber(ms,de,c_a,c_b,e0,xj))


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
    # if len(info)==0:
    #     # then this is probably a composite
    #     info = composite(name)
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


def composite_murho(name,energies):
    # 	;input:		name:  		    composite symbol
    # 	;			energies:	    photon energy [KeV]
    # 	;return:	mu_total:	    mass absorption coefficient - total [cm^2/g]
    # 	;output:	mu_energy:	    mass absorption coefficient - energy absorption [cm^2/g]
    # 	;			density:	    density [g/cm^3]
    # 	;			atom_weight:	standard atomic weight [g/mol]

    with open(r'C:\Users\qcyus\Dropbox (X-ray Imaging Group)'
              r'\IDL procedures\Dean Procedures\local\MU\COMPOSIT.DAT','r') as file:
        content = file.read().upper()
    # Use pattern1 to find the composite we want, and get the number of types of element it has
    pattern1 = ' '+name.upper()+'.+\n'
    thereitis = re.search(pattern1,content)
    if thereitis:
        n_elements = thereitis.group().split()[-1]
    else:
        raise Exception(name.upper()+' is not a legitimate composite name!')
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
    for element in composite_info[1:]:
        element_name = element[0]
        molecular_weight = float(element[1])
        atom_weight = element_info(element_name).atom_weight
        ele_murho = element_murho(element_name,energies)
        composite_mu.append(molecular_weight*atom_weight*ele_murho)
        mol_weight.append(molecular_weight*atom_weight)
    composite_mu = np.array(composite_mu).sum(axis=0)
    mol_weight   = np.array(mol_weight).sum(axis=0)
    composite_mu = composite_mu/mol_weight

    return(composite_mu)


def murho(name,energies,Verbose=False):
    elements = Constants.elements
    if name.upper() in elements:
        mu_total = element_murho(name,energies)
    else:
        mu_total = composite_murho(name,energies)
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
    from scipy.interpolate import interp1d
    path = r'C:\Users\qcyus\Dropbox (X-ray Imaging Group)\IDL procedures\Spectra Selenium compounds'
    name = name.upper()
    if name == 'SE-METH':
        file = path + r'\semet-solid.CRS'
        estop = -1
        estart = -1
    elif name == 'K2SEO3':
        file = path + r'\seo3-ph7.CRS'
        estop = -1
        estart = -1
    elif name == 'K2SEO4':
        file = path + r'\seo4-ph7-1.CRS'
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
    a = a - 0.5 * (e2 - e2.min())  # what is this used for???
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

