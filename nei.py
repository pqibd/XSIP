from near_edge_imaging import *
import numpy as np
import tkinter
import tkinter.filedialog as filedialog
from beam_near_edge_imaging import beam_near_edge_imaging
from nei_beam_parameters import nei_beam_parameters


def nei(materials='', path='', left=0, right=0, n_proj=900,
        slice=0, multislice=False, projection=False,
        display=True, pop_up_image=False, setup_type='FILE',
        order_files=False, e_range=0,
        fix_vertical_motion=False,  # maybe change default to True
        clip=False, flip=False, fix_cross_over=False,
        use_se_data=False, use_measured_standard=False,
        use_weights=False, energy_weights=0, flat_gamma=1.0,
        put_dark_back=False, fix_detector_factor=0,
        Verbose=False):
    """
    Get beam_parameters, then do the beam_near_edge_imaging. Get $\mu t$ for all
    projection images.
    :param materials:
    :param path:
    :param left:
    :param right:
    :param n_proj:
    :param slice:
    :param multislice:
    :param projection:
    :param e_range:
    :param fix_vertical_motion:
    :param fix_cross_over:
    :param flat_gamma:
    :param Verbose:
    :return:
    """

    ###############   define materials       ######################
    if materials == '':
        names = ['K2SeO4', 'K2SeO3', 'Se-Meth', 'Water']
        sources = ['FILE', 'FILE', 'FILE', 'SYSTEM']
        materials = {}
        for i in range(len(names)):
            materials[names[i]] = sources[i]

    ##############    get  Path for experiment data file  #######################

    if path == '':
        def choosePath():
            root = tkinter.Tk()
            root.withdraw()
            path = filedialog.askdirectory(title='Please select NEI directory:')
            if path == '': choosePath()
            return (path)

        path = choosePath()
    print(path)

    #############  get system setup info from arrangement.dat   ##########
    setup = nei_get_arrangement(setup_type, path)
    detector = setup.detector
    # redefine energy_range if needed
    if e_range != 0: setup.energy_range = e_range

    # process like a projection unless left and right are set
    if left == 0 & right == 0: projection = True

    ########  get beam files: averaged flat, dark, and edge  ############
    print('\n(nei) Calling "get_beam_files"')
    beam_files = get_beam_files(path=path, clip=clip, flip=flip, Verbose=Verbose)

    ######### get tomo data  ########################################
    print('(nei) Calling "get_tomo_files"')
    tomo_files = get_tomo_files(path)
    print('(nei) tomo_files loaded')

    # get the tomo for ONE slice
    if multislice == True:
        if slice < 1:
            print('------------------------'
                  '\nWarning: "Slice" starts from 1. "slice=' + str(slice) + ' is entered\n'
                  '------------------------')
        i_begin = n_proj * (slice - 1)
        i_end = i_begin + n_proj
        tomo_files = tomo_files[i_begin:i_end]
        print('(nei) Tomo files in 1 slice loaded')
    n_tomo_files = len(tomo_files)
    print('(nei) Number of Tomo files: ', n_tomo_files)  # equal to n_projections
    # tomo files to data array
    from PIL import Image
    tomo_data = []
    for i in range(n_tomo_files):
        tomo_data.append(np.array(Image.open(tomo_files[i])))
    tomo_data = np.array(tomo_data)

    ########### Get beam_parameters ##############################
    print('\n(nei) Calling "nei_beam_parameters"\n')

    beam_parameters = nei_beam_parameters(display=display, beam_files=beam_files,
                                              setup=setup, detector=detector,
                                              fix_vertical_motion=fix_vertical_motion,
                                              clip=clip, Verbose=Verbose)

    ###################  Display beam parameters  #####################
    # if display beam parameters:
    #   display beam parameters

    ###################  beam_near_edge_imaging  ######################
    print('\n(nei) Calling "beam_near_edge_imaging"\n')
    mu_rhos, rts = beam_near_edge_imaging(tomo_data, beam_parameters, materials,
                                          projection=projection, fix_vertical_motion=fix_vertical_motion,
                                          fix_cross_over=fix_cross_over)

    # Available keywords: tomo_data, beam_parameters, materials, projection=False, flat_gamma=1.0,
    #                            fix_vertical_motion=False, fix_cross_over=False, left=0, right=0,
    #                            flip=False, width_factor=1,use_sm_data=False, use_weights=False,
    #                            use_measured_standard=False, *xover
    class ProjResult:
        def __init__(self, beam_parameters, mu_rhos, rts):
            self.beam_parameters = beam_parameters
            self.mu_rhos = mu_rhos
            self.rts = rts

    return ProjResult(beam_parameters, mu_rhos, rts)
