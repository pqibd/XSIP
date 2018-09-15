# Work in progress

# Energy Dispersive Xray Absorption Spectroscopy
A python version of data analysis programs for my phd project: selenium speciation with spectral KES imaging method.
Most of the analysis programs had been written with IDL by **Dr. L. Dean Chapman**. This repository aims to 
- (1) reproduce things have been done with IDL 
- (2) optimize the computation and extend it for furthur analysis.

# How to
## Run EDXAS
### Function "nei()"
- `nei()` is the main function that goes through everything.<br>
Returns:
    - **beam_parameters**:
        - .beam_files
            - .flat
            - .dark
            - .edge
        - .beam
        - .edges
            - .top
            - .bottom
            - .peak
            - .edge
        - .mu_t : `-np.log(edge-dark/flat-dark)` of the edge_image
        - .pixel_edge_width: gaussian edge width in pixel
        - .e_width: gaussian edge width in energy
        - .edge_slope
        - .peak_slope
        - .exy : the energies at [y,x] locations
    - **mu_rhos**: A dictionary, with keys() are the names of materials, values() are the $\mu/\rho$ for every material at every [y,x] location.
    - **mu_t**: ndarray [n_projections, n_energies, n_horizontal_position]. $-\ln{[\frac{tomo-dark}{flat-dark}]}$. `nan` values from any illegal $\ln$ or divisive operation(s) are replaced by 0. Optional corrections for **mu_t** with following keywords:
        - "lowpass": Right now it is a Gaussian filter on the axis of energy. Default $\sigma$ = beam_parameters.pixel_edge_width
        - "ct": Remove air absoroption, by using the left and right sections where there is no sample in the image.
        - More keywords coming soon... 
    - **rts**: $\rho\cdot t$ for every material at every horizontal position in every projection (A sinogram ready for CT reconstruction). Dictionary: `keys()` are the names of the materials, `values()` are ndarrays [n_projections,n_horizontal_position].
    - **snrs**: "signal to noise raito". Numpy array, in shape of [n_materials,n_projections,n_horizontal_positions]
    - **recons**: Reconstruction images. Numpy array, in shape of [n_material,n_horizontal,n_horizontal]
    - **mean_rhos**: The mean values of $\rho$ in the target area in recon image. Unit: $mg/cm^3$.

### Arguments in "nei()"
- **materials**: {names:sources}. Sources mean the way in which we will get the $\mu/\rho$ for that material
    - names: Name of each material, i.e. 'K2SeO4'
    - sources: how we get the $\mu/\rho$ for that element or compound
        - **system**   : get it from calculation with known chemical information. `element_murho` and `composite_murho`
        - **file** : get it from a file
        - **standard**: get murho from experiment standard. Stardard data are collected with current experiment setting, and standard selenium compound solution ,etc.
        
    For example:
    
    >names = ['Water', 'Bone', 'Selenite', 'U'   ]<br> 
    sources = [  'system',  'system',     'file', 'standard' ]
 
- **path**: The main directory containing Flat, Dark, Edge, Tomo, etc...
- **algorithm**: The algorithm used to calculate $\rho t$. Options are:
                      'sKES_equation': Default option. A equation derived with least-square approach is used.
                                       Much faster than 'nnls'. [Ref: Ying Zhu,2012 Dissertation]
                      'nnls'         : A Non-Negative Linear Regression will be performed with `scipy.optimize.nnls`.
- **ct**: If True, a piece of left and right side of projection image will be used to correct the air absorption from sample to detector.
- **side_width**: Used with param "ct". Define the width in pixel for air absorption correction
- **n_proj**: The number of projection images for one slice of CT imaging.
- **multislice**: If True, meaning the images in the "tomo" folder contain more than one slice of CT. The 'n_proj' and 'slice' needs to be specified.
- **slice**: Which slice do we want to do the reconstruction.
- **e_range**: The energy range we want to use. Default 0, meaning the "energy_range" in "arrangement.dat" file will be used as the energy range. If not 0, this will overwrite the energy range from "arrangement.dat".
- **lowpass**: Use a lowpass filter(gaussian) on the $\mu t$ from experiment. Default is False for now(20180905)
- **use_torch**: use Pytorch.tensor instead of numpy.array for matrix operations. Default True.
- **snr**: Calculate the signal to noise ratio. Default False.
- **reconstruction**: str (default=None). Routine used for CT reconstruction after having the sinograms. Routines available: 
    - 'idl'
    - 'skimage': skimage.transform.iradon. An edited version.
- **ct_center**: Specify the rotation center for CT reconstruction if needed. Default is 0.
- **fix_vertical_motion**: Todo.
- **fix_cross_over**: Todo. May be not needed.
- **flat_gamma**: Todo. May be not needed.
- **Verbose**: If True, some detail will show up when run the program. And some matplotlib plot window might pause the program.

## Get beam properties
- `nei_beam_parameters.get_beam_parameters`. It returns the same "beam_parameters" as the one returned by `nei()`.

## Spectral K-Edge Subtraction Imaging



    
