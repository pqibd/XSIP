# Energy Dispersive X-ray Absorption Spectroscopy
A python version of data analysis software for my Ph.D. project: selenium speciation with spectral KES imaging method.
The basic algorithms were from the IDL programs written by **Dr. L. Dean Chapman**.

# How to

## Run EDXAS

### Function `nei()`

`nei()` is the main function that goes through everything.

#### Parameters
- **materials**: {name1, name2, ...}.Name of each material. eg.: `materials = ['Water', 'Bone', 'K2SeO4', 'U'   ]` 
- **path**: The main directory containing Flat, Dark, Edge, Tomo, etc...
- **algorithm**: The algorithm used to calculate $\rho t$. Options are:
  - 'sKES_equation': Default option. A equation derived with least-square approach is used. Much faster than 'nnls'. [Ref: Ying Zhu,2012 Dissertation]
  - 'nnls': A Non-Negative Linear Regression will be performed with `scipy.optimize.nnls`.
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
#### Returns
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
    - .exy : the energies at [y,x] locations on the detector.
- **mu_rhos**: A dictionary, with `keys()` is the names of materials, `values()` is the $\mu/\rho$ for every material at every [y,x] location.
- **mu_t**: ndarray [n_projections, n_energies, n_horizontal_position]. $-\ln{[\frac{tomo-dark}{flat-dark}]}$. `nan` values from any illegal $\ln$ or divisive operation(s) are replaced by 0. Optional corrections for **mu_t** with following keywords:
- **rts**: $\rho\cdot t$ for every material at every horizontal position in every projection (A sinogram ready for CT reconstruction). Dictionary: `keys()` are the names of the materials, `values()` are ndarrays [n_projections,n_horizontal_position].
- **snrs**: "signal to noise raito". Numpy array, in shape of [n_materials,n_projections,n_horizontal_positions]
- **recons**: Reconstruction images. Numpy array, in shape of [n_material,n_horizontal,n_horizontal]
- **mean_rhos**: The mean values of $\rho$ in the target area in recon image. Unit: $mg/cm^3$.


## Spectral K-Edge Subtraction Imaging

The same as 'EDXAS'.

# A Graphic User Interface

Built with `tkinter` and `tkinter.ttk`.

It has two `ttk.Notebook` tabs. 
- Tab1 is made for 'spectral KES' or 'Energy Dispersive XAS'. First, setup the parameters. Second, choose directorys (1)containing the data and (2) to save the results. Results is generated from `nei()`. Saved as `.pkl` file and images.
- Tab2 is made for CT reconstruction. The source should be sinogram(s)-like array or image. 