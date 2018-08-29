# Energy Dispersive Xray Absorption Spectroscopy
A python version of data analysis for my phd project: selenium speciation with spectral kes imaging method.
Most of the analysis programs had been written with IDL. This repository aims to 
- (1) reproduce things have been done with IDL, and 
- (2) extend it for furthur analysis.

# How to
## Run EDXAS
- `nei()` is the main function that goes through everything.
Returns:
    - beam_parameters:
        - .beam_files
            - .flat
            - .dark
            - .edge
        - .theta_b : bragg angle
        - .beam
        - .edges
            - .top
            - .bottom
            - .peak
            - .edge
        - .mu_t : `-np.log(edge-dark/flat-dark)` of the edge_image
        - .pixel_edge_width: gaussian edge width in pixel
        - .e_width: gaussian edge width in energy
        - edge_slope
        - peak_slope
        - exy : the energies at [y,x] locations
    - mu_rhos: A dictionary, with keys() are the names of materials, values() are the $\mu/\rho$ for every material at every [y,x] location.
    - mu_t: ndarray [n_projections, n_energies, n_horizontal_position]. $-\ln{\frac{tomo-dark}{flat-dark}}$. `nan` values from any illegal $ln$ or any divisive operation(s) are replaced by 0s. Optional corrections by following keywords:
        - "lowpass": Right now it is a Gaussian filter on the axis of energy. Default $\sigma = beam_parameters.pixel_edge_width$
        - "ct": Remvoe absoroption from air, by using the left and right sections where there is no sample in the image.
        - More keywords coming soon... 
    - rts: $\rho\cdot t$ for every material at every horizontal position in every projection. Dictionary: keys() are the names of the materials, values() are ndarrays [n_projections,n_horizontal_position].
    
