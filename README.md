# X-ray Spectral Imaging Program (XSIP)

A Python-based data analysis software for two imaging techniques: **spectral K-Edge Subtraction imaging (spectral KES)** and **Wide-Field Energy Dispersive X-ray Absorption Spectroscopy (Wide-Field EDXAS)**.

The core algorithm of this Python version software package is based on programs developed with **IDL** by **Prof. L. Dean Chapman**. 

## Installation

### 1. Prepare the Python environment

#### 1.1 Install Python

(If you already have Python3 (>3.6) installed on your machine, you can skip this step.)

##### Option 1. `Miniconda` (recommended)

Install the `Miniconda` Python distribution by following the official instructions [Miniconda — Conda documentation](https://docs.conda.io/en/latest/miniconda.html). 

##### Option 2. Python.org

[BeginnersGuide/Download - Python Wiki](https://wiki.python.org/moin/BeginnersGuide/Download). 

#### 1.2 Prepare a virtual environment

Create the virtual environment

```bash
conda create -n xsip python=3.6
```

Note that the `xsip` in the command above can be replaced with a name you give. The `xsip` will be used in the instruction as an example.

### 2. Install XSIP

#### 2.1 Download the repository and unzip the downloaded file

<img src="https://github.com/darwinqii/XSIP/blob/master/.figure/README/image-20210216151210487.png?raw=true" alt="Download" style="zoom:80%;" />

#### 2.2 Prepare the virtual environment

##### Activate the created virtual environment.

```bash
conda activate xsip
```

##### Install dependencies.

```
pip install -r path/to/setup/xsip_requirements.txt
```

The `xsip_requirements.txt` file is located in the `setup` folder in the downloaded repository. Use the full path to `xsip_requirements.txt` to replace the path in the command above. 

##### Add the path to the `XSIP` package to the virtual environment.

- Find the path to the `Conda` virtual environment folder by using the following command.

  Windows: `where python`

  Linux & MacOS: `which python`

  It will show you an output that looks like this on Windows

  ```
  C:\Users\user\miniconda3\envs\xsip\python.exe
  ```

  or this on Linux

  ```
  /home/user/miniconda3/envs/xsip/bin/python
  ```

- Use the part of the path to `xsip` (e.g. `..\user\miniconda3\envs\xsip`) and find the folder `site-packages` for available Python packages

  On Windows:

  ```
  C:\Users\user\miniconda3\envs\xsip\Lib\site-packages
  ```

  On Linux:

  ```bash
  /home/user/miniconda3/envs/xsip/lib/python3.6/site-packages
  ```

- Create a file with the file name `xsip.pth`, where the `xsip` can be actually anything, `.pth` is mandatory. The content of the file should be the path to the downloaded `XSIP` repository. For example, 

  ```bash
  C:\Users\user\XSIP
  ```


#### 2.3 Test the installation

Navigate to the `XSIP` folder in your terminal and enter the following command.

```bash
python gui.py
```

You will see the GUI of the XSIP if the installation is successful.

<img src="https://github.com/darwinqii/XSIP/blob/master/.figure/README/image-20210216160837492.png?raw=true" alt="GUI" style="zoom:50%;" />

## How to use the program

### The near edge imaging function: `nei()`

`nei()` is the main function that goes through everything. 

#### Example

```python
import xsip
result = xsip.nei(materials=['Na2SeO4', 'Na2SeO3', 'Se-Meth', 'Water'],
                  data_path='/directory/of/the/imaging/data',
                  multislice=True,  # Whether the imaging data is from a multislice scan
                  slice=0,  # If `multislice==True`, provide the number of slice to analyze (starting from 0)
                  n_proj=900,  # The number of projection images per slice
                  ct=True,  # Whether this is a CT scan. If `True`, `side_width` will be used.
                  side_width=20,  # The number of pixels used on the side.
                  e_range=0,  # The interested energy range for analysis. `0` for all available energies.
                  lowpass=False,  # If `True`, apply a lowpass filter to reduce high frequency noise.
                  save=False,  # Save the result or not (because the returned result is usually large)?
                  Verbose=False)  # If `True` (not suggested for general user), the program generates some figure during the data processing.
```



#### Parameters

- **materials**: [name1, name2, ...].Name of each material. eg.: `materials = ['Water', 'Bone', 'K2SeO4', 'U'   ]` 
- **path**: The main directory containing Flat, Dark, Edge, Tomo, etc..
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

### Note

#### Data folder architecture

[Todo]



## Spectral K-Edge Subtraction Imaging

The same as 'EDXAS'.

# A Graphic User Interface

Built with `tkinter` and `tkinter.ttk`.

It has two `ttk.Notebook` tabs. 
- Tab1 is made for 'spectral KES' or 'Energy Dispersive XAS'. First, setup the parameters. Second, choose directorys (1)containing the data and (2) to save the results. Results is generated from `nei()`. Saved as `.pkl` file and images.
- Tab2 is made for CT reconstruction. The source should be sinogram(s)-like array or image. 
