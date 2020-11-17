import pickle
import tkinter
import tkinter.filedialog as filedialog
import os
from pathlib import Path
import fnmatch
import imageio
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from scipy import misc


def plot(x):
    plt.figure()
    plt.plot(x)
    plt.show()

def imshow(x):
    plt.figure()
    plt.imshow(x)
    plt.show()


def choose_path(title = 'Please select data directory:'):
    """Pop-up a dialog window (using `tkinter` package) to select the needed folder, and return the path as a `string`.
    Parameters
    ----------
    None
    Returns
    -------
    path : str
        Path to the selected **folder**
    """
    root = tkinter.Tk()
    root.withdraw()
    path = filedialog.askdirectory(title=title)
    root.destroy() # this command is very important. If not destroyed,
                    # "tk" windows in the following code will be affected
    return path


def choose_file():
    """Pop-up a dialog window (using `tkinter` package) to select the needed file, and return the path as a `string`.
    Parameters
    ----------
    None
    Returns
    -------
    file : str
        Path to the selected **file**.
    """
    root = tkinter.Tk()
    root.withdraw()
    file = filedialog.askopenfilename(title='Please select file')
    return file


def file_search(path,filter='*'):
    '''This function searchs for files match provided filter in provided folder, and returns the absolute path to these files with a list of `pathlib.Path` objects.
    Parameters
    ----------
    - path : pathlib object or str
        Search files in this folder
    - filter : str
        Pattern for `fnmatch.filter` to look for target files. Default '*' returns all files regardless of the file type.
    Returns
    -------
    files : List of `pathlib` objects to files in pathlib format
    '''
    pathP = Path(path) # making sure it is pathlib class object, and use it later.
    files = fnmatch.filter(os.listdir(pathP),filter) # find files that matches the filter
    if len(files)==0:
        print('------------------------'
              '\nWarning: No files found\n'
              '------------------------')
        return(files)
    # sort file names in alphabetic order. In Linux, the `fnmatch` does not return
        # results in order by itself
    files.sort()
    # attach filenames with full path, and convert back to pathlib object
    files = [pathP/file for file in files]
    return files


def save_object(obj, fname):
    """This function is used for saving multi-level class object. `pickle` by itself cannot save class object. `class_to_dict` function is called to convert **Class object** to nested **dictionary**. 
    Parameters
    ----------
    obj : class object. 
        Class object after data analysis.
    fname : str
        File name to save the class object.
    Returns
    -------
    None
    """
    fname = str(fname) # In case fname is an pathlib class, convert to string
    # print(fname)
    with open(fname, "wb") as output:
        obj = class_to_dict(obj)
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
    return


def class_to_dict(obj):
    """Convert class object to dictionary with the same structure. 

    If the child of the parent class is still a class object, it will be converted to dictionary as well until there is no class object.
    The limit of this function is that, if there are class objects hide in non-class object, this function is not going to dig them out and do the convertion.
    Parameters
    ----------
    obj : Be careful. This function modifies the variable passed in. Please make sure to use `copy.deepcopy()` before sending the object in to avoid messing up the original object.
    Returns
    -------
    obj : Nested dictionary
    """
    try:
        obj = obj.__dict__
        for key, value in obj.items():
            obj[key] = class_to_dict(value)
        return obj
    except:
        return obj


def load_object(fname=''):
    """This function is used for loading in '.pkl' file or image file for analysis.
    Parameters
    ----------
    fname : str, default ''
        The absolute path to the saved 'pickle' file you want to load for analysis.
    Returns
    -------
    obj : class object or image arrays 
    """
    if fname=='': # if fname is not assigned, pop up a window to pick file.
        fname = choose_file()
    
    with open(fname, 'rb') as input:
        try: # If it is .pkl
            obj = pickle.load(input)
            if isinstance(obj, dict):
                obj = dict_to_class(obj)
        except: # if it is an image
            obj = imageio.imread(input)
    return obj


def dict_to_class(adict):
    """This function is used for reloading previously saved '.pkl' file to class object for the convenience of analysis. It is called by `load_object` function.
    Parameters
    ----------
    adict : dict
        Usually a nested dictionary previously saved to store results after data processing.
    Returns
    -------
    Func(adict) : class object
        class object with the same structure of the input class object
    """
    class Func(object):
        def __init__(self, d):
            for a, b in d.items():
                if isinstance(b, (list, tuple)):
                    setattr(self, a, [Func(x) if isinstance(x, dict) else x for x in b])
                else:
                    setattr(self, a, Func(b) if isinstance(b, dict) else b)

    return Func(adict)


def save_result(save_path, result, args='', values=''):
    """This function is used to save result generated in `nei()`. 
    
    It saves all the results in '.pkl' file first, then save image files for the sinogram-like rho_t arrays. 
    If reconstructions are contained in the result object, save these as images as well.
    Parameters
    ----------
    save_path : `pathlib` object or str
        Absolute path to the directory where you want to save the results.
    result : class object
        Contains all the result.
    args : list of str, default ''
        A list of names of all the arguments in the `nei()` function.
    values : list of various types of data
        A list of values of all the arguments in the `nei()` function. 
    Returns
    -------
    None
    """
    # save the parameter setup used when running the program.
    save_path= Path(save_path)# making sure save_path is 'pathlib' class object
    log_file = save_path/'log.txt'
    with open(log_file, 'w') as file:
        for i in args:
            file.write("    %s = %s\n" % (i, values[i]))

    if result==None:
        return
    # save a `.pkl` file containing all the result.
    # save 'rho_t' into images.
    # save 'recon' into images if there is 'recon'.
    else: 
        pkl_file = 'save.pkl'
        save_object(result, save_path/pkl_file)
        if result.rho_t.ndim == 3:
            for i in range(result.rho_t.shape[0]):
                # imageio.imwrite(save_path/('rho_t' + str(i) + '.png'), result.rho_t[i].astype(np.uint8))
                # misc.imsave(save_path/('rho_t' + str(i) + '.tif'), result.rho_t[i])
                imageio.imwrite(save_path/('rho_t' + str(i) + '.tif'), result.rho_t[i].astype(np.float32))

                if not isinstance(result.recons, str): # if there is reconstruction image
                    imageio.imwrite(save_path/('recon' + str(i) + '.tif'), result.recons[i].astype(np.float32))
                    plt.figure(figsize=(9, 9))
                    ax = plt.gca()
                    im = ax.imshow(result.recons[i] * 1000, cmap='gray_r')
                    divider = make_axes_locatable(ax)
                    cax = divider.append_axes('right', size='5%', pad=0.05)
                    plt.colorbar(im, cax=cax).set_label('$mg/cm^3$', rotation=0, position=(1, 1), labelpad=-5)
                    plt.tight_layout()
                    plt.savefig(save_path / ('recon_bar_{}.tif'.format(i)))
        return


def save_recon(save_path, recon):
    """This function is made for the GUI reconstruction function to save reconstruction data and make images.

    A folder named with current date, and subfolders with current time will be generated automatically for saving these results.
    Parameters
    ----------
    save_path : pathlib object or str
        Absolute path to the directory where you want to save the results.
    recon : array
        The result from CT reconstruction process.
    Returns
    -------
    None

    """
    
    path = Path(save_path)
    dt = str(datetime.today())[0:19].replace(':', '-')
    pkl_file = dt + '_recon.pkl'
    save_object(recon, (path/pkl_file))
    if recon.ndim == 3:
        for i in range(recon.shape[0]): # todo
            imageio.imwrite(path /('recon_{}_{}.tif'.format(i, dt)), recon[i].astype(np.float32))
            plt.figure(figsize=(9, 9))
            ax=plt.gca()
            im=ax.imshow(recon[i] * 1000, cmap='gray_r')
            divider=make_axes_locatable(ax)
            cax=divider.append_axes('right', size='5%', pad=0.05)
            plt.colorbar(im, cax=cax).set_label('$mg/cm^3$', rotation=0, position=(1, 1), labelpad=-5)
            plt.tight_layout()
            plt.savefig(path/ ('recon_bar_{}_{}.tif'.format(i, dt)))

    if recon.ndim == 2: # todo
        imageio.imwrite(path /('recon_{}.tif'.format(dt), recon.astype(np.float32)))

    return


def draw_square(center_yx, width, color='k'):
    """Use the center location and width in pixel to draw a square on EXISTING figure. 
    
    This function is used only for the purpose of visualization for CT reconstruction image.
    Parameters
    ----------
    center_yx : list or tuple of 2 elements
        [vertical_position(y0), horizontal_position(x0)]. The pixel location for the center of the square you want to draw.
    width : int, positive
        The edge length for the square you want to draw.
    color : str, optional, default 'k'
        The color code used to draw the square. Check out `matplotlib` online for more details.
    Returns
    -------
    None
    """
    x0 = center_yx[1]
    y0 = center_yx[0]
    plt.plot([x0 - 0.5 * width, x0 + 0.5 * width], [y0 - 0.5 * width, y0 - 0.5 * width], color=color)
    plt.plot([x0 - 0.5 * width, x0 + 0.5 * width], [y0 + 0.5 * width, y0 + 0.5 * width], color=color)
    plt.plot([x0 - 0.5 * width, x0 - 0.5 * width], [y0 - 0.5 * width, y0 + 0.5 * width], color=color)
    plt.plot([x0 + 0.5 * width, x0 + 0.5 * width], [y0 - 0.5 * width, y0 + 0.5 * width], color=color)

    return