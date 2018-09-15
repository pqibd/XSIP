import pickle
import tkinter
import tkinter.filedialog as filedialog
import os
from PIL import Image
import fnmatch
import imageio
from datetime import datetime
import matplotlib.pyplot as plt

def choose_path():
    root = tkinter.Tk()
    root.withdraw()
    path = filedialog.askdirectory(title='Please select data directory:')
    # if path == '': choose_path()
    return path


def choose_file():
    root = tkinter.Tk()
    root.withdraw()
    file = filedialog.askopenfilename(title='Please select file')
    return file


def file_search(path,filter='*'):
    '''
    :param path: search files in this path
    :param filter: change the pattern of filter to match the wanted files.
    :return:
    '''

    files = fnmatch.filter(os.listdir(path),filter)
    if len(files)==0:
        print('------------------------'
              '\nWarning: No files found\n'
              '------------------------')
        return(files)
    files = [path+r'/'+file for file in files]
    return files


def save_object(obj, fname):
    print(fname)
    with open(fname, "wb") as output:
        # try:
        #     pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
        # except AttributeError:
        obj = class_to_dict(obj)
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
    return


def load_object(fname):
    """
    Accepted file formats are ".pkl" and image files
    :param fname:
    :return:
    """
    with open(fname, 'rb') as input:
        try:
            obj = pickle.load(input)
            if isinstance(obj, dict):
                obj = dict_to_class(obj)
        except:
            obj = imageio.imread(input)
    return obj


def draw_square(center_yx, width, color='k'):
    """
    Use the center location and width in pixel to draw a square on EXISTING plot
    :param center_yx: [y0,x0].
    :param width: width of the square, unit 1.
    :param color: choose the color for square edge.
    :return: Returns nothing. But use plt.show() after this function if nothing was shown.
    """
    import matplotlib.pyplot as plt
    x0 = center_yx[1]
    y0 = center_yx[0]
    plt.plot([x0 - 0.5 * width, x0 + 0.5 * width], [y0 - 0.5 * width, y0 - 0.5 * width], color=color)
    plt.plot([x0 - 0.5 * width, x0 + 0.5 * width], [y0 + 0.5 * width, y0 + 0.5 * width], color=color)
    plt.plot([x0 - 0.5 * width, x0 - 0.5 * width], [y0 - 0.5 * width, y0 + 0.5 * width], color=color)
    plt.plot([x0 + 0.5 * width, x0 + 0.5 * width], [y0 - 0.5 * width, y0 + 0.5 * width], color=color)

    return


def class_to_dict(obj):
    """
    Convert class object to dictionary. It does the conversion for all nested child class objects as well.
    :param obj: The original input will be modified. So make sure use copy.deepcopy(), and avoid using the
                original object.
    :return: Nested dictionary.
    """
    try:
        obj = obj.__dict__
        for key, value in obj.items():
            obj[key] = class_to_dict(value)
        return obj
    except:
        return obj


def dict_to_class(adict):
    """
    Convert dictionary to object.
    :param dict:
    :return:
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
    log_file = save_path + 'log.txt'
    with open(log_file, 'w') as file:
        for i in args:
            file.write("    %s = %s\n" % (i, values[i]))
        # file.write(result.beam_parameters.setup.__dict__)

    pkl_file = r'nei.pkl'
    save_object(result, save_path + pkl_file)
    if result.rho_t.ndim == 3:
        for i in range(result.rho_t.shape[0]):
            imageio.imwrite(save_path + 'rho_t' + str(i) + '.png', result.rho_t[i])
            if not isinstance(result.recons, str):
                imageio.imsave(save_path + 'recon' + str(i) + '.png', result.recons[i])
    return


def save_recon(save_path, recon):
    """
    This save function is used to save reconstruction data and make images. Needed in the GUI reconstruction part.
    """
    path = save_path + '\\'
    dt = str(datetime.today())[0:19].replace(':', '-')
    pkl_file = dt + r'_recon.pkl'
    save_object(recon, path + pkl_file)
    if recon.ndim == 3:
        for i in range(recon.shape[0]):
            imageio.imwrite(path + dt + '_recon' + str(i) + '.png', recon[i])

    if recon.ndim == 2:
        imageio.imwrite(path + dt + '_recon' + '.png', recon)
        plt.figure()


    return
