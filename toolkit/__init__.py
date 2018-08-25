import pickle
import tkinter
import tkinter.filedialog as filedialog


def choose_path():
    root = tkinter.Tk()
    root.withdraw()
    path = filedialog.askdirectory(title='Please select NEI directory:')
    # if path == '': choose_path()
    return path


def file_search(path,filter='*'):
    '''
    :param path: search files in this path
    :param filter: change the pattern of filter to match the wanted files
    :return:
    '''
    import os
    import fnmatch
    files = fnmatch.filter(os.listdir(path),filter)
    if len(files)==0:
        print('------------------------'
              '\nWarning: No files found\n'
              '------------------------')
        return(files)
    files = [path+r'/'+file for file in files]
    return files


def save_object(obj, fname):
    with open(fname, "wb") as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
    return


def reload_object(fname):
    with open(fname, 'rb') as input:
        object = pickle.load(input)
    return object


def draw_square(center_yx, width, color='k'):
    """
    Use the center location and width in index to draw a square shape on existing plot
    :param center_yx: [y0,x0].
    :param width: width of the square, unit 1.
    :param color: choose the color for square edge.
    :return: Nothing is returned. But use plt.show() after this function if nothing was shown.
    """
    import matplotlib.pyplot as plt
    x0 = center_yx[1]
    y0 = center_yx[0]
    plt.plot([x0 - 0.5 * width, x0 + 0.5 * width], [y0 - 0.5 * width, y0 - 0.5 * width], color=color)
    plt.plot([x0 - 0.5 * width, x0 + 0.5 * width], [y0 + 0.5 * width, y0 + 0.5 * width], color=color)
    plt.plot([x0 - 0.5 * width, x0 - 0.5 * width], [y0 - 0.5 * width, y0 + 0.5 * width], color=color)
    plt.plot([x0 + 0.5 * width, x0 + 0.5 * width], [y0 - 0.5 * width, y0 + 0.5 * width], color=color)

    return
