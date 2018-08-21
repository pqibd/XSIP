import pickle

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
