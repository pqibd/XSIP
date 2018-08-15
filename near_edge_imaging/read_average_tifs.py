def read_average_tifs(files,flip=False,xlow=0,xhigh=0,
                      rotate_90=False,twelve_bit=0):
    # read all the image files into a 3D array[n_images,rows,columns],
    # take the average along the images, so that we get an average image
    from PIL import Image
    import numpy as np
    # twelve_bit=
    n_files = len(files)
    image_array = []
    for i in range(n_files):
        image_array.append(np.array(Image.open(files[i])))
    image_array = np.array(image_array)
    average = image_array.mean(axis=0)
    return(average)


