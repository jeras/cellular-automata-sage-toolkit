"""
PYCA vizualization tools
"""

import numpy
import Image

def array2image (array, k, name='') :
    """
    computes t time steps and writes the resulting configurations into an image or video
    """
    if (array.ndim == 1) :
        img = Image.fromstring('P', (1, array.shape[0]), array.tostring());
    if (array.ndim == 2) :
        img = Image.fromstring('P', (array.shape[1], array.shape[0]), array.tostring());
    else :
        print "Error: unsupported array size, from a 3D array a 2D video should be created, not yet implemented"
    if (k == 2) :
        img.putpalette([0,0,0,255,255,255]);
    else :
        print "Error: images with more than two colors not yet implemented"
    if (name != '') :
        img.save(name+".png", "PNG")
    return img

