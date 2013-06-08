from read_img import binread2D
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from os.path import splitext

def plot_img(filename):
    img, comment, t, dx, dy = binread2D(filename, False)
    extent = [0, img.shape[0]*dx, 0, img.shape[1]*dy]
    plt.imshow(img, extent=extent)
    plt.title(comment)
    plt.savefig(splitext(filename)[0]+".png",bbox_inches=0,)

if __name__=="__main__":
    import sys
    img=plot_img(sys.argv[1])
    
