from read_img import binread2D
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from os.path import splitext

def plot_img(filename):
    img, comment, t, dx, dy = binread2D(filename, False)
    fig = plt.figure()
    fig.suptitle(comment, fontsize=14, fontweight='bold')
    ax=fig.add_subplot(111)
    extent = [0, img.shape[0]*dx, 0, img.shape[1]*dy]
    ax.imshow(img, extent=extent, interpolation="nearest")
    ax.set_title("Thickness = %.3fA"%t)
    ax.set_xlabel("Angstroms")
    ax.set_ylabel("Angstroms")
    plt.savefig(splitext(filename)[0]+".png",bbox_inches=0,)

if __name__=="__main__":
    import sys
    img=plot_img(sys.argv[1])
    
