"""
/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
    Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
"""

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
    
