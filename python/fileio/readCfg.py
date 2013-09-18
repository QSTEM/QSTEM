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

import re
import numpy as np

# for more information on this regex, see the cfg_reading_regexes.txt file.
array_re = re.compile("(\d+)\s*\n([A-Za-z]{1,2})\s*\n(.*?)(?=\d+\s*\n\s*[A-Za-z]{1,2}\s*\n|\Z)", 
               flags=re.M|re.S)
box_re = re.compile("H0\((\d),(\d)\) = (\d+\.\d+)\s*")

def get_atom_arrays(filename):
    """
    Returns a dictionary with an entry per element.  The key to the dictionary 
    is the atomic number, while the value is an array containing coordinates,
    Debye-Waller factors, occupancies, and charges on a site-specific basis.
    """
    output_dict={}
    f=open(filename)
    txt=f.read()
    f.close()
    stripped_sections = array_re.findall(txt)
    for element in stripped_sections:
        # interprets the extracted string as a numpy array, and reshapes it
        #   so that we end up with one atomic coordinate per 
        output_dict[element[0]]=np.fromstring(element[2], sep=' ').reshape((-1,6))
    return output_dict

def get_cell_box(filename):
    """
    Returns a 3-element 1D numpy array with the cell dimensions, as X,Y,Z.
    """
    celldims=np.zeros(3)
    f=open(filename)
    txt=f.read()
    f.close()
    stripped_box_dims = box_re.findall(txt)
    for dim in stripped_box_dims:
        # keep it iff dimensions match - these are where
        #   the cell parameters are stored.
        if dim[0]==dim[1]:
            celldims[int(dim[0])-1]=float(dim[2])
    return celldims
            