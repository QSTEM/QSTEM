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

import numpy as np
from traits.api import HasTraits, Dict, Float, Tuple, Int, Property, \
     cached_property, Array, Bool, on_trait_change
from fileio import readCfg
import copy

class SampleModel(HasTraits):
    # elements is a dictionary representing the structure.
    #   Each entry has 4 things:
    #   - The element symbol (this is the dictionary key)
    #   - The atomic number
    #   - An array of fractional coordinates (x, y, z), 
    #          a fourth column for Debye-Waller factor,
    #          a fifth column for occupancy,
    #          a sixth column for site-specific charge
    elements = Dict
    scaled_elements = Property(depends_on="elements",cached=True)
    duplicated_elements = Property(depends_on="scaled_elements",cached=True)
    transformed_elements = Property(depends_on="duplicated_elements",cached=True)
    nCellsX = Int(4)
    nCellsY = Int(4)
    nCellsZ = Int(4)
    tiltX = Float(0)
    tiltY = Float(0)
    tiltZ = Float(0)
    boxSizeX = Float(10)
    boxSizeY = Float(10)
    boxSizeZ = Float(10)
    cellSize = Array
    tiltInDegrees = Bool(True) # when false, tilt is in radians.
    usingNCells = Bool(True) # when true, computes supercell as multiple of 
    #                           cells.  When false, fills boxSize with 
    #                           duplicated cells.
        
    def loadCfg(self,filename):
        self.cellSize=readCfg.get_cell_box(filename)
        self.elements=readCfg.get_atom_arrays(filename)
    
    def _get_scaled_elements(self):
        """
        Translate fractional coordinates into real-space coordinates
        """
        scaled_elements={}
        for key in self.elements.keys():
            scaled_elements[key]=copy.copy(self.elements[key])
            scaled_elements[key][:,:3]*=self.cellSize
        return scaled_elements    
    
    def _get_duplicated_elements(self):
        """
        Duplicates elements according to your desired number of cells.
        """
        return self._duplicateCell(self.nCellsX, self.nCellsY, self.nCellsZ)

    def _get_transformed_elements(self):
        if not self.tiltInDegrees:
            pass
        return self._tiltCell(self.tiltX, self.tiltY, self.tiltZ)
        
    def _duplicateCell(self, nx, ny, nz):
        duplicated_elements={}
        for key in self.elements.keys():
            num_repeats = nx*ny*nz
            dup_array = np.tile(self.scaled_elements[key], [num_repeats,1])
            num_entries = self.scaled_elements[key].shape[0]
            ctr=0
            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        offset=ctr*num_entries
                        cell_offset = np.array([x*self.cellSize[0], 
                                                y*self.cellSize[1], 
                                                z*self.cellSize[2],
                                                0,0,0])
                        dup_array[offset:offset+num_entries]+=cell_offset
                        ctr += 1
            duplicated_elements[key]=dup_array
        return duplicated_elements
                
    def _tiltCell(self, x, y, z):
        from numpy import cos, sin
        transformed_elements=copy.copy(self.duplicated_elements)
        
        Rx = np.array([[1,0,0],
                       [0, cos(x), -sin(x)],
                       [0, sin(x), cos(x)]])
        Ry = np.array([[cos(y),0,sin(y)],
                       [0, 1, 0],
                       [-sin(y), 0, cos(y)]])
        Rz = np.array([[cos(z),-sin(z),0],
                       [sin(z), cos(z),0],
                       [0, 0, 1]])
        for key in self.duplicated_elements.keys():
            transformed_elements[key][:,:3]=self.duplicated_elements[key][:,:3].dot(Rx).dot(Ry).dot(Rz)
        return transformed_elements
            
    def _get_element_arrays(self):
        """
        This is an alias for the transformed arrays, which represent the structure that has
        been scaled, duplicated, and rotated as specified.
        """
        return self.transformed_elements