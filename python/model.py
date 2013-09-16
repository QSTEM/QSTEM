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
     cached_property, Array, Bool

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
    nCellsX = Int(4)
    nCellsY = Int(4)
    nCellsZ = Int(4)
    tiltX = Float(0)
    tiltY = Float(0)
    tiltZ = Float(0)
    boxSizeX = Float(10)
    boxSizeY = Float(10)
    boxSizeZ = Float(10)
    cellSizeX = Float
    cellSizeY = Float
    cellSizeZ = Float
    tiltInDegrees = Bool(True) # when false, tilt is in radians.
    usingNCells = Bool(True) # when true, computes supercell as multiple of cells
    #                            When false, fills boxSize with duplicated cells.
    
    def duplicateCell(self):
        
    
    