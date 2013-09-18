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

"""
% function [img,t,dx,dy] = binread2D(fileName,printFlag,flag)
% file for reading binary data of different formats
% input: 
% fileName - name of file (string)
% printFlag - 0: non-verbose, 1: verbose (default)
% output:
% img      - image (complex or real
% t        - thickness
% dx, dy   - pixel size (for calibration)
function [img,t,dx,dy] = binread2D(fileName,printFlag,flag)

file format from:
http://elim.physik.uni-ulm.de/?page_id=873
"""

import numpy as np
import struct

headerlength = 4;

def binread2D(filename, printFlag=True):
    comment=""
    # open the file and define the file ID (fid):
    f=open(filename,"rb")
    #with open(filename,'rb') as f:
    if True:
        # there are 8 ints followed by 3 doubles
        header = struct.unpack("iiiiiiiiddd",f.read(56))

        print header

        headerSize = header[0]

        # read additional parameters from file, if any exist:
        paramSize = header[1]

        commentSize = header[2]

        Nx = header[3]
        Ny = header[4]

        complexFlag = header[5]
        dataSize = header[6]
        doubleFlag = (dataSize==8*(complexFlag+1))
        complexFlag = bool(complexFlag)

        version = header[7]
        
        thicknessOrDefocus=header[8]

        dx = header[9]
        dy = header[10]

        if (paramSize > 0):
            params = np.fromfile(file=f, dtype=np.float64, count=paramSize);
            if printFlag:
                print '%d Parameters:'%paramSize
                print params

        # read comments from file, if any exist:
        if (commentSize > 0):
            comment = struct.unpack("%ds"%commentSize,f.read(commentSize))[0]

    if printFlag:
        print('binread2D %s: %d x %d pixels'%(filename,Nx,Ny))

    if complexFlag:
        if doubleFlag:
            if printFlag:
                print '64-bit complex data, %.3fMB)\n'%(Nx*Ny*16/1048576)
            img = np.fromfile(file=f, dtype=np.complex128, count=Nx*Ny)
        else:
            if printFlag:
                fprintf('32-bit complex data, %.3fMB)\n',Nx*Ny*8/1048576);
            img = np.fromfile(file=f, dtype=np.complex64, count = Nx*Ny)
    else:
        if doubleFlag:
            if printFlag:
                fprintf('64-bit real data, %.3fMB)\n',Nx*Ny*8/1048576);
            img = np.fromfile(file=f, dtype=np.float64, count=Nx*Ny)
        else:
            if printFlag:
                print '32-bit real data, %.3fMB)\n'%(Nx*Ny*4/1048576)
            img = np.fromfile(file=f, dtype=np.float32, count=Nx*Ny)
    img=img.reshape(Ny,Nx)
    
    return img, comment, thicknessOrDefocus, dx, dy

if __name__=="__main__":
    import sys
    filename = sys.argv[1]
    img, comment, t, dx, dy = binread2D(filename,False)
    print img
    print comment
    print "Thickness/defocus: %.3f"%t
    print "X pixel size (A): %.3f"%dx
    print "Y pixel size (A): %.3f"%dy
