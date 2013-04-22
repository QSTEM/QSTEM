#!/usr/bin/env python

import numpy as np

def get_dtype(hdr):
    element_size=hdr['element_size']
    if element_size==4: # 32-bit float
        dtype=np.float32
    elif element_size==8:
        if hdr['is_complex']==1: # 32-bit complex
            dtype=np.complex64
        else: #double
            dtype=np.float64
    elif element_size==16:
        dtype=np.complex128
    return dtype

def readIMG (filename, debug=False):
    # total header size: 56 bytes
    hdr_dtype=np.dtype([
        ('hdr_size','<i4'),
        ('param_size','<i4'),
        ('comment_size','<i4'),
        ('nx','<i4'),
        ('ny','<i4'),
        ('is_complex','<i4'),
        ('element_size','<i4'),
        ('qstem_version','<i4'),
        ('thickness','<f8'),
        ('x_px_size','<f8'),
        ('y_px_size','<f8')])
    hdr=np.fromfile(filename,dtype=hdr_dtype,count=1)
    if debug:
        print "header, param, comment sizes: "
        print hdr['hdr_size'], hdr['param_size'], hdr['comment_size']
        print "Image size: "
        print hdr['nx'], hdr['ny']
        if hdr['is_complex']==1: print "Data is complex."
        print "Data element size: "
        print hdr['element_size']
        print "Thickness at this image: "
        print hdr['thickness']
        print "Pixel size: "
        print hdr['x_px_size'], hdr['y_px_size']
    data_dtype=get_dtype(hdr)
    f=open(filename,"rb")
    # skip the header bytes
    f.read(56)
    aux_data = f.read(8*hdr['param_size'])
    comments = str(f.read(hdr['comment_size']))
    data = f.read(hdr['nx']*hdr['ny']*hdr['element_size'])
    data = np.frombuffer(data,count=hdr['nx']*hdr['ny'],dtype=data_dtype)
    data = data.reshape((hdr['nx'],hdr['ny']))
    return data

def save_to_image(data,filename,format='.tif'):
    import Image
    if np.iscomplexobj(data):
        realimg=Image.fromarray(data.real)
        imagimg=Image.fromarray(data.imag)
        realimg.save(filename+"_amp"+format)
        imagimg.save(filename+"_phase"+format)
    else:
        img=Image.fromarray(data)
        img.save(filename)
