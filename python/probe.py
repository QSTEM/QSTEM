# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 23:32:26 2010

@author: Mike Sarahan
Code ported from Christoph Koch's MATLAB original.
"""

import numpy as np
from traits.api import HasTraits, Dict, Float, Tuple, Int, Property, \
     cached_property, Array, Bool
from matplotlib.pyplot import imshow, figure, title
from numpy.fft import fft2, ifft2, fftshift, ifftshift

import sys

class probe(HasTraits):
    # TODO: might these be better as numpy arrays (or recarrays?)
    a = Array()
    phi = Dict()
    HT = Float(200)
    alpha = Float(15)
    tilt_x = Float(0)
    tilt_y = Float(0)
    N=Int(200) # the realspace number of pixels
    d=Float(0.05) # the sampling in real space
    """
    dk = Property(depends_on = ['N','d'], cached=True)
    wavelength = Property(depends_on = 'HT', cached=True)
    kx = Property(depends_on = ['dk','tilt_x'], 
                  cached=True) # implicitly depends on N, d
    ky = Property(depends_on = ['dk','tilt_y'], 
                  cached=True) # implicitly depends on N, d
    kx2 = Property(depends_on = 'kx', cached=True)
    ky2 = Property(depends_on = 'ky', cached=True)
    k2 = Property(depends_on = ['kx2','ky2'], cached=True)
    qmax = Property(depends_on = ['alpha','wavelength'], cached=True)
    ktheta = Property(depends_on = ['kx2','ky2','wavelength'], cached=True)
    kphi = Property(depends_on = ['ky', 'kx'], cached=True)
    ktm = Property(depends_on = ['qmax', 'wavelength'], cached=True)
    # the aberration function value
    chi = Property(depends_on = ['a','phi','kphi','ktheta'], cached=True)
    array = Property(depends_on = ['chi', 'N', 'ktheta', 'ktm', 'qmax', 'dk'], 
                     cached=True)  
    """
    dk = Property
    wavelength = Property
    kx = Property
    ky = Property
    kx2 = Property
    ky2 = Property
    k2 = Property
    qmax = Property
    ktheta = Property
    kphi = Property
    ktm = Property
    # the aberration function value
    chi = Property
    array = Property
    debug = Bool(False)
    
    def __init__(self, HT=200, N=400, d=0.05, alpha=15, debug=False):
        self.HT=HT
        # Amplitudes of aberration coefficients
        self.a=self.get_template_a()
        # Angles of aberration coefficients   
        self.phi=self.get_template_phi()
        self.N = N #number of real space pixels
        self.d = d #sampling in real space
        self.alpha=alpha
        self.debug=debug

    def print_description(self):
        print "Voltage: %d"%self.HT
        print "Array size (pixels): %d"%self.N
        print "Real space pixel size (A): %.5f"%self.d
        print "Convergence angle (mrad): %.2f"%self.alpha
        print "Aberration amplitudes: "
        # self.a is a numpy structured array.  Turn it into
        #   a dictionary so the keys get printed with the values.
        names = self.a.dtype.names 
        a=[dict(zip(names, record)) for record in self.a]
        print a
        print "Aberration angles: "
        print self.phi
        
    #@function_name_decorator
    def set_coefficient(self, a_array=None, phi_dict=None, name=None, amplitude=0, rotation=0):
        """
        Use this to set one or more aberrations.  For setting a_array
        or phi_dict, you'll need to pass complete arrays/dicts.  You can
        get these by calling the get_template_a() and get_template_phi() 
        methods.
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        if a_array is not None:
            self.a = a_array
        if phi_dict is not None:
            self.phi = phi_dict
        # TODO: handle name as list?
        if name is not None:
            self.a[name]=amplitude
            # only set rotation for non-centrosymmetric aberrations
            if name in self.phi:
                self.phi[name]=rotation

    def set_scherzer_defocus(self):
        self.a['a20'] = -np.sign(self.a['a40'])*np.sqrt(1.5*np.abs(
            self.a['a40'])*self.wavelength(self.HT))

    def propagate(self, potential_slice):
        """
        Propagate this probe wavefunction through the given potential slice

        Input: a complex numpy array that is the same size as the probe's
               wavefunction array

        Output: none, but alters this probe's wavefunction array in place.
        """
        pass
        
    # this is a numpy recarray because we want to broadcast
    #   the scaling step.
    def get_template_a(self):
        """
        Get the numpy structured array representing the amplitudes
        of aberrations.
        """
        a=np.zeros(1, self._get_a_dtype())
        a['20']=-60  # defocus: -60 nm
        a['40']=1000 # C3/spherical aberration: 1000 um
        a['60']=1    # C5/Chromatic aberration: 1 mm
        return a
    
    def get_template_phi(self):
        """
        Get the python dictionary representing the angles of non-centro-
        symmetric abberations
        """
        return {
            '22':0.0,
            '33':0.0,
            '44':0.0,
            '55':0.0,
            '66':0.0,
            '31':0.0,
            '42':0.0,
            '53':0.0,
            '64':0.0,
            '51':0.0,
            '62':0.0
        }

    def get_physical_extent(self):
        return [-self.dk/2*self.N,self.dk/2*self.N,
                 -self.dk/2*self.N,self.dk/2*self.N]

    def get_real_space_amplitude(self):
        return fftshift(np.abs(ifft2(ifftshift(self.array))))

    def get_real_space_phase(self):
        return fftshift(np.angle(ifft2(ifftshift(self.array))))
        

    def plot_real_space_amplitude(self):
        imshow(self.get_real_space_amplitude(),
               extent=self.get_physical_extent())

    def plot_real_space_phase(self):
        imshow(self.get_real_space_phase(),
               extent=self.get_physical_extent())
                
    def get_reciprocal_space_amplitude(self):
        return abs(self.array)

    def get_reciprocal_space_phase(self):
        return np.angle(self.array)

    def plot_reciprocal_space_amplitude(self):
        #if get(handles.radiobutton_Amplitude,'Value')
        imshow(self.get_reciprocal_space_amplitude(),
               extent = self.get_physical_extent());

    def plot_reciprocal_space_phase(self):
        imshow(self.get_reciprocal_space_phase(),
               extent = self.get_physical_extent());
            
    ########  From here down, we calculate all the interrelated properties.
    ########    The user shouldn't need to directly use any of these methods.
    
    def _get_a_dtype(self):
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return {'names':['20', '22',
                         '31','33',
                         '40', '42', '44',
                         '51', '53', '55',
                         '60', '62', '64', '66'], 
                'formats':['f8']*14}    
        
    def _get_scaled_amplitudes(self):
        """
        translate from "reasonable" values from, say, a UI,
        into uniform units for the computation of chi
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        scaled_a = self.a.copy().view(np.float64)
        scales =np.array([10, 10,           #a_2x, nm->A
                          10, 10,             #a_3x, nm->A
                          1E4, 1E4, 1E4,      #a_4x, um->A
                          1E4, 1E4, 1E4,      #a_5x, um->A
                          1E7, 1E7, 1E7, 1E7,  #a_6x, mm->A
                          ],dtype=np.float64)
        scaled_a *= scales
        return scaled_a.view(self._get_a_dtype())    
    
    #@cached_property
    def _get_qmax(self):
        """
        TODO: define this
        qmax is 
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return np.sin(self.alpha/1000)/self.wavelength

    #@cached_property
    def _get_ktheta(self):
        """
        TODO: define this comment
        ktheta is 
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return np.arcsin(np.sqrt(self.kx2+self.ky2)*self.wavelength)
    
    #@cached_property
    def _get_kphi(self):
        """
        kphi is
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return np.arctan2(self.ky,self.kx)
    
    #@cached_property
    def _get_ktm(self):
        """
        ktm is
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return np.arcsin(self.qmax*self.wavelength)

    #@cached_property
    def _get_kx(self):
        """
        kx is
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        N = self.N
        tilt_x = self.tilt_x
        tilt_y = self.tilt_y
        dk = self.dk
        return np.meshgrid(dk*(-N/2.+np.arange(N))+tilt_x,dk*
                           (-N/2.+np.arange(N))+tilt_y)[0]
    
    #@cached_property
    def _get_ky(self):
        """
        ky is 
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        N = self.N
        tilt_x = self.tilt_x
        tilt_y = self.tilt_y
        dk = self.dk        
        return np.meshgrid(dk*(-N/2.+np.arange(N))+tilt_x,dk*
                               (-N/2.+np.arange(N))+tilt_y)[1]    

    #@cached_property
    def _get_dk(self):
        """
        dk is resolution in reciprocal space
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return 1.0/(self.N*self.d)

    #@cached_property
    def _get_kx2(self):
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return self.kx**2
    
    #@cached_property
    def _get_ky2(self):
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return self.ky**2    

    #@cached_property
    def _get_dk(self):
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return 1.0/(self.N*self.d)

    #@cached_property
    def _get_k2(self):
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        return self.kx2+self.ky2

    # in Angstroem
    #@cached_property
    def _get_wavelength ( self ):
        """
        Computes relativistic wavelength according to stored HT
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        emass = 510.99906;   # electron rest mass in keV
        hc = 12.3984244;     # h*c
        ht = self.HT
        return hc/np.sqrt(ht*(2*emass+ht));
    
    #@cached_property
    def _get_chi(self):
        """
        returns the aberration function
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        a=self._get_scaled_amplitudes()
        print a
        phi=self.phi
        kphi=self.kphi
        ktheta = self.ktheta
        pi=np.pi
        cos=np.cos
        wavelength=self.wavelength
        chi = 2.0*pi/wavelength*(1.0/2*(a['22']*cos(2*(kphi-phi['22']))+a['20'])*ktheta**2 +
                1.0/3*(a['33']*cos(3*(kphi-phi['33']))+a['31']*cos(1*(kphi-phi['31'])))*ktheta**3 +
                1.0/4*(a['44']*cos(4*(kphi-phi['44']))+a['42']*cos(2*(kphi-phi['42']))+a['40'])*ktheta**4 +
                1.0/5*(a['55']*cos(5*(kphi-phi['55']))+a['53']*cos(3*(kphi-phi['53']))+a['51']*cos(1*(kphi-phi['51'])))*ktheta**5 +
                1.0/6*(a['66']*cos(6*(kphi-phi['66']))+a['64']*cos(4*(kphi-phi['64']))+a['62']*cos(2*(kphi-phi['62']))+a['60'])*ktheta**6)
        return chi
        
    ##@cached_property
    def _get_array(self):
        """
        returns the probe wavefunction array.  This function ultimately
        depends on every other function/value.
        """
        if self.debug:
            print "executing " + sys._getframe().f_code.co_name
        arr = np.zeros((self.N,self.N),dtype=np.complex);
        # MATLAB: probe(find(ktheta < ktm)) = 1;
        arr[self.ktheta<self.ktm] = 1+1j
        Nedge = 2;
        dEdge = Nedge/(self.qmax/self.dk);  # fraction of aperture radius that will be smoothed
        # some fancy indexing: pull out array elements that are within
        #    our smoothing edges
        ind = np.bitwise_and((self.ktheta/self.ktm > (1-dEdge)),
                             (self.ktheta/self.ktm < (1+dEdge)))
        arr[ind] = 0.5*(1-np.sin(np.pi/(2*dEdge)*(self.ktheta[ind]/self.ktm-1)))
        # add in the complex part
        # MATLAB: probe = probe.*exp(i*chi);
        arr*=np.exp(1j*self.chi);
        return arr

if __name__ == "__main__":
    from matplotlib import pyplot as plt
    p = probe()
    plt.plot_real_space_amplitude()
