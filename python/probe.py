# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 23:32:26 2010

@author: Mike Sarahan
Code ported from Christoph Koch's MATLAB original.
"""

import numpy as np
from traits.api import HasTraits, Dict, Float, Tuple, Int, Property, \
     cached_property, Array
from matplotlib.pyplot import imshow

class probe(HasTraits):
    # TODO: might these be better as numpy arrays (or recarrays?)
    a = Array()
    phi = Dict()
    HT = Float(400)
    alpha = Float(15)
    tilt_x = Float(0)
    tilt_y = Float(0)
    N=Int(200) # the realspace number of pixels
    d=Float(0.05) # the sampling in real space
    dk = Property(depends_on = ['N','d'])
    wavelength = Property(depends_on = 'HT' )
    kx = Property(depends_on = ['dk','tilt_x']) # implicitly depends on N, d
    ky = Property(depends_on = ['dk','tilt_y']) # implicitly depends on N, d
    kx2 = Property(depends_on = 'kx')
    ky2 = Property(depends_on = 'ky')
    k2 = Property(depends_on = ['kx2','ky2'])
    qmax = Property(depends_on = ['alpha','wavelength'])
    ktheta = Property(depends_on = ['kx2','ky2','wavelength'])
    kphi = Property(depends_on = ['ky', 'kx'])
    ktm = Property(depends_on = ['qmax', 'wavelength'])
    # the aberration function value
    chi = Property(depends_on = ['a','phi','kphi','ktheta'])
    array = Property(depends_on = ['chi', 'N', 'ktheta', 'ktm', 'qmax', 'dk'])  
    
    def __init__(self, HT=200, N=400, d=0.05, alpha=15):
        self.HT=HT
        # Amplitudes of aberration coefficients
        self.a=self.get_template_a()
        # Angles of aberration coefficients   
        self.phi=self.get_template_phi()
        self.N = N #number of real space pixels
        self.d = d #sampling in real space
        self.alpha=alpha

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

    def set_coefficient(self, a_dict=None, phi_dict=None, name=None, amplitude=0, rotation=0):
        if a_dict is not None:
            self.a = a_dict
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
        Get the numpy structured array for 
        """
        a=np.zeros(1, self._get_a_dtype())
        a['20']=-60  # defocus: -60 nm
        a['40']=1000 # C3/spherical aberration: 1000 um
        a['60']=1    # C5/Chromatic aberration: 1 mm
        return a
    
    def get_template_phi(self):
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

    def plot_real_space_amplitude(self):
        from numpy.fft import ifft2, fftshift, ifftshift
        
        imshow(fftshift(np.abs(ifft2(ifftshift(self.array)))),
                   extent=[-self.dk/2*self.N,self.dk/2*self.N,
                             -self.dk/2*self.N,self.dk/2*self.N]);
            #else
                #imagesc(d*(-N/2+[0:N-1]),d*(-N/2+[0:N-1]),fftshift(angle(ifft2(ifftshift(probe)))));
                
    def plot_reciprocal_space_amplitude(self):
        #if get(handles.radiobutton_Amplitude,'Value')
        imshow(abs(self.array),
                   extent = [-self.dk/2*self.N,self.dk/2*self.N,
                             -self.dk/2*self.N,self.dk/2*self.N]);
            #else
                #imagesc(dkx*(-Nx/2+[0:Nx-1]),dky*(-Ny/2+[0:Ny-1]),angle(probe));
            
    ########  From here down, we calculate all the interrelated properties.
    ########    The user shouldn't need to directly use any of these methods.
    
    def _get_a_dtype(self):
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
        scaled_a = self.a.copy().view(np.float64)
        scales =np.array([10, 10,           #a_2x
                          10, 10,             #a_3x
                          1E4, 1E4, 1E4,      #a_4x
                          1E4, 1E4, 1E4,      #a_5x
                          1E7, 1E7, 1E7, 1E7,  #a_6x
                          ],dtype=np.float64)
        scaled_a *= scales
        return scaled_a.view(self._get_a_dtype())    
    
    @cached_property
    def _get_qmax(self):
        """
        TODO: define this
        qmax is 
        """
        return np.sin(self.alpha/1000)/self.wavelength

    @cached_property
    def _get_ktheta(self):
        """
        TODO: define this comment
        ktheta is 
        """
        return np.arcsin(np.sqrt(self.kx2+self.ky2)*self.wavelength)
    
    @cached_property
    def _get_kphi(self):
        """
        kphi is
        """
        return np.arctan2(self.ky,self.kx)
    
    @cached_property
    def _get_ktm(self):
        """
        ktm is
        """
        return np.arcsin(self.qmax*self.wavelength)

    @cached_property
    def _get_kx(self):
        """
        kx is
        """
        N = self.N
        tilt_x = self.tilt_x
        tilt_y = self.tilt_y
        dk = self.dk
        return np.meshgrid(dk*(-N/2.+np.arange(N))+tilt_x,dk*
                           (-N/2.+np.arange(N))+tilt_y)[0]
    
    @cached_property
    def _get_ky(self):
        """
        ky is 
        """
        N = self.N
        tilt_x = self.tilt_x
        tilt_y = self.tilt_y
        dk = self.dk        
        return np.meshgrid(dk*(-N/2.+np.arange(N))+tilt_x,dk*
                               (-N/2.+np.arange(N))+tilt_y)[1]    

    @cached_property
    def _get_dk(self):
        """
        dk is resolution in reciprocal space
        """
        return 1.0/(self.N*self.d)

    @cached_property
    def _get_kx2(self):
        return self.kx**2
    
    @cached_property
    def _get_ky2(self):
        return self.ky**2    

    @cached_property
    def _get_dk(self):
        return 1.0/(self.N*self.d)

    @cached_property
    def _get_k2(self):
        return self.kx2+self.ky2

    # in Angstroem
    @cached_property
    def _get_wavelength ( self ):
        """
        Computes relativistic wavelength according to stored HT
        """
        emass = 510.99906;   # electron rest mass in keV
        hc = 12.3984244;     # h*c
        ht = self.HT
        return hc/np.sqrt(ht*(2*emass+ht));
    
    @cached_property
    def _get_chi(self):
        """
        returns the aberration function
        """
        a=self._get_scaled_amplitudes()
        phi=self.phi
        kphi=self.kphi
        ktheta = self.ktheta
        pi=np.pi
        cos=np.cos
        chi = 2*pi/self.wavelength*(1/2*(a['22']*cos(2*(kphi-phi['22']))+a['20'])*ktheta**2 +
                1/3*(a['33']*cos(3*(kphi-phi['33']))+a['31']*cos(1*(kphi-phi['31'])))*ktheta**3 +
                1/4*(a['44']*cos(4*(kphi-phi['44']))+a['42']*cos(2*(kphi-phi['42']))+a['40'])*ktheta**4 +
                1/5*(a['55']*cos(5*(kphi-phi['55']))+a['53']*cos(3*(kphi-phi['53']))+a['51']*cos(1*(kphi-phi['51'])))*ktheta**5 +
                1/6*(a['66']*cos(6*(kphi-phi['66']))+a['64']*cos(4*(kphi-phi['64']))+a['62']*cos(2*(kphi-phi['62']))+a['60'])*ktheta**6)
        return chi
        
    @cached_property
    def _get_array(self):
        """
        returns the probe wavefunction array.  This function ultimately
        depends on every other function/value.
        """
        arr = np.zeros((self.N,self.N),dtype=np.complex);
        # MATLAB: probe(find(ktheta < ktm)) = 1;
        arr[self.ktheta<self.ktm] = 1+1j
        Nedge = 2;
        dEdge = Nedge/(self.qmax/self.dk);  # fraction of aperture radius that will be smoothed
        # some fancy indexing: pull out array elements that are within
        #    our smoothing edges
        ind = np.bitwise_and((self.ktheta/self.ktm > (1-dEdge)),
                             (self.ktheta/self.ktm < (1+dEdge)))
        arr[ind] = 0.5*(1-np.sin(np.pi/(2*dEdge)*(self.ktheta[ind]/self.ktm-1)));
        # add in the complex part
        # MATLAB: probe = probe.*exp(i*chi);
        arr*=np.exp(1j*self.chi);
        return arr

if __name__ == "__main__":
    from matplotlib import pyplot as plt
    p = probe()
    plt.plot_real_space_amplitude()