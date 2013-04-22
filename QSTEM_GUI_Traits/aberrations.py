# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 15:55:58 2010

@author: Mike
"""

from enthought.traits.api import *
from enthought.traits.ui.api import \
    View, Item, VGroup, HGroup, spring, Handler, Group, TextEditor
from enthought.chaco.chaco_plot_editor import ChacoPlotItem
import enthought.traits.ui
import numpy as np

from enthought.chaco.api import Plot, ArrayPlotData, jet
from enthought.enable.component_editor import ComponentEditor

from numpy import cos, sqrt, arctan2, arcsin, pi, arange

class ImagePlot(HasTraits):
    plot = Instance(Plot)
    traits_view = View(
        Item('plot', editor=ComponentEditor(), show_label=False),
        width=500, height=500, resizable=True, title="Chaco Plot")

    def __init__(self):
        super(ImagePlot, self).__init__()
        x = np.linspace(0, 10, 50)
        y = np.linspace(0, 5, 50)
        xgrid, ygrid = np.meshgrid(x, y)
        z = np.exp(-(xgrid*xgrid+ygrid*ygrid)/100)
        plotdata = ArrayPlotData(imagedata = z)
        plot = Plot(plotdata)
        plot.img_plot("imagedata", xbounds=x, ybounds=y, colormap=jet)
        self.plot = plot

class Aberrations(HasTraits):
    """
    A generic class for aberrations.  This does not work by itself - classes
    below are defined for several notations, and each of these define how the
    a (aberration coefficients) and phi (angle) arrays are populated.
    """
    #alpha=Float(10.)
    #wl=Float(0.0000000002)
    N=Int(256)
    d=Float(0.5)
    k=Property(Array, depends_on=['N','dk','tiltX','tiltY'])
    dk=Property(Float, depends_on=['N','d'])
    qmax=Property(Float, depends_on=['alpha','wl'])
    ktheta=Property(Array, depends_on=(['k','wl']))
    kphi=Property(Array, depends_on=['k'])
    tiltX=Float
    tiltY=Float
    a=Array(dtype=np.float32, shape=(7,7))
    phi=Array(dtype=np.float32, shape=(7,7))
    abFunc=Property(depends_on=['a','phi','kphi','ktheta','wl'])
    
    def _get_dk(self):
        return 1.0/(self.N*self.d)
        
    def _get_qmax(self):
        return np.sin(self.alpha)/self.wl    
    
    def _get_k(self):
        N=self.N
        dk=self.dk
        tiltX=self.tiltX
        tiltY=self.tiltY
        return np.meshgrid(dk*(-N/2.+arange(N))+tiltX,dk*(-N/2.+arange(N))+tiltY)
    
    def _get_ktheta(self):
       return arcsin(sqrt(self.k[0]**2+self.k[1]**2)*self.wl)
    
    def _get_kphi(self):
        return arctan2(self.k[1],self.k[0])
    
    def _get_abFunc(self):
        """Returns a 2D array representing the aberration function"""
        a=self.a
        phi=self.phi
        kphi=self.kphi
        ktheta=self.ktheta
        wl=self.wl
        a[2:4]=a[2:4]*10**-9
        a[4:6]=a[4:6]*10**-6
        a[6]=a[6]*10**-3
        return 2*pi/wl*(1/2*(a[2,2]*cos(2*(kphi-phi[2,2]))+a[2,0])*ktheta**2+
                     1/3*(a[3,3]*cos(3*(kphi-phi[3,3]))+a[3,1]*cos(1*(kphi-phi[3,1])))*ktheta**3+
                     1/4*(a[4,4]*cos(4*(kphi-phi[4,4]))+a[4,2]*cos(2*(kphi-phi[4,2]))+a[4,0])*ktheta**4+
                     1/5*(a[5,5]*cos(5*(kphi-phi[5,5]))+a[5,3]*cos(3*(kphi-phi[5,3]))+a[5,1]*cos(1*(kphi-phi[5,1])))*ktheta**5+
                     1/6*(a[6,6]*cos(6*(kphi-phi[6,6]))+a[6,4]*cos(4*(kphi-phi[6,4]))+a[6,2]*cos(2*(kphi-phi[6,2]))+a[6,0])*ktheta**6)

    def setScherzer(self):
        pass
        
            
class AbKrivanek(Aberrations):
    C1_0=Float(desc='Defocus')
    C3_0=Float(desc = 'Spherical Aberration')
    C5_0=Float(desc='5th order spherical aberration')
    C1_2a=Float(desc='2-Fold astigmatism')
    C1_2b=Float(desc='2-Fold astigmatism')
    C2_3a=Float(desc='3-Fold astigmatism')
    C2_3b=Float(desc='3-Fold astigmatism')
    C3_4a=Float(desc='4-Fold astigmatism')
    C3_4b=Float(desc='4-Fold astigmatism')
    C4_5a=Float(desc='5-Fold astigmatism')
    C4_5b=Float(desc='5-Fold astigmatism')
    C5_6a=Float(desc='6-Fold astigmatism')
    C5_6b=Float(desc='6-Fold astigmatism')
    C2_1a=Float(desc='Axial coma')
    C2_1b=Float(desc='Axial coma')
    C4_1a=Float(desc='Axial coma')
    C4_1b=Float(desc='Axial coma')
    C3_2a=Float(desc='Axial star')
    C3_2b=Float(desc='Axial star')
    C5_2a=Float(desc='Axial star')
    C5_2b=Float(desc='Axial star')
    C4_3a=Float(desc='3-lobe aberration')
    C4_3b=Float(desc='3-lobe aberration')
    C5_4a=Float(desc='4-lobe aberration')
    C5_4b=Float(desc='4-lobe aberration')
    
    a = Property(Array, depends_on=['C1_0', 'C3_0', 'C5_0', 'C1_2a', 'C1_2b', 
                                    'C2_3a', 'C2_3b', 'C3_4a', 'C3_4b',
                                    'C4_5a', 'C4_5b', 'C5_6a', 'C5_6b',
                                    'C2_1a', 'C2_1b', 'C4_1a', 'C4_1b',
                                    'C3_2a', 'C3_2b', 'C5_2a', 'C5_2b',
                                    'C4_3a', 'C4_3b', 'C5_4a', 'C5_4b'])
                                    
    phi = Property(Array, depends_on=['C1_0', 'C3_0', 'C5_0', 'C1_2a', 'C1_2b',
                                      'C2_3a', 'C2_3b', 'C3_4a', 'C3_4b',
                                      'C4_5a', 'C4_5b', 'C5_6a', 'C5_6b',
                                      'C2_1a', 'C2_1b', 'C4_1a', 'C4_1b',
                                      'C3_2a', 'C3_2b', 'C5_2a', 'C5_2b',
                                      'C4_3a', 'C4_3b', 'C5_4a', 'C5_4b'])
                                    
    def _calc(self):
        aa = np.array([
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [self.C1_0, 0, self.C1_2a, 0, 0, 0, 0],
            [0, self.C2_1a, 0, self.C2_3a, 0, 0, 0],
            [self.C3_0, 0, self.C3_2a, 0, self.C3_4a, 0, 0],
            [0, self.C4_1a, 0, self.C4_3a, 0, self.C4_5a, 0],
            [self.C5_0, 0, self.C5_2a, 0, self.C5_4a, 0, self.C5_6a]])

        ab = np.array([
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, self.C1_2b, 0, 0, 0, 0],
            [0, self.C2_1b, 0, self.C2_3b, 0, 0, 0],
            [0, 0, self.C3_2b, 0, self.C3_4b, 0, 0],
            [0, self.C4_1b, 0, self.C4_3b, 0, self.C4_5b, 0],
            [0, 0, self.C5_2b, 0, self.C5_4b, 0, self.C5_6b]])
        
        cparr=np.zeros((7,7),dtype=np.complex)        
        cparr.real=aa; cparr.imag=ab
        return cparr
        
    def _get_a(self):
        return np.abs(self._calc())
        
    def _get_phi(self):
        return np.angle(self._calc())
    
    def _set_a(self,a=a,phi=phi):
        print type(a)
        print a, phi
        hv = a*np.cos(phi)
        vv = a*np.sin(phi)
        C1_0 = a[2,0]
        C3_0 = a[4,0]
        C5_0 = a[6,0]
        C1_2a = hv[2,2]
        C1_2b = vv[2,2]
        C2_3a = hv[3,3]
        C2_3b = vv[3,3]
        C3_4a = hv[4,4]
        C3_4b = vv[4,4]
        C4_5a = hv[5,5]
        C4_5b = vv[5,5]
        C5_6a = hv[6,6]
        C5_6b = vv[6,6]
        C2_1a = hv[3,1]
        C2_1b = vv[3,1]
        C4_1a = hv[5,1]
        C4_1b = vv[5,1]
        C3_2a = hv[4,2]
        C3_2b = vv[4,2]
        C5_2a = hv[6,2]
        C5_2b = vv[6,2]
        C4_3a = hv[5,3]
        C4_3b = vv[5,3]
        C5_4a = hv[6,4]
        C5_4b = vv[6,4]
        
    def _set_phi(self,a=a,phi=phi):
        self._set_a(a,phi)
        
    traits_view=View(
        VGroup(
            HGroup(
                Item('C1_2a',label='C1,2a  (nm)'),
                Item('C1_2b',label='C1,2b  (nm)'),
                Item('C1_0',label='C1       (nm)'),
                spring),
            HGroup(
                Item('C2_3a',label='C2,3a  (nm)'),
                Item('C2_3b',label='C2,3b  (nm)'),
                Item('C2_1a',label='C2,1a  (nm)'),
                Item('C2_1b',label='C2,1b  (nm)'),
                     spring),
            HGroup(
                Item('C3_4a',label='C3,4a  (um)'),
                Item('C3_4b',label='C3,4b  (um)'),
                Item('C3_2a',label='C3,2a  (um)'),
                Item('C3_2b',label='C3,2b  (um)'),
                Item('C3_0',label='C3        (um)'),
                spring),
            HGroup(
                Item('C4_5a',label='C4,5a  (um)'),
                Item('C4_5b',label='C4,5b  (um)'),
                Item('C4_3a',label='C4,3a  (um)'),
                Item('C4_3b',label='C4,3b  (um)'),
                Item('C4_1a',label='C4,1a   (um)'),
                Item('C4_1b',label='C4,1b  (um)'),
                spring),
            HGroup(
                Item('C5_6a',label='C5,6a (mm)'),
                Item('C5_6b',label='C5,6b (mm)'),
                Item('C5_4a',label='C5,4a (mm)'),
                Item('C5_4b',label='C5,4b (mm)'),
                Item('C5_2a',label='C5,2a (mm)'),
                Item('C5_2b',label='C5,2b (mm)'),
                Item('C5_0',label='C5      (mm)'),
                spring)
            )
        )

class AbRose(Aberrations):
    C1=Float(desc='Defocus')
    C3=Float(desc = 'Spherical Aberration')
    C5=Float(desc='5th order spherical aberration')
    A1=Float(desc='2-Fold astigmatism')
    A2=Float(desc='3-Fold astigmatism')
    A3=Float(desc='4-Fold astigmatism')
    A4=Float(desc='5-Fold astigmatism')
    A5=Float(desc='6-Fold astigmatism')
    B2=Float(desc='Axial coma')
    B4=Float(desc='Axial coma')
    S3=Float(desc='Axial star')
    S5=Float(desc='Axial star')
    D4=Float(desc='3-lobe aberration')
    R5=Float(desc='4-lobe aberration')
    PhiA1=Float(desc='2-Fold astigmatism angle')
    PhiA2=Float(desc='3-Fold astigmatism angle')
    PhiB2=Float(desc='Axial coma angle')
    PhiA3=Float(desc='4-Fold astigmatism angle')
    PhiS3=Float(desc='Axial star angle')
    PhiA4=Float(desc='5-Fold astigmatism angle')
    PhiD4=Float(desc='3-lobe aberration angle')
    PhiB4=Float(desc='Axial coma angle')
    PhiA5=Float(desc='6-Fold astigmatism angle')
    PhiR5=Float(desc='4-lobe aberration angle')
    PhiS5=Float(desc='Axial star angle')
    ang = Enum('mrad','degrees')
    
    a = Property(Array, depends_on=['C1', 'C3', 'C5', 'A1', 'A2', 'A3', 'A4',
                                    'A5', 'B2', 'B4', 'S3', 'S5', 'D4', 'B4',
                                    'R5'])
    phi = Property(Array, depends_on=['PhiA1', 'PhiA2', 'PhiA3', 'PhiA4',
                                      'PhiA5', 'PhiB2', 'PhiB4', 'PhiS3',
                                      'PhiS5', 'PhiD4', 'PhiR5','ang'])
    
    def _get_a(self):
        return np.array([
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [self.C1, 0, self.A1, 0, 0, 0, 0],
            [0, self.B2, 0, self.A2, 0, 0, 0],
            [self.C3, 0, self.S3, 0, self.A3, 0, 0],
            [0, self.B4, 0, self.D4, 0, self.A4, 0],
            [self.C5, 0, self.S5, 0, self.R5, 0, self.A5]])
    
    def _get_phi(self):
        return np.array([
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, self.PhiA1, 0, 0, 0, 0],
            [0, self.PhiB2, 0, self.PhiA2, 0, 0, 0],
            [0, 0, self.PhiS3, 0, self.PhiA3, 0, 0],
            [0, self.PhiB4, 0, self.PhiD4, 0, self.PhiA4, 0],
            [0, 0, self.PhiS5, 0, self.PhiR5, 0, self.PhiA5]])
            
    def _set_a(self,a=a):
        C1=a[2,0]
        C3=a[4,0]
        C5=a[6,0]
        A1=a[2,2]
        B2=a[3,1]
        A2=a[3,3]
        S3=a[4,2]
        A3=a[4,4]
        B4=a[5,1]
        D4=a[5,3]
        A4=a[5,5]
        S5=a[6,2]
        R5=a[6,4]
        A5=a[6,6]
        
    def _set_phi(self,phi=phi):
        PhiA1=phi[2,2]
        PhiB2=phi[3,1]
        PhiA2=phi[3,3]
        PhiS3=phi[4,2]
        PhiA3=phi[4,4]
        PhiB4=phi[5,1]
        PhiD4=phi[5,3]
        PhiA4=phi[5,5]
        PhiS5=phi[6,2]
        PhiR5=phi[6,4]
        PhiA5=phi[6,6]
    
    @on_trait_change('ang')
    def angConvert(self,trait,value):
        if value is 'mrad': self.phi=self.phi*17.453293
        if value is 'degrees':self.phi=self.phi/17.453293
        
    #@on_trait_change()
    def mrad2deg():
        return        
    
    traits_view=View(
        VGroup(
            HGroup(
                Item('A1',label='|A1| (nm)'),
                Item('PhiA1',label='PhiA1   '),
                Item('C1',label='|C1| (nm)'),
                spring,
                Item('ang',label='Angle notation',style='custom')
                ),
            HGroup(
                Item('A2',label='|A2| (nm)'),
                Item('PhiA2',label='PhiA2   '),
                Item('B2',label='|B2| (nm)'),
                Item('PhiB2',label='PhiB2   '),
                     spring),
            HGroup(
                Item('A3',label='|A3| (um)'),
                Item('PhiA3',label='PhiA3   '),
                Item('S3',label='|S3| (um)'),
                Item('PhiS3',label='PhiS3   '),
                Item('C3',label='|C3| (um)'),
                spring),
            HGroup(
                Item('A4',label='|A4| (um)'),
                Item('PhiA4',label='PhiA4   '),
                Item('D4',label='|D4| (um)'),
                Item('PhiD4',label='PhiD4   '),
                Item('B4',label='|B4| (um)'),
                Item('PhiB4',label='PhiB4   '),
                spring),
            HGroup(
                Item('A5',label='|A5| (mm)'),
                Item('PhiA5',label='PhiA5   '),
                Item('R5',label='|R5| (mm)'),
                Item('PhiR5',label='PhiR5   '),
                Item('S5',label='|S5| (mm)'),
                Item('PhiS5',label='PhiS5   '),
                Item('C5',label='|C5| (mm)'),
                spring)
            )
        )
  
class Ab3(Aberrations):
    a20=Range(low=-10000,high=10000,value=0)
    a40=Range(low=-10000,high=10000,value=1000)
    a60=Range(low=-100,high=100,value=5)
    a22=Range(low=0,high=1000,value=10)
    a31=Range(low=0,high=1000,value=10)
    a33=Range(low=0,high=1000,value=10)
    a42=Range(low=0,high=1000,value=10)
    a44=Range(low=0,high=1000,value=10)
    a51=Range(low=0,high=1000,value=10)
    a53=Range(low=0,high=1000,value=10)
    a55=Range(low=0,high=1000,value=10)
    a62=Range(low=0,high=1000,value=10)
    a64=Range(low=0,high=1000,value=10)
    a66=Range(low=0,high=1000,value=10)
    Phi22=Range(low=-pi,high=pi,value=0)
    Phi31=Range(low=-pi,high=pi,value=0)
    Phi33=Range(low=-pi,high=pi,value=0)
    Phi42=Range(low=-pi,high=pi,value=0)
    Phi44=Range(low=-pi,high=pi,value=0)
    Phi51=Range(low=-pi,high=pi,value=0)
    Phi53=Range(low=-pi,high=pi,value=0)
    Phi55=Range(low=-pi,high=pi,value=0)
    Phi62=Range(low=-pi,high=pi,value=0)
    Phi64=Range(low=-pi,high=pi,value=0)
    Phi66=Range(low=-pi,high=pi,value=0)

    a = Property( Array, depends_on=['a20', 'a40', 'a60', 'a22', 'a31', 'a33',
                                     'a42', 'a44', 'a51', 'a53', 'a55', 'a62',
                                     'a64','a66'])
                                    
    phi = Property( Array, depends_on=['Phi22', 'Phi31', 'Phi33', 'Phi42',
                                       'Phi44', 'Phi51', 'Phi53', 'Phi55',
                                       'Phi62', 'Phi64', 'Phi66'])
                                       
    def _get_a(self):
        return np.array([
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],
            [self.a20, 0, self.a22, 0, 0, 0, 0],
            [0, self.a31, 0, self.a33, 0, 0, 0],
            [self.a40, 0, self.a42, 0, self.a44, 0, 0],
            [0, self.a51, 0, self.a53, 0, self.a55, 0],
            [self.a60, 0, self.a62, 0, self.a64, 0, self.a66]])

    def _get_phi(self):
        return np.array([
            [0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0],                    
            [0, 0, self.Phi22, 0, 0, 0, 0],
            [0, self.Phi31, 0, self.Phi33, 0, 0, 0],
            [0, 0, self.Phi42, 0, self.Phi44, 0, 0],
            [0, self.Phi51, 0, self.Phi53, 0, self.Phi55, 0],
            [0, 0, self.Phi62, 0, self.Phi64, 0, self.Phi66]])
            
    def _set_a(self,a=a):
        a20=a[2,0]
        a40=a[4,0]
        a60=a[6,0]
        a22=a[2,2]
        a31=a[3,1]
        a33=a[3,3]
        a42=a[4,2]
        a44=a[4,4]
        a51=a[5,1]
        a53=a[5,3]
        a55=a[5,5]
        a62=a[6,2]
        a64=a[6,4]
        a66=a[6,6]
        
    def _set_phi(self,phi=phi):
        Phi22=phi[2,2]
        Phi31=phi[3,1]
        Phi33=phi[3,3]
        Phi42=phi[4,2]
        Phi44=phi[4,4]
        Phi51=phi[5,1]
        Phi53=phi[5,3]
        Phi55=phi[5,5]
        Phi62=phi[6,2]
        Phi64=phi[6,4]
        Phi66=phi[6,6]
        
    traits_view=View(
        VGroup(
            HGroup(
                Item('a22',label='|a2,2| (nm)'),
                Item('Phi22',label='Phi22'),
                Item('a20',label='|a2,0| (nm)'),
                spring),
            HGroup(
                Item('a33',label='|a3,3| (nm)'),
                Item('Phi33',label='Phi3,3'),
                Item('a31',label='|a3,1| (nm)'),
                Item('PhiB2',label='Phi3,1'),
                     spring),
            HGroup(
                Item('a44',label='|a4,4| (um)'),
                Item('Phi44',label='Phi44'),
                Item('a42',label='|a4,2| (um)'),
                Item('Phi42',label='Phi42 (um)'),
                Item('a40',label='|a4,0| (um)'),
                spring),
            HGroup(
                Item('a55',label='|a5,5| (um)'),
                Item('Phi55',label='Phi55'),
                Item('a53',label='|a5,3| (um)'),
                Item('Phi53',label='Phi53'),
                Item('a51',label='|a5,1| (um)'),
                Item('Phi51',label='Phi51'),
                spring),
            HGroup(
                Item('a66',label='|a6,6| (mm)'),
                Item('Phi66',label='Phi66'),
                Item('a64',label='|a6,4| (mm)'),
                Item('Phi64',label='Phi64'),
                Item('a62',label='|a6,2| (mm)'),
                Item('Phi62',label='Phi62'),
                Item('a60',label='|a6,0| (mm)'),
                spring)
            )
        )

class ProbeHandler( Handler ):
    def object_nomenclature_changed(self, info):
        if (info.object.nomenclature is 'Krivanek' and 
            (not isinstance( info.object.ab, AbKrivanek ))):
            info.object.ab = AbKrivanek(kw={'a':info.object.ab.a, 
                                            'phi':info.object.ab.phi,
                                            'wl':info.object.wl,
                                            'alpha':info.object.alpha})
        if (info.object.nomenclature is 'Rose' and 
            (not isinstance( info.object.ab, AbRose ))):
            info.object.ab = AbRose(kw={'a':info.object.ab.a,
                                        'phi':info.object.ab.phi,
                                        'wl':info.object.wl,
                                        'alpha':info.object.alpha})
        if (info.object.nomenclature is 'Number3' and
            (not isinstance( info.object.ab, Ab3 ))):
            info.object.ab = Ab3(kw={'a':info.object.ab.a, 'phi':info.object.ab.phi,
                                    'wl':info.object.wl, 'alpha':info.object.alpha})
                                              
class Probe(HasTraits):
    HT=Range(low=40,high=3000.0,value=200.0)
    alpha=Range(low=0.0,high=80.0,value=15.0)
    wl=Property(Float, depends_on=['HT'])
    nomenclature=Enum('Krivanek', 'Rose', 'Number3')
    ab=Instance(Aberrations)
    
    def _get_wl(self):
        h=6.626*10**-34
        m0=9.109*10**-31
        eV=1.602*10**-19*self.HT*1000
        C=2.998*10**8
        return h/np.sqrt(2*m0*eV*(1+eV/(2*m0*C**2)))*10**12
    
    gen_group = Group(
        HGroup(
        Item(name='nomenclature', label='Nomenclature'),
        spring,
        Item(name='HT',label="High Tension, kV",
             help='The microscope accelerating voltage'),
        Item('wl', label="Wavelength, pm ",style = 'readonly',
             editor=TextEditor(format_str='%3.2f')),
        spring,
        Item('alpha', label="Conv. Angle")),
        )
    ab_group = Group(
        Group(
            Item(name='ab',style='custom'),
            show_labels=False
        ),
        show_border = True
        )
        
    view=View(
        Group(gen_group,ab_group),
        title     = 'Higher-order Aberrations',
        buttons   = [ 'OK', 'Cancel' ],
        resizable = True,
        handler   = ProbeHandler()
        )

class ProbePlot(HasTraits):
    
    def __init__(self):
        super(Probe, self).__init__()
        x = linspace(0, self.N, self.N/self.d)
        y = linspace(0, self.N, self.N/self.d)
        xgrid, ygrid = meshgrid(x[1:], y[1:])
        z = exp(-(xgrid*xgrid+ygrid*ygrid)/10000)
        plotdata = ArrayPlotData(imagedata = z)
        plot = Plot(plotdata)
        self.renderer=plot.img_plot("imagedata", xbounds=x, ybounds=y, colormap=bone)
        #self.renderer = plot.plot(("x", "y"), type="scatter", color="blue")[0]
        self.plot = plot
    
    traits_view = View(
        VGroup(
            HGroup(Item('HT', label="High Tension, kV",help='The microscope accelerating voltage'),spring,Item('wl', label="Wavelength, nm",style = 'readonly')),
            HGroup(spring,Item('alpha', label="Conv. Angle")),
            HGroup(
                VGroup(
                    Item('noms',label="Nomenclature")),
                    Item('notations[self.noms]',label='Labels')),
            HGroup(Item('plot', editor=ComponentEditor(), show_label=False))),
        width=800, height=600, resizable=True, title="Chaco Plot"
        )
    
if __name__=='__main__':
    p=Probe(ab=AbKrivanek())
    p.configure_traits()