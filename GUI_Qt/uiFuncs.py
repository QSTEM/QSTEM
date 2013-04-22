# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 23:32:26 2010

@author: Mike Sarahan
Code ported from Christoph Koch's MATLAB original.
"""

import numpy as np
    
class probe(object):
    def __init__(self, mode='STEM', HT=200, a=None, phi=None, alpha=15,
                 tilt=None):
        self.mode=mode
        self.HT=HT
        # Amplitudes of aberration coefficients
        self.a={
            'C1':0.0,
            'C3':0.0,
            'C5':0.0,
            'A1':0.0,
            'A2':0.0,
            'A3':0.0,
            'A4':0.0,
            'A5':0.0,
            'B2':0.0,
            'B4':0.0,
            'S3':0.0,
            'S5':0.0,
            'D4':0.0,
            'R5':0.0
            }
        # Angles of aberration coefficients   
        self.phi={
            'A1':0.0,
            'A2':0.0,
            'A3':0.0,
            'A4':0.0,
            'A5':0.0,
            'B2':0.0,
            'B4':0.0,
            'S3':0.0,
            'S5':0.0,
            'D4':0.0,
            'R5':0.0
            }
        self.N #number of real space pixels
        self.d #sampling in real space
        self.dk=1.0/(self.N*self.d)
        self.kx,self.ky=np.meshgrid(dk[0]*(-N[0]/2.+arange(N[0]))+tilt[0],dk[0]*(-N[1]/2.+arange(N[1]))+tilt[1])
        # kx2 = kx.^2;
        # ky2 = ky.^2;
        # k2 = kx2+ky2;
        self.l=wavelength(self.HT)
        self.alpha=alpha
        self.qmax  = np.sin(self.alpha)/self.l 
        self.tilt=tilt
        self.probeArray=np.zeros((self.N))
        self.ktheta = arcsin(sqrt(self.kx**2+self.ky**2)*self.l)
        self.kphi = np.arctan2(self.ky,self.kx)
        self.ktm = arcsin(qmax*l)

    def scherzer(self):
        self.a['C1']=-np.sign(self.a['C3'])*np.sqrt(1.5*np.abs(self.a['C3'])*wavelength(self.HT))
        
    def aberrationFunc(self):
        phi=self.phi
        kphi=self.kphi
        pi=np.pi
        cos=np.cos
        self.chi = 2*pi/self.l*(1/2*(a['A1']*cos(2*(kphi-phi['A1']))+a['C1'])*ktheta**2+
                     1/3*(a['A2']*cos(3*(kphi-phi['A2']))+a['B2']*cos(1*(kphi-phi['B2'])))*ktheta**3+
                     1/4*(a['A3']*cos(4*(kphi-phi['A3']))+a['S3']*cos(2*(kphi-phi['S3']))+a['C3'])*ktheta**4+
                     1/5*(a['A4']*cos(5*(kphi-phi['A4']))+a['D4']*cos(3*(kphi-phi['D4']))+a['B4']*cos(1*(kphi-phi['B4'])))*ktheta**5+
                     1/6*(a['A5']*cos(6*(kphi-phi['A5']))+a['R5']*cos(4*(kphi-phi['R5']))+a['S5']*cos(2*(kphi-phi['S5']))+a['C5'])*ktheta**6)

    def 
    probe = zeros(Ny,Nx);
    probe(find(ktheta < ktm)) = 1;
    Nedge = 2;
    dEdge = Nedge/(qmax/dkx);  % fraction of aperture radius that will be smoothed
    ind = find((ktheta/ktm > 1-dEdge) & (ktheta/ktm < 1+dEdge));
    probe(ind) = 0.5*(1-sin(pi/(2*dEdge)*(ktheta(ind)/ktm-1)));


handles = readParams(handles);
guidata(hObject, handles);
a = handles.a;          % Array of amplitudes of aberrations (according to nomenclature 3)  
phi = handles.phi;      % Array of angles of aberrations   (according to nomenclature 3)    
c = handles.c;          % List of symmetric aberrations  (c(2) = def, c(4) = Cs, c(6) = C5) 




chi   = 2*np.pi/l*(1/2*(a(2,2).*cos(2*(kphi-phi(2,2)))+c(2)).*ktheta.^2+...
                     1/3*(a(3,3).*cos(3*(kphi-phi(3,3)))+a(3,1).*cos(1*(kphi-phi(3,1)))).*ktheta.^3+...
                     1/4*(a(4,4).*cos(4*(kphi-phi(4,4)))+a(4,2).*cos(2*(kphi-phi(4,2)))+c(4)).*ktheta.^4+...
                     1/5*(a(5,5).*cos(5*(kphi-phi(5,5)))+a(5,3).*cos(3*(kphi-phi(5,3)))+a(5,1).*cos(1*(kphi-phi(5,1)))).*ktheta.^5+...
                     1/6*(a(6,6).*cos(6*(kphi-phi(6,6)))+a(6,4).*cos(4*(kphi-phi(6,4)))+a(6,2).*cos(2*(kphi-phi(6,2)))+c(6)).*ktheta.^6);
% compute the probe
probe = probe.*exp(i*chi);

if (0)
    figure('Name','probe intensity');
    if get(handles.radiobutton_RealSpace,'Value')
        if get(handles.radiobutton_Amplitude,'Value')
            imagesc(dx*(-Nx/2+[0:Nx-1]),dy*(-Ny/2+[0:Ny-1]),fftshift(abs(ifft2(ifftshift(probe)))));
        else
            imagesc(dx*(-Nx/2+[0:Nx-1]),dy*(-Ny/2+[0:Ny-1]),fftshift(angle(ifft2(ifftshift(probe)))));
        end
    else
        if get(handles.radiobutton_Amplitude,'Value')
            imagesc(dkx*(-Nx/2+[0:Nx-1]),dky*(-Ny/2+[0:Ny-1]),abs(probe));
        else
            imagesc(dkx*(-Nx/2+[0:Nx-1]),dky*(-Ny/2+[0:Ny-1]),angle(probe));
        end
    end
else
    axes(handles.axes_Probe);
    if get(handles.checkbox_Surf,'Value')
        surf(dx*(-Nx/2+[0:Nx-1]),dy*(-Ny/2+[0:Ny-1]),fftshift(abs(ifft2(ifftshift(probe))).^2));   
        shading interp;
    else
        imagesc(dx*(-Nx/2+[0:Nx-1]),dy*(-Ny/2+[0:Ny-1]),fftshift(abs(ifft2(ifftshift(probe))).^2));
        set(gca,'YDir','normal');
    end
    if get(handles.checkbox_Color,'Value')
        colormap('default');
    else
        colormap('gray');
    end
    axes(handles.axes_Phasemap);
    imagesc(dkx*(-Nx/2+[0:Nx-1]),dky*(-Ny/2+[0:Ny-1]),angle(probe));    
    set(gca,'YDir','normal');
    if get(handles.checkbox_Color,'Value')
        colormap('default');
    else
        colormap('gray');
    end
end