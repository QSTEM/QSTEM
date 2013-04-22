% Through focus series reconstruction
% file throughFocus.m

format compact;
saveMemory = 2;

directory = 'sphere1';
directory = 'STO_gb';
directory = 'F:\Karleen4';
lambda = wavelength(300);  % electron wavelength for 200kV beam energy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Microscope parameters:
kmax = 1/0.1;     % objective aperture in 1/A set very large (0.1A)
Cs = 0.6e7;       % Cs in A  (1mm=1e7A) 1.2mm SESAM
Cc = 1.1e7;         % Cc in A 1mm = 1e7A
deltaE_E = 0.8/300e3;  % relative energy spread (0.4eV at 200kV)
alpha = 0.1e-3;  % half-angle convergence
Astig = 0;    % astigmatism in A
AstigPhi = 30*pi/180;  % angle of astigmatism


Df  = Cc*deltaE_E;
fprintf('Focal spread: %g A\n',Df);

def0 = -sqrt(1.5*Cs*lambda)       % basis of defocus (extended Scherzer)
% def0 = -sqrt(Cs*lambda)       % basis of defocus (Scherzer)
def0 = 200;       % basis of defocus (Scherzer)
defStep = -100;
Nimg = 20;

plotFlag = 1;

waveName = sprintf('%s/wave.img',directory);
[wave,t,dX,dY] = binread2D(waveName);
fprintf('dx: %g, dy: %g, thickness: %g\n',dX,dY,t);
[Ny,Nx] = size(wave);
if (0),    imagesc(dY*[0:Ny-1],dX*[0:Nx-1],abs(wave)); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add some carbon to the specimen surface:
carbon = exp(i*0.3*rand(size(wave)));
wave = fft2(wave.*carbon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the reciprocal space sampling:
dkx = 1/(Nx*dX);
dky = 1/(Ny*dY);

[kx,ky] = meshgrid(dkx*[-Nx/2:Nx/2-1],dky*[-Ny/2:Ny/2-1]);
k2 = fftshift(kx.^2+ky.^2);
phi = atan2(ky,kx);
clear kx ky


aperture = zeros(size(wave));
aperture(k2 <kmax*kmax) = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% temporal coherence envelope:
Etemp = exp(-2*(Df*lambda*k2).^2);   % all units in A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% def = def0+(idf-1)*defStep;
defArray = def0 + defStep*[0:Nimg-1];  % regular defocus array
defArray = [def0+40 def0+20 defArray]; % 2 additional images fro alignment
for idf=1:length(defArray)
    def = defArray(idf);
    Espat = exp(-(pi*alpha*def+pi*alpha*lambda^2*Cs.*k2).^2.*k2);
    chi = exp(-i*pi*(lambda*(def+Astig.*cos(2*(phi-AstigPhi))).*k2+0.5*lambda^3*Cs.*k2.^2));
    % chi = exp(-i*pi*(lambda*(def).*k2+0.5*lambda^3*Cs*k2.^2));
    
    img = abs(ifft2(aperture.*wave.*chi.*Etemp.*Espat)).^2+0.02*randn(size(wave));

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Save image
  if (1)
      fileName = sprintf('F:/Karleen4/STOimg_%d.img',idf);
      fprintf('Image %d (def=%gnm):  ',idf,def/10);
      binwrite2D(img,fileName,0,t);
  else
  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Display a few figures:
      subplot(2,2,1)
      imagesc(dY*[0:Ny-1],dX*[0:Nx-1],img);
      axis equal; axis tight;
      colormap('gray');
      title(sprintf('def=%gnm',def/10));
      
      subplot(2,2,2)
      img = abs(fft2(img));
      img(1,1) = 0;
      imagesc(dky*[-Ny/2:Ny/2-1],dkx*[-Nx/2:Nx/2-1],log(1+fftshift(img)));
      title('Diffractogram');
      xlim([-1.5 1.5])
      ylim([-1.5 1.5])
      
      
      subplot(2,2,3)
      plot(dkx*[-Nx/2:Nx/2-1],[fftshift(Etemp(:,1)) fftshift(Espat(:,1)) ...
          fftshift(-imag(chi(:,1))).*fftshift(Etemp(:,1)).*fftshift(Espat(:,1))]);
      legend('temporal coherence','spatial coherence', 'CTF')
      xlim([0 2]);
      
      subplot(2,2,4)
      imagesc(dky*[-Ny/2:Ny/2-1],dkx*[-Nx/2:Nx/2-1],- ...
          fftshift(-imag(chi).*Etemp.*Espat));
      title('2D transfer function (Im[chi(q)])');
      xlim([-1.5 1.5])
      ylim([-1.5 1.5])
      
      
      pause(0.1)
  end
end





