function varargout = imageSim(varargin)
% IMAGESIM M-file for imageSim.fig
%      IMAGESIM, by itself, creates a new IMAGESIM or raises the existing
%      singleton*.
%
%      H = IMAGESIM returns the handle to a new IMAGESIM or the handle to
%      the existing singleton*.
%
%      IMAGESIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGESIM.M with the given input arguments.
%
%      IMAGESIM('Property','Value',...) creates a new IMAGESIM or raises
%      the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imageSim_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to imageSim_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imageSim

% Last Modified by GUIDE v2.5 10-Nov-2011 20:42:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imageSim_OpeningFcn, ...
                   'gui_OutputFcn',  @imageSim_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before imageSim is made visible.
function imageSim_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imageSim (see VARARGIN)

% Choose default command line output for imageSim
handles.output = hObject;
handles.initialized = 0;
handles.diffractogramMode = 0;
handles.applyVibration = 0;
handles.vibration_Axis1 = 0;
handles.vibration_Axis2 = 0;
handles.vibration_Angle = 0;
handles.c = zeros(6,1);
handles.a = zeros(6);
handles.phi = zeros(6);

% Read the wave function from disk
handles = readWave(handles);
% initialize sampling distance.
% handles.dx = 10*str2double(get(handles.editSampling,'String'));
handles = updateKArrays(handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);
% handles = updateImages(handles);
% Update handles structure

guidata(hObject, handles);

% UIWAIT makes imageSim wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imageSim_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main drawing function:
%% This is the most relevant function.  It contains all the maths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function handles = updateImages(handles)
[Ny,Nx] = size(handles.wave0);
handles.Amorph   = str2double(get(handles.editAmorph,'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the wave function from the different weightings and amorphous
% signal:
if (0)
    if handles.MinAmplitude < 1
        ampl = handles.MinAmplitude+abs(handles.wave0)*(1-handles.MinAmplitude);
    else
        ampl = ones(size(handles.wave0));
    end
    phase = handles.MaxPhase*angle(handles.wave0);
    amorph = handles.Amorph*(rand(Ny,Nx)-0.5);
    handles.wave = ampl.*exp(i*(phase+amorph));
else
    amorph = handles.Amorph*(rand(Ny,Nx)-0.5);
    handles.wave = handles.wave0+i*amorph;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the CTF (this is where the Math is):

% apply the beam tilt as an offset to the optic axis:
kx2 = handles.kx2;
ky2 = handles.ky2;
k2  = (kx2.^2+ky2.^2);
theta2  = k2*(handles.lambda.^2);
theta   = sqrt(theta2);
phi     = atan2(ky2,kx2);

pl      = pi*handles.lambda;
p_l     = pi/handles.lambda;
fa      = handles.Defocus+abs(handles.Astig)*cos(2*(phi-angle(handles.Astig)));
coma3a  = abs(handles.Coma)*cos((phi-angle(handles.Coma)))+...
          abs(handles.Astig3)*cos(3*(phi-angle(handles.Astig3)));
% ctf =   (defocus+astigmatism) + spherical aberration + coma)
handles.c(2)     = handles.Defocus;
handles.c(4)     = handles.Cs;
handles.phi(2,2) = angle(handles.Astig);
handles.a(2,2)   = abs(handles.Astig);
handles.a(3,1)   = abs(handles.Coma);
handles.phi(3,1) = angle(handles.Coma);
handles.a(3,3)   = abs(handles.Astig3);
handles.phi(3,3) = angle(handles.Astig3);
if get(handles.checkbox_IncludeHigherOrders,'Value')

    chi=2*p_l*(1/2*(handles.a(2,2).*cos(2*(phi-handles.phi(2,2)))+handles.c(2)).*theta2+...
               1/3*(handles.a(3,3).*cos(3*(phi-handles.phi(3,3)))+handles.a(3,1).*cos(1*(phi-handles.phi(3,1)))).*theta2.*theta+...
               1/4*(handles.a(4,4).*cos(4*(phi-handles.phi(4,4)))+handles.a(4,2).*cos(2*(phi-handles.phi(4,2)))+handles.c(4)).*(theta2.^2)+...
               1/5*(handles.a(5,5).*cos(5*(phi-handles.phi(5,5)))+handles.a(5,3).*cos(3*(phi-handles.phi(5,3)))+handles.a(5,1).*cos(1*(phi-handles.phi(5,1)))).*(theta.^5)+...
               1/6*(handles.a(6,6).*cos(6*(phi-handles.phi(6,6)))+handles.a(6,4).*cos(4*(phi-handles.phi(6,4)))+handles.a(6,2).*cos(2*(phi-handles.phi(6,2)))+handles.c(6)).*(theta2.^3));
    ctf     = exp(-i*chi);    
    
    if (get(handles.checkboxConvAngle,'Value'))
        % I will use: [dchi(q)/dq]^2 = [dchi/d|q|]^2+[dchi/|q|dphi]^2
        % first let's compute the q-derivative (the division by lambda is removed in both cases):
        dchi_dq = ((handles.a(2,2).*cos(2*(phi-handles.phi(2,2)))+handles.c(2)).*theta+...
                   (handles.a(3,3).*cos(3*(phi-handles.phi(3,3)))+handles.a(3,1).*cos(1*(phi-handles.phi(3,1)))).*theta2+...
                   (handles.a(4,4).*cos(4*(phi-handles.phi(4,4)))+handles.a(4,2).*cos(2*(phi-handles.phi(4,2)))+handles.c(4)).*theta.^3+...
                   (handles.a(5,5).*cos(5*(phi-handles.phi(5,5)))+handles.a(5,3).*cos(3*(phi-handles.phi(5,3)))+handles.a(5,1).*cos(1*(phi-handles.phi(5,1)))).*theta2.^2+...
                   (handles.a(6,6).*cos(6*(phi-handles.phi(6,6)))+handles.a(6,4).*cos(4*(phi-handles.phi(6,4)))+handles.a(6,2).*cos(2*(phi-handles.phi(6,2)))+handles.c(6)).*theta.^5)*2*pi;
        dchi_dphi=(1/2*(2*handles.a(2,2).*sin(2*(phi-handles.phi(2,2)))).*theta+...
                   1/3*(3*handles.a(3,3).*sin(3*(phi-handles.phi(3,3)))+1*handles.a(3,1).*sin(1*(phi-handles.phi(3,1)))).*theta2+...
                   1/4*(4*handles.a(4,4).*sin(4*(phi-handles.phi(4,4)))+2*handles.a(4,2).*sin(2*(phi-handles.phi(4,2)))).*theta.^3+...
                   1/5*(5*handles.a(5,5).*sin(5*(phi-handles.phi(5,5)))+3*handles.a(5,3).*sin(3*(phi-handles.phi(5,3)))+1*handles.a(5,1).*sin(1*(phi-handles.phi(5,1)))).*theta2.^2+...
                   1/6*(6*handles.a(6,6).*sin(6*(phi-handles.phi(6,6)))+4*handles.a(6,4).*sin(4*(phi-handles.phi(6,4)))+2*handles.a(6,2).*sin(2*(phi-handles.phi(6,2)))).*theta.^5)*(-2*pi);
        Espat = exp(-sign(handles.Convergence)*(handles.Convergence/(2*handles.lambda)).^2.*(dchi_dq.^2+dchi_dphi.^2));
        Espat(find(Espat > 5)) = 5;
        ctf = ctf.*Espat;
    end
else
    chi=2*p_l*(1/2*(handles.a(2,2).*cos(2*(phi-handles.phi(2,2)))+handles.c(2)).*theta2+...
               1/3*(handles.a(3,3).*cos(3*(phi-handles.phi(3,3)))+handles.a(3,1).*cos(1*(phi-handles.phi(3,1)))).*theta2.*theta+...
               1/4*(handles.c(4)).*(theta2.^2));
    ctf     = exp(-i*chi);    

    % ctf = exp(-i*2*p_l*theta2.*(1/2*fa+1/3*coma3a.*theta+1/4*handles.Cs*theta2));
    if (get(handles.checkboxConvAngle,'Value'))
       dchi_dq = ((handles.a(2,2).*cos(2*(phi-handles.phi(2,2)))+handles.c(2)).*theta+...
                   (handles.a(3,3).*cos(3*(phi-handles.phi(3,3)))+handles.a(3,1).*cos(1*(phi-handles.phi(3,1)))).*theta2+...
                   (handles.c(4)).*theta.^3)*2*pi;
        dchi_dphi=(1/2*(2*handles.a(2,2).*sin(2*(phi-handles.phi(2,2)))).*theta+...
                   1/3*(3*handles.a(3,3).*sin(3*(phi-handles.phi(3,3)))+1*handles.a(3,1).*sin(1*(phi-handles.phi(3,1)))).*theta2)*(-2*pi);
        Espat = exp(-sign(handles.Convergence)*(handles.Convergence/(2*handles.lambda)).^2.*(dchi_dq.^2+dchi_dphi.^2));
        % Espat = exp(-sign(handles.Convergence)*(pi*handles.Convergence*(fa+handles.lambda^2*handles.Cs.*k2)).^2.*k2);
        Espat(find(Espat > 5)) = 5;
        ctf = ctf.*Espat;

    end
end

if (handles.ObjApertureTheta > 0)
    handles.ObjApertureEdge = 1e-3*str2num(get(handles.edit_ObjApertureEdge,'String'));
    ctf(find(theta >= (handles.ObjApertureTheta+handles.ObjApertureEdge))) = 0;
    ind = find((theta < (handles.ObjApertureTheta+handles.ObjApertureEdge)) & (theta > handles.ObjApertureTheta));
    ctf(ind) = 0.5*ctf(ind).*(1+cos(pi*(theta(ind)-handles.ObjApertureTheta)/handles.ObjApertureEdge));
    clear ind
end

% Compute the damping envelopes:
% This is what I used in simImageFromWave - check this!!!:
% Etemp = exp(-2.0*(1.0/handles.lambda*handles.Delta*theta2).^2);  



if (get(handles.checkboxDelta,'Value'))
    Etemp = exp(-sign(handles.Delta)*(0.5*pi/handles.lambda*handles.Delta*theta2).^2);
    Etemp(find(Etemp > 3)) = 3;
    
    % figure; plot(handles.dkx*10*[-Nx/2:Nx/2-1],[fftshift(Espat(:,1)) fftshift(Etemp(:,1))]); pause(0.1);
    % Apply the damping envelopes as well:
    ctf = Etemp.*ctf;
end

if (handles.TiltX ~= 0) || (handles.TiltY ~= 0)
    tx = round(handles.TiltX/(handles.lambda*handles.dkx));
    ty = round(handles.TiltY/(handles.lambda*handles.dky));
    handles.Psi = ifft2(circshift(fft2(handles.wave),[ty tx]).*ctf);
else
    tx = 0;
    ty = 0;
    handles.Psi = ifft2(fft2(handles.wave).*ctf);
end
% handles.Cs*max(max(theta2)).^2

% The radiobutton 'Wave' decides whether to plot the wave, or the CTF
if get(handles.radiobuttonWave,'Value')  
    % show the wave function:
    axes(handles.axesAmplitude);
    xl = get(handles.axesAmplitude,'XLim');  yl = get(handles.axesAmplitude,'YLim');    cla;
    imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],abs(handles.wave));
    set(gca,'YDir','normal');
    title('amplitide');
    xlabel('x in nm');
    ylabel('y in nm');    
    axis equal; axis tight;
    if (handles.initialized)
       % xlim(xl); ylim(yl); 
       set(handles.axesAmplitude,'XLim',xl); set(handles.axesAmplitude,'YLim',yl);
    end
    
    axes(handles.axesPhase);
    xl = get(handles.axesPhase,'XLim');  yl = get(handles.axesPhase,'YLim');    cla;
    imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],angle(handles.wave));
    set(gca,'YDir','normal');
    title('phase');
    axis equal; axis tight;
    if (handles.initialized)
       % xlim(xl); ylim(yl); 
       set(handles.axesPhase,'XLim',xl); set(handles.axesPhase,'YLim',yl);
    end
else
    % show the CTF:
    axes(handles.axesAmplitude);
    xl = get(handles.axesAmplitude,'XLim');  yl = get(handles.axesAmplitude,'YLim');    cla;
    % handles.kx
    % imagesc(handles.kx,handles.ky,real(ctf));
    imagesc(handles.dkx*10*[-Nx/2:Nx/2-1],handles.dky*10*[-Ny/2:Ny/2-1],fftshift(real(ctf)));
    set(gca,'YDir','normal');
    title('Real(CTF)');
    xlabel('kx in 1/nm');
    ylabel('ky in 1/nm');    

    axis equal; axis tight;
    if (handles.initialized)
       %xlim(xl); ylim(yl); 
       set(handles.axesAmplitude,'XLim',xl); set(handles.axesAmplitude,'YLim',yl);
    end

    
    axes(handles.axesPhase);
    xl = get(handles.axesPhase,'XLim');  yl = get(handles.axesPhase,'YLim');    cla;
    imagesc(handles.dkx*10*[-Nx/2:Nx/2-1],handles.dky*10*[-Ny/2:Ny/2-1],fftshift(imag(ctf)));
    set(gca,'YDir','normal');
    hold on
    if (get(handles.checkboxDelta,'Value') && get(handles.checkboxConvAngle,'Value'))
        plot(handles.dky*10*[-Ny/2:Ny/2-1],handles.dky*5*Ny*(0.9/max([Espat(:,1); Etemp(:,1)])*[fftshift(Espat(:,1)) fftshift(Etemp(:,1))]),'Linewidth',2);
        legend('E_{spat}','E_{temp}');
    end
    if (get(handles.checkboxDelta,'Value') && (~get(handles.checkboxConvAngle,'Value')))
        plot(handles.dky*10*[-Ny/2:Ny/2-1],handles.dky*5*Ny*(0.9/max(Etemp(:,1))*[fftshift(Etemp(:,1))]),'Linewidth',2);
        legend('E_{temp}');
    end
    if (~get(handles.checkboxDelta,'Value') && get(handles.checkboxConvAngle,'Value'))
        plot(handles.dky*10*[-Ny/2:Ny/2-1],handles.dky*5*Ny*(0.9/max(Espat(:,1))*[fftshift(Espat(:,1))]),'Linewidth',2);
        legend('E_{spat}');
    end
    
    hold off
    title('Imag(CTF)');  
    axis equal; axis tight;
    if (handles.initialized)
       % xlim(xl); ylim(yl); 
       set(handles.axesPhase,'XLim',xl); set(handles.axesPhase,'YLim',yl);
    end
end

axes(handles.axesImage);
xl = get(handles.axesImage,'XLim');  yl = get(handles.axesImage,'YLim');    cla;
if get(handles.radiobuttonImage,'Value')  
    % show the image
    simIm = abs(handles.Psi).^2;
    if get(handles.checkbox_Vibration,'Value')
        vib1 = (pi*handles.vibration_Axis1/(-log(0.5))).^2;
        vib2 = (pi*handles.vibration_Axis2/(-log(0.5))).^2;
        simIm = real(ifft2(fft2(simIm).*exp(-vib1*k2.*abs(cos(phi-handles.vibration_Angle))-vib2*k2.*abs(sin(phi-handles.vibration_Angle)))));
    end
    imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],simIm);
    title(sprintf('Image t=%.1f nm [%.4f .. %.4f]',0.1*handles.thickness,min(min(simIm)),max(max(simIm))));
elseif get(handles.radiobutton_WaveAmplitude,'Value')
    % simIm = abs(handles.Psi);
    % imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],simIm);
    % title(sprintf('Amplitude t=%.1f nm [%.4f .. %.4f]',0.1*handles.thickness,min(min(simIm)),max(max(simIm))));
    diffg = log(1+abs(fft2(handles.Psi)));
    diffg(mod(ty+Ny,Ny)+1,1+mod(tx+Nx,Nx)) = 0;
    imagesc(handles.dkx*10*[-Nx/2:Nx/2-1],handles.dky*10*[-Ny/2:Ny/2-1],fftshift(diffg));
    title('Diffraction Pattern (log scale)');
elseif get(handles.radiobutton_WavePhase,'Value')
    simIm = angle(handles.Psi);
    imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],simIm);
    title(sprintf('Phase t=%.1f nm [%.4f .. %.4f]',0.1*handles.thickness,min(min(simIm)),max(max(simIm))));
else
    % show the diffractogram
    diffg = abs(fft2(abs(handles.Psi).^2));
    if get(handles.checkbox_Vibration,'Value')
        vib1 = (pi*handles.vibration_Axis1/(-log(0.5))).^2;
        vib2 = (pi*handles.vibration_Axis2/(-log(0.5))).^2;
        diffg = diffg.*exp(-vib1*k2.*abs(cos(phi-handles.vibration_Angle))-vib2*k2.*abs(sin(phi-handles.vibration_Angle)));
    end

    diffg = log(1+diffg);
    diffg(1,1) = 0;
    imagesc(handles.dkx*10*[-Nx/2:Nx/2-1],handles.dky*10*[-Ny/2:Ny/2-1],fftshift(diffg));
    title('Diffractogram (log scale)');
end
colormap('gray');
set(gca,'YDir','normal');
if (get(handles.radiobuttonDiffractogram,'Value') == 0) && (get(handles.radiobutton_WaveAmplitude,'Value') == 0)
    xlabel('x in nm');
    ylabel('y in nm');
else
    xlabel('x in 1/nm');
    ylabel('y in 1/nm');    
end
axis equal; axis tight;
if (handles.initialized)
    % xlim(xl); ylim(yl);
    set(handles.axesImage,'XLim',xl); set(handles.axesImage,'YLim',yl);
end


handles.initialized = 1;

% guidata(gcbo, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A few utility functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = readWave(handles)
if (0)
    % img1 = sum(imread('amplitude.jpg'),3);
    % img2 = sum(imread('phase.jpg'),3);
    % normalize the image intensities:
    img1 = img1/max(max(img1));
    img2 = img2/max(max(img2))-0.5;
    
    handles.wave0 = img1.*exp(i*img2);
    clear img1 img2
else
    [filename, pathname] = uigetfile({'*.img','Wave function';'*.*','All files'},'Select a complex wave function');
    if isequal(filename,0) || isequal(pathname,0)
        % disp('User pressed cancel')
    else
        filename =  fullfile(pathname, filename);
        
        [img,t,dx,dy,params,comment] = binread2D(filename);
        [Ny,Nx] = size(img);
        realFlag = isreal(img);
        comment = char(comment);
        handles.wave0 = img;
        handles.dx = dx;
        handles.dy = dy;
        handles.thickness = t;
        
        if ~isempty(params)
            if length(params) > 8
                handles.highVoltage = params(1);
                handles.Cs = 1e-7*params(2);
                handles.Defocus = params(3);
                handles.Astig = 0.01*params(4)*exp(i*params(5));		% astigmatism
                handles.Delta = params(6);           % Delta = Cc * dE_E = focal spread
				handles.Convergence = params(7);		% illumination convergence angle
				handles.TiltX = 1e-3*180/pi*params(8);			% beam tilt in deg -> rad
				handles.TiltY = 1e-3*180/pi*params(9);			% beam tilt in deg -> rad

                set(handles.editAccVoltage,'String',handles.highVoltage);
                set(handles.editAstig1,'String',10*abs(handles.Astig));
                set(handles.editAstig2,'String',180/pi*angle(handles.Astig));
                set(handles.editDefocus,'String',0.1*handles.Defocus);
                set(handles.editCs,'String',sprintf('%.4f',1e-7*handles.Cs));
                set(handles.editConvAngle,'String',handles.Convergence);
                set(handles.editDelta,'String',0.1*handles.Delta);
                set(handles.editTiltX,'String',1e3*handles.TiltX);
                set(handles.editTiltY,'String',1e3*handles.TiltY);
            end
        end
        % fprintf('Pixel size: %g x %g\n',dx,dy);
        set(handles.editSampling,'String',0.1*min(dx,dy));
    end    
end
handles.numWaves = 1;             
% guidata(gcbo, handles);



function handles = updateKArrays(handles)
[h,w] = size(handles.wave0);
handles.dkx       = 1/(w*handles.dx);
handles.dky       = 1/(h*handles.dy);
kx  = fftshift(handles.dkx*[-w/2:w/2-1]);
ky  = fftshift(handles.dky*[-h/2:h/2-1].');
handles.kx2 = repmat(kx,h,1);
handles.ky2 = repmat(ky,1,w);
% handles.phi = atan2(handles.ky2,handles.kx2);
% guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Handler for the different checkbox- and edit-controls:
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbuttonUpdate.
function handles = pushbuttonUpdate_Callback(hObject, eventdata, handles)

% First, we need to retreive the relevant imaging parameters
handles.MinAmplitude  = str2double(get(handles.editMinAmplitude,'String'));
handles.MaxPhase = str2double(get(handles.editMaxPhase,'String'));
handles.Amorph   = str2double(get(handles.editAmorph,'String'));
handles.highVoltage = str2double(get(handles.editAccVoltage,'String'));
handles.lambda   = wavelength(handles.highVoltage);
set(handles.textAccVoltageUnit,'String',sprintf('kV, lambda=%.4f A',handles.lambda));
% Read Defocus, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxDefocus,'Value')
    % convert Defocus from nm to A
    handles.Defocus = 10*str2double(get(handles.editDefocus,'String'));
else, handles.Defocus=0; end
% Read Cs, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxCs,'Value')
    % convert Cs from mm to A
    handles.Cs = 1e7*str2double(get(handles.editCs,'String'));
else, handles.Cs=0; end
% Read Tilt, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxTilt,'Value')
    % convert Tilt from mrad to rad
    handles.TiltX = 1e-3*str2double(get(handles.editTiltX,'String'));
    handles.TiltY = 1e-3*str2double(get(handles.editTiltY,'String'));
else, handles.TiltX=0; handles.TiltY=0; end
if get(handles.checkboxComa,'Value')
    % convert Coma from nm and deg to A and rad
    handles.Coma = 1e4*abs(str2double(get(handles.editComa,'String')))*...
        exp(i*pi/180*str2double(get(handles.editComaAngle,'String')));
else, handles.Coma=0;  end
if get(handles.checkbox3foldAstig,'Value')
    % convert Coma from nm and deg to A and rad
    handles.Astig3 = 1e4*abs(str2double(get(handles.edit3foldAstig,'String')))*...
        exp(i*pi/180*str2double(get(handles.edit3foldAstigAngle,'String')));
else, handles.Astig3=0; end
% Read Astigmatism, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxAstig,'Value')
    % convert Tilt from mrad to rad
    handles.Astig = 10*str2double(get(handles.editAstig1,'String'))*...
        exp(i*pi/180*str2double(get(handles.editAstig2,'String')));
else, handles.Astig=0; end

% Read Delta, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxDelta,'Value')
    % convert Cs from mm to A
    handles.Delta = 10*str2double(get(handles.editDelta,'String'));
else, handles.Delta=0; end
% Read Convergence angle, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxConvAngle,'Value')
    % convert Cs from mm to A
    handles.Convergence = 1e-3*str2double(get(handles.editConvAngle,'String'));
else, handles.Convergence=0; end
% Read objective aperture angle:
if get(handles.checkbox_ObjAperture,'Value')
    % convert angle from mrad to rad
    handles.ObjApertureTheta = 1e-3*str2double(get(handles.edit_ObjApertureTheta,'String'));
else, handles.ObjApertureTheta=0; end
% Read sample vibration
if get(handles.checkbox_Vibration,'Value')
    % convert vibration from pm to A:
    handles.vibration_Axis1 = 0.1*str2double(get(handles.edit_Vibration_Axis1,'String'));
    handles.vibration_Axis2 = 0.1*str2double(get(handles.edit_Vibration_Axis2,'String'));
    handles.vibration_Angle = pi/180*str2double(get(handles.edit_Vibration_Angle,'String'));
    handles.applyVibration = 1;
else
    handles.vibration_Axis1 = 0;
    handles.vibration_Axis2 = 0;
    handles.vibration_Angle = 0;
    handles.applyVibration = 0;
end


%% Redraw Images
handles = updateImages(handles);
guidata(hObject, handles);


% --- Executes on button press in checkboxDefocus.
function checkboxDefocus_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkboxTilt.
function checkboxTilt_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkboxComa.
function checkboxComa_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkboxImage.
function checkboxImage_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkboxDiffractogram.
function checkboxDiffractogram_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkboxRealCTF.
function checkboxRealCTF_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkboxImagCTF.
function checkboxImagCTF_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkboxDelta.
function checkboxDelta_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkboxConvAngle.
function checkboxConvAngle_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkbox_ObjAperture.
function checkbox_ObjAperture_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edit fields:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_ObjApertureEdge_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_ObjApertureEdge_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ObjApertTheta_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_ObjApertTheta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDefocus_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editDefocus_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCs_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editCs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editAstig1_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editAstig1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editAstig2_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editAstig2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCurv1_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editCurv1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editDelta_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editDelta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editConvAngle_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editConvAngle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editPincushion_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editPincushion_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSpiral_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editSpiral_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSampling_Callback(hObject, eventdata, handles)
%handles.dx = 10*str2double(get(hObject,'String'));
%handles = updateKArrays(handles);
%pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function editSampling_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editAccVoltage_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editAccVoltage_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editTiltX_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editTiltX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTiltY_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editTiltY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function editComa_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function editComa_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in radiobuttonWave.
function radiobuttonWave_Callback(hObject, eventdata, handles)
set(handles.axesAmplitude,'XLimMode','auto');
set(handles.axesAmplitude,'XLimMode','auto');
set(handles.axesPhase,'XLimMode','auto');
set(handles.axesPhase,'XLimMode','auto');
handles.initialized = 0;
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in radiobuttonCTF.
function radiobuttonCTF_Callback(hObject, eventdata, handles)
set(handles.axesAmplitude,'XLimMode','auto');
set(handles.axesAmplitude,'XLimMode','auto');
set(handles.axesPhase,'XLimMode','auto');
set(handles.axesPhase,'XLimMode','auto');
handles.initialized = 0;
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in radiobuttonImage.
function radiobuttonImage_Callback(hObject, eventdata, handles)
if (handles.diffractogramMode)
    set(handles.axesImage,'XLimMode','auto');
    set(handles.axesImage,'XLimMode','auto');
    handles.initialized = 0;
end
handles.diffractogramMode = 0;
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in radiobuttonDiffractogram.
function radiobuttonDiffractogram_Callback(hObject, eventdata, handles)
if (handles.diffractogramMode == 0)
    set(handles.axesImage,'XLimMode','auto');
    set(handles.axesImage,'XLimMode','auto');
    handles.initialized = 0;
end
handles.diffractogramMode = 1;
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in radiobutton_WaveAmplitude.
function radiobutton_WaveAmplitude_Callback(hObject, eventdata, handles)
%if (handles.diffractogramMode)
%    set(handles.axesImage,'XLimMode','auto');
%    set(handles.axesImage,'XLimMode','auto');
%    handles.initialized = 0;
%end
%handles.diffractogramMode = 0;
if (handles.diffractogramMode == 0)
    set(handles.axesImage,'XLimMode','auto');
    set(handles.axesImage,'XLimMode','auto');
    handles.initialized = 0;
end
handles.diffractogramMode = 1;
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on button press in radiobutton_WavePhase.
function radiobutton_WavePhase_Callback(hObject, eventdata, handles)
if (handles.diffractogramMode)
    set(handles.axesImage,'XLimMode','auto');
    set(handles.axesImage,'XLimMode','auto');
    handles.initialized = 0;
end
handles.diffractogramMode = 0;
handles.initialized = 0;
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);





function editMinAmplitude_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function editMinAmplitude_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMaxPhase_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function editMaxPhase_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editAmorph_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function editAmorph_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function editComaAngle_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function editComaAngle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox3foldAstig.
function checkbox3foldAstig_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


function edit3foldAstig_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit3foldAstig_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3foldAstigAngle_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit3foldAstigAngle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_DefStep_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_DefStep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_DefCount_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_DefCount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
[Ny,Nx] = size(handles.wave0);
switch get(hObject,'Value')
    case 1
        img1 = sum(imread('amplitude.jpg'),3);
        img1 = img1/max(max(img1));
        handles.wave0 = img1.*exp(i*angle(handles.wave0));        
    case 2
        img1 = zeros(Ny,Nx);
        img1(floor(Ny/2)+1,floor(Nx/2)+1) = 1;
        img2 = angle(handles.wave0);
        handles.wave0 = img1.*exp(i*img2);       
    case 3
        img1 = zeros(Ny,Nx);
        img1(floor(Ny/2)+1,floor(Nx/2)+1+floor([-Nx/4 Nx/4])) = 1;
        img2 = angle(handles.wave0);
        handles.wave0 = img1.*exp(i*img2);       
end
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)

[Ny,Nx] = size(handles.wave0);
switch get(hObject,'Value')
    case 1
        img2 = sum(imread('phase.jpg'),3);
        % normalize the image intensities:
        img2 = img2/max(max(img2))-0.5;
        handles.wave0 = abs(handles.wave0).*exp(i*img2);
    case 2
        img2 = zeros(Ny,Nx);
        img2(floor(Ny/2)+1,floor(Nx/2)+1) = 1;
        img2 = img2/max(max(img2))-0.5;
        handles.wave0 = abs(handles.wave0).*exp(i*img2);
    case 3 % 2 delta functions
        img2 = zeros(Ny,Nx);
        img2(floor(Ny/2)+1,floor(Nx/2)+1+floor([-Nx/4 Nx/4])) = 1;
        img2 = img2/max(max(img2))-0.5;
        handles.wave0 = abs(handles.wave0).*exp(i*img2);
end
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ThicknessDefocusSeries.
function pushbutton_ThicknessDefocusSeries_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile('*.img', 'Select set of wave functions','MultiSelect','on');

filename

% find the number of files that has been selected
if iscell(filename)
    handles.numWaves = size(filename,2);
    handles.waveNames = filename;
    handles.wavePath = pathname;
    fprintf('Will process %d wave functions\n',handles.numWaves);
else
    if (isempty(filename))
        fprintf('No file selected!\n');
        return
    end
    handles.numWaves = 1;
    if (strcmp(filename,'wave_0_0.img') == 1)
        handles.numWaves = -1;  % let the user define number of thicknesses and TDS runs
    end
    [img,t,dx,dy,params,comment] = binread2D(fullfile(pathname,filename));
    handles.wave0 = img;
    handles.dx = dx;
    handles.dy = dy;
    handles.thickness = t;
    handles.wavePath = pathname;
end
    
pushbutton_FocSeries_Callback(hObject, eventdata, handles);

% --- Executes on button press in pushbutton_saveImg.
function pushbutton_saveImg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uiputfile('*.img', 'Save image as');
if isequal(filename,0) return; end
 
[Ny,Nx] = size(handles.wave0);
sampling = 10*str2double(get(handles.editSampling,'String'));
pushbuttonUpdate_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
dx = handles.dx;
dy = handles.dy;

handles = updateImages(handles);
img = abs(handles.Psi).^2;
if get(handles.checkbox_Vibration,'Value')
    kx2 = handles.kx2;
    ky2 = handles.ky2;
    k2  = (kx2.^2+ky2.^2);
    phi     = atan2(ky2,kx2);

    vib1 = (pi*handles.vibration_Axis1/(-log(0.5))).^2;
    vib2 = (pi*handles.vibration_Axis2/(-log(0.5))).^2;
    img = real(ifft2(fft2(img).*exp(-vib1*k2.*abs(cos(phi-handles.vibration_Angle))-vib2*k2.*abs(sin(phi-handles.vibration_Angle)))));
end


if get(handles.checkbox_Sampling,'Value')
    Nx2 = round((Nx-1)*dx/sampling)+1;
    Ny2 = round((Ny-1)*dy/sampling)+1;
    xi = sampling*[0:Nx2-1];
    yi = sampling*[0:Ny2-1];
    if xi(end) > dx*(Nx-1), xi(end) = dx*(Nx-1); end;
    if yi(end) > dy*(Ny-1), yi(end) = dy*(Ny-1); end;
    
    img2 = interp2(dx*[0:Nx-1],dy*[0:Ny-1].',img,xi,yi.','linear',0);
        
    binwrite2D(img2,fullfile(pathname, filename),sampling,sampling,handles.thickness,0,0);
else
    binwrite2D(img,fullfile(pathname, filename),dx,dy,handles.thickness,0,0);
end







% --- Executes on button press in pushbutton_FocSeries.
function pushbutton_FocSeries_Callback(hObject, eventdata, handles)
folder = uigetdir('', 'Select the target folder');
if isequal(folder,0) return; end
 
[Ny,Nx] = size(handles.wave0);
set(handles.checkboxDefocus,'Value',1);
% convert Defocus from nm to A
defocus0 = str2double(get(handles.editDefocus,'String'));
defCount = round(str2double(get(handles.edit_DefCount,'String')));
defStep  = str2double(get(handles.edit_DefStep,'String'));
sampling = 10*str2double(get(handles.editSampling,'String'));
pushbuttonUpdate_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
dx = handles.dx;
dy = handles.dy;
k2  = (handles.kx2.^2+handles.ky2.^2);
phi = atan2(handles.ky2,handles.kx2);

Ntds = 0;
Nthickness = handles.numWaves;
if (handles.numWaves < 0)  % means, we have read the first wave of a TDS series
    prompt={'Number of thicknesses (>=1):','Number of TDS runs (>=0):'};
    name='Select range of frames for Video output';
    numlines=1;
    defaultanswer={'1','0'};
 
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    Nthickness = str2num(cell2mat(answer(1)));
    Ntds = str2num(cell2mat(answer(2)));
end
    
for it = 0:Nthickness-1
    if (handles.numWaves > 1)
        filenameT = fullfile(handles.wavePath,handles.waveNames{it+1});
       % read the wave function:
       [newWave,t,dx,dy,params,comment] = binread2D(filenameT);
       handles.wave0 = newWave;
       handles.dx = dx;
       handles.dy = dy;
       handles.thickness = t;        
       fprintf('Processing wave %s at thickness %.1fnm\n',filenameT,0.1*t);
    elseif handles.numWaves < 0
        fprintf('Thickness index: %d\n',it);       
    end


    for idf=0:defCount-1
        def = defocus0+idf*defStep;
        set(handles.editDefocus,'String',def);
        handles.Defocus = 10*def;
        for itds=0:max(Ntds-1,0)
            if Ntds > 0
                filenameT = fullfile(handles.wavePath,sprintf('wave_%d_%d.img',itds,it));
                % read the wave function:
                [newWave,t,dx,dy,params,comment] = binread2D(filenameT);
                handles.wave0 = newWave;
                handles.dx = dx;
                handles.dy = dy;
                handles.thickness = t;
                % fprintf('Processing wave %s at thickness %.1fnm\n',filenameT,0.1*t);
            end                  
            handles = updateImages(handles);
            % handles = guidata(hObject);
            img = abs(handles.Psi).^2;
            if get(handles.checkbox_Vibration,'Value')
                vib1 = (pi*handles.vibration_Axis1/(-log(0.5))).^2;
                vib2 = (pi*handles.vibration_Axis2/(-log(0.5))).^2;
                img = real(ifft2(fft2(img).*exp(-vib1*k2.*abs(cos(phi-handles.vibration_Angle))-vib2*k2.*abs(sin(phi-handles.vibration_Angle)))));
            end

            % resample the image, if requested:
            if get(handles.checkbox_Sampling,'Value')
                Nx2 = round((Nx-1)*dx/sampling)+1;
                Ny2 = round((Ny-1)*dy/sampling)+1;
                xi = sampling*[0:Nx2-1];
                yi = sampling*[0:Ny2-1];
                if xi(end) > dx*(Nx-1), xi(end) = dx*(Nx-1); end;
                if yi(end) > dy*(Ny-1), yi(end) = dy*(Ny-1); end;
                
                img2 = interp2(dx*[0:Nx-1],dy*[0:Ny-1].',img,xi,yi.','linear',0);
                img = img2;
                
                if (handles.numWaves > 1)
                    % extract thickness index from wave name:
                    idString = handles.waveNames{it+1};
                    idString = idString(6:end-4);
                    if defCount == 1
                        fileName = sprintf('%s\\img_%s.img',folder,idString);                        
                    else
                        fileName = sprintf('%s\\img_%s_%d.img',folder,idString,idf);
                    end
                    fprintf('Image %s (thickness=%.1f nm, def=%g nm):  ',fileName,0.1*handles.thickness,def);                    
                else
                    if (Ntds == 0)
                        if defCount == 1
                            fileName = sprintf('%s\\img.img',folder);
                        else
                            fileName = sprintf('%s\\img_%d.img',folder,idf);
                        end
                        fprintf('Image %d (def=%gnm):  ',idf,def);
                    else
                        if (itds == 0)
                            fprintf('Thickness = %.1f nm, def = %g nm) ... \n    ... ',0.1*handles.thickness,def);
                        end
                    end
                end
                if (Ntds == 0)
                    binwrite2D(img2,fileName,sampling,sampling,handles.thickness,0,0);
                    
                else
                   if (itds == 0)
                       imgTDS = img2;
                   else
                       imgTDS = imgTDS+img2;
                   end
                end
            else
                if (handles.numWaves > 1)
                    idString = handles.waveNames{it+1};
                    idString = idString(6:end-4);
                    if defCount == 1
                        fileName = sprintf('%s\\img_%s.img',folder,idString);
                    else
                        fileName = sprintf('%s\\img_%s_%d.img',folder,idString,idf);
                    end
                    fprintf('Image %s (thickness=%.1fnm, def=%gnm):  ',fileName,0.1*handles.thickness,def);
                else
                    if defCount == 1
                        fileName = sprintf('%s\\img.img',folder);
                    else
                        fileName = sprintf('%s\\img_%d.img',folder,idf);
                    end                    
                    fprintf('Image %d (def=%gnm):  ',idf,def);
                end
                if (Ntds == 0)
                    binwrite2D(img,fileName,dx,dy,handles.thickness,0,0);
                else
                   if (itds == 0)
                       imgTDS = img2;
                   else
                       imgTDS = imgTDS+img2;
                   end
                end
            end            
        end % end of TDS loop
        if (Ntds > 0)
           imgTDS = imgTDS/Ntds;
           if defCount == 1
               fileName = sprintf('%s\\img_%d.img',folder,it);
           else
               fileName = sprintf('%s\\img_%d_%d.img',folder,it,idf);
           end
           fprintf('TDS-Image img_%d_%d.img (thickness=%.1fnm, def=%gnm):  ',it,idf,0.1*handles.thickness,def);
           if get(handles.checkbox_Sampling,'Value')
               binwrite2D(imgTDS,fileName,sampling,sampling,handles.thickness,0,0);
           else
               binwrite2D(imgTDS,fileName,dx,dy,handles.thickness,0,0);
           end
        end
    end % end of defocus loop
end % end of thickness loop
set(handles.editDefocus,'String',defocus0);
guidata(hObject, handles);
clear img
    
    
% --- Executes on button press in pushbutton_NewWindow.
function pushbutton_NewWindow_Callback(hObject, eventdata, handles)

pushbuttonUpdate_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
%    handles = updateImages(handles);
[Ny,Nx] = size(handles.wave0);

figure;
% ax = gca;
if get(handles.radiobuttonImage,'Value')  
    % show the image
    simIm = abs(handles.Psi).^2;
    imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],simIm);
    title(sprintf('Image t=%.1f nm [%.2f .. %.2f]',0.1*handles.thickness,min(min(simIm)),max(max(simIm))));
elseif get(handles.radiobutton_WaveAmplitude,'Value')
    simIm = abs(handles.Psi);
    imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],simIm);
    title(sprintf('Amplitude t=%.1f nm [%.2f .. %.2f]',0.1*handles.thickness,min(min(simIm)),max(max(simIm))));
elseif get(handles.radiobutton_WavePhase,'Value')
    simIm = angle(handles.Psi);
    imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],simIm);
    title(sprintf('Phase t=%.1f nm [%.2f .. %.2f]',0.1*handles.thickness,min(min(simIm)),max(max(simIm))));
else
    % show the diffractogram
    diffg = log(1+abs((fft2(abs(handles.Psi).^2))));
    diffg(1,1) = 0;
    imagesc(handles.dkx*10*[-Nx/2:Nx/2-1],handles.dky*10*[-Ny/2:Ny/2-1],fftshift(diffg));
    title('Diffractogram (log scale)');
end
% imagesc(handles.dx/10*[0:Nx-1],handles.dy/10*[0:Ny-1],abs(handles.Psi).^2,'Parent',gca);
colormap('gray');
set(gca,'YDir','normal');
if (get(handles.radiobuttonDiffractogram,'Value') == 0)
    xlabel('x in nm');
    ylabel('y in nm');
else
    xlabel('x in 1/nm');
    ylabel('y in 1/nm');    
end
colorbar;
axis equal; axis tight;

% --- Executes on button press in pushbutton_Scherzer.
function pushbutton_Scherzer_Callback(hObject, eventdata, handles)
handles.Cs = str2num(get(handles.editCs,'String'));
handles.AccVoltage = str2num(get(handles.editAccVoltage,'String'));
handles.Defocus = -0.1*sign(handles.Cs)*sqrt(1.5*abs(handles.Cs)*1e7*wavelength(handles.AccVoltage));

if abs(handles.Defocus) > 1
    set(handles.editDefocus,'String',sprintf('%.1f',handles.Defocus));    
else
    set(handles.editDefocus,'String',sprintf('%.4f',handles.Defocus));
end
guidata(hObject, handles);
pushbuttonUpdate_Callback(hObject, eventdata, handles);




% function [img,t,dx,dy] = binread2D(fileName,flag)
% file for reading binary data of different formats
% input: 
% fileName - name of file (string)
% printFlag - 0: non-verbose, 1: verbose (default)
% flag: for overriding file flags.
% output:
% img      - image (complex or real
% t        - thickness
% dx, dy   - pixel size (for calibration)
function [img,t,dx,dy,params,comment] = binread2D(fileName,printFlag,flag)

headerlength = 4;
if nargin < 2
    printFlag = 0;
end
% open the file and define the file ID (fid):
fid=fopen(fileName,'r','ieee-le');

% header = [headersize(bytes) paramSize commentSize Nx Ny complFlag doubleFlag dataSize version]

header = fread(fid,8,'int32');
Ny = header(5);
Nx = header(4);
t = fread(fid,1,'float64');
dx = fread(fid,1,'float64');
dy = fread(fid,1,'float64');


% read additional parameters from file, if any exist:
paramSize = header(2);
if (paramSize > 0)
    params = fread(fid,paramSize,'float64');
    % params
else
   params = []; 
end

% read comments from file, if any exist:
commentSize = header(3);
if (commentSize > 0)
    comment = fread(fid,commentSize,'char').';
    % comment(commentSize+1) = 0;
    % fprintf('Comment: %s\n',comment);
else 
    comment = '';
end


    % flag = 0;

integerFlag = 0;
complexFlag = header(6);
doubleFlag = (header(7) == 8*(complexFlag+1));
if (nargin <3)
    flag = integerFlag+doubleFlag*4+complexFlag*2;    
end

% fprintf('Header Size: %dbytes, t: %g, dx: %g, dy: %g\n',header(1),t,dx,dy);


% bit 2: double
% bit 1: complex
% bit 0: integer
if printFlag
    fprintf('binread2D %s: %d x %d pixels flag = %d (',fileName,Nx,Ny,flag);
end
switch bitand(flag,7)
case 0
    complexFlag = 0;
    doubleFlag = 0;
    if printFlag
        fprintf('32-bit real data, %.3fMB)\n',Nx*Ny*4/1048576);
    end
    img = fread(fid,[Ny,Nx],'float32');
case 4
    complexFlag = 0;
    doubleFlag = 1;
    if printFlag
    fprintf('64-bit real data, %.3fMB)\n',Nx*Ny*8/1048576);
    end
    img = fread(fid,[Ny,Nx],'float64');
case 2
    complexFlag = 1;
    doubleFlag = 0;
    if printFlag
    fprintf('32-bit complex data, %.3fMB)\n',Nx*Ny*8/1048576);
    end
    % img = fread(fid,[Nx,2*Ny],'float32');
    % img = img(1:Nx,1:2:2*Ny-1)+i*img(1:Nx,2:2:2*Ny);
    img = fread(fid,[2*Ny,Nx],'float32');
    img = img(1:2:2*Ny-1,1:Nx)+i*img(2:2:2*Ny,1:Nx);
case 6
    complexFlag = 1;
    doubleFlag = 1;
    if printFlag
    fprintf('64-bit complex data, %.3fMB)\n',Nx*Ny*16/1048576);
    end
    img = fread(fid,[2*Ny,Nx],'float64');
    img = img(1:2:2*Ny-1,1:Nx)+i*img(2:2:2*Ny,1:Nx);
case 1
    complexFlag = 0;
    doubleFlag = 0;
    if printFlag
        fprintf('16-bit integer data, %.3fMB)\n',Nx*Ny*8/1048576);
    end
    img = fread(fid,[Ny,Nx],'int16');    
case 5
    complexFlag = 0;
    doubleFlag = 0;
    if printFlag
    fprintf('32-bit integer data)\n');
    end
    img = fread(fid,[Ny,Nx],'int32');    
end
% imagesc(img); colormap('gray'); axis equal; axis tight;

fclose(fid);














% --- Executes on button press in checkbox_Sampling.
function checkbox_Sampling_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Sampling








% --- Executes on button press in pushbutton_SaveWave.
function pushbutton_SaveWave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SaveWave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uiputfile('*.img', 'Save wave function as');
if isequal(filename,0) return; end
 
[Ny,Nx] = size(handles.wave0);
set(handles.checkboxDefocus,'Value',1);
% convert Defocus from nm to A
defocus0 = str2double(get(handles.editDefocus,'String'));
defCount = round(str2double(get(handles.edit_DefCount,'String')));
defStep  = str2double(get(handles.edit_DefStep,'String'));
sampling = 10*str2double(get(handles.editSampling,'String'));
pushbuttonUpdate_Callback(hObject, eventdata, handles);
handles = guidata(hObject);
dx = handles.dx;
dy = handles.dy;


if get(handles.checkbox_Sampling,'Value')
    Nx2 = round((Nx-1)*dx/sampling)+1;
    Ny2 = round((Ny-1)*dy/sampling)+1;
    xi = sampling*[0:Nx2-1];
    yi = sampling*[0:Ny2-1];
    if xi(end) > dx*(Nx-1), xi(end) = dx*(Nx-1); end;
    if yi(end) > dy*(Ny-1), yi(end) = dy*(Ny-1); end;
    
    Psi2 = interp2(dx*[0:Nx-1],dy*[0:Ny-1].',Psi,xi,yi.','linear',0);
        
    binwrite2D(Psi2,fullfile(pathname, filename),sampling,sampling,handles.thickness,0,0);
else
    binwrite2D(handles.Psi,fileName,dx,dy,handles.thickness,0,0);
end






function handles = readAllFields(handles)
handles.MinAmplitude  = str2double(get(handles.editMinAmplitude,'String'));
handles.MaxPhase = str2double(get(handles.editMaxPhase,'String'));
handles.Amorph   = str2double(get(handles.editAmorph,'String'));
handles.highVoltage = str2double(get(handles.editAccVoltage,'String'));
handles.lambda   = wavelength(handles.highVoltage);
set(handles.textAccVoltageUnit,'String',sprintf('kV, lambda=%.4f A',handles.lambda));
% Read Defocus, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxDefocus,'Value')
    % convert Defocus from nm to A
    handles.Defocus = 10*str2double(get(handles.editDefocus,'String'));
else, handles.Defocus=0; end
% Read Cs, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxCs,'Value')
    % convert Cs from mm to A
    handles.Cs = 1e7*str2double(get(handles.editCs,'String'));
else, handles.Cs=0; end
% Read Tilt, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxTilt,'Value')
    % convert Tilt from mrad to rad
    handles.TiltX = 1e-3*str2double(get(handles.editTiltX,'String'));
    handles.TiltY = 1e-3*str2double(get(handles.editTiltY,'String'));
else, handles.TiltX=0; handles.TiltY=0; end
if get(handles.checkboxComa,'Value')
    % convert Coma from nm and deg to A and rad
    handles.Coma = 1e4*abs(str2double(get(handles.editComa,'String')))*...
        exp(i*pi/180*str2double(get(handles.editComaAngle,'String')));
else, handles.Coma=0;  end
if get(handles.checkbox3foldAstig,'Value')
    % convert Coma from nm and deg to A and rad
    handles.Astig3 = 1e4*abs(str2double(get(handles.edit3foldAstig,'String')))*...
        exp(i*pi/180*str2double(get(handles.edit3foldAstigAngle,'String')));
else, handles.Astig3=0; end
% Read Astigmatism, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxAstig,'Value')
    % convert Tilt from mrad to rad
    handles.Astig = 10*str2double(get(handles.editAstig1,'String'))*...
        exp(i*pi/180*str2double(get(handles.editAstig2,'String')));
else, handles.Astig=0; end

% Read Delta, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxDelta,'Value')
    % convert Cs from mm to A
    handles.Delta = 10*str2double(get(handles.editDelta,'String'));
else, handles.Delta=0; end
% Read Convergence angle, if the corresponding checkbox is "on", set it to 0 otherwise:
if get(handles.checkboxConvAngle,'Value')
    % convert Cs from mm to A
    handles.Convergence = 1e-3*str2double(get(handles.editConvAngle,'String'));
else, handles.Convergence=0; end
% Read objective aperture angle:
if get(handles.checkbox_ObjAperture,'Value')
    % convert angle from mrad to rad
    handles.ObjApertureTheta = 1e-3*str2double(get(handles.edit_ObjApertureTheta,'String'));
else, handles.ObjApertureTheta=0; end
% Read sample vibration
if get(handles.checkbox_Vibration,'Value')
    % convert vibration from pm to A:
    handles.vibration_Axis1 = 0.1*str2double(get(handles.edit_Vibration_Axis1,'String'));
    handles.vibration_Axis2 = 0.1*str2double(get(handles.edit_Vibration_Axis2,'String'));
    handles.vibration_Angle = pi/180*str2double(get(handles.edit_Vibration_Angle,'String'));
    handles.applyVibration = 1;
else
    handles.vibration_Axis1 = 0;
    handles.vibration_Axis2 = 0;
    handles.vibration_Angle = 0;
    handles.applyVibration = 0;
end




% --- Executes on button press in pushbutton_moreAberrations.
function pushbutton_moreAberrations_Callback(hObject, eventdata, handles)
% handles = readAllFields(hObject, eventdata, handles);
% handles = readAllFields(handles);
handles = pushbuttonUpdate_Callback(hObject, eventdata, handles);

handles.c(2)     = handles.Defocus;
handles.c(4)     = handles.Cs;
handles.phi(2,2) = angle(handles.Astig);
handles.a(2,2)   = abs(handles.Astig);
handles.a(3,1)   = abs(handles.Coma);
handles.phi(3,1) = angle(handles.Coma);
handles.a(3,3)   = abs(handles.Astig3);
handles.phi(3,3) = angle(handles.Astig3);

[a,phi,c,params] = aberrations_TEM(handles,fft2(handles.wave));
if length(params) > 1
    handles.highVoltage = params(5);
    handles.Convergence    = params(6);
    handles.Delta    = params(7);
    handles.a        = a;
    handles.phi      = phi;
    handles.Defocus  = 1e-1*c(2);
    handles.Cs       = c(4);
    handles.Astig    = 1e-1*a(2,2)*exp(i*phi(2,2));
    handles.Coma     = a(3,1)*exp(i*phi(3,1));
    handles.Astig3   = a(3,3)*exp(i*phi(3,3));
    % Coma
    set(handles.editComa,'String',1e-4*abs(handles.Coma))
    set(handles.editComaAngle,'String',180/pi*angle(handles.Coma))
    if abs(handles.Coma) > 0,  set(handles.checkboxComa,'Value',1); end
    % Astig3
    set(handles.edit3foldAstig,'String',1e-4*abs(handles.Astig3));
    set(handles.edit3foldAstigAngle,'String',180/pi*angle(handles.Astig3));
    if abs(handles.Astig3) > 0,  set(handles.checkbox3foldAstig,'Value',1); end
    % Defocus, 2-fold Astig, etc.
    set(handles.editDefocus,'String',sprintf('%.1f',handles.Defocus));
    if abs(handles.Defocus) > 0,  set(handles.checkboxDefocus,'Value',1); end
    
    set(handles.editAstig1,'String',sprintf('%.1f',abs(handles.Astig)));
    set(handles.editAstig2,'String',sprintf('%.1f',angle(handles.Astig)*180/pi));
    if abs(handles.Astig) > 0,  set(handles.checkboxAstig,'Value',1); end
    
    set(handles.editCs,'String',sprintf('%.4f',1e-7*handles.Cs));
    if abs(handles.Cs) > 0,  set(handles.checkboxCs,'Value',1); end
    
    set(handles.editConvAngle,'String',handles.Convergence);
    if abs(handles.Convergence) > 0,  set(handles.checkboxConvAngle,'Value',1); end
    
    set(handles.editDelta,'String',0.1*handles.Delta);
    if abs(handles.Delta) > 0,  set(handles.checkboxDelta,'Value',1); end

    set(handles.editAccVoltage,'String',handles.highVoltage);
    guidata(hObject, handles);
    handles = pushbuttonUpdate_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
end



% --- Executes on button press in checkbox_Vibration.
function checkbox_Vibration_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);



function edit_Vibration_Axis1_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_Vibration_Axis1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Vibration_Axis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Vibration_Axis2_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_Vibration_Axis2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Vibration_Axis2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Vibration_Angle_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_Vibration_Angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Vibration_Angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_IncludeHigherOrders.
function checkbox_IncludeHigherOrders_Callback(hObject, eventdata, handles)
pushbuttonUpdate_Callback(hObject, eventdata, handles);


