function varargout = qstem(varargin)
% QSTEM M-file for qstem.fig
%      QSTEM, by itself, creates a new QSTEM or raises the existing
%      singleton*.
%
%      H = QSTEM returns the handle to a new QSTEM or the handle to
%      the existing singleton*.
%
%      QSTEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QSTEM.M with the given input arguments.
%
%      QSTEM('Property','Value',...) creates a new QSTEM or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before qstem_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to qstem_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qstem

% Last Modified by GUIDE v2.5 14-Apr-2011 11:25:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @qstem_OpeningFcn, ...
                   'gui_OutputFcn',  @qstem_OutputFcn, ...
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


% --- Executes just before qstem is made visible.
function qstem_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to qstem (see VARARGIN)

% Choose default command line output for qstem
handles.modeNum = 1;
handles.mode    = 'STEM';



handles.output = hObject;

handles.binPath = '';
% Root: HKLM; Subkey: "SOFTWARE\TEM Solutions\QSTEM"; ValueType: string; ValueName: "Installation Directory"; ValueData: "{app}"
% try, handles.binPath  = winqueryreg('HKEY_LOCAL_MACHINE', 'SOFTWARE\TEM Solutions\QSTEM', 'Installation Directory'); catch, handles.binPath  = pwd; end;


handles           = readAllFields(hObject, eventdata, handles);
set(gcf,'Name','qstem V2.15');
handles.Mm        = eye(3);
handles.Detectors = [70.0 200.0 0 0; 0.0 40.0 0 0];  % initialize 1: ADF, 2: BF 
handles.atomPos   = [];
handles.RbboxHandle     = -1;
handles.saveLevel     = 0;
handles.printLevel    = 2;
handles.atomRadius    = 5.0;
handles.savePotential = 0;
handles.saveTotalPotential = 0;
handles.propProgInterval = 10;
handles.potProgInterval  = 1000;
handles.a = zeros(6);
handles.phi = zeros(6);
cla(handles.axes_model);
set(handles.axes_model,'YDir','reverse');
handles.atomSize = 15;
set(handles.edit_atomSize,'String',handles.atomSize);

handles.configPath = pwd();
handles.configFile = 'qstem.qsc';

if (length(varargin) > 0)
    inputFile = varargin{1};
    if isempty(findstr(inputFile,':\'))
        inputFile = fullfile(pwd(),inputFile);
    end
    fprintf('Starting configuration: %s\n',inputFile);

    [filepath,filename,ext] = fileparts(inputFile);
    handles.configPath = filepath;
    % cd(filepath);
    
    if (ext == '.cfg') % load a model
        handles.posFileName = varargin{1};
        handles = pushbuttonLoadModel_Callback(hObject, eventdata, handles,0);
    else % load a simulation config
        handles.configFile = [filename, ext];
        handles = pushbutton_LoadConfigFile_Callback(hObject, eventdata, handles,0);        
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes qstem wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = qstem_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if nargout > 0
    varargout{1} = handles.output;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton_LoadConfigFile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = pushbutton_LoadConfigFile_Callback(hObject, eventdata, handles,askUser)
if (nargin < 4)
    askUser = 1;
end



configFile = fullfile(handles.configPath, handles.configFile);
% handles = readConfigFile(configFile,handles,askUser)
% if we use the above function we only need to keep those lines below which
% set the controls values.

pathname = handles.configPath;
filename = handles.configFile;

if (askUser)
    [filename, pathname] = uigetfile({'*.qsc;*.dat', 'QSTEM simulation config file (*.qsc, *.dat)';...
        '*.*', 'All files (*.*)'},...
        'Select a file',configFile);
    if isequal(filename,0) || isequal(pathname,0)
        % disp('User pressed cancel')
        return;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if a valid file has been selected:
    configFile = fullfile(pathname, filename);
    handles.configPath = pathname;
    handles.configFile = filename;
end
oldPath = pwd();
cd(pathname)
fid = fopen(filename);
if fid == -1
    % Display message box only if the user has selected the file:
    if (askUser)
        msgbox(sprintf('Could not open file %s\n',configFile));
    end
    cd(oldPath);
    return
end


% mode:  STEM
C = textscan(fid,'%s%q','delimiter',':','commentStyle','%');
fclose(fid);
pNames  = C{1};
pValues = C{2};
Nparams = length(pNames);
nameLen = 0;
valueLen  = 0;
for j=1:Nparams
    name  = cell2mat(pNames(j));
    value = sscanf(cell2mat(pValues(j)),'%s',1);
    if nameLen < length(name), nameLen = length(name); end
    if valueLen < length(value), valueLen = length(value); end
    % nval  = str2num(value);
    % fprintf('%d: %s: %s (%f)\n',j,name,value,nval);
end
nameArray = zeros(Nparams,nameLen);
for j=1:Nparams
    name  = cell2mat(pNames(j));
    nameArray(j,1:length(name)) = name;
end
% nameArray
% determine the index for the parameter that defines the mode STEM:
% seekStr = ['mode' zeros(1,nameLen-length('mode'))];
% indPar = find(sum(nameArray-repmat(seekStr,Nparams,1),2) == 0);
% value = sscanf(cell2mat(pValues(indPar)),'%s',1);
value = extractParam('mode',nameArray,pValues);
if strcmp(value,'STEM') == 1
    handles.mode = 'STEM';
    handles.modeNum = 1;
    set(handles.uipanel_ScanningWindow,'Title','Scan Window');
    set(handles.radiobutton_STEM,'Value',1);
elseif strcmp(value,'TEM') == 1
    handles.mode = 'TEM';
    handles.modeNum = 0;
    % set(handles.uipanel_ScanningWindow,'Title','Center');
    set(handles.radiobutton_TEM,'Value',1);
elseif strcmp(value,'CBED') == 1
    handles.mode = 'CBED';
    handles.modeNum = 2;
    % set(handles.uipanel_ScanningWindow,'Title','Center');
    set(handles.radiobutton_CBED,'Value',1);
else
    msgbox('This file does not seem to be a valid input file');
    cd(oldPath);
    return
end
% fprintf('Verified configuration for STEM mode!\n');


handles.posFileName = extractParam('filename',nameArray,pValues);
handles.NcellX = str2num(extractParam('NCELLX',nameArray,pValues));
handles.NcellY = str2num(extractParam('NCELLY',nameArray,pValues));
value = extractParam('NCELLZ',nameArray,pValues);
if isempty(findstr(value,'/'))
    handles.NcellZ = str2num(value);
    handles.NsubSlabs = 1;
else
    ind = findstr(value,'/');
    handles.NcellZ = str2num(value(1:ind(1)-1));
    handles.NsubSlabs = str2num(value(ind(1)+1:end));
end
handles.Nslice = str2num(extractParam('slices',nameArray,pValues));
handles.PotentialOffsetX = str2num(extractParam('xOffset',nameArray,pValues));
handles.PotentialOffsetY = str2num(extractParam('yOffset',nameArray,pValues));
handles.PotentialOffsetZ = str2num(extractParam('zOffset',nameArray,pValues));

handles.SliceThickness = str2num(extractParam('slice-thickness',nameArray,pValues));
handles.CenterSlices = 0;
value = extractParam('center slices',nameArray,pValues);
if value(1) == 'y'
    handles.CenterSlices = 1;
end
set(handles.checkbox_CenterSlices,'Value',handles.CenterSlices);

% Periodicity of super-cell:
handles.PeriodicZ = 0;
handles.PeriodicXY = 0;
value = extractParam('periodicXY',nameArray,pValues);
if value(1) == 'y'
    handles.PeriodicXY = 1;
end
set(handles.checkbox_PeriodicXY,'Value',handles.PeriodicXY);
value = extractParam('periodicZ',nameArray,pValues);
if value(1) == 'y'
    handles.PeriodicZ = 1;
end
set(handles.checkbox_PeriodicZ,'Value',handles.PeriodicZ);


% number of cells:
set(handles.edit_NcellX,'String',sprintf('%d',handles.NcellX));
set(handles.edit_NcellY,'String',sprintf('%d',handles.NcellY));
set(handles.edit_NcellZ,'String',sprintf('%d',handles.NcellZ));
set(handles.edit_NsubSlabs,'String',sprintf('%d',handles.NsubSlabs));
set(handles.edit_Nslice,'String',sprintf('%d',handles.Nslice));
set(handles.edit_PotentialOffsetZ,'String',sprintf('%.3f',handles.PotentialOffsetZ));
set(handles.edit_PotentialOffsetY,'String',sprintf('%.3f',handles.PotentialOffsetY));
set(handles.edit_PotentialOffsetX,'String',sprintf('%.3f',handles.PotentialOffsetX));
handles.NsliceTot      = handles.NsubSlabs*handles.Nslice;
set(handles.text_NsliceTot,'String',sprintf('(Total: %d)',handles.NsliceTot));


% TDS settings
handles.TDS = 0;
value = extractParam('tds',nameArray,pValues);
if value(1) == 'y'
    handles.TDS = 1;
end
set(handles.checkbox_TDS,'Value',handles.TDS);
handles.Temperature = str2num(extractParam('temperature',nameArray,pValues));
handles.TDSruns = str2num(extractParam('Runs for averaging',nameArray,pValues));

% extract Box definition:
set(handles.radiobutton_Ncells,'Value',1);
handles.BoxMode = 0;
value = extractParam('Cube',nameArray,pValues);
if ~isempty(value)
    handles.BoxX = value(1);
    handles.BoxY = value(2);
    handles.BoxZ = value(3);
    set(handles.radiobutton_Box,'Value',1);
    set(handles.edit_NcellX,'String',handles.BoxX);
    set(handles.edit_NcellY,'String',handles.BoxY);
    set(handles.edit_NcellZ,'String',handles.BoxZ);
    if handles.SliceThickness == 0
        handles.SliceThickness = handles.BoxZ/handles.NsliceTot;
    end
    handles.BoxMode = 1;
    set(handles.text_Nx,'String','Box:         ax:');
    set(handles.text_Ny,'String','by:');
    set(handles.text_Nz,'String','cz:');
else
    if handles.SliceThickness == 0
        handles.SliceThickness = handles.Mm(3,3)*handles.NcellZ/handles.NsliceTot;
    end
end
set(handles.edit_SliceThickness,'String',sprintf('%.4f',handles.SliceThickness));


handles.TiltX = 180/pi*str2num(extractParam('Crystal tilt X',nameArray,pValues));
handles.TiltY = 180/pi*str2num(extractParam('Crystal tilt Y',nameArray,pValues));
handles.TiltZ = 180/pi*str2num(extractParam('Crystal tilt Z',nameArray,pValues));

handles.BeamTiltX = str2num(extractParam('Beam tilt X',nameArray,pValues));
handles.BeamTiltY = str2num(extractParam('Beam tilt Y',nameArray,pValues));
value = extractParam('Tilt back',nameArray,pValues);
if value(1) == 'y'
    handles.TiltBack = 1;
end

handles.AstigMag = str2num(extractParam('astigmatism',nameArray,pValues));
handles.AstigAngle = str2num(extractParam('astigmatism angle',nameArray,pValues));

set(handles.edit_TiltX,'String',handles.TiltX);
set(handles.edit_TiltY,'String',handles.TiltY);
set(handles.edit_TiltZ,'String',handles.TiltZ);
set(handles.edit_BtiltX,'String',handles.BeamTiltX);
set(handles.edit_BtiltY,'String',handles.BeamTiltY);
set(handles.checkbox_TiltBack,'Value',handles.TiltBack);
set(handles.edit_AstigMag,'String',handles.AstigMag);
set(handles.edit_AstigAngle,'String',handles.AstigAngle);


% microscope parameters:
handles.HighVoltage = str2num(extractParam('v0',nameArray,pValues));
set(handles.edit_HighVoltage,'String',handles.HighVoltage);
handles.Wavelength = wavelength(handles.HighVoltage);
set(handles.text_wavelength,'String',sprintf('kV    (wavelength = %.2fpm)',100*handles.Wavelength));

value = extractParam('defocus',nameArray,pValues);
handles.Defocus = str2num(value);
if abs(handles.Defocus) < 1
    set(handles.edit_Defocus,'String',sprintf('%.4f',handles.Defocus));
else
    set(handles.edit_Defocus,'String',sprintf('%.1f',handles.Defocus));
end

handles.alpha = str2num(extractParam('alpha',nameArray,pValues));
set(handles.edit_alpha,'String',handles.alpha);
handles.C3 = str2num(extractParam('Cs',nameArray,pValues));
set(handles.edit_C3,'String',handles.C3);

set(handles.edit_AstigMag,'String',handles.AstigMag);
set(handles.edit_AstigAngle,'String',handles.AstigAngle);
set(handles.edit_C3,'String',handles.C3);
set(handles.edit_alpha,'String',handles.alpha);


handles.C5 = str2num(extractParam('C5',nameArray,pValues));  % C5 in mm
handles.c(6) = 1e7*handles.C5;
% set(handles.edit_C5,'String',handles.C5);

% other aberrations:
% Write all the other aberrations as well:
for ix=3:6
    for iy=1:6
        % str = sprintf('a[%d,%d]',ix,iy);
        % str = sprintf('a_%d%d',ix,iy);
        value = extractParam(sprintf('a_%d%d',ix,iy),nameArray,pValues);
        % if strcmp(value,'0') == 0
        if str2num(value) > 0
            handles.a(ix,iy) = str2num(value);
            % handles.phi(ix,iy) = pi/180*str2num(extractParam(sprintf('phi[%d,%d]',ix,iy),nameArray,pValues));
            handles.phi(ix,iy) = pi/180*str2num(extractParam(sprintf('phi_%d%d',ix,iy),nameArray,pValues));
            % fprintf('Found: a(%d,%d) = %g, phi = %g\n',ix,iy,handles.a(ix,iy),handles.phi(ix,iy));
        end
    end
end


% Cc in mm:
handles.Cc   = str2num(extractParam('Cc',nameArray,pValues));
handles.dE_E = str2num(extractParam('dV/V',nameArray,pValues));
handles.dE    = handles.dE_E*handles.HighVoltage*1e3;
set(handles.edit_dE,'String',handles.dE);
% Delta in nm:
handles.Delta   = handles.Cc*1e3*handles.dE/handles.HighVoltage;

% folder and output options:
handles.SlicesBetweenOutputs = str2num(extractParam('slices between outputs',nameArray,pValues));
handles.OutputFolder = extractParam('Folder',nameArray,pValues);
set(handles.edit_OutputFolder,'String',handles.OutputFolder);
set(handles.edit_SlicesBetweenOutputs,'String',handles.SlicesBetweenOutputs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the coordinates:
% Load the structure model:
% check to see whether the file needs extension by the path of the
% config file:
filename = fullfile(handles.configPath,handles.posFileName);
fid = fopen(filename);
if fid ~= -1
    handles.posFileName = filename;
    fclose(fid);
end
cd(oldPath);

if (0)
    [coords, aType, Mm,DW,occ,charge] = readCFG_qstem(handles.posFileName,[0 0 0],0);
    if isempty(coords)
        msgbox(sprintf('Could not open .cfg file %s!',handles.posFileName));
    else
        set(handles.uipanel_model,'Title',sprintf('Model: %s',handles.posFileName));
        handles.atomPos = [aType coords DW occ charge];
        clear coords aType
        handles.Mm = Mm;
    end
else
    handles = pushbuttonLoadModel_Callback(hObject, eventdata, handles,0,0);
end

% Read the probe array parameters:
handles.ProbeNx = str2num(extractParam('nx',nameArray,pValues));
handles.ProbeNy = str2num(extractParam('ny',nameArray,pValues));
if handles.ProbeNy <= 0
    handles.ProbeNy = handles.ProbeNx;
end
handles.ProbeResolutionX = str2num(extractParam('resolutionX',nameArray,pValues));
handles.ProbeResolutionY = str2num(extractParam('resolutionY',nameArray,pValues));
set(handles.edit_ProbeNx,'String',handles.ProbeNx);
set(handles.edit_ProbeNy,'String',handles.ProbeNy);

% fprintf('Resol: %g %g %d %d %g %g\n',handles.ProbeResolutionX,handles.ProbeResolutionY,handles.ProbeNx,handles.ProbeNy,handles.Mm(1,1),handles.Mm(2,2));
if handles.ProbeResolutionX <= 1e-4
    handles.ProbeResolutionX = handles.Mm(1,1)/handles.ProbeNx;
end
if handles.ProbeResolutionY <= 1e-4
    handles.ProbeResolutionY = handles.Mm(2,2)/handles.ProbeNy;
end
set(handles.edit_ProbeResolutionX,'String',handles.ProbeResolutionX);
set(handles.edit_ProbeResolutionY,'String',handles.ProbeResolutionY);

% Define the additional parameters:
handles.ProbeWindowX = handles.ProbeResolutionX*handles.ProbeNx;
handles.ProbeWindowY = handles.ProbeResolutionY*handles.ProbeNy;
set(handles.edit_ProbeWindowX,'String',handles.ProbeWindowX);
set(handles.edit_ProbeWindowY,'String',handles.ProbeWindowX);

handles.ProbeMaxAngleX = 2/6*handles.Wavelength * (handles.ProbeNx/handles.ProbeWindowX);
handles.ProbeMaxAngleY = 2/6*handles.Wavelength * (handles.ProbeNy/handles.ProbeWindowY);
set(handles.edit_ProbeMaxAngleX,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleX));
set(handles.edit_ProbeMaxAngleY,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleY));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scan parameters:
switch handles.modeNum
    case 1
        handles.Xstart = str2num(extractParam('scan_x_start',nameArray,pValues));
        handles.Xstop  = str2num(extractParam('scan_x_stop',nameArray,pValues));
        handles.Xpixels = str2num(extractParam('scan_x_pixels',nameArray,pValues));
        
        handles.Ystart = str2num(extractParam('scan_y_start',nameArray,pValues));
        handles.Ystop  = str2num(extractParam('scan_y_stop',nameArray,pValues));
        handles.Ypixels = str2num(extractParam('scan_y_pixels',nameArray,pValues));
        
        % Read the detector configuration:
        % handles.Detectors = [70.0 200.0 0 0; 0.0 40.0 0 0];  % initialize 1: ADF, 2: BF
        handles.Detectors = extractParam('detector',nameArray,pValues);
        set(handles.edit_DetectorNumber,'String',1);
        if size(handles.Detectors,1) > 0
            j = 1;
            set(handles.edit_InnerAngle,'String',handles.Detectors(j,1));
            set(handles.edit_OuterAngle,'String',handles.Detectors(j,2));
            set(handles.edit_OffsetX,'String',handles.Detectors(j,3));
            set(handles.edit_OffsetY,'String',handles.Detectors(j,4));
        else
            set(handles.edit_InnerAngle,'String',' ');
            set(handles.edit_OuterAngle,'String',' ');
            set(handles.edit_OffsetX,'String',' ');
            set(handles.edit_OffsetY,'String',' ');
        end
    case 2  % CBED
        handles.Xstart = str2num(extractParam('scan_x_start',nameArray,pValues));
        handles.Xstop  = handles.Xstart;
        handles.Xpixels = 1;
        
        handles.Ystart = str2num(extractParam('scan_y_start',nameArray,pValues));
        handles.Ystop  = handles.Ystart;
        handles.Ypixels = 1;
        set(handles.edit_InnerAngle,'String',' ');
        set(handles.edit_OuterAngle,'String',' ');
        set(handles.edit_OffsetX,'String',' ');
        set(handles.edit_OffsetY,'String',' ');
        
    case 0  % TEM mode:
        handles.Xstart = 0.5*handles.ProbeWindowX;
        handles.Xstop = handles.Xstart;
        handles.Xpixels = 1;
        handles.Ystart = 0.5*handles.ProbeWindowY;
        handles.Ystop = handles.Ystart;
        handles.Ypixels = 1;
        set(handles.edit_InnerAngle,'String',' ');
        set(handles.edit_OuterAngle,'String',' ');
        set(handles.edit_OffsetX,'String',' ');
        set(handles.edit_OffsetY,'String',' ');
end
set(handles.edit_Xstart,'String',handles.Xstart);
set(handles.edit_Xstop,'String',handles.Xstop);
set(handles.edit_Xpixels,'String',handles.Xpixels);
set(handles.edit_Ystart,'String',handles.Ystart);
set(handles.edit_Ystop,'String',handles.Ystop);
set(handles.edit_Ypixels,'String',handles.Ypixels);




handles.saveLevel     = str2num(extractParam('save level',nameArray,pValues));
handles.printLevel    = str2num(extractParam('print level',nameArray,pValues));
handles.atomRadius    = str2num(extractParam('atom radius',nameArray,pValues));
handles.savePotential   = strcmp(extractParam('save potential',nameArray,pValues),'yes');
handles.saveTotalPotential   = strcmp(extractParam('save projected potential',nameArray,pValues),'yes');
handles.potProgInterval = str2num(extractParam('potential progress interval',nameArray,pValues));
handles.propProgInterval = str2num(extractParam('propagation progress interval',nameArray,pValues));

guidata(hObject, handles);
switch handles.modeNum
    case 2
        handles = radiobutton_CBED_Callback(hObject, eventdata, handles);        
    case 1
        handles = radiobutton_STEM_Callback(hObject, eventdata, handles);
    case 0
        handles = radiobutton_TEM_Callback(hObject, eventdata, handles);
end
guidata(hObject, handles);
% end of LoadConfigFile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = extractParam(name,nameArray,pValues);
[Nparams,nameLen] = size(nameArray);
if strcmp(name,'detector') == 0
    if strcmp(name,'Cube') == 0
        indPar = find(sum(abs(nameArray-repmat([name zeros(1,nameLen-length(name))],Nparams,1)),2) == 0);
        if isempty(indPar)
            value = '0';
            return
        else
            % pValues(indPar(1))
            value = sscanf(cell2mat(pValues(indPar(1))),'%s',1);
            % fprintf('%s (%d): %s\n',name,length(name),value);
        end
    else % extracting Box information:
        indPar = find(sum(abs(nameArray-repmat([name zeros(1,nameLen-length(name))],Nparams,1)),2) == 0);
        if isempty(indPar)
            value = [];
        else
            value = sscanf(cell2mat(pValues(indPar(1))),'%f',3);
        end
    end
else 
   % fprintf('Looking for detectors\n')
   % value = [70.0 200.0 0 0; 0.0 40.0 0 0];    
   indPar = find(sum(abs(nameArray-repmat([name zeros(1,nameLen-length(name))],Nparams,1)),2) == 0);
   value = [];
   for j=1:length(indPar)
       det = textscan(cell2mat(pValues(indPar(j))),'%f %f %s %f %f');
       if prod(size(det{4})) > 0
           value(j,1:4) = [det{1} det{2} det{4} det{5}];
       else
           value(j,1:4) = [det{1} det{2} 0 0];
       end
   end
end


% --- Executes on button press in pushbutton_SaveConfigFile.
function pushbutton_SaveConfigFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SaveConfigFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check for TDS flag:
if handles.modeNum == 0
    handles.TDS         = get(handles.checkbox_TDS,'Value');
    if handles.TDS
        ButtonName=questdlg({'Are you shure about including TDS?','This will greatly increase the computational effort.'}, ...
            'Verify TDS flag', ...
            'Keep TDS on','Turn TDS off','Turn TDS off');
        switch ButtonName
            case 'Turn TDS off',
                handles.TDS = 0;
                set(handles.checkbox_TDS,'Value',handles.TDS);
                pause(0.2);
        end % switch
    end
end


configFile = fullfile(handles.configPath, handles.configFile);
[filename, pathname] = uiputfile({'*.qsc;*.dat', 'QSTEM simulation config files (*.qsc, *.dat)';...
    '*.*', 'All files (*.*)'},...
    'Save file',configFile);
if isequal(filename,0) || isequal(pathname,0)
    % disp('User pressed cancel')
else
    % Attach proper ending if necessary:
    if isempty(length(filename) - findstr(filename,'.qsc'))
        filename = [filename '.qsc'];
    end
    handles.configPath = pathname;
    handles.configFile = filename;
    configFile = fullfile(pathname, filename);
    handles = readAllFields(hObject, eventdata, handles);
    guidata(hObject, handles);
    writeFields2File(configFile,handles);
end

% --- Executes on button press in pushbutton_RunStem3.
function pushbutton_RunStem3_Callback(hObject, eventdata, handles)
configFile = fullfile(handles.configPath, handles.configFile);
[filename,folder] = uigetfile({'*.qsc;*.dat','QSTEM simulation config files (*.qsc, *.dat)'; ...
        '*.*',  'All Files (*.*)'},'Select config file to execute',configFile);
if ~(isequal(filename,0) || isequal(folder,0))
    oldFolder = pwd();
    % cmd = sprintf('cd "%s" & "%s" %s &',folder,fullfile(handles.binPath,'stem3.exe'),filename);
    systemStr = computer();
    switch systemStr
        case 'PCWIN'
            cd(folder);
            cmd = sprintf('cd "%s" & stem3 %s &',folder,filename);
            system(cmd);
            cd(oldFolder);
        case 'GLNX86'           
            warndlg('PC Linux is not yet supported\nPlease run stem3 from the command line','Running stem3','modal');
            
        otherwise
            warndlg(sprintf('This system (%s) is not supported\nPlease run stem3 from the command line',systemStr),'Running stem3','modal')
    end
end


% --- Executes on button press in pushbutton_DisplayResults.
function pushbutton_DisplayResults_Callback(hObject, eventdata, handles)

if isfield(handles,'OutputFolder')
    if isempty(findstr(handles.OutputFolder,':\'))
        % system(sprintf('showimage %s',fullfile(handles.configPath,handles.OutputFolder)));
        % system(sprintf('"%s" %s',fullfile(handles.binPath,'showimage'),fullfile(handles.configPath,handles.OutputFolder)));
        system(sprintf('showimage %s',fullfile(handles.configPath,handles.OutputFolder)));
    else
        % system(sprintf('showimage %s',handles.OutputFolder));
        % system(sprintf('"%s" %s',fullfile(handles.binPath,'showimage'),handles.OutputFolder));
        system(sprintf('"%s" %s','showimage',handles.OutputFolder));
    end
else    
    % system(sprintf('showimage %s',handles.configPath));    
    % system(sprintf('"%s" %s',fullfile(handles.binPath,'showimage'),handles.configPath));    
    system(sprintf('"%s" %s',showimage,handles.configPath));    
end

% --- Executes on button press in pushbutton_UpdateView.
function handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_UpdateView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.pushbutton_ScanAccept,'Visible','off');
set(handles.pushbutton_ScanReject,'Visible','off');
handles.RbboxHandle = -1;

handles = readAllFields(hObject, eventdata, handles);
if isempty(handles.atomPos)
   return 
end
coords = handles.atomPos;

modelSize = get(handles.radiobutton_SuperCell,'Value')*handles.NcellZ;
if get(handles.radiobutton_Slab,'Value') == 1
   modelSize = ceil(handles.NcellZ / handles.NsubSlabs); 
end
    


if get(handles.radiobutton_Ncells,'Value')  % in NCELLS mode
    center = 0.5*[handles.NcellX handles.NcellY handles.NcellZ];
    box = [0 0 0;handles.NcellX 0 0; handles.NcellX handles.NcellY 0; ...
        0 handles.NcellY 0;0 0 handles.NcellZ;handles.NcellX 0 handles.NcellZ; ...
        handles.NcellX handles.NcellY handles.NcellZ; 0 handles.NcellY handles.NcellZ];
    % box = (handles.Mm.'*(box-repmat(center,size(box,1),1)).').';
    box = (handles.Mm*(box-repmat(center,size(box,1),1)).').';
    rotationCenter = (handles.Mm*(center.')).';
    % rotationCenter = center*handles.Mm;
    box = repmat(rotationCenter,size(box,1),1)+(handles.Mrot*(box.')).';
    % This moves the most left, front and bottom point of the super cell to
    % the origin.
    NcellPotOffset = min(box,[],1); % +[handles.PotentialOffsetX handles.PotentialOffsetY 0];
    box = box-repmat(NcellPotOffset,size(box,1),1);
    if (1)
        fprintf('Super cell size: %.3f x %.3f x %.3f A\n',max(box(:,1)),max(box(:,2)),max(box(:,3)));
    end
    handles.NcellPotOffset = NcellPotOffset;
else  % in box mode
    box = [0 0 0;1 0 0; 1 1 0; 0 1 0;0 0 1;1 0 1; 1 1 1; 0 1 1]*[handles.BoxX 0 0;0 handles.BoxY 0; 0 0 handles.BoxZ];
    % box = box+repmat([handles.PotentialOffsetX handles.PotentialOffsetY 0],size(box,1),1);
    rotationCenter = 0.5*[handles.BoxX handles.BoxY handles.BoxZ]; % +[handles.PotentialOffsetX handles.PotentialOffsetY 0];
    handles.NcellPotOffset = [0 0 0];
end
% This bounding box contains the limits of the rotated structure:
handles.boundingBox = box;

% Adjust the slice Thickness:
handles.SliceThickness = (max(handles.boundingBox(:,3))-min(handles.boundingBox(:,3)))/handles.NsliceTot;
set(handles.edit_SliceThickness,'String',sprintf('%.4f',handles.SliceThickness));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the probe window parameters:
if handles.modeNum == 0  % if in TEM mode
    handles.ProbeWindowX = max(handles.boundingBox(:,1))-min(handles.boundingBox(:,1));
    handles.ProbeWindowY = max(handles.boundingBox(:,2))-min(handles.boundingBox(:,2));
    % fprintf('probe: %d x %d\n',handles.ProbeNx,handles.ProbeNy);
    handles.ProbeResolutionX = handles.ProbeWindowX/handles.ProbeNx;
    handles.ProbeResolutionY = handles.ProbeWindowY/handles.ProbeNy;
    set(handles.edit_ProbeWindowX,'String',handles.ProbeWindowX);
    set(handles.edit_ProbeWindowY,'String',handles.ProbeWindowY);
    set(handles.edit_ProbeResolutionX,'String',handles.ProbeResolutionX);
    set(handles.edit_ProbeResolutionY,'String',handles.ProbeResolutionY);
    set(handles.pushbutton_ImageSim,'Visible','on');
    set(handles.pushbutton_MoreAberrations,'Visible','off');
    % fprintf('in TEM mode ...\n');
else
    set(handles.pushbutton_ImageSim,'Visible','off');
    set(handles.pushbutton_MoreAberrations,'Visible','on');
    % fprintf('in STEM mode ...\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the TOP VIEW:
if  get(handles.radiobutton_viewFront,'Value') == 0
    % show top view:
    % If the user selected Ncell mode:
    if get(handles.radiobutton_Ncells,'Value')
        % replicate in X- and Y-direction, if requested
        if modelSize > 0
            coords0 = coords;
            Natom = size(coords,1);
            for iz = 0:handles.NcellZ-1
                for ix = 0:handles.NcellX-1
                    if (ix == 0) && (iz == 0)
                        for iy = 1:handles.NcellY-1
                            % [aType coords DW occ charge]
                            coords = [coords; coords0+repmat([0 ix iy iz 0 0 0],Natom,1)];
                        end
                    else
                        for iy = 0:handles.NcellY-1
                            coords = [coords; coords0+repmat([0 ix iy iz 0 0 0],Natom,1)];
                        end
                    end
                end
            end
        end
        % coords(:,2:4) = (handles.Mrot*((coords(:,2:4)*handles.Mm).')).';
        % coords(:,2:4) = (handles.Mm.'*(coords(:,2:4).')).';
        coords(:,2:4) = (handles.Mm*(coords(:,2:4).')).';
        % coords(:,2:4) = coords(:,2:4)*handles.Mm;
        Ncoords = size(coords,1);
        % Do the rotation:
        if sum(sum(handles.Mrot == eye(3))) < 9
            % coords(:,2:4) = repmat(rotationCenter,Ncoords,1)+(coords(:,2:4)-repmat(rotationCenter,Ncoords,1))*(handles.Mrot.');
            coords(:,2:4) = repmat(rotationCenter,Ncoords,1)+(handles.Mrot*(coords(:,2:4)-repmat(rotationCenter,Ncoords,1)).').';
            % coords(:,2:4) = coords(:,2:4)-repmat(min(coords(:,2:4),[],1),Ncoords,1);
            % coords(:,2:4) = coords(:,2:4)-repmat(NcellPotOffset,Ncoords,1);
        end        
        
        coords(:,2:4) = coords(:,2:4)-repmat(NcellPotOffset-[handles.PotentialOffsetX,handles.PotentialOffsetY,0],Ncoords,1);
       
    else  % end of NCELLS mode - start BOX mode
       % MmOrigInv = inv(handles.Mm);
       % Mm        = handles.Mm*(handles.Mrot.')  % this is even different than handles.Mm*handles.Mrot       
       % ((handles.Mm.')*(handles.Mrot.')).'  is the same as below
       Mm        = handles.Mrot*handles.Mm;
       
       Mminv     = inv(Mm);
       % find the integer positions of the corners of the big cube:
       Mbox   = [handles.BoxX 0 0;0 handles.BoxY 0; 0 0 handles.BoxZ];
       box    = [0 0 0;1 0 0; 1 1 0; 0 1 0;0 0 1;1 0 1; 1 1 1; 0 1 1]*Mbox;
       % boxred = box*Mminv;
       boxred = (Mminv*(box.')).';
       NxMin = min(floor(boxred(:,1)))-1;         NxMax =  max(ceil(boxred(:,1)));
       NyMin = min(floor(boxred(:,2)))-1;         NyMax =  max(ceil(boxred(:,2)));
       NzMin = min(floor(boxred(:,3)))-1;         NzMax =  max(ceil(boxred(:,3)));
       % fprintf('%d .. %d | %d .. %d | %d .. %d\n',NxMin,NxMax,NyMin,NyMax,NzMin,NzMax);
       coords0 = coords;
       Natom = size(coords,1);
       coords = [];
       for ix = NxMin:NxMax
           for iy = NyMin:NyMax
               for iz = NzMin:NzMax               
                   % newcoords = (coords0+repmat([0 ix iy iz 0 0 0],Natom,1))*Mm;
                   % newcoords = [coords0(:,1) (coords0(:,2:4)+repmat([ix iy iz],Natom,1))*Mm+repmat([handles.PotentialOffsetX,handles.PotentialOffsetY,0],Natom,1)  coords0(:,5:7)];
                   newcoords = [coords0(:,1) (Mm*(coords0(:,2:4).'+repmat([ix iy iz].',1,Natom))).'+repmat([handles.PotentialOffsetX,handles.PotentialOffsetY,0],Natom,1)  coords0(:,5:7)];
                   ind = find((sum(newcoords(:,2:4)<0,2)==0) & (newcoords(:,2)<=handles.BoxX) & (newcoords(:,3)<=handles.BoxY)& (newcoords(:,4)<=handles.BoxZ));
                   if ~isempty(ind)
                       coords = [coords; newcoords(ind,:)];
                   end
               end
           end
       end
       if (0)
           ind = find((sum(coords(:,2:4)<0,2)==0) & (coords(:,2)<=handles.BoxX) & (coords(:,3)<=handles.BoxY)& (coords(:,4)<=handles.BoxZ));
           coords = coords(ind,:);
       end
    end
    hold off
    if (0)
        figure;
        h = scatter(handles.PotentialOffsetX+coords(:,2),handles.PotentialOffsetY+coords(:,3),handles.atomSize,[0 0 0],'filled');
        set(h,'ButtonDownFcn',{@scatter_ButtonDownFcn,handles,hObject});
    else
        % This is the top view:
        if (0)
            h = scatter(coords(:,2),coords(:,3),handles.atomSize,coords(:,1),'filled');
            set(h,'ButtonDownFcn',{@scatter_ButtonDownFcn,handles,hObject});
        else
            Zmin = min(coords(:,1));
            Zmax = max(coords(:,1));
            if Zmin == Zmax
                Zmax = Zmin+1;
            end
            cmap = colormap();
            % sortrows(coords,[1 4]);
            coordsNext = coords;
            while ~isempty(coordsNext)
                Z = coordsNext(1,1);
                ind = find(coordsNext(:,1) == Z);
                col = 1+round((Z-Zmin)/(Zmax-Zmin)*(size(cmap,1)-1));
                h = plot(coordsNext(ind,2),coordsNext(ind,3),'.','MarkerSize',handles.atomSize,'MarkerEdgeColor',cmap(col,:));
                hold on
                set(h,'ButtonDownFcn',{@scatter_ButtonDownFcn,handles,hObject});                
                coordsNext = coordsNext(find(coordsNext(:,1) ~= Z),:);                
            end
            set(handles.axes_model,'YDir','normal');
            hold off
        end
    end
    
    % show the potential window as well as the scan positions:
    if (handles.modeNum == 1) || (handles.modeNum == 2)  % if we are in scan (STEM) mode

        hold on
        if (handles.Xpixels > 0) && (handles.Ypixels > 0)
            dx = (handles.Xstop-handles.Xstart)/handles.Xpixels;
            dy = (handles.Ystop-handles.Ystart)/handles.Ypixels;
        end
        x1 = handles.Xstart-0.5*handles.ProbeWindowX;
        x2 = handles.Xstop+0.5*handles.ProbeWindowX;
        y1 = handles.Ystart-0.5*handles.ProbeWindowY;
        y2 = handles.Ystop+0.5*handles.ProbeWindowY;
        h = patch([x1 x2 x2 x1],[y1 y1 y2 y2],[1 0 0]); alpha(0.2);
        set(h,'ButtonDownFcn',{@scatter_ButtonDownFcn,handles,hObject});
        x1 = handles.Xstart;
        x2 = handles.Xstop;
        y1 = handles.Ystart;
        y2 = handles.Ystop;
        h = patch([x1 x2 x2 x1],[y1 y1 y2 y2],[0 1 0]); alpha(0.2);
        set(h,'ButtonDownFcn',{@scatter_ButtonDownFcn,handles,hObject});
        [X,Y] = meshgrid(handles.Xstart+dx*[0:handles.Xpixels],handles.Ystart+dy*[0:handles.Ypixels]);
        h = scatter(X(1:end),Y(1:end),50,'k+');
        set(h,'ButtonDownFcn',{@scatter_ButtonDownFcn,handles,hObject});
        
        hold off
    else  % if we are in TEM mode:

        x1 = min([handles.boundingBox(:,1); 0]);
        x2 = max([handles.boundingBox(:,1); handles.ProbeWindowX]);
        y1 = min([handles.boundingBox(:,2); 0]);
        y2 = max([handles.boundingBox(:,2); handles.ProbeWindowY]);
        h = patch([0 handles.ProbeWindowX handles.ProbeWindowX 0],[0 0 handles.ProbeWindowY handles.ProbeWindowY],[1 0 0]); alpha(0.2);
        set(h,'ButtonDownFcn',{@scatter_ButtonDownFcn,handles,hObject});
        xlim([x1 x2]);
        ylim([y1 y2]);
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
    % show front view:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
    
    % replicate in X- and Z-direction, if requested
    if get(handles.radiobutton_Ncells,'Value')
        if modelSize > 0
            coords0 = coords;
            Natom = size(coords,1);
        
            % If the user selected Ncell mode:
            for iy=0:handles.NcellY-1
                for ix = 0:handles.NcellX-1
                    if (ix == 0) && (iy == 0)
                        for iz = 1:modelSize-1
                            coords = [coords; coords0+repmat([0 ix iy iz 0 0 0],Natom,1)];
                        end
                    else
                        for iz = 0:modelSize-1
                            coords = [coords; coords0+repmat([0 ix iy iz 0 0 0],Natom,1)];
                        end
                    end
                end
            end
        end
        % coords(:,2:4) = (handles.Mrot*((coords(:,2:4)*handles.Mm).')).';
        coords(:,2:4) = (handles.Mm*(coords(:,2:4).')).';
        Ncoords = size(coords,1);
        % Do the rotation:
        if sum(sum(handles.Mrot == eye(3))) < 9
            coords(:,2:4) = repmat(rotationCenter,Ncoords,1)+(handles.Mrot*(coords(:,2:4)-repmat(rotationCenter,Ncoords,1)).').';
            % coords(:,2:4) = repmat(rotationCenter,Ncoords,1)+(coords(:,2:4)-repmat(rotationCenter,Ncoords,1))*(handles.Mrot.');
            % coords(:,2:4) = coords(:,2:4)-repmat(NcellPotOffset,Ncoords,1);
        end
        coords(:,2:4) = coords(:,2:4)-repmat(NcellPotOffset-[handles.PotentialOffsetX,handles.PotentialOffsetY,0],Ncoords,1);
    else  % Box mode !!!  front view
        % msgbox('Box mode not implemented yet!');
        Mm        = handles.Mrot*handles.Mm;
        
        Mminv     = inv(Mm);
        % find the integer positions of the corners of the big cube:
        Mbox   = [handles.BoxX 0 0;0 handles.BoxY 0; 0 0 handles.BoxZ];
        box    = [0 0 0;1 0 0; 1 1 0; 0 1 0;0 0 1;1 0 1; 1 1 1; 0 1 1]*Mbox;
        % boxred = box*Mminv;
        boxred = (Mminv*(box.')).';
        NxMin = min(floor(boxred(:,1)))-1;         NxMax =  max(ceil(boxred(:,1)));
        NyMin = min(floor(boxred(:,2)))-1;         NyMax =  max(ceil(boxred(:,2)));
        NzMin = min(floor(boxred(:,3)))-1;         NzMax =  max(ceil(boxred(:,3)));

        % fprintf('%d .. %d | %d .. %d | %d .. %d\n',NxMin,NxMax,NyMin,NyMax,NzMin,NzMax);
        coords0 = coords;
        Natom = size(coords,1);
        coords = [];
        for ix = NxMin:NxMax
            for iy = NyMin:NyMax
                for iz = NzMin:NzMax
                    % fprintf('%d %d %d: (%f %f %f) %f %f %f\n',ix,iy,iz,coords0(1,2:4),(Mm*(coords0(1,2:4).'+[ix iy iz].')).');
                    % coords = [coords; coords0+repmat([0 ix iy iz 0 0 0],Natom,1)];
                    % newcoords = [coords0(:,1) (coords0(:,2:4)+repmat([ix iy iz],Natom,1))*Mm+repmat([handles.PotentialOffsetX,handles.PotentialOffsetY,0],Natom,1) coords0(:,5:7)];                      
                    newcoords = [coords0(:,1) (Mm*(coords0(:,2:4).'+repmat([ix iy iz].',1,Natom))).'+repmat([handles.PotentialOffsetX,handles.PotentialOffsetY,0],Natom,1)  coords0(:,5:7)];
                    ind = find((sum(newcoords(:,2:4)<0,2)==0) & (newcoords(:,2)<=handles.BoxX) & (newcoords(:,3)<=handles.BoxY)& (newcoords(:,4)<=handles.BoxZ));
                    if ~isempty(ind)
                        coords = [coords; newcoords(ind,:)];
                    end
                end                
            end
        end
        if (0)
            % coords(:,2:4) = (handles.Mrot*((coords(:,2:4)*handles.Mm).')).';
            %coords(:,2:4) = coords(:,2:4)*Mm;
            ind = find((sum(coords(:,2:4)<0,2)==0) & (coords(:,2)<=handles.BoxX) & (coords(:,3)<=handles.BoxY)& (coords(:,4)<=handles.BoxZ));
            coords = coords(ind,:);
        end
        if (get(handles.radiobutton_Slab,'Value') == 1)
            ind = find(coords(:,4)<=handles.BoxZ/handles.NsubSlabs);
            coords = coords(ind,:);
        elseif (get(handles.radiobutton_UnitCell,'Value') == 1)
            ind = find(coords(:,4)<=handles.Mm(3,3));
            coords = coords(ind,:);            
        end
        
    end  % end of box mode in front view
    handles.CenterSlices = get(handles.checkbox_CenterSlices,'Value') == 1;
    if (handles.modeNum == 1) || (handles.modeNum == 2)  % for STEM mode
        x1 = handles.Xstart-0.5*handles.ProbeWindowX;
        x2 = handles.Xstop+0.5*handles.ProbeWindowX;
        dz = 0; -min(coords(:,4));
        dz = -min(handles.boundingBox(:,3));
        % fprintf('dz= %f, min(coords) = %f\n',dz,min(coords(:,4)));
    else                     % for TEM mode
       x1 = min(handles.boundingBox(:,1));
       x2 = max(handles.boundingBox(:,1));
       sliceThickness = max(handles.boundingBox(:,3))-min(handles.boundingBox(:,3));
       dz = -min(handles.boundingBox(:,3));
    end
    % fprintf('%d atoms\n',size(coords,1))
    hold off    
    if (0)
        scatter(coords(:,2),handles.PotentialOffsetZ+dz+coords(:,4),handles.atomSize,coords(:,1),'filled');
    else
        Zmin = min(coords(:,1));
        Zmax = max(coords(:,1));
        if Zmin == Zmax
            Zmax = Zmin+1;
        end
        cmap = colormap();
        % sortrows(coords,[1 4]);
        coordsNext = coords;
        while ~isempty(coordsNext)
            Z = coordsNext(1,1);
            ind = find(coordsNext(:,1) == Z);
            col = 1+round((Z-Zmin)/(Zmax-Zmin)*(size(cmap,1)-1));
            plot(coordsNext(ind,2),handles.PotentialOffsetZ+dz+coordsNext(ind,4),'.','MarkerSize',handles.atomSize,'MarkerEdgeColor',cmap(col,:));
            hold on
            coordsNext = coordsNext(find(coordsNext(:,1) ~= Z),:);
        end
        hold off        
        set(handles.axes_model,'YDir','reverse');
    end
    hold on
    
    % x1 = min(coords(:,2));
    % x2 = max(coords(:,2));
    if get(handles.radiobutton_Ncells,'Value')  % in NCELLS mode (front view)
        sliceThickness = handles.SliceThickness; % handles.NcellZ*handles.Mm(3,3)/handles.NsliceTot;
        % if only a single sub-slab is shown:
        if (get(handles.radiobutton_UnitCell,'Value') == 1)
            for iz=0:ceil(handles.NsliceTot/handles.NcellZ)-1+handles.CenterSlices
                patch([x1 x2 x2 x1],sliceThickness*([iz iz iz+1 iz+1]-0.5*handles.CenterSlices),[(mod(iz,3)==0),(mod(iz,3)==1),(mod(iz,3)==2)]); alpha(0.2);
            end
        elseif (get(handles.radiobutton_Slab,'Value') == 1)
            for iz=0:handles.Nslice-1+handles.CenterSlices
                patch([x1 x2 x2 x1],sliceThickness*([iz iz iz+1 iz+1]-0.5*handles.CenterSlices),[(mod(iz,3)==0),(mod(iz,3)==1),(mod(iz,3)==2)]); alpha(0.2);
            end
        elseif (get(handles.radiobutton_SuperCell,'Value') == 1)
            for iz=0:handles.NsubSlabs-1
                patch([x1 x2 x2 x1],handles.Nslice*sliceThickness*[iz iz iz+1 iz+1],[(mod(iz,3)==0),(mod(iz,3)==1),(mod(iz,3)==2)]); alpha(0.2);
            end
        end
    else  % in box mode (front view)
        sliceThickness = handles.BoxZ/handles.NsliceTot;
        if (get(handles.radiobutton_UnitCell,'Value') == 1)
            for iz=0:ceil(handles.NsliceTot/(handles.BoxZ/handles.Mm(3,3)))-1+handles.CenterSlices
                patch([x1 x2 x2 x1],sliceThickness*([iz iz iz+1 iz+1]-0.5*handles.CenterSlices),[(mod(iz,3)==0),(mod(iz,3)==1),(mod(iz,3)==2)]); alpha(0.2);
            end
        elseif (get(handles.radiobutton_Slab,'Value') == 1)
            for iz=0:ceil(handles.BoxZ/(handles.NsubSlabs*sliceThickness))-1+handles.CenterSlices
                patch([x1 x2 x2 x1],sliceThickness*([iz iz iz+1 iz+1]-0.5*handles.CenterSlices),[(mod(iz,3)==0),(mod(iz,3)==1),(mod(iz,3)==2)]); alpha(0.2);
            end
        elseif (get(handles.radiobutton_SuperCell,'Value') == 1)
            for iz=0:handles.NsubSlabs-1
                patch([x1 x2 x2 x1],handles.Nslice*sliceThickness*[iz iz iz+1 iz+1],[(mod(iz,3)==0),(mod(iz,3)==1),(mod(iz,3)==2)]); alpha(0.2);
            end
        end
        
    end  % end of plotting the patch (front view)

    
    hold off
end
axis equal; axis tight

handles.coords = coords;
clear coords
pause(0.2)
if handles.modeNum == 0
    set(handles.pushbutton_ImageSim,'Visible','on');
    set(handles.pushbutton_MoreAberrations,'Visible','off');
    % fprintf('in TEM mode ...\n');
else
    set(handles.pushbutton_ImageSim,'Visible','off');
    set(handles.pushbutton_MoreAberrations,'Visible','on');
    % fprintf('in STEM mode ...\n');
end

pause(0.2)
guidata(hObject, handles);


% --- Executes on button press in pushbuttonLoadModel.
function handles = pushbuttonLoadModel_Callback(hObject, eventdata, handles,askUser,showStruct)
% hObject    handle to pushbuttonLoadModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if nargin < 4
    showStruct = 1;
end
if (nargin < 4)
    askUser = 1;
end
[pathname, filename, ext] = fileparts(handles.posFileName);
filename = [filename ext];
if (askUser)
    if isfield(handles,'posFileName')
        [filename, pathname] = uigetfile({'*.cfg', 'CFG file (*.cfg)';...
            '*.*', 'All files (*.*)'},...
            'Select a file',handles.posFileName);
    else
        [filename, pathname] = uigetfile({'*.cfg', 'CFG file (*.cfg)';...
            '*.*', 'All files (*.*)'},...
            'Select a file');
    end
    if isequal(filename,0) || isequal(pathname,0)
        % disp('User pressed cancel')
        return
    end
end

if (0)
    oldPath = pwd()
    cd(pathname);
    
    [coords, aType, Mm, DW, occ, charge] = readCFG_qstem(filename,[0 0 0],2);
    cd(oldPath);
else
    [coords, aType, Mm, DW, occ, charge] = readCFG_qstem(fullfile(pathname,filename),[0 0 0],2);
end

if ~isempty(coords)
    handles.posFileName = fullfile(pathname,filename);
    set(handles.uipanel_model,'Title',sprintf('Model: %s',filename));
    handles.atomPos = [aType coords  DW occ charge];
    % [aType coords  DW occ charge]
    clear coords aType DW occ charge
    handles.Mm = Mm;
    guidata(hObject, handles);

    if showStruct
        handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
        guidata(hObject, handles);
    end
end


% --- Executes on button press in pushbutton_show3Dmodel.
function pushbutton_show3Dmodel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_show3Dmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% coords = handles.atomPos(:,2:4)*(handles.Mm);
% drawXtal([handles.atomPos(:,1) coords],handles.Mm);
% clear coords
if isfield(handles,'coords')
    if isfield(handles,'NcellPotOffset')
        drawXtal([handles.coords(:,1) handles.coords(:,2:4)],handles.Mm,-handles.NcellPotOffset);
    else
        drawXtal([handles.coords(:,1) handles.coords(:,2:4)],handles.Mm);
    end
else
    drawCFG(handles.posFileName,1);
end

function edit_Xstart_Callback(hObject, eventdata, handles)
if handles.modeNum == 2 % CBED
    handles.Xstart = str2double(get(handles.edit_Xstart,'String'));
    handles.Xstop = handles.Xstart;
    handles.Ystart = str2double(get(handles.edit_Ystart,'String'));
    handles.Ystop = handles.Ystart;
    set(handles.edit_Xstop,'String',handles.Xstop);
    set(handles.edit_Ystop,'String',handles.Ystop);

    handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_Xstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Xstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Xstop_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Xstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Xstop as text
%        str2double(get(hObject,'String')) returns contents of edit_Xstop as a double


% --- Executes during object creation, after setting all properties.
function edit_Xstop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ystart_Callback(hObject, eventdata, handles)
if handles.modeNum == 2 % CBED
    handles.Xstart = str2double(get(handles.edit_Xstart,'String'));
    handles.Xstop = handles.Xstart;
    handles.Ystart = str2double(get(handles.edit_Ystart,'String'));
    handles.Ystop = handles.Ystart;
    set(handles.edit_Xstop,'String',handles.Xstop);
    set(handles.edit_Ystop,'String',handles.Ystop);

    handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function edit_Ystart_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ystop_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_Ystop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ystop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Xpixels_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Xpixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Xpixels as text
%        str2double(get(hObject,'String')) returns contents of edit_Xpixels as a double


% --- Executes during object creation, after setting all properties.
function edit_Xpixels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Xpixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ypixels_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Ypixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Ypixels as text
%        str2double(get(hObject,'String')) returns contents of edit_Ypixels as a double


% --- Executes during object creation, after setting all properties.
function edit_Ypixels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ypixels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NcellX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NcellX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NcellX as text
%        str2double(get(hObject,'String')) returns contents of edit_NcellX as a double


% --- Executes during object creation, after setting all properties.
function edit_NcellX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NcellX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NcellY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NcellY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NcellY as text
%        str2double(get(hObject,'String')) returns contents of edit_NcellY as a double


% --- Executes during object creation, after setting all properties.
function edit_NcellY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NcellY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tiltX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tiltX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tiltX as text
%        str2double(get(hObject,'String')) returns contents of edit_tiltX as a double


% --- Executes during object creation, after setting all properties.
function edit_tiltX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tiltX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tiltY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tiltY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tiltY as text
%        str2double(get(hObject,'String')) returns contents of edit_tiltY as a double


% --- Executes during object creation, after setting all properties.
function edit_tiltY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tiltY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tiltZ_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tiltZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tiltZ as text
%        str2double(get(hObject,'String')) returns contents of edit_tiltZ as a double


% --- Executes during object creation, after setting all properties.
function edit_tiltZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tiltZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ProbeNx_Callback(hObject, eventdata, handles)
if handles.modeNum > 0  % if STEM mode
    handles.ProbeNx = str2num(get(handles.edit_ProbeNx,'String'));
    handles.ProbeResolutionX = str2num(get(handles.edit_ProbeResolutionX,'String'));
    if handles.ProbeNx < 8
        handles.ProbeNx = 8;
        set(handles.edit_ProbeNx,'String',handles.ProbeNx);
    end    
    %set X-parameters:
    handles.ProbeWindowX = handles.ProbeResolutionX*handles.ProbeNx;
    handles.ProbeMaxAngleX = 2/6*handles.Wavelength * (handles.ProbeNx/handles.ProbeWindowX);
else  % if in TEM mode
    handles.ProbeNx = str2num(get(handles.edit_ProbeNx,'String'));
    handles.ProbeWindowX = str2num(get(handles.edit_ProbeWindowX,'String'));
    handles.ProbeResolutionX = handles.ProbeWindowX/handles.ProbeNx;    
    handles.ProbeMaxAngleX = 2/6*handles.Wavelength * (handles.ProbeNx/handles.ProbeWindowX);
end    
set(handles.edit_ProbeWindowX,'String',handles.ProbeWindowX);
set(handles.edit_ProbeResolutionX,'String',handles.ProbeResolutionX);
set(handles.edit_ProbeMaxAngleX,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleX));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_ProbeNx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ProbeNy_Callback(hObject, eventdata, handles)
if handles.modeNum > 0
    handles.ProbeNy = str2num(get(handles.edit_ProbeNy,'String'));
    handles.ProbeResolutionY = str2num(get(handles.edit_ProbeResolutionY,'String'));
    if handles.ProbeNy < 8
        handles.ProbeNy = 8;
        set(handles.edit_ProbeNy,'String',handles.ProbeNy);
    end
    
    % set Y-parameters:
    handles.ProbeWindowY = handles.ProbeResolutionY*handles.ProbeNy;
    
    handles.ProbeMaxAngleY = 2/6*handles.Wavelength * (handles.ProbeNy/handles.ProbeWindowY);
else  % if in TEM mode
    handles.ProbeNy = str2num(get(handles.edit_ProbeNy,'String'));
    handles.ProbeWindowY = str2num(get(handles.edit_ProbeWindowY,'String'));
    handles.ProbeResolutionY = handles.ProbeWindowY/handles.ProbeNy;    
    handles.ProbeMaxAngleY = 2/6*handles.Wavelength * (handles.ProbeNy/handles.ProbeWindowY);
end
set(handles.edit_ProbeWindowY,'String',handles.ProbeWindowY);
set(handles.edit_ProbeResolutionY,'String',handles.ProbeResolutionY);
set(handles.edit_ProbeMaxAngleY,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleY));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_ProbeNy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ProbeNy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ProbeResolutionX_Callback(hObject, eventdata, handles)
handles.ProbeNx = str2num(get(handles.edit_ProbeNx,'String'));
handles.ProbeResolutionX = str2num(get(handles.edit_ProbeResolutionX,'String'));
if handles.ProbeResolutionX < 1e-4
    handles.ProbeResolutionX = 1e-4;
    set(handles.edit_ProbeResolutionX,'String',handles.ProbeResolutionX);
end

%set X-parameters:
handles.ProbeWindowX = handles.ProbeResolutionX*handles.ProbeNx;
handles.ProbeMaxAngleX = 2/6*handles.Wavelength * (handles.ProbeNx/handles.ProbeWindowX);
set(handles.edit_ProbeWindowX,'String',handles.ProbeWindowX);
set(handles.edit_ProbeMaxAngleX,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleX));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_ProbeResolutionX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ProbeResolutionX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ProbeResolutionY_Callback(hObject, eventdata, handles)
handles.ProbeNy = str2num(get(handles.edit_ProbeNy,'String'));
handles.ProbeResolutionY = str2num(get(handles.edit_ProbeResolutionY,'String'));
if handles.ProbeResolutionY < 1e-4
    handles.ProbeResolutionY = 1e-4;
    set(handles.edit_ProbeResolutionY,'String',handles.ProbeResolutionY);
end

% set Y-parameters:
handles.ProbeWindowY = handles.ProbeResolutionY*handles.ProbeNy;
handles.ProbeMaxAngleY = 2/6*handles.Wavelength * (handles.ProbeNy/handles.ProbeWindowY);
set(handles.edit_ProbeWindowY,'String',handles.ProbeWindowX);
set(handles.edit_ProbeMaxAngleY,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleY));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_ProbeResolutionY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ProbeWindowX_Callback(hObject, eventdata, handles)
handles.ProbeNx = str2num(get(handles.edit_ProbeNx,'String'));
handles.ProbeWindowX = str2num(get(handles.edit_ProbeWindowX,'String'));
if handles.ProbeWindowX < 1
    handles.ProbeWindowX = 1;
    set(handles.edit_ProbeWindowX,'String',handles.ProbeWindowX);
end

%set X-parameters:
handles.ProbeResolutionX = handles.ProbeWindowX/handles.ProbeNx;
handles.ProbeMaxAngleX = 2/6*handles.Wavelength * (handles.ProbeNx/handles.ProbeWindowX);
set(handles.edit_ProbeResolutionX,'String',handles.ProbeResolutionX);
set(handles.edit_ProbeMaxAngleX,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleX));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_ProbeWindowX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ProbeWindowY_Callback(hObject, eventdata, handles)
handles.ProbeNy = str2num(get(handles.edit_ProbeNy,'String'));
handles.ProbeWindowY = str2num(get(handles.edit_ProbeWindowY,'String'));
if handles.ProbeWindowY < 1
    handles.ProbeWindowY = 1;
    set(handles.edit_ProbeWindowY,'String',handles.ProbeWindowY);
end

%set X-parameters:
handles.ProbeResolutionY = handles.ProbeWindowY/handles.ProbeNy;
handles.ProbeMaxAngleY = 2/6*handles.Wavelength * (handles.ProbeNy/handles.ProbeWindowY);
set(handles.edit_ProbeResolutionY,'String',handles.ProbeResolutionY);
set(handles.edit_ProbeMaxAngleY,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleY));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_ProbeWindowY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ProbeMaxAngleX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ProbeMaxAngleX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ProbeMaxAngleX as text
%        str2double(get(hObject,'String')) returns contents of edit_ProbeMaxAngleX as a double


% --- Executes during object creation, after setting all properties.
function edit_ProbeMaxAngleX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ProbeMaxAngleX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ProbeMaxAngleY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ProbeMaxAngleY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ProbeMaxAngleY as text
%        str2double(get(hObject,'String')) returns contents of edit_ProbeMaxAngleY as a double


% --- Executes during object creation, after setting all properties.
function edit_ProbeMaxAngleY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ProbeMaxAngleY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NsubSlabs_Callback(hObject, eventdata, handles)

handles.NcellZ         = floor(str2num(get(handles.edit_NcellZ,'String')));
handles.NsubSlabs      = floor(str2num(get(handles.edit_NsubSlabs,'String')));
handles.Nslice         = floor(str2num(get(handles.edit_Nslice,'String')));

handles.NsliceTot      = handles.NsubSlabs*handles.Nslice;
if get(handles.radiobutton_Box,'Value')
    handles.SliceThickness = handles.BoxZ/handles.NsliceTot;
else
    handles.SliceThickness = handles.Mm(3,3)*handles.NcellZ/handles.NsliceTot;
end

set(handles.text_NsliceTot,'String',sprintf('(Total: %d)',handles.NsliceTot));
set(handles.edit_SliceThickness,'String',sprintf('%.4f',handles.SliceThickness));

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_NsubSlabs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NsubSlabs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Nslice_Callback(hObject, eventdata, handles)
handles.NcellZ         = floor(str2num(get(handles.edit_NcellZ,'String')));
handles.NsubSlabs      = floor(str2num(get(handles.edit_NsubSlabs,'String')));
handles.Nslice         = floor(str2num(get(handles.edit_Nslice,'String')));

handles.NsliceTot      = handles.NsubSlabs*handles.Nslice;
if get(handles.radiobutton_Box,'Value')
    handles.SliceThickness = handles.BoxZ/handles.NsliceTot;
else
    if isfield(handles,'boundingBox')
        handles.SliceThickness = (max(handles.boundingBox(:,3))-min(handles.boundingBox(:,3)))/handles.NsliceTot;
    else
        handles.SliceThickness = handles.Mm(3,3)*handles.NcellZ/handles.NsliceTot;
    end
end

set(handles.text_NsliceTot,'String',sprintf('(Total: %d)',handles.NsliceTot));
set(handles.edit_SliceThickness,'String',sprintf('%.4f',handles.SliceThickness));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_Nslice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Nslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_SliceThickness_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SliceThickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SliceThickness as text
%        str2double(get(hObject,'String')) returns contents of edit_SliceThickness as a double


% --- Executes during object creation, after setting all properties.
function edit_SliceThickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SliceThickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit_C3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C3 as text
%        str2double(get(hObject,'String')) returns contents of edit_C3 as a double


% --- Executes during object creation, after setting all properties.
function edit_C3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_HighVoltage_Callback(hObject, eventdata, handles)

handles.highVoltage = str2double(get(hObject,'String'));

% update the wavelength:
handles.Wavelength = wavelength(handles.HighVoltage);
set(handles.text_wavelength,'String',sprintf('kV    (wavelength = %.2fpm)',100*handles.Wavelength));

% update the scattering angles:
handles.ProbeMaxAngleX = 2/6*handles.Wavelength * (handles.ProbeNx/handles.ProbeWindowX);
handles.ProbeMaxAngleY = 2/6*handles.Wavelength * (handles.ProbeNy/handles.ProbeWindowY);
set(handles.edit_ProbeMaxAngleX,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleX));
set(handles.edit_ProbeMaxAngleY,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleY));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_HighVoltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_HighVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Defocus_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Defocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Defocus as text
%        str2double(get(hObject,'String')) returns contents of edit_Defocus as a double


% --- Executes during object creation, after setting all properties.
function edit_Defocus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Defocus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AstigMag_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AstigMag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AstigMag as text
%        str2double(get(hObject,'String')) returns contents of edit_AstigMag as a double


% --- Executes during object creation, after setting all properties.
function edit_AstigMag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AstigMag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_AstigAngle_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AstigAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AstigAngle as text
%        str2double(get(hObject,'String')) returns contents of edit_AstigAngle as a double


% --- Executes during object creation, after setting all properties.
function edit_AstigAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AstigAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Scherzer.
function pushbutton_Scherzer_Callback(hObject, eventdata, handles)
handles.C3 = str2num(get(handles.edit_C3,'String'));
handles.HighVoltage = str2num(get(handles.edit_HighVoltage,'String'));
handles.Defocus = -0.1*sign(handles.C3)*sqrt(1.5*abs(handles.C3)*1e7*wavelength(handles.HighVoltage));

if abs(handles.Defocus) > 1
    set(handles.edit_Defocus,'String',sprintf('%.1f',handles.Defocus));    
else
    set(handles.edit_Defocus,'String',sprintf('%.4f',handles.Defocus));
end
guidata(hObject, handles);


% --- Executes on button press in checkbox_CenterSlices.
function checkbox_CenterSlices_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CenterSlices (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CenterSlices


% --- Executes on button press in checkbox_PeriodicXY.
function checkbox_PeriodicXY_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PeriodicXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PeriodicXY


% --- Executes on button press in checkbox_PeriodicZ.
function checkbox_PeriodicZ_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PeriodicZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_PeriodicZ



function edit_PotentialOffsetX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PotentialOffsetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PotentialOffsetX as text
%        str2double(get(hObject,'String')) returns contents of edit_PotentialOffsetX as a double


% --- Executes during object creation, after setting all properties.
function edit_PotentialOffsetX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PotentialOffsetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PotentialOffsetY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PotentialOffsetY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PotentialOffsetY as text
%        str2double(get(hObject,'String')) returns contents of edit_PotentialOffsetY as a double


% --- Executes during object creation, after setting all properties.
function edit_PotentialOffsetY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PotentialOffsetY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PotentialOffsetZ_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PotentialOffsetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PotentialOffsetZ as text
%        str2double(get(hObject,'String')) returns contents of edit_PotentialOffsetZ as a double


% --- Executes during object creation, after setting all properties.
function edit_PotentialOffsetZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PotentialOffsetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








function edit_atomSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_atomSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_atomSize as text
%        str2double(get(hObject,'String')) returns contents of edit_atomSize as a double


% --- Executes during object creation, after setting all properties.
function edit_atomSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_atomSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton_IncreaseDetector.
function pushbutton_IncreaseDetector_Callback(hObject, eventdata, handles)

if isempty(handles.Detectors)
    set(handles.edit_DetectorNumber,'String',1);
    set(handles.edit_InnerAngle,'String',' ');
    set(handles.edit_OuterAngle,'String',' ');
    set(handles.edit_OffsetX,'String',' ');
    set(handles.edit_OffsetY,'String',' ');    
else

    handles.DetectorNumber = floor(str2num(get(handles.edit_DetectorNumber,'String')));
    if handles.DetectorNumber < 1
        handles.DetectorNumber = 1;
    end
    if handles.DetectorNumber > size(handles.Detectors,1)
        handles.DetectorNumber = size(handles.Detectors,1);
    end
    if handles.DetectorNumber < size(handles.Detectors,1)
        handles.DetectorNumber = handles.DetectorNumber+1;
    end
    
    set(handles.edit_DetectorNumber,'String',handles.DetectorNumber);
    j = handles.DetectorNumber;
    set(handles.edit_InnerAngle,'String',handles.Detectors(j,1));
    set(handles.edit_OuterAngle,'String',handles.Detectors(j,2));
    set(handles.edit_OffsetX,'String',handles.Detectors(j,3));
    set(handles.edit_OffsetY,'String',handles.Detectors(j,4));
end    


% --- Executes on button press in pushbutton_DecreaseDetector.
function pushbutton_DecreaseDetector_Callback(hObject, eventdata, handles)
if isempty(handles.Detectors)
    set(handles.edit_DetectorNumber,'String',1);
    set(handles.edit_InnerAngle,'String',' ');
    set(handles.edit_OuterAngle,'String',' ');
    set(handles.edit_OffsetX,'String',' ');
    set(handles.edit_OffsetY,'String',' ');    
else

    handles.DetectorNumber = floor(str2num(get(handles.edit_DetectorNumber,'String')));
    if handles.DetectorNumber < 1
        handles.DetectorNumber = 1;
    end
    if handles.DetectorNumber > size(handles.Detectors,1)
        handles.DetectorNumber = size(handles.Detectors,1);
    end
    if handles.DetectorNumber > 1
        handles.DetectorNumber = handles.DetectorNumber-1;
    end
    
    set(handles.edit_DetectorNumber,'String',handles.DetectorNumber);
    j = handles.DetectorNumber;
    set(handles.edit_InnerAngle,'String',handles.Detectors(j,1));
    set(handles.edit_OuterAngle,'String',handles.Detectors(j,2));
    set(handles.edit_OffsetX,'String',handles.Detectors(j,3));
    set(handles.edit_OffsetY,'String',handles.Detectors(j,4));
end    





function edit_DetectorNumber_Callback(hObject, eventdata, handles)
if isempty(handles.Detectors)
    set(handles.edit_DetectorNumber,'String',1);
    set(handles.edit_InnerAngle,'String',' ');
    set(handles.edit_OuterAngle,'String',' ');
    set(handles.edit_OffsetX,'String',' ');
    set(handles.edit_OffsetY,'String',' ');    
else

    handles.DetectorNumber = floor(str2num(get(handles.edit_DetectorNumber,'String')));
    if handles.DetectorNumber < 1
        handles.DetectorNumber = 1;
    end
    if handles.DetectorNumber > size(handles.Detectors,1)
        handles.DetectorNumber = size(handles.Detectors,1);
    end
    
    set(handles.edit_DetectorNumber,'String',handles.DetectorNumber);
    j = handles.DetectorNumber;
    set(handles.edit_InnerAngle,'String',handles.Detectors(j,1));
    set(handles.edit_OuterAngle,'String',handles.Detectors(j,2));
    set(handles.edit_OffsetX,'String',handles.Detectors(j,3));
    set(handles.edit_OffsetY,'String',handles.Detectors(j,4));
end    



% --- Executes during object creation, after setting all properties.
function edit_DetectorNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DetectorNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_InnerAngle_Callback(hObject, eventdata, handles)
handles.DetectorNumber = floor(str2num(get(handles.edit_DetectorNumber,'String')));
handles.Detectors(handles.DetectorNumber,1) = abs(str2num(get(handles.edit_InnerAngle,'String')));
set(handles.edit_InnerAngle,'String',handles.Detectors(handles.DetectorNumber,1));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_InnerAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_InnerAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_OuterAngle_Callback(hObject, eventdata, handles)
handles.DetectorNumber = floor(str2num(get(handles.edit_DetectorNumber,'String')));
handles.Detectors(handles.DetectorNumber,2) = abs(str2num(get(handles.edit_OuterAngle,'String')));
set(handles.edit_OuterAngle,'String',handles.Detectors(handles.DetectorNumber,2));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_OuterAngle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OuterAngle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_OffsetX_Callback(hObject, eventdata, handles)
handles.DetectorNumber = floor(str2num(get(handles.edit_DetectorNumber,'String')));
handles.Detectors(handles.DetectorNumber,3) = (str2num(get(handles.edit_OffsetX,'String')));
set(handles.edit_OffsetX,'String',handles.Detectors(handles.DetectorNumber,3));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_OffsetX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OffsetX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_OffsetY_Callback(hObject, eventdata, handles)
handles.DetectorNumber = floor(str2num(get(handles.edit_DetectorNumber,'String')));
handles.Detectors(handles.DetectorNumber,4) = (str2num(get(handles.edit_OffsetY,'String')));
set(handles.edit_OffsetY,'String',handles.Detectors(handles.DetectorNumber,4));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_OffsetY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OffsetY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_AddDetector.
function pushbutton_AddDetector_Callback(hObject, eventdata, handles)
handles.DetectorNumber = size(handles.Detectors,1)+1;
set(handles.edit_DetectorNumber,'String',handles.DetectorNumber);
j = handles.DetectorNumber;
handles.Detectors(j,1:4) = 0;
set(handles.edit_InnerAngle,'String',handles.Detectors(j,1));
set(handles.edit_OuterAngle,'String',handles.Detectors(j,2));
set(handles.edit_OffsetX,'String',handles.Detectors(j,3));
set(handles.edit_OffsetY,'String',handles.Detectors(j,4));
guidata(hObject, handles);


% --- Executes on button press in pushbutton_DeleteDetector.
function pushbutton_DeleteDetector_Callback(hObject, eventdata, handles)
handles.DetectorNumber = floor(str2num(get(handles.edit_DetectorNumber,'String')));
if ~isempty(handles.Detectors)
    j = handles.DetectorNumber;
    handles.Detectors = [handles.Detectors(1:j-1,:);handles.Detectors(j+1:end,:)];
    if handles.DetectorNumber > size(handles.Detectors,1)
        handles.DetectorNumber = size(handles.Detectors,1);
    end
    set(handles.edit_DetectorNumber,'String',handles.DetectorNumber);
    j = handles.DetectorNumber;
    set(handles.edit_InnerAngle,'String',handles.Detectors(j,1));
    set(handles.edit_OuterAngle,'String',handles.Detectors(j,2));
    set(handles.edit_OffsetX,'String',handles.Detectors(j,3));
    set(handles.edit_OffsetY,'String',handles.Detectors(j,4));
end
guidata(hObject, handles);




function edit_SlicesBetweenOutputs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SlicesBetweenOutputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SlicesBetweenOutputs as text
%        str2double(get(hObject,'String')) returns contents of edit_SlicesBetweenOutputs as a double


% --- Executes during object creation, after setting all properties.
function edit_SlicesBetweenOutputs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SlicesBetweenOutputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_Temperature_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Temperature as text
%        str2double(get(hObject,'String')) returns contents of edit_Temperature as a double


% --- Executes during object creation, after setting all properties.
function edit_Temperature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Temperature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_TDSruns_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TDSruns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TDSruns as text
%        str2double(get(hObject,'String')) returns contents of edit_TDSruns as a double


% --- Executes during object creation, after setting all properties.
function edit_TDSruns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TDSruns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_TDS.
function checkbox_TDS_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_TDS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_TDS





function edit_Cc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Cc as text
%        str2double(get(hObject,'String')) returns contents of edit_Cc as a double


% --- Executes during object creation, after setting all properties.
function edit_Cc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Cc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha as a double


% --- Executes during object creation, after setting all properties.
function edit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_OutputFolder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OutputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OutputFolder as text
%        str2double(get(hObject,'String')) returns contents of edit_OutputFolder as a double


% --- Executes during object creation, after setting all properties.
function edit_OutputFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OutputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_SelectFolder.
function pushbutton_SelectFolder_Callback(hObject, eventdata, handles)
folder = get(handles.edit_OutputFolder,'String');

if isstr(folder)
    if isempty(findstr(folder,':\'))
        folder = fullfile(handles.configPath,folder);
    end
else
    folder = handles.configPath;
end
s = dir(folder);
if size(s,1) == 0
    folder = handles.configPath;
end
folder = uigetdir(folder, 'Select target folder for simulation results');
if folder ~= 0
    % First let's try to remove the config folder, if we can:
    if ~isempty(strmatch(handles.configPath,folder))
       folder = folder(length(handles.configPath)+1:end); 
       if folder(1) == '\', folder = folder(2:end); end;
    end
    handles.OutputFolder = folder;
    set(handles.edit_OutputFolder,'String',handles.OutputFolder);
    guidata(hObject,handles);
end

% --- Executes on button press in radiobutton_SuperCell.
function radiobutton_SuperCell_Callback(hObject, eventdata, handles)
handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on button press in radiobutton_Slab.
function radiobutton_Slab_Callback(hObject, eventdata, handles)
handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in radiobutton_UnitCell.
function radiobutton_UnitCell_Callback(hObject, eventdata, handles)
handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in radiobutton_viewTop.
function radiobutton_viewTop_Callback(hObject, eventdata, handles)
handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in radiobutton_viewFront.
function radiobutton_viewFront_Callback(hObject, eventdata, handles)
handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function reads all relevant fields and should be executed before
% procesing the user-entered data, especially before displaying or saving
% the config file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = readAllFields(hObject, eventdata, handles)

% Microscope parameters:
handles.HighVoltage = str2num(get(handles.edit_HighVoltage,'String'));
handles.Wavelength = wavelength(handles.HighVoltage);
set(handles.text_wavelength,'String',sprintf('kV    (wavelength = %.2fpm)',100*handles.Wavelength));

handles.Defocus     = str2num(get(handles.edit_Defocus,'String'));
handles.AstigMag    = str2num(get(handles.edit_AstigMag,'String'));
handles.AstigAngle  = str2num(get(handles.edit_AstigAngle,'String'));
handles.alpha       = str2num(get(handles.edit_alpha,'String'));
handles.C3          = str2num(get(handles.edit_C3,'String'));
% handles.C5          = str2num(get(handles.edit_C5,'String'));

handles.Cc          = str2num(get(handles.edit_Cc,'String'));
handles.dE          = str2num(get(handles.edit_dE,'String'));
handles.Delta       = handles.Cc*1e3*handles.dE/handles.HighVoltage;
handles.TDS         = get(handles.checkbox_TDS,'Value');
handles.Temperature = str2num(get(handles.edit_Temperature,'String'));
handles.TDSruns     = floor(str2num(get(handles.edit_TDSruns,'String')));
handles.BeamTiltX   = str2num(get(handles.edit_BtiltX,'String'));
handles.BeamTiltY   = str2num(get(handles.edit_BtiltY,'String'));
handles.TiltBack    = get(handles.checkbox_TiltBack,'Value');

% handles.Delta = handles.Cc *
% str2num(extractParam('dV/V',nameArray,pValues));

% supercell and slicing data:
handles.atomSize = floor(str2num(get(handles.edit_atomSize,'String')));

if get(handles.radiobutton_Ncells,'Value');
    handles.NcellX = floor(str2num(get(handles.edit_NcellX,'String')));
    handles.NcellY = floor(str2num(get(handles.edit_NcellY,'String')));
    handles.NcellZ = floor(str2num(get(handles.edit_NcellZ,'String')));
else
    handles.BoxX = str2num(get(handles.edit_NcellX,'String'));
    handles.BoxY = str2num(get(handles.edit_NcellY,'String'));
    handles.BoxZ = str2num(get(handles.edit_NcellZ,'String'));    
end

% Read the specimen tilt and compute the rotation matrix:
handles.degree = get(handles.radiobutton_degree,'Value');

handles.TiltX  = str2num(get(handles.edit_TiltX,'String'));
handles.TiltY  = str2num(get(handles.edit_TiltY,'String'));
handles.TiltZ  = str2num(get(handles.edit_TiltZ,'String'));
if handles.degree
   handles.TiltX = handles.TiltX*pi/180; 
   handles.TiltY = handles.TiltY*pi/180; 
   handles.TiltZ = handles.TiltZ*pi/180; 
end
cx = cos(handles.TiltX);
sx = sin(handles.TiltX);
cy = cos(handles.TiltY);
sy = sin(handles.TiltY);
cz = cos(handles.TiltZ);
sz = sin(handles.TiltZ);
if (1)
    % Mx = [1 0 0;0 cx sx;0 -sx cx];
    % My = [cy 0 -sy;0 1 0; sy 0 cy];
    % Mz = [cz sz 0;-sz cz 0; 0 0 1];
    Mx = [1 0 0;0 cx -sx;0 sx cx];
    My = [cy 0 sy;0 1 0; -sy 0 cy];
    Mz = [cz -sz 0;sz cz 0; 0 0 1];
    handles.Mrot = Mz*My*Mx;
else
    handles.Mrot(1,1) = cz*cy;
    handles.Mrot(1,2) = cz*sy*sx-sz*cx;
    handles.Mrot(1,3) = cz*sy*cx+sz*sx;

    handles.Mrot(2,1) = sz*cy;
    handles.Mrot(2,2) = sz*sy*sx+cz*cx;
    handles.Mrot(2,3) = sz*sy*cx-cz*sx;

    handles.Mrot(3,1) = -sy;
    handles.Mrot(3,2) = cy*sx;
    handles.Mrot(3,3) = cy*cx;
end

% Read the slicing information:
handles.NsubSlabs = floor(str2num(get(handles.edit_NsubSlabs,'String')));
handles.Nslice = floor(str2num(get(handles.edit_Nslice,'String')));
handles.NsliceTot = handles.NsubSlabs*handles.Nslice;
handles.PotentialOffsetZ = str2num(get(handles.edit_PotentialOffsetZ,'String'));
handles.CenterSlices = get(handles.checkbox_CenterSlices,'Value');

% scan parameters:
handles.Xstart = str2num(get(handles.edit_Xstart,'String'));
handles.Xstop  = str2num(get(handles.edit_Xstop,'String'));
handles.Ystart = str2num(get(handles.edit_Ystart,'String'));
handles.Ystop  = str2num(get(handles.edit_Ystop,'String'));
handles.Xpixels = floor(str2num(get(handles.edit_Xpixels,'String')));
handles.Ypixels = floor(str2num(get(handles.edit_Ypixels,'String')));


% Probe window parameters:
handles.ProbeNx = floor(str2num(get(handles.edit_ProbeNx,'String')));
handles.ProbeNy = floor(str2num(get(handles.edit_ProbeNy,'String')));
handles.ProbeResolutionX = str2num(get(handles.edit_ProbeResolutionX,'String'));
handles.ProbeResolutionY = str2num(get(handles.edit_ProbeResolutionY,'String'));
handles.ProbeWindowX = str2num(get(handles.edit_ProbeWindowX,'String'));
handles.ProbeWindowY = str2num(get(handles.edit_ProbeWindowY,'String'));
handles.PotentialOffsetX = str2num(get(handles.edit_PotentialOffsetX,'String'));
handles.PotentialOffsetY = str2num(get(handles.edit_PotentialOffsetY,'String'));
handles.PotentialOffsetZ = str2num(get(handles.edit_PotentialOffsetZ,'String'));
handles.PeriodicZ  = get(handles.checkbox_PeriodicZ,'Value') == 1;
handles.PeriodicXY  = get(handles.checkbox_PeriodicXY,'Value') == 1;


% update the scattering angles:
handles.ProbeMaxAngleX = 2/6*handles.Wavelength * (handles.ProbeNx/handles.ProbeWindowX);
handles.ProbeMaxAngleY = 2/6*handles.Wavelength * (handles.ProbeNy/handles.ProbeWindowY);
set(handles.edit_ProbeMaxAngleX,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleX));
set(handles.edit_ProbeMaxAngleY,'String',sprintf('%.1f',1e3*handles.ProbeMaxAngleY));
        
% Output options:
handles.OutputFolder = get(handles.edit_OutputFolder,'String');
handles.SlicesBetweenOutputs = str2num(get(handles.edit_SlicesBetweenOutputs,'String'));
 
guidata(hObject, handles);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of readAllFields()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes on button press in pushbutton_AutoScanX.
function pushbutton_AutoScanX_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_AutoScanX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_AutoScanY.
function pushbutton_AutoScanY_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_AutoScanY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_AutoZOffset.
function pushbutton_AutoZOffset_Callback(hObject, eventdata, handles)

handles.sliceThickness = str2num(get(handles.edit_SliceThickness,'String'));
handles.PotentialOffsetZ = 0.5*handles.sliceThickness;
set(handles.edit_PotentialOffsetZ,'String',sprintf('%.3f',handles.PotentialOffsetZ));
set(handles.checkbox_CenterSlices,'Value',0);
guidata(hObject, handles);



% --- Executes on button press in radiobutton_Box.
function radiobutton_Box_Callback(hObject, eventdata, handles)
handles.BoxMode = 1;
set(handles.text_Nx,'String','Box:         ax:');
set(handles.text_Ny,'String','by:');
set(handles.text_Nz,'String','cz:');

if ~isfield(handles,'BoxX');
    if isfield(handles,'Mm')
        handles.BoxX = handles.Mm(1,1)*handles.NcellX;
        handles.BoxY = handles.Mm(2,2)*handles.NcellY;
        handles.BoxZ = handles.Mm(3,3)*handles.NcellZ;
    else
       handles.BoxX = 50;
       handles.BoxY = 50;
       handles.BoxZ = 100;
    end
end

set(handles.edit_NcellX,'String',handles.BoxX);
set(handles.edit_NcellY,'String',handles.BoxY);
set(handles.edit_NcellZ,'String',handles.BoxZ);

if get(handles.radiobutton_Box,'Value')
    handles.SliceThickness = handles.BoxZ/handles.NsliceTot;
else
    handles.SliceThickness = handles.Mm(3,3)*handles.NcellZ/handles.NsliceTot;
end
set(handles.edit_SliceThickness,'String',sprintf('%.4f',handles.SliceThickness));

guidata(hObject, handles);

% --- Executes on button press in radiobutton_Ncells.
function radiobutton_Ncells_Callback(hObject, eventdata, handles)

handles.BoxMode = 0;
set(handles.text_Nx,'String','Unit cells:  Nx:');
set(handles.text_Ny,'String','Ny:');
set(handles.text_Nz,'String','Nz:');

set(handles.edit_NcellX,'String',handles.NcellX);
set(handles.edit_NcellY,'String',handles.NcellY);
set(handles.edit_NcellZ,'String',handles.NcellZ);

if get(handles.radiobutton_Box,'Value')
    handles.SliceThickness = handles.BoxZ/handles.NsliceTot;
else
    handles.SliceThickness = handles.Mm(3,3)*handles.NcellZ/handles.NsliceTot;
end
set(handles.edit_SliceThickness,'String',sprintf('%.4f',handles.SliceThickness));


guidata(hObject, handles);




% --- Executes on button press in radiobutton_degree.
function radiobutton_degree_Callback(hObject, eventdata, handles)

% convert, if the unit was rad previously:
if handles.degree == 0
    handles.TiltX  = str2num(get(handles.edit_TiltX,'String'));
    handles.TiltY  = str2num(get(handles.edit_TiltY,'String'));
    handles.TiltZ  = str2num(get(handles.edit_TiltZ,'String'));
   
    handles.TiltX = handles.TiltX*180.0/pi;
    handles.TiltY = handles.TiltY*180.0/pi;
    handles.TiltZ = handles.TiltZ*180.0/pi;
       
    set(handles.edit_TiltX,'String',handles.TiltX);
    set(handles.edit_TiltY,'String',handles.TiltY);
    set(handles.edit_TiltZ,'String',handles.TiltZ);
 
end
handles.degree = 1;
guidata(hObject, handles);




% --- Executes on button press in radiobutton_radians.
function radiobutton_radians_Callback(hObject, eventdata, handles)

% convert, if the unit was rad previously:
if handles.degree == 1
    handles.TiltX  = str2num(get(handles.edit_TiltX,'String'));
    handles.TiltY  = str2num(get(handles.edit_TiltY,'String'));
    handles.TiltZ  = str2num(get(handles.edit_TiltZ,'String'));

    handles.TiltX = handles.TiltX*pi/180;
    handles.TiltY = handles.TiltY*pi/180;
    handles.TiltZ = handles.TiltZ*pi/180;
    
    set(handles.edit_TiltX,'String',handles.TiltX);
    set(handles.edit_TiltY,'String',handles.TiltY);
    set(handles.edit_TiltZ,'String',handles.TiltZ);
 
end
handles.degree = 0;
guidata(hObject, handles);





function edit_TiltX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TiltX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TiltX as text
%        str2double(get(hObject,'String')) returns contents of edit_TiltX as a double


% --- Executes during object creation, after setting all properties.
function edit_TiltX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TiltX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_TiltY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TiltY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TiltY as text
%        str2double(get(hObject,'String')) returns contents of edit_TiltY as a double


% --- Executes during object creation, after setting all properties.
function edit_TiltY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TiltY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_TiltZ_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TiltZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TiltZ as text
%        str2double(get(hObject,'String')) returns contents of edit_TiltZ as a double


% --- Executes during object creation, after setting all properties.
function edit_TiltZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TiltZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in pushbutton_Advanced.
function pushbutton_Advanced_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data.saveLevel     = handles.saveLevel;
data.printLevel    = handles.printLevel;
data.atomRadius    = handles.atomRadius;
data.savePotential = handles.savePotential;
data.saveTotalPotential = handles.saveTotalPotential;
data.potProgInterval = handles.potProgInterval;
data.propProgInterval = handles.propProgInterval;
data = advancedSettings(data);

handles.saveLevel     = data.saveLevel;
handles.printLevel    = data.printLevel;
handles.atomRadius    = data.atomRadius;
handles.savePotential = data.savePotential;
handles.saveTotalPotential = data.saveTotalPotential;
handles.potProgInterval = data.potProgInterval;
handles.propProgInterval = data.propProgInterval;
% fprintf('saveLevel: %d, printLevel: %d, atomRadius: %.1f, savePot: %d (%d, %d)\n',handles.saveLevel,...
%    handles.printLevel,handles.atomRadius,handles.savePotential,handles.propProgInterval,handles.potProgInterval); 
guidata(hObject, handles);

function edit_dE_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dE as text
%        str2double(get(hObject,'String')) returns contents of edit_dE as a double


% --- Executes during object creation, after setting all properties.
function edit_dE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_NcellZ_Callback(hObject, eventdata, handles)
% hObject    handle to edit_NcellZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NcellZ as text
%        str2double(get(hObject,'String')) returns contents of edit_NcellZ as a double

if get(handles.radiobutton_Box,'Value')
    handles.BoxZ = str2double(get(hObject,'String'));
    handles.SliceThickness = handles.BoxZ/handles.NsliceTot;
else
    handles.NcellZ = str2double(get(hObject,'String'));
    handles.SliceThickness = handles.Mm(3,3)*handles.NcellZ/handles.NsliceTot;
end
set(handles.edit_SliceThickness,'String',sprintf('%.4f',handles.SliceThickness));
guidata(hObject, handles);



% --- Executes on mouse press over axes background.
function axes_model_ButtonDownFcn(hObject, eventdata, handles)


function scatter_ButtonDownFcn(src,eventdata,handles,hObject)
persistent RbboxHandle
if nargin < 4
    if ~isempty(RbboxHandle)
        set(RbboxHandle,'Visible','off');
    end
    RbboxHandle = [];
    return;
end

point1 = get(gca,'CurrentPoint');    % button down detected
rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
hold on
% If there was already a previous box drawn:
set(handles.pushbutton_ScanAccept,'Visible','on');
set(handles.pushbutton_ScanReject,'Visible','on');
set(handles.edit_Xstart,'String',sprintf('%.2f',p1(1)));
set(handles.edit_Xstop,'String',sprintf('%.2f',p1(1)+offset(1)));
set(handles.edit_Ystart,'String',sprintf('%.2f',p1(2)));
set(handles.edit_Ystop,'String',sprintf('%.2f',p1(2)+offset(2)));

if ~isempty(RbboxHandle)
    set(RbboxHandle,'Visible','off');
    RbboxHandle = [];
end

axis manual
RbboxHandle = plot(x,y);
%    answer = msgbox(sprintf('Dou you want to define the following scan window:\nX: %.2f .. %.2f, Y: %.2f .. %.2fA?',...
%    p1(1),p1(1)+offset(1),p1(2),p1(2)+offset(2)));
hold off    
guidata(hObject, handles);


% --- Executes on button press in pushbutton_ScanAccept.
function pushbutton_ScanAccept_Callback(hObject, eventdata, handles)

scatter_ButtonDownFcn();
set(handles.pushbutton_ScanAccept,'Visible','off');
set(handles.pushbutton_ScanReject,'Visible','off');
handles.Xstart = str2num(get(handles.edit_Xstart,'String'));
handles.Xstop  = str2num(get(handles.edit_Xstop,'String'));
handles.Ystart = str2num(get(handles.edit_Ystart,'String'));
handles.Ystop  = str2num(get(handles.edit_Ystop,'String'));

guidata(hObject, handles);
handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in pushbutton_ScanReject.
function pushbutton_ScanReject_Callback(hObject, eventdata, handles)

set(handles.edit_Xstart,'String',handles.Xstart);
set(handles.edit_Xstop,'String',handles.Xstop);
set(handles.edit_Ystart,'String',handles.Ystart);
set(handles.edit_Ystop,'String',handles.Ystop);

scatter_ButtonDownFcn();
set(handles.pushbutton_ScanAccept,'Visible','off');
set(handles.pushbutton_ScanReject,'Visible','off');

guidata(hObject, handles);





% --- Executes on button press in pushbutton_MoreAberrations.
function pushbutton_MoreAberrations_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_MoreAberrations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = readAllFields(hObject, eventdata, handles);
params = [handles.ProbeNx,handles.ProbeNy,handles.ProbeResolutionX,handles.ProbeResolutionY,handles.HighVoltage,handles.alpha];

c = zeros(6,1);
c(2) = 10*handles.Defocus;
c(4) = 1e7*handles.C3;
c(6) = 1e7*handles.C5;
handles.phi(2,2) = handles.AstigAngle*pi/180;
handles.a(2,2)   = 10*handles.AstigMag;

[a,phi,c,params] = aberrations(handles.a,handles.phi,c,params);
if length(params) > 1
    handles.HighVoltage = params(5);
    handles.alpha    = params(6);
    handles.a        = a;
    handles.phi      = phi;
    handles.Defocus  = 1e-1*c(2);
    handles.C3       = 1e-7*c(4);
    handles.C5       = 1e-7*c(6);
    handles.AstigAngle = phi(2,2)*180/pi;
    handles.AstigMag   = 1e-1*a(2,2);
    
    set(handles.edit_Defocus,'String',sprintf('%.1f',handles.Defocus));
    set(handles.edit_AstigMag,'String',sprintf('%.1f',handles.AstigMag));
    set(handles.edit_AstigAngle,'String',sprintf('%.1f',handles.AstigAngle));
    set(handles.edit_C3,'String',sprintf('%.3f',handles.C3));
    set(handles.edit_alpha,'String',handles.alpha);
    set(handles.edit_HighVoltage,'String',handles.HighVoltage);
    
    guidata(hObject, handles);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEM Radiobutton Callback
% --- Executes on button press in radiobutton_STEM.
function handles = radiobutton_STEM_Callback(hObject, eventdata, handles)
set(handles.uipanel_ScanningWindow,'Title','Scan Window');
set(handles.uipanel_ScanningWindow,'Visible','on');
handles.modeNum = 1;
handles.mode = 'STEM';
% set(handles.uipanel_TEMcover1,'Visible','off');

% reload the old setting, if they have been saved previously
if isfield(handles,'XstartOld')
    handles.Xstart = handles.XstartOld;
    handles.Ystart = handles.YstartOld;
    handles.Xstop = handles.XstopOld;
    handles.Ystop = handles.YstopOld;
    handles.Xpixels = handles.XpixelsOld;
    handles.Ypixels = handles.YpixelsOld;
    
    handles.ProbeWindowX = handles.ProbeWindowXOld;
    handles.ProbeWindowY = handles.ProbeWindowYOld;
    handles.ProbeResolutionX = handles.ProbeResolutionXOld;
    handles.ProbeResolutionY = handles.ProbeResolutionYOld;    
end
set(handles.edit_ProbeWindowX,'Enable','on');
set(handles.edit_ProbeWindowY,'Enable','on');
set(handles.edit_ProbeResolutionX,'Enable','on');
set(handles.edit_ProbeResolutionY,'Enable','on');

set(handles.edit_Xstart,'Visible','on');
set(handles.edit_Xstop,'Visible','on');
set(handles.edit_Xpixels,'Visible','on');
set(handles.edit_Ystart,'Visible','on');
set(handles.edit_Ystop,'Visible','on');
set(handles.edit_Ypixels,'Visible','on');

set(handles.text3,'String','X     start:');  % Xstart
set(handles.text6,'String','Y     start:');  % Ystart
set(handles.text3,'Visible','on');
set(handles.text6,'Visible','on');

set(handles.text4,'Visible','on');
set(handles.text5,'Visible','on');
set(handles.text7,'Visible','on');
set(handles.text8,'Visible','on');
set(handles.text9,'Visible','on');
set(handles.text10,'Visible','on');
set(handles.pushbutton_ImageSim,'Visible','off');
set(handles.pushbutton_MoreAberrations,'Visible','on');

set(handles.edit_Xstart,'String',handles.Xstart);
set(handles.edit_Xstop,'String',handles.Xstop);
set(handles.edit_Xpixels,'String',handles.Xpixels);
set(handles.edit_Ystart,'String',handles.Ystart);
set(handles.edit_Ystop,'String',handles.Ystop);
set(handles.edit_Ypixels,'String',handles.Ypixels);
    
set(handles.edit_ProbeWindowX,'String',handles.ProbeWindowX);
set(handles.edit_ProbeWindowY,'String',handles.ProbeWindowX);
set(handles.edit_ProbeResolutionX,'String',handles.ProbeResolutionX);
set(handles.edit_ProbeResolutionY,'String',handles.ProbeResolutionX);


handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
% End of STEM radiobutton Callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CBED Radiobutton Callback
% --- Executes on button press in radiobutton_CBED.
function handles = radiobutton_CBED_Callback(hObject, eventdata, handles)
set(handles.uipanel_ScanningWindow,'Title','Probe Position');
set(handles.uipanel_ScanningWindow,'Visible','on');
handles.modeNum = 2;
handles.mode = 'CBED';
% set(handles.uipanel_TEMcover1,'Visible','off');


% reload the old setting, if they have been saved previously
if isfield(handles,'XstartOld')
    handles.Xstart = 0.5*(handles.XstartOld+handles.XstopOld);
    handles.Ystart = 0.5*(handles.YstartOld+handles.YstopOld);
    
    handles.ProbeWindowX = handles.ProbeWindowXOld;
    handles.ProbeWindowY = handles.ProbeWindowYOld;
    handles.ProbeResolutionX = handles.ProbeResolutionXOld;
    handles.ProbeResolutionY = handles.ProbeResolutionYOld;    
end
handles.Xstop = handles.Xstart;
handles.Ystop = handles.Ystart;
handles.Xpixels = 1;
handles.Ypixels = 1;

set(handles.edit_ProbeWindowX,'Enable','on');
set(handles.edit_ProbeWindowY,'Enable','on');
set(handles.edit_ProbeResolutionX,'Enable','on');
set(handles.edit_ProbeResolutionY,'Enable','on');

set(handles.edit_Xstart,'Visible','on');
set(handles.edit_Ystart,'Visible','on');
set(handles.text3,'String','X  position:');  % Xstart
set(handles.text6,'String','Y  position:');  % Ystart
set(handles.text3,'Visible','on');  % Xstart
set(handles.text6,'Visible','on');  % Ystart
    
set(handles.edit_Xstop,'Visible','off');
set(handles.edit_Xpixels,'Visible','off');
set(handles.edit_Ystop,'Visible','off');
set(handles.edit_Ypixels,'Visible','off');

set(handles.text4,'Visible','off');
set(handles.text5,'Visible','off');
set(handles.text7,'Visible','off');
set(handles.text8,'Visible','off');
set(handles.text9,'Visible','off');
set(handles.text10,'Visible','off');
set(handles.pushbutton_ImageSim,'Visible','off');
set(handles.pushbutton_MoreAberrations,'Visible','on');


set(handles.edit_ProbeWindowX,'String',handles.ProbeWindowX);
set(handles.edit_ProbeWindowY,'String',handles.ProbeWindowX);
set(handles.edit_ProbeResolutionX,'String',handles.ProbeResolutionX);
set(handles.edit_ProbeResolutionY,'String',handles.ProbeResolutionX);

set(handles.edit_Xstart,'String',handles.Xstart);
set(handles.edit_Xstop,'String',handles.Xstop);
set(handles.edit_Xpixels,'String',handles.Xpixels);
set(handles.edit_Ystart,'String',handles.Ystart);
set(handles.edit_Ystop,'String',handles.Ystop);
set(handles.edit_Ypixels,'String',handles.Ypixels);


guidata(hObject, handles);
handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEM Radiobutton Callback
% --- Executes on button press in radiobutton_TEM.
function handles = radiobutton_TEM_Callback(hObject, eventdata, handles)

if isempty(handles.atomPos) 
    uiwait(msgbox('Please load a model structure first','Error','modal'));
    switch handles.modeNum
        case 1
            set(handles.radiobutton_STEM,'Value',1);
            handles = radiobutton_STEM_Callback(hObject, eventdata, handles);
        case 2
            set(handles.radiobutton_CBED,'Value',1);
            handles = radiobutton_CBED_Callback(hObject, eventdata, handles);
    end
    guidata(hObject, handles);
    return;
end



set(handles.uipanel_ScanningWindow,'Title','');
set(handles.uipanel_ScanningWindow,'Visible','off');

% Save the STEM scan settings:
if handles.modeNum > 0
    handles = readAllFields(hObject, eventdata, handles);
    handles.XstartOld = handles.Xstart;
    handles.YstartOld = handles.Ystart;
    handles.XstopOld = handles.Xstop;
    handles.YstopOld = handles.Ystop;
    handles.XpixelsOld = handles.Xpixels;
    handles.YpixelsOld = handles.Ypixels;

    handles.ProbeWindowXOld = handles.ProbeWindowX;
    handles.ProbeWindowYOld = handles.ProbeWindowY;
    handles.ProbeResolutionXOld = handles.ProbeResolutionX;
    handles.ProbeResolutionYOld = handles.ProbeResolutionY;    
end

handles.modeNum = 0;
handles.mode = 'TEM';

if isfield(handles,'boundingBox')
    xmax = max(handles.boundingBox(:,1));
    xmin = min(handles.boundingBox(:,1));
    ymax = max(handles.boundingBox(:,2));
    ymin = min(handles.boundingBox(:,2));    
else
    coords = handles.atomPos(:,2:4);
    coords = [coords; [1 0 0]; [0 1 0]];
    coords = coords*handles.Mm;
    xmin = min(coords(:,1));
    xmax = max(coords(:,1));
    ymin = min(coords(:,2));
    ymax = max(coords(:,2));
end

handles.ProbeWindowX = xmax -xmin;
handles.ProbeWindowY = ymax -ymin;
handles.ProbeResolutionX = handles.ProbeWindowX/handles.ProbeNx;
handles.ProbeResolutionY = handles.ProbeWindowY/handles.ProbeNy;
set(handles.edit_ProbeWindowX,'String',handles.ProbeWindowX);
set(handles.edit_ProbeWindowY,'String',handles.ProbeWindowY);
set(handles.edit_ProbeResolutionX,'String',handles.ProbeResolutionX);
set(handles.edit_ProbeResolutionY,'String',handles.ProbeResolutionY);
set(handles.edit_ProbeWindowX,'Enable','off');
set(handles.edit_ProbeWindowY,'Enable','off');
set(handles.edit_ProbeResolutionX,'Enable','off');
set(handles.edit_ProbeResolutionY,'Enable','off');

% Adjust the Scan setting to match the TEM requirements:
if (0) 
    handles.Xstart = 0.5*handles.ProbeWindowX;
    handles.Xstop = handles.Xstart;
    handles.Xpixels = 1;
    handles.Ystart = 0.5*handles.ProbeWindowY;
    handles.Ystop = handles.Ystart;
    handles.Ypixels = 1;
    
    set(handles.edit_Xstart,'String',handles.Xstart);
    set(handles.edit_Xstop,'String',handles.Xstop);
    set(handles.edit_Xpixels,'String',handles.Xpixels);
    set(handles.edit_Ystart,'String',handles.Ystart);
    set(handles.edit_Ystop,'String',handles.Ystop);
    set(handles.edit_Ypixels,'String',handles.Ypixels);
else
    set(handles.edit_Xstart,'Visible','off');
    set(handles.edit_Xstop,'Visible','off');
    set(handles.edit_Xpixels,'Visible','off');
    set(handles.edit_Ystart,'Visible','off');
    set(handles.edit_Ystop,'Visible','off');
    set(handles.edit_Ypixels,'Visible','off');
    set(handles.text3,'Visible','off');
    set(handles.text4,'Visible','off');
    set(handles.text5,'Visible','off');
    set(handles.text6,'Visible','off');
    set(handles.text7,'Visible','off');
    set(handles.text8,'Visible','off');
    set(handles.text9,'Visible','off');
    set(handles.text10,'Visible','off');
    set(handles.pushbutton_ImageSim,'Visible','on');
    set(handles.pushbutton_MoreAberrations,'Visible','off');
end
% set(handles.uipanel_TEMcover1,'Visible','on');
set(handles.pushbutton_ImageSim,'Visible','on');

handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
% End of TEM Radiobutton Callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







    
    % --- Executes on button press in pushbutton_ModelProperties.
function pushbutton_ModelProperties_Callback(hObject, eventdata, handles)
if (0)
    set(handles.radiobutton_viewTop,'Value',1);
    handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
end
DisplayModelProperties(handles);



% --- Executes on button press in pushbutton_ImageSim.
function pushbutton_ImageSim_Callback(hObject, eventdata, handles)
% system(fullfile(handles.binPath,'imageSim.exe'));
system('imageSim');





function edit_BtiltX_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_BtiltX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_BtiltY_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_BtiltY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_TiltBack.
function checkbox_TiltBack_Callback(hObject, eventdata, handles)




function edit_DwellTime_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_DwellTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Brightness_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_Brightness_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton_Zone.
function pushbutton_Zone_Callback(hObject, eventdata, handles)
if isfield(handles,'zone')
    defaultanswer = handles.zone;
else
    defaultanswer = {'0','0','1'};
end
prompt = {'h:','k:','l:'};
numlines = 1;
answer=inputdlg(prompt,'Define zone axis',numlines,defaultanswer);
if length(answer) < 3
    return
end
zone = [str2num(answer{1}); str2num(answer{2}); str2num(answer{3})];

refZone = [0; 0; 1];
[rotAngles,Mrot] = rotateZoneAxis(zone,handles.Mm,refZone);

handles.Mrot = Mrot;
handles.zone  = answer;
handles.TiltX = rotAngles(1)*pi/180;
handles.TiltY = rotAngles(2)*pi/180;
handles.TiltZ = rotAngles(3)*pi/180;

handles.degree = get(handles.radiobutton_degree,'Value');

if handles.degree
    set(handles.edit_TiltX,'String',sprintf('%.3f',rotAngles(1)));
    set(handles.edit_TiltY,'String',sprintf('%.3f',rotAngles(2)));
    set(handles.edit_TiltZ,'String',sprintf('%.3f',rotAngles(3)));
else
    set(handles.edit_TiltX,'String',sprintf('%.3f',rotAngles(1)*pi/180));
    set(handles.edit_TiltY,'String',sprintf('%.3f',rotAngles(2)*pi/180));  
    set(handles.edit_TiltZ,'String',sprintf('%.3f',rotAngles(3)*pi/180));  
end
guidata(hObject, handles);
handles = pushbutton_UpdateView_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
