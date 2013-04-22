function handles = readConfigFile(fileName,handles,askUser)

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
elseif strcmp(value,'TEM') == 1
    handles.mode = 'TEM';
    handles.modeNum = 0;
elseif strcmp(value,'CBED') == 1
    handles.mode = 'CBED';
    handles.modeNum = 2;
else
    msgbox('This file does not seem to be a valid input file');
    cd(oldPath);
    return
end
% fprintf('Verified configuration for STEM mode!\n');


handles.BeamTiltX = 0;
handles.BeamTiltY = 0;
handles.TiltBack  = 0;
handles.TDSruns = 30;



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
if value(1) == 'y', handles.CenterSlices = 1; end

% Periodicity of super-cell:
handles.PeriodicZ = 0;
handles.PeriodicXY = 0;
value = extractParam('periodicXY',nameArray,pValues);
if value(1) == 'y', handles.PeriodicXY = 1; end
value = extractParam('periodicZ',nameArray,pValues);
if value(1) == 'y', handles.PeriodicZ = 1; end

% number of cells:
handles.NsliceTot      = handles.NsubSlabs*handles.Nslice;

% TDS settings
handles.TDS = 0;
value = extractParam('tds',nameArray,pValues);
if value(1) == 'y'
    handles.TDS = 1;
end
handles.Temperature = str2num(extractParam('temperature',nameArray,pValues));
handles.TDSruns = str2num(extractParam('Runs for averaging',nameArray,pValues));

% extract Box definition:
handles.BoxMode = 0;
value = extractParam('Cube',nameArray,pValues);
if ~isempty(value)
    handles.BoxX = value(1);
    handles.BoxY = value(2);
    handles.BoxZ = value(3);
    if handles.SliceThickness == 0
        handles.SliceThickness = handles.BoxZ/handles.NsliceTot;
    end
    handles.BoxMode = 1;
else
    if handles.SliceThickness == 0
        handles.SliceThickness = handles.Mm(3,3)*handles.NcellZ/handles.NsliceTot;
    end
end


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


% microscope parameters:
handles.HighVoltage = str2num(extractParam('v0',nameArray,pValues));
handles.Wavelength = wavelength(handles.HighVoltage);

value = extractParam('defocus',nameArray,pValues);
handles.Defocus = str2num(value);

handles.alpha = str2num(extractParam('alpha',nameArray,pValues));
handles.C3 = str2num(extractParam('Cs',nameArray,pValues));
handles.C5 = str2num(extractParam('C5',nameArray,pValues));  % C5 in mm
handles.c(6) = 1e7*handles.C5;
% set(handles.edit_C5,'String',handles.C5);
handles.a = zeros(6);
handles.phi = zeros(6);

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
% Delta in nm:
handles.Delta   = handles.Cc*1e3*handles.dE/handles.HighVoltage;

% folder and output options:
handles.SlicesBetweenOutputs = str2num(extractParam('slices between outputs',nameArray,pValues));
handles.OutputFolder = extractParam('Folder',nameArray,pValues);



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

handles = LoadModel(handles);

% Read the probe array parameters:
handles.ProbeNx = str2num(extractParam('nx',nameArray,pValues));
handles.ProbeNy = str2num(extractParam('ny',nameArray,pValues));
if handles.ProbeNy <= 0
    handles.ProbeNy = handles.ProbeNx;
end
handles.ProbeResolutionX = str2num(extractParam('resolutionX',nameArray,pValues));
handles.ProbeResolutionY = str2num(extractParam('resolutionY',nameArray,pValues));

% fprintf('Resol: %g %g %d %d %g %g\n',handles.ProbeResolutionX,handles.ProbeResolutionY,handles.ProbeNx,handles.ProbeNy,handles.Mm(1,1),handles.Mm(2,2));
if handles.ProbeResolutionX <= 1e-4
    handles.ProbeResolutionX = handles.Mm(1,1)/handles.ProbeNx;
end
if handles.ProbeResolutionY <= 1e-4
    handles.ProbeResolutionY = handles.Mm(2,2)/handles.ProbeNy;
end

% Define the additional parameters:
handles.ProbeWindowX = handles.ProbeResolutionX*handles.ProbeNx;
handles.ProbeWindowY = handles.ProbeResolutionY*handles.ProbeNy;

handles.ProbeMaxAngleX = 2/6*handles.Wavelength * (handles.ProbeNx/handles.ProbeWindowX);
handles.ProbeMaxAngleY = 2/6*handles.Wavelength * (handles.ProbeNy/handles.ProbeWindowY);

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
    case 2  % CBED
        handles.Xstart = str2num(extractParam('scan_x_start',nameArray,pValues));
        handles.Xstop  = handles.Xstart;
        handles.Xpixels = 1;
        
        handles.Ystart = str2num(extractParam('scan_y_start',nameArray,pValues));
        handles.Ystop  = handles.Ystart;
        handles.Ypixels = 1;
        
    case 0  % TEM mode:
        handles.Xstart = 0.5*handles.ProbeWindowX;
        handles.Xstop = handles.Xstart;
        handles.Xpixels = 1;
        handles.Ystart = 0.5*handles.ProbeWindowY;
        handles.Ystop = handles.Ystart;
        handles.Ypixels = 1;
end




handles.saveLevel     = str2num(extractParam('save level',nameArray,pValues));
handles.printLevel    = str2num(extractParam('print level',nameArray,pValues));
handles.atomRadius    = str2num(extractParam('atom radius',nameArray,pValues));
handles.savePotential   = strcmp(extractParam('save potential',nameArray,pValues),'yes');
handles.saveTotalPotential   = strcmp(extractParam('save projected potential',nameArray,pValues),'yes');
handles.potProgInterval = str2num(extractParam('potential progress interval',nameArray,pValues));
handles.propProgInterval = str2num(extractParam('propagation progress interval',nameArray,pValues));



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


function handles = LoadModel(handles)
[pathname, filename, ext] = fileparts(handles.posFileName);
filename = [filename ext];

[coords, aType, Mm, DW, occ, charge] = readCFG_qstem(fullfile(pathname,filename),[0 0 0],2);

if ~isempty(coords)
    handles.posFileName = fullfile(pathname,filename);
    handles.atomPos = [aType coords  DW occ charge];
    % [aType coords  DW occ charge]
    clear coords aType DW occ charge
    handles.Mm = Mm;
end
