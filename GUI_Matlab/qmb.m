function varargout = qmb(varargin)
% QMB M-file for qmb.fig
%      QMB, by itself, creates a new QMB or raises the existing
%      singleton*.
%
%      H = QMB returns the handle to a new QMB or the handle to
%      the existing singleton*.
%
%      QMB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QMB.M with the given input arguments.
%
%      QMB('Property','Value',...) creates a new QMB or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before qmb_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to qmb_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help qmb

% Last Modified by GUIDE v2.5 21-Nov-2014 12:52:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @qmb_OpeningFcn, ...
                   'gui_OutputFcn',  @qmb_OutputFcn, ...
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


% --- Executes just before qmb is made visible.
function qmb_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to qmb (see VARARGIN)

% Choose default command line output for qmb
handles.output = hObject;
set(gcf,'Name','qstem model builder V2.3');
handles.imageRotation = 0;
handles.img = [];
handles.imgRot = [];
handles.indSel = [];
handles.model_path = [];
handles.grainCount = 0;
handles.BoxX = 100;
handles.BoxY = 100;
handles.offsWp = [];
handles.NwarpX = 0;
handles.NwarpY = 0;


grainIndex = 1;
handles.grain{grainIndex}.TiltX = 0;
handles.grain{grainIndex}.TiltY = 0;
handles.grain{grainIndex}.TiltZ = 0;
handles.grain{grainIndex}.OffsetX = 0;
handles.grain{grainIndex}.OffsetY = 0;
handles.grain{grainIndex}.OffsetZ = 0;
handles.grain{grainIndex}.minZ = str2num(get(handles.edit_minZ,'String'));
handles.grain{grainIndex}.maxZ = str2num(get(handles.edit_maxZ,'String'));


axes(handles.axes1);
xlim([0 handles.BoxX]);
ylim([0 handles.BoxY]);
set(handles.edit_BoxX,'String',handles.BoxX);
set(handles.edit_BoxY,'String',handles.BoxY);
handles.indSel = [];

colormap('default');
handles.cmap = colormap();

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes qmb wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = qmb_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in pushbutton_LoadImage.
function pushbutton_LoadImage_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LoadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (0)
    img_path = 'C:\Users\koch\Christoph\FRWR3\Roesner\recon_10images\phase.img';
else
   [filename, pathname] = uigetfile({'*.dm3';'*.img';'*.tif';'*.jpg'},'Load an image');
   if isequal(filename,0) || isequal(pathname,0)
       % disp('User pressed cancel')
       return
   else
       img_path = fullfile(pathname, filename);
   end 
end
[p,name,ext] = fileparts(img_path);

if strcmp(ext,'.img')
    [img,t,dx,dy] = binread2D(img_path);
else
    if strcmp(ext,'.dm3')
        [img,dx,dy] = ReadDM3_Matlab(img_path);
        dx = dx*10;
        dy = dy*10;
        img = single(img.');
        % [img params tilt tags] = readDM3(img_path);
        % dx = str2num(get(handles.edit_ImageDx,'String'));
        % dy = str2num(get(handles.edit_ImageDy,'String'));
    else
        img = imread(img_path);
        if size(img,3) > 1
            img = single(sum(img,3));
        else
            img = single(img);
        end
        dx = str2num(get(handles.edit_ImageDx,'String'));
        dy = str2num(get(handles.edit_ImageDy,'String'));
    end
end
limLow = min(min(img));
limHigh = max(max(img));
img = 100*(img-limLow)/(limHigh-limLow);
limLow = 0;
limHigh = 100;

handles.imgRot = [];
handles.img = img;
handles.imageDx = dx;
handles.imageDy = dy;
set(handles.edit_ImageDx,'String',sprintf('%.5f',handles.imageDx));
set(handles.edit_ImageDy,'String',sprintf('%.5f',handles.imageDy));
set(handles.edit_LimLow,'String',sprintf('%.0f',limLow));
set(handles.edit_LimHigh,'String',sprintf('%.0f',limHigh));

[Ny,Nx]= size(img);
set(handles.uipanel1,'Title',sprintf('Image: %s (%d x %d pixels)',name,Nx,Ny));

handles.BoxX = handles.imageDx*Nx;
handles.BoxY = handles.imageDy*Ny;
set(handles.edit_BoxX,'String',handles.BoxX );
set(handles.edit_BoxY,'String',handles.BoxY);
set(handles.edit_ImageOffsetX,'String',0);
set(handles.edit_ImageOffsetY,'String',0);


guidata(hObject, handles);
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_AutoSize.
function pushbutton_AutoSize_Callback(hObject, eventdata, handles)
handles.imageDx = str2double(get(handles.edit_ImageDx,'String'));
handles.imageDy = str2double(get(handles.edit_ImageDy,'String'));
[Ny,Nx] = size(handles.img);
offsetX = str2double(get(handles.edit_ImageOffsetX,'String'));
offsetY = str2double(get(handles.edit_ImageOffsetY,'String'));

handles.BoxX = handles.imageDx*Nx-offsetX;
handles.BoxY = handles.imageDy*Ny-offsetY;
set(handles.edit_BoxX,'String',sprintf('%.1f',handles.BoxX));
set(handles.edit_BoxY,'String',sprintf('%.1f',handles.BoxY));
guidata(hObject, handles);
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_Redraw.
function pushbutton_Redraw_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Redraw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
offsetX = str2double(get(handles.edit_ImageOffsetX,'String'));
offsetY = str2double(get(handles.edit_ImageOffsetY,'String'));
handles.BoxX = str2double(get(handles.edit_BoxX,'String'));
handles.BoxY = str2double(get(handles.edit_BoxY,'String'));
handles.imageDx = str2double(get(handles.edit_ImageDx,'String'));
handles.imageDy = str2double(get(handles.edit_ImageDy,'String'));
limLow = str2double(get(handles.edit_LimLow,'String'));
limHigh = str2double(get(handles.edit_LimHigh,'String'));


% if a background image has already been loaded:
if ~isempty(handles.img)
    rot = str2double(get(handles.edit_ImageRotation,'String'));
    % if the rotation has changed from last time:
    [Ny,Nx] = size(handles.img);
    [x,y] = meshgrid(-Nx/2+[0:Nx-1],-Ny/2+[0:Ny-1]);
    if (rot ~= handles.imageRotation) || isempty(handles.imgRot)
        handles.imageRotation = rot;
        set(handles.edit_ImageRotation,'String',rot);
        if (rot == 0)
            handles.imgRot = handles.img;
            guidata(hObject, handles);
        else
            xp = cos(-rot*pi/180.0)*x+sin(-rot*pi/180.0)*y;
            yp = cos(-rot*pi/180.0)*y-sin(-rot*pi/180.0)*x;
            handles.imgRot = interp2(x,y,handles.img,xp,yp,'linear',0);
            guidata(hObject, handles);
        end        
    end        
    set(gca,'ButtonDownFcn',{@assign_ButtonDownFcn,handles,hObject});
    h = imagesc(handles.imageDx*[0:Nx-1]-offsetX,handles.imageDy*[Ny-1:-1:0]-offsetY,handles.imgRot);
    % imagesc(handles.imageDx*[0:Nx-1]-offsetX,handles.imageDy*[0:Ny-1]-offsetY,handles.imgRot);
    set(h,'ButtonDownFcn',{@assign_ButtonDownFcn,handles,hObject});
    guidata(hObject, handles);

    set(gca,'YDir','normal');
    colormap('gray');
    xlabel('x in Å');
    ylabel('y in Å');
    axis equal; axis tight;
    
    %    if isfield(handles,'maxX')
    %       xlim([0 handles.maxX]);
    %       ylim([0 handles.maxY]);
    %    end
    caxis([limLow limHigh]);
else
    cla(handles.axes1);
end
xlim([0 handles.BoxX]);
ylim([0 handles.BoxY]);



% --- Executes on button press in pushbutton_LoadModel.
function pushbutton_LoadModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LoadModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (0)
    % handles.model_path = 'C:\Users\koch\Christoph\Matlab\QSTEM\Examples\Pbcorner_inAl_new.cfg';
    handles.model_path = 'C:\Users\koch\Christoph\Matlab\QSTEM\Examples\Pbcorner_inAl3.cfg';
else
   [filename, pathname] = uigetfile('*.cfg','Load model');
   if isequal(filename,0) || isequal(pathname,0)
       % disp('User pressed cancel')
       return
   else
       handles.model_path = fullfile(pathname, filename);
   end 
end
% Zmax = str2double(get(handles.edit_Zmax,'String'));


[coords, aType, Mm, DW, occ, charge] = readCFG_qstem(handles.model_path,[0 0 0],1);
handles.Mm = Mm;
% handles.indDisp = find(coords(:,3) < Zmax).';
handles.atomPos = [coords aType DW occ charge];

% sort atom positions:
sortrows(handles.atomPos,[1 2 3 4]);

handles.maxX = Mm(1,1);
handles.maxY = Mm(2,2);
handles.Natoms = size(handles.atomPos,1);
set(handles.uipanel2,'Title',sprintf('Model (%d atoms)',handles.Natoms));

clear coords aType DW charge occ

guidata(hObject, handles);
% Update plot:
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);




% --- Executes on button press in pushbutton_UpdateModel.
function pushbutton_UpdateModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_UpdateModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.atomSize = str2double(get(handles.edit_AtomSize,'String'));
cmap = colormap('jet');
cmap(1,:) = [1 1 1];
% cmap(81:85,:) = repmat([1 0 0],5,1);
xlimits = xlim();
ylimits = ylim();


pushbutton_Redraw_Callback(hObject, eventdata, handles);
% xlim([0 handles.maxX]);
% ylim([0 handles.maxY]);
% handles.atomPos(1:10,:)

hold on
Zmin = min(handles.atomPos(:,4));
Zmax = max(handles.atomPos(:,4));
if Zmin == Zmax
    Zmax = Zmin+1;
end
cmap = handles.cmap;

Zscale = mean(handles.atomPos(:,4));
% sortrows(coords,[1 4]);
coordsNext = handles.atomPos(:,1:4);
set(gca,'ButtonDownFcn',{@assign_ButtonDownFcn,handles,hObject});
while ~isempty(coordsNext)
    Z = coordsNext(1,4);
    ind = find(coordsNext(:,4) == Z);
    col = 1+round((Z-Zmin)/(Zmax-Zmin)*(size(cmap,1)-1));
    h = plot(coordsNext(ind,1),coordsNext(ind,2),'.','MarkerSize',handles.atomSize*sqrt(Z/Zscale),'MarkerEdgeColor',cmap(col,:));
    hold on
    set(h,'ButtonDownFcn',{@assign_ButtonDownFcn,handles,hObject});
    coordsNext = coordsNext(find(coordsNext(:,4) ~= Z),:);    
end
if max(handles.indSel) > size(handles.atomPos,1)
    handles.indSel = [];
end
if ~isempty(handles.indSel)
    h = plot(handles.atomPos(handles.indSel,1),handles.atomPos(handles.indSel,2),'.','MarkerSize',handles.atomSize,'MarkerEdgeColor',[0.1 1 0.1]);
    hold on
    set(h,'ButtonDownFcn',{@assign_ButtonDownFcn,handles,hObject});    
end
if get(handles.checkbox_showLines,'Value')
    if isfield(handles,'grainCount')
        for grainIndex = 1:handles.grainCount
            if (isfield(handles.grain{grainIndex},'boundX')) && (isfield(handles.grain{grainIndex},'atomPos'))
                l = line(handles.grain{grainIndex}.boundX,handles.grain{grainIndex}.boundY);
                set(l,'LineWidth',1,'Color',[1 0 0]);
            end
        end
    end
end
% fprintf('Will display WP\n');
displayWarpPoints(handles);
hold off
xlim(xlimits);
ylim(ylimits);

guidata(hObject, handles);






% --- Executes on button press in pushbutton_WarpPositions.
function pushbutton_WarpPositions_Callback(hObject, eventdata, handles)
scale = 5e3;
% This step produces first a mesh of warping vertices by which the atom
% positions are warped.
NwarpX = round(str2double(get(handles.edit_WarpNx,'String')));
NwarpY = round(str2double(get(handles.edit_warpNy,'String')));

% define warping points:
[wpX,wpY] = meshgrid(handles.BoxX/(NwarpX-1)*[0:NwarpX-1],handles.BoxY/(NwarpY-1)*[0:NwarpY-1]);
Nwp = prod(size(wpX));

offsWp = zeros(2*Nwp,1);

options=optimset('fminsearch');
options.TolFun = 1e-10;
options.TolX   = 1e-8;
option.Display = 'iter';

fprintf('Start refinement - please wait!\n');
offsWp = fminsearch(@(offs) warpPositions(offs,wpX,wpY,handles,scale),offsWp,options);
% offsX = scale*reshape(offsWp(1:Nwp),[NwarpX,NwarpY])
% offsY = scale*reshape(offsWp(Nwp+1:2*Nwp),[NwarpX,NwarpY])
[chi2,p] = warpPositions(offsWp,wpX,wpY,handles,scale);
fprintf('chi2=%g: %g .. %g, %g .. %g (%g, %g)\n',chi2,min(p(:,1)),max(p(:,1)),min(p(:,2)),max(p(:,2)),handles.imageDx,handles.imageDy);

% update atom positions:

handles.atomPos(:,1:2) = p;
guidata(hObject, handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);


function [chi2,p] = warpPositions(offsWp,wpX,wpY,handles,scale)

p0  = handles.atomPos(handles.indSel,1:2);        % original atom positions
img = handles.imgRot;
img = img-min(min(img));
[Ny,Nx] = size(img);
offsetX = str2double(get(handles.edit_ImageOffsetX,'String'));
offsetY = str2double(get(handles.edit_ImageOffsetY,'String'));
NwarpX = round(str2double(get(handles.edit_WarpNx,'String')));
NwarpY = round(str2double(get(handles.edit_warpNy,'String')));


% create offsets at atom positions
Nwp = length(offsWp)/2;
cOffsX = scale*reshape(offsWp(1:Nwp),[NwarpY,NwarpX]);
cOffsY = scale*reshape(offsWp(Nwp+1:2*Nwp),[NwarpY,NwarpX]);

if (Nwp < 9)
    atomOffsX = interp2(wpX,wpY,cOffsX,p0(:,1),p0(:,2),'linear',0);
    atomOffsY = interp2(wpX,wpY,cOffsY,p0(:,1),p0(:,2),'linear',0);
else
    atomOffsX = interp2(wpX,wpY,cOffsX,p0(:,1),p0(:,2),'cubic',0);
    atomOffsY = interp2(wpX,wpY,cOffsY,p0(:,1),p0(:,2),'cubic',0);
end
% atomOffsY(find(isnan(atomOffsY))) = 0;
% atomOffsX(find(isnan(atomOffsX))) = 0;

% create image index pixels:
p  = p0+[offsetX+atomOffsX, offsetY+atomOffsY];        % warped atom positions
p(:,1) = round(p(:,1)/handles.imageDx);
p(:,2) = round(p(:,2)/handles.imageDy);

% fprintf('%g .. %g, %g .. %g (%g, %g)\n',min(p(:,1)),max(p(:,1)),min(p(:,2)),max(p(:,2)),handles.imageDx,handles.imageDy);
% [p p0 atomOffsX atomOffsY]
% ind2sub:  first y-coordinate, then x-coordinate
ind = sub2ind([Ny,Nx],p(:,2),p(:,1));
% [p ind p0 atomOffsX atomOffsY]
chi2 = -sum(img(ind));
% fprintf('Chi2: %g (%d, %d non-zero entries)\n',-chi2,length(find(atomOffsX ~= 0)),length(find(atomOffsY ~= 0)));
if nargout > 1
    p(:,1) = p(:,1)*handles.imageDx-offsetX;
    p(:,2) = p(:,2)*handles.imageDy-offsetY;
end


function p = warpPositions_noCompare(offsWp,wpX,wpY,handles)

p0  = handles.atomPos(handles.indSel,1:2);        % original atom positions
NwarpX = round(str2double(get(handles.edit_WarpNx,'String')));
NwarpY = round(str2double(get(handles.edit_warpNy,'String')));


% create offsets at atom positions
Nwp = length(offsWp)/2;
cOffsX = reshape(offsWp(1:Nwp),[NwarpY,NwarpX]);
cOffsY = reshape(offsWp(Nwp+1:2*Nwp),[NwarpY,NwarpX]);

if (Nwp < 9)
    atomOffsX = interp2(wpX,wpY,cOffsX,p0(:,1),p0(:,2),'linear',0);
    atomOffsY = interp2(wpX,wpY,cOffsY,p0(:,1),p0(:,2),'linear',0);
else
    atomOffsX = interp2(wpX,wpY,cOffsX,p0(:,1),p0(:,2),'cubic',0);
    atomOffsY = interp2(wpX,wpY,cOffsY,p0(:,1),p0(:,2),'cubic',0);
end
% atomOffsY(find(isnan(atomOffsY))) = 0;
% atomOffsX(find(isnan(atomOffsX))) = 0;

% create image index pixels:
p  = [atomOffsX, atomOffsY];        % warped atom position offsets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function which displays the reference points used for warping
function displayWarpPoints(handles)
if get(handles.radiobutton_Warp,'Value') && get(handles.checkbox_ShowWarpPoints,'Value') && ~isempty(handles.indSel)
		hold on
		NwarpX = round(str2double(get(handles.edit_WarpNx,'String')));
		NwarpY = round(str2double(get(handles.edit_warpNy,'String')));
		Nwp = NwarpX*NwarpY;
		if  (NwarpX ~= handles.NwarpX) || (NwarpY ~= handles.NwarpY) 
			handles.NwarpX = NwarpX;
			handles.NwarpY = NwarpY;
			handles.offsWp = zeros(2*Nwp,1);
			handles.atomPosUnWarped = handles.atomPos(handles.indSel,:);
			% define the range
			maxX = max(handles.atomPosUnWarped(:,1));
			minX = min(handles.atomPosUnWarped(:,1));
			maxY = max(handles.atomPosUnWarped(:,2));
			minY = min(handles.atomPosUnWarped(:,2));
			% define warping points:
			[wpX,wpY] = meshgrid(minX+(maxX-minX)/(NwarpX-1)*[0:NwarpX-1],minY+(maxY-minY)/(NwarpY-1)*[0:NwarpY-1]);
			handles.wpX = wpX;
			handles.wpY = wpY;
		end

		
		wpX = handles.wpX;
		wpY = handles.wpY;
		handles.NwarpX = NwarpX;
		handles.NwarpY = NwarpY;
		
		ipointX = round(str2double(get(handles.edit_PointX,'String')));
		ipointY = round(str2double(get(handles.edit_PointY,'String')));
		if (ipointX < 1), ipointX = 1; end;
		if (ipointX > NwarpX), ipointX = NwarpX; end;
		if (ipointY < 1), ipointY = 1; end;
		if (ipointY > NwarpY), ipointY = NwarpY; end;
		set(handles.edit_PointX,'String',ipointX);
		set(handles.edit_PointY,'String',ipointY);
		
		if ~isempty(handles.offsWp)
			wpX = wpX+reshape(handles.offsWp(1:Nwp),[NwarpY,NwarpX]);
			wpY = wpY+reshape(handles.offsWp(Nwp+1:2*Nwp),[NwarpY,NwarpX]);
		end
		plot(reshape(wpX,NwarpX*NwarpY,1),reshape(wpY,NwarpX*NwarpY,1),'+','MarkerSize',handles.atomSize,'Color',[0 0 1])
		plot(wpX(1+NwarpY-ipointY,ipointX),wpY(1+NwarpY-ipointY,ipointX),'+','MarkerSize',handles.atomSize*1.3,'Color',[1 0 1])
		hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function which translates the relative positions of the warping mesh to
% the atom positions
function handles = applyWarping(handles,deltaX,deltaY)
NwarpX = round(str2double(get(handles.edit_WarpNx,'String')));
NwarpY = round(str2double(get(handles.edit_warpNy,'String')));
Nwp = NwarpX*NwarpY;
if (isempty(handles.offsWp)) || (NwarpX ~= handles.NwarpX) || (NwarpY ~= handles.NwarpY) || isempty(handles.atomPosUnWarped)
	handles.NwarpX = NwarpX;
	handles.NwarpY = NwarpY;
	handles.offsWp = zeros(2*Nwp,1);
	handles.atomPosUnWarped = handles.atomPos(handles.indSel,:);
	% define the range
	maxX = max(handles.atomPosUnWarped(:,1));
	minX = min(handles.atomPosUnWarped(:,1));
	maxY = max(handles.atomPosUnWarped(:,2));
	minY = min(handles.atomPosUnWarped(:,2));
	% define warping points:
	[wpX,wpY] = meshgrid(minX+(maxX-minX)/(NwarpX-1)*[0:NwarpX-1],minY+(maxY-minY)/(NwarpY-1)*[0:NwarpY-1]);
	handles.wpX = wpX;
	handles.wpY = wpY;
end
wpX = handles.wpX;
wpY = handles.wpY;
	
ipointX = round(str2double(get(handles.edit_PointX,'String')));
ipointY = round(str2double(get(handles.edit_PointY,'String')));
if (ipointX < 1), ipointX = 1; end; 
if (ipointX > NwarpX), ipointX = NwarpX; end; 
if (ipointY < 1), ipointY = 1; end; 
if (ipointY > NwarpY), ipointY = NwarpY; end; 
ipoint = sub2ind([NwarpX,NwarpY],1+NwarpY-ipointY,ipointX);
	
% Apply the shift of the currently selected warping node:	
handles.offsWp(ipoint) = handles.offsWp(ipoint) +deltaX;
handles.offsWp(Nwp+ipoint) = handles.offsWp(Nwp+ipoint) + deltaY;
p = warpPositions_noCompare(handles.offsWp,wpX,wpY,handles);
    
% update atom positions:
handles.atomPos(handles.indSel,1:2) = handles.atomPosUnWarped(:,1:2)+p;



% --- Executes on button press in pushbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)

distance = str2double(get(handles.edit_Distance,'String'));
if isempty(handles.indSel)
	handles.atomPosUnWarped = [];
	return;
end

if get(handles.radiobutton_Warp,'Value')
	% the following function shifts the current point to the left, i.e.
	% by the vector (-distance,0)
	handles = applyWarping(handles,-distance,0);
else
    handles.atomPos(handles.indSel,1) = handles.atomPos(handles.indSel,1) - distance;
end
guidata(hObject, handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);
% displayWarpPoints(handles);


% --- Executes on button press in pushbutton_up.
function pushbutton_up_Callback(hObject, eventdata, handles)

distance = str2double(get(handles.edit_Distance,'String'));
if isempty(handles.indSel)
	handles.atomPosUnWarped = [];
	return;
end

if get(handles.radiobutton_Warp,'Value')
	% the following function shifts the current point up, i.e.
	% by the vector (0,distance)
	handles = applyWarping(handles,0,distance);
else
    handles.atomPos(handles.indSel,2) = handles.atomPos(handles.indSel,2) + distance;
end
guidata(hObject, handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);
% displayWarpPoints(handles);



% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)

distance = str2double(get(handles.edit_Distance,'String'));
if isempty(handles.indSel)
	handles.atomPosUnWarped = [];
	return;
end

if get(handles.radiobutton_Warp,'Value')
	% the following function shifts the current point up, i.e.
	% by the vector (0,distance)
	handles = applyWarping(handles,distance,0);
else
    handles.atomPos(handles.indSel,1) = handles.atomPos(handles.indSel,1) + distance;
end
guidata(hObject, handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);
% displayWarpPoints(handles);

% --- Executes on button press in pushbutton_down.
function pushbutton_down_Callback(hObject, eventdata, handles)

distance = str2double(get(handles.edit_Distance,'String'));
if isempty(handles.indSel)
	handles.atomPosUnWarped = [];
	return;
end

if get(handles.radiobutton_Warp,'Value')
	% the following function shifts the current point up, i.e.
	% by the vector (0,distance)
	handles = applyWarping(handles,0,-distance);
else
    handles.atomPos(handles.indSel,2) = handles.atomPos(handles.indSel,2) - distance;
end
guidata(hObject, handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);
% displayWarpPoints(handles);





% --- Executes on button press in pushbutton_DeleteColumn.
function pushbutton_DeleteColumn_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_DeleteColumn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.indSel)
	handles.atomPosUnWarped = [];
   return
end
handles.atomPos(handles.indSel,4) = 0;
handles.indSel = [];
handles.offsWp = [];
ind = find(handles.atomPos(:,4) > 0);
handles.atomPos = handles.atomPos(ind,:);

% Zmax = str2double(get(handles.edit_Zmax,'String'));
% handles.indDisp = find(handles.atomPos(:,3) < Zmax).';

handles.Natoms = size(handles.atomPos,1);
set(handles.uipanel2,'Title',sprintf('Model (%d atoms)',handles.Natoms));


guidata(hObject, handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);



    



function assign_ButtonDownFcn(src,eventdata,handles,hObject)
persistent RbboxHandle

handles.atomPosUnWarped = [];
handles.offsWp = [];

point1 = get(gca,'CurrentPoint');    % button down detected
rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
if (offset(1) == 0) && (offset(2) == 0)
    handles.indSel = [];
    guidata(hObject, handles);
    pushbutton_UpdateModel_Callback(hObject, eventdata, handles);
    return
end
x1 = p1(1);
x2 = p1(1)+offset(1);
y1 = p1(2);
y2 = p1(2)+offset(2);

% determine selected atoms:
% handles.indSel = [];
if get(handles.checkbox_Sublattice,'Value')
	Zsel = str2num(get(handles.edit_Zsublattice,'String'));
	handles.indSel = [handles.indSel; find((handles.atomPos(:,1) > x1) & (handles.atomPos(:,1) < x2) & (handles.atomPos(:,2) > y1) & (handles.atomPos(:,2) < y2) & (handles.atomPos(:,4)==Zsel))];	
else
	handles.indSel = [handles.indSel; find((handles.atomPos(:,1) > x1) & (handles.atomPos(:,1) < x2) & (handles.atomPos(:,2) > y1) & (handles.atomPos(:,2) < y2))];
end
indSel = handles.indSel;
% coords = handles.atomPos(indSel,:)
guidata(hObject, handles);

% Update plot:
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);




% --- Executes on button press in pushbutton_SaveModel.
function pushbutton_SaveModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_SaveModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (0)
    fileName = [handles.model_path];
    if strcmp(fileName(end-7:end),'_new.cfg') == 0
        fileName = [handles.model_path(1:end-4), '_new.cfg'];
    end
else
    [name, pathname] = uiputfile('*.cfg', 'Write structure file');
    if isequal(name,0) || isequal(pathname,0)
        return
    else
        fileName = fullfile(pathname, name);
    end
    
end
% Mm = diag(max(handles.atomPos(:,1:3)))
handles.BoxX = str2double(get(handles.edit_BoxX,'String'));
handles.BoxY = str2double(get(handles.edit_BoxY,'String'));

Mm = eye(3);
 
Mm(3,3) = max(handles.atomPos(:,3));
% find the thickest grain:
if handles.grainCount > 0
    maxZ = 0;
    for grainIndex = 1:handles.grainCount
        if isfield(handles.grain{grainIndex},'maxZ')
            if handles.grain{grainIndex}.maxZ > maxZ
                maxZ = handles.grain{grainIndex}.maxZ;
            end
        end
    end
    if maxZ > 0
        Mm(3,3) = maxZ;
    end
end
    


Mm(1,1) = handles.BoxX;
Mm(2,2) = handles.BoxY;
writeCFG(fileName,Mm,handles.atomPos(:,4),handles.atomPos(:,1:3),[0 0 0],1,handles.atomPos(:,5),handles.atomPos(:,6));
fprintf('Saved atom coordinates to %s\n',fileName);



% --- Executes on button press in pushbutton_ForgetSelection.
function pushbutton_ForgetSelection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ForgetSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.indSel = [];
handles.offsWp = [];
fprintf('Selection deleted\n');
guidata(hObject, handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);



% --- Executes on button press in pushbutton_Zone.
function pushbutton_Zone_Callback(hObject, eventdata, handles)
if handles.grainCount < 1
    errordlg('Please load a crystal structure first !');
    return
end
grainIndex = str2double(get(handles.edit_GrainIndex,'String'));
if ~isfield(handles.grain{grainIndex},'Mm')
    errordlg(sprintf('Please load a crystal structure for grain %d first!',grainIndex));
    return
end

if isfield(handles.grain{grainIndex},'zone')
    defaultanswer = handles.grain{grainIndex}.zone;
else
    defaultanswer = {'0','0','1'};
end
prompt = {'h:','k:','l:'};
numlines = 1;
answer=inputdlg(prompt,'Define zone axis',numlines,defaultanswer);
if length(answer) < 3
    return
end

if (strcmp(answer{1},'0') && strcmp(answer{2},'0') && strcmp(answer{3},'0'))
    waitfor(errordlg({'[0 0 0] is an invalid zone axis! ','Setting to [0 0 1] instead.'}));
    answer = {'0';'0';'1'};
end
    
zone = [str2double(answer{1}); str2double(answer{2}); str2double(answer{3})];

refZone = [0; 0; 1];
[rotAngles,Mrot] = rotateZoneAxis(zone,handles.grain{grainIndex}.Mm,refZone);

handles.grain{grainIndex}.zone  = answer;
handles.grain{grainIndex}.TiltX = rotAngles(1);
handles.grain{grainIndex}.TiltY = rotAngles(2);
handles.grain{grainIndex}.TiltZ = rotAngles(3);
set(handles.edit_TiltX,'String',sprintf('%.3f',rotAngles(1)));
set(handles.edit_TiltY,'String',sprintf('%.3f',rotAngles(2)));
set(handles.edit_TiltZ,'String',sprintf('%.3f',rotAngles(3)));
set(handles.info_h,'String',handles.grain{grainIndex}.zone(1));
set(handles.info_k,'String',handles.grain{grainIndex}.zone(2));
set(handles.info_l,'String',handles.grain{grainIndex}.zone(3));
guidata(hObject, handles);
pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function builts the actual model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton_constructSuperCell(hObject, eventdata, handles)
clearedAtomPos = 0;
useMex = 1;
wb = -1;
% leave this function if we have not collected all the data yet.
if handles.grainCount < 1
    return
end

handles.BoxX = str2double(get(handles.edit_BoxX,'String'));
handles.BoxY = str2double(get(handles.edit_BoxY,'String'));


for grainIndex = 1:handles.grainCount
    if (useMex == 0)
        if (wb < 0)
            wb = waitbar((grainIndex-1)/handles.grainCount,'Constructing super cell');
        else
            waitbar((grainIndex-1)/handles.grainCount,wb,'Constructing super cell');
        end
    end
    if (isfield(handles.grain{grainIndex},'boundX')) && (isfield(handles.grain{grainIndex},'atomPos'))
        
        cx = 1; sx = 0;
        cy = 1; sy = 0;
        cz = 1; sz = 0;
        % offsetX = str2num(get(handles.edit_OffsetX,'String'));
        % offsetY = str2num(get(handles.edit_OffsetY,'String'));;
        % offsetZ = str2num(get(handles.edit_OffsetZ,'String'));;
        offsetX = 0;
        offsetY = 0;
        offsetZ = 0;
        minZ = str2num(get(handles.edit_minZ,'String'));
        maxZ = str2num(get(handles.edit_maxZ,'String'));
        if (isfield(handles.grain{grainIndex},'TiltX'))
            cx = cos(handles.grain{grainIndex}.TiltX*pi/180);
            sx = sin(handles.grain{grainIndex}.TiltX*pi/180);
            cy = cos(handles.grain{grainIndex}.TiltY*pi/180);
            sy = sin(handles.grain{grainIndex}.TiltY*pi/180);
            if (isfield(handles.grain{grainIndex},'TiltZ'))
                cz = cos(handles.grain{grainIndex}.TiltZ*pi/180);
                sz = sin(handles.grain{grainIndex}.TiltZ*pi/180);
            end
            if (isfield(handles.grain{grainIndex},'OffsetX'))
                offsetX = handles.grain{grainIndex}.OffsetX;
                offsetY = handles.grain{grainIndex}.OffsetY;
                offsetZ = handles.grain{grainIndex}.OffsetZ;
            end
            if (isfield(handles.grain{grainIndex},'minZ'))
                minZ = handles.grain{grainIndex}.minZ;
                maxZ = handles.grain{grainIndex}.maxZ;
            end
        end        
        Mx = [1 0 0;0 cx -sx;0 sx cx];
        My = [cy 0 sy;0 1 0; -sy 0 cy];
        Mz = [cz -sz 0;sz cz 0; 0 0 1];
        Mrot = Mz*My*Mx;

        % Now we have the rotation and offset for this unit cell.  All we
        % have left to do is to fill the whole box with this content and
        % cut out what we want to see.
        
        Mm        = Mrot*handles.grain{grainIndex}.Mm;
        Mminv     = inv(Mm);
        coords = handles.grain{grainIndex}.atomPos;
        % find the integer positions of the corners of the big cube:
        % Mbox   = [handles.BoxX 0 0;0 handles.BoxY 0; 0 0 maxZ];
        maxX = max(handles.grain{grainIndex}.boundX);
        maxY = max(handles.grain{grainIndex}.boundY);
        minX = min(handles.grain{grainIndex}.boundX);
        minY = min(handles.grain{grainIndex}.boundY);
        Mbox = [maxX 0 0;0 maxY 0; 0 0 maxZ];
        box    = [0 0 0;1 0 0; 1 1 0; 0 1 0;0 0 1;1 0 1; 1 1 1; 0 1 1]*Mbox;
        boxred = (Mminv*(box.')).';
        NxMin = min(floor(boxred(:,1)));         NxMax =  max(ceil(boxred(:,1)));
        NyMin = min(floor(boxred(:,2)));         NyMax =  max(ceil(boxred(:,2)));
        NzMin = min(floor(boxred(:,3)));         NzMax =  max(ceil(boxred(:,3)));
        % fprintf('%d .. %d | %d .. %d | %d .. %d\n',NxMin,NxMax,NyMin,NyMax,NzMin,NzMax);
        coords0 = coords;
        Natom = size(coords,1);
        coords = [];
        % handles.grain{grainIndex}.Mm
        Mm        = Mrot*handles.grain{grainIndex}.Mm;
        if (useMex == 0)
            for ix = NxMin:NxMax
                for iy = NyMin:NyMax
                    for iz = NzMin:NzMax
                        % newcoords = (coords0+repmat([0 ix iy iz 0 0 0],Natom,1))*Mm;
                        newcoords = [coords0(:,1) (Mm*(coords0(:,2:4).'+repmat([ix iy iz].',1,Natom))).'+repmat([offsetX,offsetY,offsetZ],Natom,1)  coords0(:,5:7)];
                        ind = find((newcoords(:,2)>=minX)  & (newcoords(:,2)<=maxX) & (newcoords(:,3)>=minY) & (newcoords(:,3)<=maxY)& (newcoords(:,4)<=maxZ) & (newcoords(:,4)>=minZ));
                        if ~isempty(ind)
                            coords = [coords; newcoords(ind,:)];
                        end
                    end
                end
                % fprintf('current range: %f..%f,%f..%f,%f..%f\n',min(newcoords(:,2)),min(newcoords(:,2)),min(newcoords(:,3)),min(newcoords(:,3)),min(newcoords(:,4)),min(newcoords(:,4)));
                % waitbar((grainIndex-1+(ix-NxMin)/(NxMax-NxMin))/handles.grainCount,wb,sprintf('Constructing super cell (%d)',size(coords,1)));
                waitbar((grainIndex-1+(ix-NxMin)/(NxMax-NxMin))/handles.grainCount,wb,sprintf('Constructing super cell'));
            end
        else
            % Call a mex file which implements the 3-foled nested loop
            % above.
            % Input: coords0, Mm,
            %        offsetX,offsetY,offsetZ,
            %        minX,maxX,minY,maxY,minZ,maxZ
            %        NxMin,NxMax,NyMin,NyMax,NzMin,NzMax
            % organize as arrays: offset, limits, countLimits
            blockOffset = [offsetX,offsetY,offsetZ];
            blockLimits = [minX,maxX,minY,maxY,minZ,maxZ];
            blockCounts = [NxMin,NxMax,NyMin,NyMax,NzMin,NzMax];
            coords = constructBlock_mex(coords0,Mm,blockOffset,blockLimits,blockCounts);
            coords = coords.';
            % size(coords)
            

        end

        if (~isempty(coords))
            % Now we need to apply the polygonal boundary constraint:
            inside = inpolygon(coords(:,2),coords(:,3),handles.grain{grainIndex}.boundX,handles.grain{grainIndex}.boundY);
            ind = find(inside == 1);
            % format of coords: handles.grain{grainIndex}.atomPos = [aType coords  DW occ charge];
            % format of handles.atomPos: handles.atomPos = [coords aType DW charge];
            if (~clearedAtomPos)
                handles.atomPos = [coords(ind,2:4) coords(ind,1) coords(ind,5:7)];
                clearedAtomPos = 1;
            else
                handles.atomPos = [handles.atomPos; coords(ind,2:4) coords(ind,1) coords(ind,5:7)];
            end
            size(handles.atomPos);
        end
    end
    if (useMex == 0)
        waitbar(grainIndex/handles.grainCount,wb,'Constructing super cell');
        if (grainIndex == handles.grainCount)
            close(wb);
        end
    end
end

axes(handles.axes1);
if isfield(handles,'atomPos')
    if ~isempty(handles.atomPos)
        guidata(hObject, handles);
        % Now that we have created the model, let's update the display of it:
        pushbutton_UpdateModel_Callback(hObject, eventdata, handles);
    end
end




% --- Executes on button press in pushbutton_UpdateGrains.
function pushbutton_UpdateGrains_Callback(hObject, eventdata, handles)
pushbutton_constructSuperCell(hObject, eventdata, handles)

% guidata(hObject, handles);





% --- Executes on button press in pushbutton_BoundaryAll.
function pushbutton_BoundaryAll_Callback(hObject, eventdata, handles)
boxX = str2double(get(handles.edit_BoxX,'String'));
boxY = str2double(get(handles.edit_BoxY,'String'));

boundX = [0 boxX boxX 0    0];
boundY = [0 0    boxY boxY 0];

grainIndex = str2double(get(handles.edit_GrainIndex,'String'));
handles.grain{grainIndex}.boundX = boundX;
handles.grain{grainIndex}.boundY = boundY;
guidata(hObject, handles);
pushbutton_constructSuperCell(hObject, eventdata, handles);


% --- Executes on button press in pushbutton_SelectPolygon.
function pushbutton_SelectPolygon_Callback(hObject, eventdata, handles)
button = 1;
count = 0;
boundX = [];
boundY = [];
while (button == 1)
    [x,y,button] = ginput(1);
    if (button == 1)
        boundX = [boundX;x];
        boundY = [boundY;y];
        if count > 0
            l = line(boundX(count:count+1),boundY(count:count+1));
            set(l,'LineWidth',2,'Color',[0 1 0]);
        end
        count = count+1;
    else
        if (count > 2)
            boundX = [boundX;boundX(1)];
            boundY = [boundY;boundY(1)];
            l = line(boundX(count:count+1),boundY(count:count+1));
            set(l,'LineWidth',2,'Color',[0 1 0]);
        else
           % if no boundary with more than 2 points has been defined, then
           % abort this procedure.
           boundX = [];
           boundY = [];
        end
	end
end
if ~isempty(boundX)
	if get(handles.checkbox_Sublattice,'Value')
		Zsel = str2num(get(handles.edit_Zsublattice,'String'));
		ind = find(handles.atomPos(:,4)==Zsel);
		ind2 = find(inpolygon(handles.atomPos(ind,1),handles.atomPos(ind,2),boundX,boundY));
		% fprintf('%d atoms with Z=%d\n',length(ind),Zsel);
		handles.indSel = [handles.indSel; ind(ind2)];
		% handles.indSel = [handles.indSel; find(inpolygon(handles.atomPos(:,1),handles.atomPos(:,2),boundX,boundY))];
	else
		handles.indSel = [handles.indSel; find(inpolygon(handles.atomPos(:,1),handles.atomPos(:,2),boundX,boundY))];
	end
end
guidata(hObject, handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);




% --- Executes on button press in pushbutton_Boundaries.
function pushbutton_Boundaries_Callback(hObject, eventdata, handles)

button = 1;
count = 0;
boundX = [];
boundY = [];
while (button == 1)
    [x,y,button] = ginput(1);
    if (button == 1)
        boundX = [boundX;x];
        boundY = [boundY;y];
        if count > 0
            l = line(boundX(count:count+1),boundY(count:count+1));
            set(l,'LineWidth',2,'Color',[1 0 0]);
        end
        count = count+1;
    else
        if (count > 2)
            boundX = [boundX;boundX(1)];
            boundY = [boundY;boundY(1)];
            l = line(boundX(count:count+1),boundY(count:count+1));
            set(l,'LineWidth',2,'Color',[1 0 0]);
        else
           % if no boundary with more than 2 points has been defined, then
           % abort this procedure.
           boundX = [];
           boundY = [];
        end
    end
end
% keep this boundary for later
if ~isempty(boundX)
    grainIndex = str2double(get(handles.edit_GrainIndex,'String'));
    handles.grain{grainIndex}.boundX = boundX;
    handles.grain{grainIndex}.boundY = boundY;
    guidata(hObject, handles);
    pushbutton_constructSuperCell(hObject, eventdata, handles);
end




% --- Executes on button press in pushbutton_LoadStructure.
function pushbutton_LoadStructure_Callback(hObject, eventdata, handles)

% [pathname, filename, ext] = fileparts(handles.posFileName);
% filename = [filename ext];
askUser = 1;
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
% get the current grain index:
grainIndex = str2double(get(handles.edit_GrainIndex,'String'));
% read the file
[coords, aType, Mm, DW, occ, charge] = readCFG_qstem(fullfile(pathname,filename),[0 0 0],2);

if ~isempty(coords)
    handles.posFileName = fullfile(pathname,filename);
    handles.grain{grainIndex}.cfgName = filename;
    handles.grain{grainIndex}.atomPos = [aType coords  DW occ charge];
    set(handles.edit_CFGName,'String',handles.grain{grainIndex}.cfgName);
    clear coords aType DW occ charge
    handles.grain{grainIndex}.Mm = Mm;
    if (handles.grainCount < grainIndex) handles.grainCount = grainIndex; end;
    guidata(hObject, handles);
    pushbutton_constructSuperCell(hObject, eventdata, handles)
end

% Set Zone Axis
pushbutton_Zone_Callback(hObject, eventdata, handles);




function edit_dx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dx as text
%        str2double(get(hObject,'String')) returns contents of edit_dx as a double


% --- Executes during object creation, after setting all properties.
function edit_dx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dy as text
%        str2double(get(hObject,'String')) returns contents of edit_dy as a double


% --- Executes during object creation, after setting all properties.
function edit_dy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [handles, needUpdate] = readParams(handles)
needUpdate = 1;
grainIndex = str2double(get(handles.edit_GrainIndex,'String'));
if (handles.grainCount < grainIndex) handles.grainCount = grainIndex; end

handles.grain{grainIndex}.TiltX  = str2num(get(handles.edit_TiltX,'String'));
handles.grain{grainIndex}.TiltY  = str2num(get(handles.edit_TiltY,'String'));
handles.grain{grainIndex}.TiltZ  = str2num(get(handles.edit_TiltZ,'String'));

handles.grain{grainIndex}.OffsetX  = str2num(get(handles.edit_OffsetX,'String'));
handles.grain{grainIndex}.OffsetY  = str2num(get(handles.edit_OffsetY,'String'));
handles.grain{grainIndex}.OffsetZ  = str2num(get(handles.edit_OffsetZ,'String'));

handles.grain{grainIndex}.minZ  = str2num(get(handles.edit_minZ,'String'));
handles.grain{grainIndex}.maxZ  = str2num(get(handles.edit_maxZ,'String'));




function edit_GrainIndex_Callback(hObject, eventdata, handles)
grainIndex = str2double(get(handles.edit_GrainIndex,'String'));

if (grainIndex > handles.grainCount) 
    handles.grain{grainIndex}.TiltX = 0;
    handles.grain{grainIndex}.TiltY = 0;
    handles.grain{grainIndex}.TiltZ = 0;
    handles.grain{grainIndex}.OffsetX = 0;
    handles.grain{grainIndex}.OffsetY = 0;
    handles.grain{grainIndex}.OffsetZ = 0;
    handles.grain{grainIndex}.minZ = str2num(get(handles.edit_minZ,'String'));
    handles.grain{grainIndex}.maxZ = str2num(get(handles.edit_maxZ,'String'));
    set(handles.edit_CFGName,'String',' ');
end
if isfield(handles.grain{grainIndex},'cfgName')
    set(handles.edit_CFGName,'String',handles.grain{grainIndex}.cfgName);
else
    set(handles.edit_CFGName,'String',' ');   
end

if isfield(handles.grain{grainIndex},'TiltX')
    set(handles.edit_TiltX,'String',handles.grain{grainIndex}.TiltX);
    set(handles.edit_TiltY,'String',handles.grain{grainIndex}.TiltY);
    set(handles.edit_TiltZ,'String',handles.grain{grainIndex}.TiltZ);
else
    set(handles.edit_TiltX,'String','0');
    set(handles.edit_TiltY,'String','0');
    set(handles.edit_TiltZ,'String','0');    
end

if isfield(handles.grain{grainIndex},'OffsetX')
    set(handles.edit_OffsetX,'String',handles.grain{grainIndex}.OffsetX);
    set(handles.edit_OffsetY,'String',handles.grain{grainIndex}.OffsetY);
    set(handles.edit_OffsetZ,'String',handles.grain{grainIndex}.OffsetZ);
else
    set(handles.edit_OffsetX,'String','0');
    set(handles.edit_OffsetY,'String','0');
    set(handles.edit_OffsetZ,'String','0');
end

if isfield(handles.grain{grainIndex},'minZ')
    set(handles.edit_minZ,'String',handles.grain{grainIndex}.minZ);
    set(handles.edit_maxZ,'String',handles.grain{grainIndex}.maxZ);
else
    set(handles.edit_minZ,'String','0');
    set(handles.edit_maxZ,'String','10');
end

if isfield(handles.grain{grainIndex},'zone')
    set(handles.info_h,'String',handles.grain{grainIndex}.zone(1));
    set(handles.info_k,'String',handles.grain{grainIndex}.zone(2));
    set(handles.info_l,'String',handles.grain{grainIndex}.zone(3));
else
    set(handles.info_h,'String','0');
    set(handles.info_k,'String','0');
    set(handles.info_l,'String','1');
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_GrainIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_GrainIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_CFGName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CFGName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CFGName as text
%        str2double(get(hObject,'String')) returns contents of edit_CFGName as a double


% --- Executes during object creation, after setting all properties.
function edit_CFGName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CFGName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_TiltX_Callback(hObject, eventdata, handles)
[handles, needUpdate] = readParams(handles);
guidata(hObject, handles);
if  (needUpdate)
    pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);
end



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
[handles, needUpdate] = readParams(handles);
guidata(hObject, handles);
if  (needUpdate)
    pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);
end


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
[handles, needUpdate] = readParams(handles);
guidata(hObject, handles);
if  (needUpdate)
    pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);
end

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



function edit_OffsetX_Callback(hObject, eventdata, handles)
[handles, needUpdate] = readParams(handles);
guidata(hObject, handles);
if  (needUpdate)
    pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);
end


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
[handles, needUpdate] = readParams(handles);
guidata(hObject, handles);
if  (needUpdate)
    pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);
end

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



function edit_OffsetZ_Callback(hObject, eventdata, handles)
[handles, needUpdate] = readParams(handles);
guidata(hObject, handles);
if  (needUpdate)
    pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);
end


% --- Executes during object creation, after setting all properties.
function edit_OffsetZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OffsetZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minZ_Callback(hObject, eventdata, handles)
[handles, needUpdate] = readParams(handles);
guidata(hObject, handles);
if  (needUpdate)
    pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_minZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_maxZ_Callback(hObject, eventdata, handles)
[handles, needUpdate] = readParams(handles);
guidata(hObject, handles);
if  (needUpdate)
    pushbutton_UpdateGrains_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function edit_maxZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ImageDx_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);
pushbutton_constructSuperCell(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_ImageDx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ImageDy_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);
pushbutton_constructSuperCell(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_ImageDy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ImageRotation_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);
pushbutton_constructSuperCell(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_ImageRotation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ImageOffsetX_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);
pushbutton_constructSuperCell(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_ImageOffsetX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ImageOffsetY_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);
pushbutton_constructSuperCell(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_ImageOffsetY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Zmax_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_Zmax_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AtomSize_Callback(hObject, eventdata, handles)
pushbutton_UpdateModel_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_AtomSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LimLow_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_LimLow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LimLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LimHigh_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_LimHigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LimHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_WarpNx_Callback(hObject, eventdata, handles)
handles.offsWp = [];
handles.atomPosUnWarped = [];
guidata(hObject,handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_WarpNx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_warpNy_Callback(hObject, eventdata, handles)
handles.offsWp = [];
handles.atomPosUnWarped = [];
guidata(hObject,handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_warpNy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_PointX_Callback(hObject, eventdata, handles)
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit_PointX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Distance_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_Distance_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function edit_BoxX_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_BoxX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_BoxX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_BoxY_Callback(hObject, eventdata, handles)
% Update plot:
pushbutton_Redraw_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_BoxY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_BoxY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton_Export.
function pushbutton_Export_Callback(hObject, eventdata, handles)

[name, pathname] = uiputfile('*.xyz', 'Export XYZ file');
if isequal(name,0) || isequal(pathname,0)
    return
else
    fileName = fullfile(pathname, name);
end


name = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl',...
    'Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
    'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
    'I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
    'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'];


Natom = size(handles.atomPos,1);

fid = fopen(fileName,'w');
fprintf(fid,'%d\n',Natom);
fprintf(fid,'exported from qstem model builder\n');

% writeCFG(fileName,Mm,handles.atomPos(:,4),handles.atomPos(:,1:3),[0 0 0],1,handles.atomPos(:,5),handles.atomPos(:,6));
for j=1:Natom
    at = handles.atomPos(j,4);
    fprintf(fid,'%s %.4f %.4f %.4f\n',name(2*at-1:2*at),handles.atomPos(j,1:3));  % specify type and mass of first atom (all others are the same)
end

fclose(fid);
fprintf('Exported atom coordinates to %s\n',fileName);




% --- Executes on button press in checkbox_showLines.
function checkbox_showLines_Callback(hObject, eventdata, handles)
pushbutton_UpdateModel_Callback(hObject, eventdata, handles)












function edit_PointY_Callback(hObject, eventdata, handles)
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_PointY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PointY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ShowWarpPoints.
function checkbox_ShowWarpPoints_Callback(hObject, eventdata, handles)
if get(handles.checkbox_ShowWarpPoints,'Value')
end
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);





% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
handles.atomPosUnWarped = [];
handles.offsWp = [];
guidata(hObject,handles);
pushbutton_UpdateModel_Callback(hObject, eventdata, handles);


% --- Executes on button press in checkbox_Sublattice.
function checkbox_Sublattice_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Sublattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Sublattice



function edit_Zsublattice_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Zsublattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Zsublattice as text
%        str2double(get(hObject,'String')) returns contents of edit_Zsublattice as a double


% --- Executes during object creation, after setting all properties.
function edit_Zsublattice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Zsublattice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function offsetEditAmount_Callback(hObject, eventdata, handles)
% hObject    handle to offsetEditAmount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offsetEditAmount as text
%        str2double(get(hObject,'String')) returns contents of offsetEditAmount as a double


% --- Executes during object creation, after setting all properties.
function offsetEditAmount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offsetEditAmount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in shift_up.
function shift_up_Callback(hObject, eventdata, handles)
% hObject    handle to shift_up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shiftAmount = str2double(get(handles.offsetEditAmount,'String'));
newOffset = str2double(get(handles.edit_OffsetY,'String')) + shiftAmount;
set(handles.edit_OffsetY,'String',num2str(newOffset));
edit_OffsetY_Callback(hObject, eventdata, handles);

% --- Executes on button press in shift_right.
function shift_right_Callback(hObject, eventdata, handles)
% hObject    handle to shift_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shiftAmount = str2double(get(handles.offsetEditAmount,'String'));
newOffset = str2double(get(handles.edit_OffsetX,'String')) + shiftAmount;
set(handles.edit_OffsetX,'String',num2str(newOffset));
edit_OffsetX_Callback(hObject, eventdata, handles);

% --- Executes on button press in shift_down.
function shift_down_Callback(hObject, eventdata, handles)
% hObject    handle to shift_down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shiftAmount = str2double(get(handles.offsetEditAmount,'String'));
newOffset = str2double(get(handles.edit_OffsetY,'String')) - shiftAmount;
set(handles.edit_OffsetY,'String',num2str(newOffset));
edit_OffsetY_Callback(hObject, eventdata, handles);

% --- Executes on button press in shift_left.
function shift_left_Callback(hObject, eventdata, handles)
% hObject    handle to shift_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
shiftAmount = str2double(get(handles.offsetEditAmount,'String'));
newOffset = str2double(get(handles.edit_OffsetX,'String')) - shiftAmount;
set(handles.edit_OffsetX,'String',num2str(newOffset));
edit_OffsetX_Callback(hObject, eventdata, handles);


function info_k_Callback(hObject, eventdata, handles)
% hObject    handle to info_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of info_k as text
%        str2double(get(hObject,'String')) returns contents of info_k as a double


% --- Executes during object creation, after setting all properties.
function info_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to info_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function info_h_Callback(hObject, eventdata, handles)
% hObject    handle to info_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of info_h as text
%        str2double(get(hObject,'String')) returns contents of info_h as a double


% --- Executes during object creation, after setting all properties.
function info_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to info_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function info_l_Callback(hObject, eventdata, handles)
% hObject    handle to info_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of info_l as text
%        str2double(get(hObject,'String')) returns contents of info_l as a double


% --- Executes during object creation, after setting all properties.
function info_l_CreateFcn(hObject, eventdata, handles)
% hObject    handle to info_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in rotate_clock.
function rotate_clock_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_clock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotateAmount = str2double(get(handles.offsetEditAmount,'String'));
newTiltZ = str2double(get(handles.edit_TiltZ,'String')) - rotateAmount;
set(handles.edit_TiltZ,'String',num2str(newTiltZ));
edit_TiltZ_Callback(hObject, eventdata, handles);


% --- Executes on button press in rotate_cclock.
function rotate_cclock_Callback(hObject, eventdata, handles)
% hObject    handle to rotate_cclock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotateAmount = str2double(get(handles.offsetEditAmount,'String'));
newTiltZ = str2double(get(handles.edit_TiltZ,'String')) + rotateAmount;
set(handles.edit_TiltZ,'String',num2str(newTiltZ));
edit_TiltZ_Callback(hObject, eventdata, handles);
