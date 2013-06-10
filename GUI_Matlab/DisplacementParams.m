function varargout = DisplacementParams(varargin)
% DISPLACEMENTPARAMS MATLAB code for DisplacementParams.fig
%      DISPLACEMENTPARAMS, by itself, creates a new DISPLACEMENTPARAMS or raises the existing
%      singleton*.
%
%      H = DISPLACEMENTPARAMS returns the handle to a new DISPLACEMENTPARAMS or the handle to
%      the existing singleton*.
%
%      DISPLACEMENTPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLACEMENTPARAMS.M with the given input arguments.
%
%      DISPLACEMENTPARAMS('Property','Value',...) creates a new DISPLACEMENTPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DisplacementParams_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DisplacementParams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DisplacementParams

% Last Modified by GUIDE v2.5 07-Sep-2012 17:29:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DisplacementParams_OpeningFcn, ...
                   'gui_OutputFcn',  @DisplacementParams_OutputFcn, ...
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


% --- Executes just before DisplacementParams is made visible.
function DisplacementParams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DisplacementParams (see VARARGIN)

% Choose default command line output for DisplacementParams
handles.output = hObject;
if size(varargin,2) > 0
    if isstruct(varargin{1})
        if isfield(varargin{1},'a')
            handles.a = varargin{1}.a;  % [coords aType DW occ charge];
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% First, we should remove replicates in Z-direction
			% sort atoms by x
			% columns in a: Z x y z DW occ charge
			a = sortrows([zeros(1,size(handles.a,2)); handles.a],[1 2 3 4]);
			% difference in:       Z       or        x         or      y
			ind = find(abs(diff(a(:,1)))+abs(diff(a(:,2)))+abs(diff(a(:,3))) > 0);
			handles.a = a(ind,:);
		end
    end
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DisplacementParams wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DisplacementParams_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_Zreference_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Zreference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Zreference as text
%        str2double(get(hObject,'String')) returns contents of edit_Zreference as a double


% --- Executes during object creation, after setting all properties.
function edit_Zreference_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Zreference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ztarget_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Ztarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Ztarget as text
%        str2double(get(hObject,'String')) returns contents of edit_Ztarget as a double


% --- Executes during object creation, after setting all properties.
function edit_Ztarget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ztarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Radius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Radius as text
%        str2double(get(hObject,'String')) returns contents of edit_Radius as a double


% --- Executes during object creation, after setting all properties.
function edit_Radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Compute.
function pushbutton_Compute_Callback(hObject, eventdata, handles)
Zref    = str2num(get(handles.edit_Zreference,'String'));
Ztarget = str2num(get(handles.edit_Ztarget,'String'));
radius  = str2num(get(handles.edit_Radius,'String'));
Nneighbor  = str2num(get(handles.edit_Nneighbors,'String'));
radius2 = radius^2;
% Warn, if reference or target Z = 0
if (Zref == 0) || (Ztarget == 0)
	errordlg('Please make sure the atomic numbers of the desired species are properly defined.', 'Create displacement map');
	return;
end

% plot(handles.a(:,2),handles.a(:,3),'o');
% find the atoms that correspond to our search criteria
indTarget = find(handles.a(:,1) == Ztarget);
xyTarget=handles.a(indTarget,[2 3]);
indRef    = find(handles.a(:,1) == Zref);
xyRef=handles.a(indRef,[2 3]);

Ntarget = length(indTarget);
% Warn, if reference or target Z wrong
if (Ntarget <= 0)
	errordlg(sprintf('Could not find atoms with Z=%d',Ztarget), 'Create displacement map');
	return;
end
if (isempty(indRef))
	errordlg(sprintf('Could not find reference atoms with Z=%d',Zref), 'Create displacement map');
	return;
end

displacementX = zeros(Ntarget,1);
displacementY = zeros(Ntarget,1);
for j=1:Ntarget
	ind = find((xyRef(:,1)-xyTarget(j,1)).^2+(xyRef(:,2)-xyTarget(j,2)).^2 < radius2);
	% sort the neighbors by distance and take only the closest ones:
	if length(ind) > Nneighbor
		[dis,inds] = sort(((xyRef(ind,1)-xyTarget(j,1)).^2+(xyRef(ind,2)-xyTarget(j,2)).^2),1);
		ind = ind(inds(1:Nneighbor));
	end
	% center = [mean(xyRef(ind,1)),mean(xyRef(ind,2))]; 
	displacementX(j) = xyTarget(j,1)-mean(xyRef(ind,1));
	displacementY(j) = xyTarget(j,2)-mean(xyRef(ind,2));
end
ind = find(sqrt(displacementX.^2+displacementY.^2) > 1.2);
displacementX(ind) = 0;
displacementY(ind) = 0;
figure
quiver(xyTarget(:,1),xyTarget(:,2),displacementX,displacementY);


% --- Executes on button press in pushbutton_Export.
function pushbutton_Export_Callback(hObject, eventdata, handles)
Zref    = str2num(get(handles.edit_Zreference,'String'));
Ztarget = str2num(get(handles.edit_Ztarget,'String'));
radius  = str2num(get(handles.edit_Radius,'String'));
Nneighbor  = str2num(get(handles.edit_Nneighbors,'String'));
radius2 = radius^2;
% Warn, if reference or target Z = 0
if (Zref == 0) || (Ztarget == 0)
	errordlg('Please make sure the atomic numbers of the desired species are properly defined.', 'Create displacement map');
	return;
end

% find the atoms that correspond to our search criteria
indTarget = find(handles.a(:,1) == Ztarget);
xyTarget=handles.a(indTarget,[2 3]);
indRef    = find(handles.a(:,1) == Zref);
xyRef=handles.a(indRef,[2 3]);

Ntarget = length(indTarget);
% Warn, if reference or target Z wrong
if (Ntarget <= 0)
	errordlg(sprintf('Could not find atoms with Z=%d',Ztarget), 'Create displacement map');
	return;
end
if (isempty(indRef))
	errordlg(sprintf('Could not find reference atoms with Z=%d',Zref), 'Create displacement map');
	return;
end

displacementX = zeros(Ntarget,1);
displacementY = zeros(Ntarget,1);
for j=1:Ntarget
	ind = find((xyRef(:,1)-xyTarget(j,1)).^2+(xyRef(:,2)-xyTarget(j,2)).^2 < radius2);
	% center = [mean(xyRef(ind,1)),mean(xyRef(ind,2))]; 
	if length(ind) > Nneighbor
		[dis,inds] = sort(((xyRef(ind,1)-xyTarget(j,1)).^2+(xyRef(ind,2)-xyTarget(j,2)).^2),1);
		ind = ind(inds(1:Nneighbor));
	end

	displacementX(j) = xyTarget(j,1)-mean(xyRef(ind,1));
	displacementY(j) = xyTarget(j,2)-mean(xyRef(ind,2));
end
ind = find(sqrt(displacementX.^2+displacementY.^2) > 1.2);
displacementX(ind) = 0;
displacementY(ind) = 0;


[filename, pathname] = uiputfile({'*.*'},'Save Vedtor data');
fn = fullfile(pathname,filename);

M = [xyTarget,displacementX,displacementY,sqrt(displacementX.^2+displacementY.^2),180/pi*atan2(displacementY,displacementX)];
if (0)
	save(fn,'M','-ascii','-tabs');
	fprintf('Saved vector plots in the following format:\n');
	fprintf('posX  posY  vectX  vectY  vectLength vectAngle\n');
	fprintf('All length units in A, angle in deg.\n');
else
	fid = fopen(fn,'w');
	fprintf(fid,'posX  \tposY  \tvectX \tvectY \tlength \tangle\n');
	for j=1:Ntarget
		fprintf(fid,'%.3f \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f\n',M(j,:));
	end
	fclose(fid);
	fprintf('Saved vector plots in the following format:\n');
	fprintf('posX  posY  vectX  vectY  vectLength vectAngle\n');
	fprintf('All length units in A, angle in deg.\n');
end



function edit_Nneighbors_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_Nneighbors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Nneighbors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
