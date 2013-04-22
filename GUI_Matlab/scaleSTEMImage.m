function varargout = scaleSTEMImage(varargin)
% SCALESTEMIMAGE M-file for scaleSTEMImage.fig
%      SCALESTEMIMAGE, by itself, creates a new SCALESTEMIMAGE or raises the existing
%      singleton*.
%
%      H = SCALESTEMIMAGE returns the handle to a new SCALESTEMIMAGE or the handle to
%      the existing singleton*.
%
%      SCALESTEMIMAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCALESTEMIMAGE.M with the given input arguments.
%
%      SCALESTEMIMAGE('Property','Value',...) creates a new SCALESTEMIMAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before scaleSTEMImage_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to scaleSTEMImage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help scaleSTEMImage

% Last Modified by GUIDE v2.5 29-Jun-2010 10:37:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @scaleSTEMImage_OpeningFcn, ...
                   'gui_OutputFcn',  @scaleSTEMImage_OutputFcn, ...
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


% --- Executes just before scaleSTEMImage is made visible.
function scaleSTEMImage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to scaleSTEMImage (see VARARGIN)

% Choose default command line output for scaleSTEMImage
handles.output = hObject;
handles.h = varargin{1};
set(handles.edit_SourceSize,'String',handles.h.sourceSize);
handles = readAllFields(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes scaleSTEMImage wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = scaleSTEMImage_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function handles = readAllFields(handles)
brightness = str2num(get(handles.edit_Brightness,'String'))*1e8;
dwellTime  = str2num(get(handles.edit_DwellTime,'String'))*1e-6;  % µs -> sec
sourceSize = str2num(get(handles.edit_SourceSize,'String'));
alpha      = str2num(get(handles.edit_alpha,'String'))*1e-3;
lambda     = wavelength(str2num(get(handles.edit_HighVoltage,'String')));
noiseFlag  = get(handles.checkbox_AddNoise,'Value');

% Re-compute image based on these parameters:
scale = brightness*(1e-8*pi*alpha*sqrt((lambda/(pi*alpha))^2+(0.5*sourceSize)^2))^2;
% unit so far: A/cm^2sr*cm^2*sr = C/sec
% Need to divide by e and multiply by dwell time to get number of
% electrons:
scale = scale*dwellTime/1.602176487e-19;
handles.scale = scale;
handles.noiseFlag = noiseFlag;
set(handles.text_scale,'String',sprintf('Incident electrons per pixel: %.1f',scale));


% --- Executes on button press in pushbutton_Display.
function pushbutton_Display_Callback(hObject, eventdata, handles)
handles = readAllFields(handles);
[Ny,Nx] = size(handles.h.imgFinal);
figure;
img = handles.scale*handles.h.imgFinal;
if handles.noiseFlag
    imagesc(handles.h.dx*[0:Nx],handles.h.dy*[0:Ny],img+sqrt(img).*randn(size(img)));    
else
    imagesc(handles.h.dx*[0:Nx],handles.h.dy*[0:Ny],img);
end
colormap('gray');
colorbar

% --- Executes on button press in pushbutton_Save.
function pushbutton_Save_Callback(hObject, eventdata, handles)
handles = readAllFields(handles);
[Ny,Nx] = size(handles.h.imgFinal);

img = handles.scale*handles.h.imgFinal;
if handles.noiseFlag
    img = img+sqrt(img).*randn(size(img));    
end
[filename, pathname] = uiputfile('*.img', 'Save STEM image in .img format');
if isequal(filename,0) || isequal(pathname,0)
else
    binwrite2D(img,fullfile(pathname, filename),handles.h.dx,handles.h.dy,handles.h.t);
end


function edit_Brightness_Callback(hObject, eventdata, handles)
handles = readAllFields(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_Brightness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_DwellTime_Callback(hObject, eventdata, handles)
handles = readAllFields(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_DwellTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DwellTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alpha_Callback(hObject, eventdata, handles)
handles = readAllFields(handles);
% Update handles structure
guidata(hObject, handles);


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



function edit_SourceSize_Callback(hObject, eventdata, handles)
handles = readAllFields(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_SourceSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SourceSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_AddNoise.
function checkbox_AddNoise_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_AddNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_AddNoise





function edit_HighVoltage_Callback(hObject, eventdata, handles)
handles = readAllFields(handles);
% Update handles structure
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


