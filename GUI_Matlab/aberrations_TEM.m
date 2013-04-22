function varargout = aberrations_TEM(varargin)
% ABERRATIONS_TEM M-file for aberrations_TEM.fig
%      ABERRATIONS_TEM by itself, creates a new ABERRATIONS_TEM or raises
%      the
%      existing singleton*.
%
%      H = ABERRATIONS_TEM returns the handle to a new ABERRATIONS_TEM or the handle to
%      the existing singleton*.
%
%      ABERRATIONS_TEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABERRATIONS_TEM.M with the given input arguments.
%
%      ABERRATIONS_TEM('Property','Value',...) creates a new ABERRATIONS_TEM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aberrations_TEM_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aberrations_TEM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aberrations_TEM

% Last Modified by GUIDE v2.5 16-Feb-2012 15:11:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aberrations_TEM_OpeningFcn, ...
                   'gui_OutputFcn',  @aberrations_TEM_OutputFcn, ...
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

function s = signpm(v)
s = 2*(v > 0)-1;

% --- Executes just before aberrations_TEM is made visible.
function aberrations_TEM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aberrations_TEM (see VARARGIN)

% Choose default command line output for aberrations_TEM
handles.oldNomenclature = 1;
handles.oldDeg = 1;
set(handles.radiobutton_Nomenclature1,'Value',1);
set(handles.radiobutton_deg,'Value',1);
handles.parentHandles = varargin{1};
handles.a      = handles.parentHandles.a;
handles.phi    = handles.parentHandles.phi;
handles.c      = handles.parentHandles.c;
handles.waveRec= varargin{2};
[Ny,Nx] = size(handles.waveRec);
handles.params = [Nx,Ny,handles.parentHandles.dx,handles.parentHandles.dy,handles.parentHandles.highVoltage,1e3*handles.parentHandles.Convergence];

set(handles.edit_Alpha,'String',1e3*handles.parentHandles.Convergence);
set(handles.edit_Delta,'String',0.1*handles.parentHandles.Delta);
if (0)
    set(handles.edit_C1,'String',sprintf('%.3f',1e-3*handles.c(2)));
    set(handles.edit_C3,'String',sprintf('%.3f',1e-3*handles.c(4)));
    set(handles.edit_A1,'String',sprintf('%.3f',1e-3*handles.a(2,2)));
    set(handles.edit_Phi22,'String',sprintf('%.3f',1e-3*handles.phi(2,2)));
end
guidata(hObject, handles);
writeParams(handles);
handles.pressedCancel = 0;


% Update handles structure
guidata(hObject, handles);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.
if(nargin > 3)
end

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Show a question icon from dialogicons.mat - variables questIconData
% and questIconMap

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal');
pushbutton_Test_Callback(hObject, eventdata, handles);

% UIWAIT makes aberrations_TEM wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = aberrations_TEM_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles = readParams(handles);
handles.params(6) = str2num(get(handles.edit_Alpha,'String'));
handles.params(7) = 10*str2num(get(handles.edit_Delta,'String'));


varargout{1} = handles.a;
varargout{2} = handles.phi;
varargout{3} = handles.c;

if handles.pressedCancel;
    varargout{4} = 0;   
else
    varargout{4} = handles.params;
end
% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in pushbutton_Ok.
function pushbutton_Ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);

% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Update handles structure
handles.pressedCancel = 1;
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end


% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.figure1);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.figure1);
end    



function edit_A1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_A1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_A1 as text
%        str2double(get(hObject,'String')) returns contents of edit_A1 as a double


% --- Executes during object creation, after setting all properties.
function edit_A1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_A1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_A2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_A2 as text
%        str2double(get(hObject,'String')) returns contents of edit_A2 as a double


% --- Executes during object creation, after setting all properties.
function edit_A2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_A2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_A3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_A3 as text
%        str2double(get(hObject,'String')) returns contents of edit_A3 as a double


% --- Executes during object creation, after setting all properties.
function edit_A3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_A3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_A4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_A4 as text
%        str2double(get(hObject,'String')) returns contents of edit_A4 as a double


% --- Executes during object creation, after setting all properties.
function edit_A4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_A4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_A5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_A5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_A5 as text
%        str2double(get(hObject,'String')) returns contents of edit_A5 as a double


% --- Executes during object creation, after setting all properties.
function edit_A5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_A5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi22_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi22 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi22 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi33_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi33 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi33 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi44_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi44 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi44 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi55_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi55 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi55 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi66_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi66 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi66 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi66_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi66 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B2 as text
%        str2double(get(hObject,'String')) returns contents of edit_B2 as a double


% --- Executes during object creation, after setting all properties.
function edit_B2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_B2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi31_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi31 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi31 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_S3 as text
%        str2double(get(hObject,'String')) returns contents of edit_S3 as a double


% --- Executes during object creation, after setting all properties.
function edit_S3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_S3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi42_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi42 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi42 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_D4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_D4 as text
%        str2double(get(hObject,'String')) returns contents of edit_D4 as a double


% --- Executes during object creation, after setting all properties.
function edit_D4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_D4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi53_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi53 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi53 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi53_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B4 as text
%        str2double(get(hObject,'String')) returns contents of edit_B4 as a double


% --- Executes during object creation, after setting all properties.
function edit_B4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_B4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi51_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi51 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi51 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_D5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_D5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_D5 as text
%        str2double(get(hObject,'String')) returns contents of edit_D5 as a double


% --- Executes during object creation, after setting all properties.
function edit_D5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_D5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi64_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi64 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi64 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi64_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi64 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_B5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_B5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_B5 as text
%        str2double(get(hObject,'String')) returns contents of edit_B5 as a double


% --- Executes during object creation, after setting all properties.
function edit_B5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_B5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Phi62_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Phi62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Phi62 as text
%        str2double(get(hObject,'String')) returns contents of edit_Phi62 as a double


% --- Executes during object creation, after setting all properties.
function edit_Phi62_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Phi62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_C1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C1 as text
%        str2double(get(hObject,'String')) returns contents of edit_C1 as a double


% --- Executes during object creation, after setting all properties.
function edit_C1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C1 (see GCBO)
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



function edit_C5_Callback(hObject, eventdata, handles)
% hObject    handle to edit_C5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_C5 as text
%        str2double(get(hObject,'String')) returns contents of edit_C5 as a double


% --- Executes during object creation, after setting all properties.
function edit_C5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_C5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in radiobutton_Nomenclature1.
function radiobutton_Nomenclature1_Callback(hObject, eventdata, handles)
handles = readParams(handles);
handles.oldNomenclature = 1; 

set(handles.text_A1,'String','|A_1| (nm):');
set(handles.text_Phi22,'String','Phi_2,2:');
set(handles.text_C1,'String','C1 (nm):');

set(handles.text_A2,'String','|A_2| (nm):');
set(handles.text_Phi33,'String','Phi_3,3:');
set(handles.text_B2,'String','|B_2| (nm):');
set(handles.text_Phi31,'String','Phi_3,1:');

set(handles.text_A3,'String','|A_3| (um):');
set(handles.text_Phi44,'String','Phi_4,4:');
set(handles.text_S3,'String','|S_3| (um):');
set(handles.text_Phi42,'String','Phi_4,2:');
set(handles.text_C3,'String','C3 (um):');

set(handles.text_A4,'String','|A_4| (um):');
set(handles.text_Phi55,'String','Phi_5,5:');
set(handles.text_D4,'String','|D_4| (um):');
set(handles.text_Phi53,'String','Phi_5,3:');
set(handles.text_B4,'String','|B_4| (um):');
set(handles.text_Phi51,'String','Phi_5,1:');

set(handles.text_A5,'String','|A_5| (mm):');
set(handles.text_Phi66,'String','Phi_6,6:');
set(handles.text_D5,'String','|R_5| (mm):');
set(handles.text_Phi64,'String','Phi_6,4:');
set(handles.text_B5,'String','|S_5| (mm):');
set(handles.text_Phi62,'String','Phi_6,2:');
set(handles.text_C5,'String','C_5 (mm):');

set(handles.radiobutton_deg,'Enable','on');
set(handles.radiobutton_rad,'Enable','on');
writeParams(handles);
guidata(hObject, handles);

% --- Executes on button press in radiobutton_Nomenclature3.
function radiobutton_Nomenclature3_Callback(hObject, eventdata, handles)
handles = readParams(handles);
handles.oldNomenclature = 3; 

set(handles.text_A1,'String','a_2,2 (nm):');
set(handles.text_Phi22,'String','Phi_2,2:');
set(handles.text_C1,'String','a_2,0 (nm):');

set(handles.text_A2,'String','a_3,3 (nm):');
set(handles.text_Phi33,'String','Phi_3,3:');
set(handles.text_B2,'String','a_3,1 (nm):');
set(handles.text_Phi31,'String','Phi_3,1:');

set(handles.text_A3,'String','a_4,4 (um):');
set(handles.text_Phi44,'String','Phi_4,4:');
set(handles.text_S3,'String','a_4,2 (um):');
set(handles.text_Phi42,'String','Phi_4,2:');
set(handles.text_C3,'String','a_4,0 (um):');

set(handles.text_A4,'String','a_5,5 (um):');
set(handles.text_Phi55,'String','Phi_5,5:');
set(handles.text_D4,'String','a_5,3 (um):');
set(handles.text_Phi53,'String','Phi_5,3:');
set(handles.text_B4,'String','a_5,1 (um):');
set(handles.text_Phi51,'String','Phi_5,1:');

set(handles.text_A5,'String','a_6,6 (mm):');
set(handles.text_Phi66,'String','Phi_6,6:');
set(handles.text_D5,'String','a_6,4 (mm):');
set(handles.text_Phi64,'String','Phi_6,4:');
set(handles.text_B5,'String','a_6,2 (mm):');
set(handles.text_Phi62,'String','Phi_6,2:');
set(handles.text_C5,'String','a_6,0 (mm):');

set(handles.radiobutton_deg,'Enable','on');
set(handles.radiobutton_rad,'Enable','on');
writeParams(handles);
guidata(hObject, handles);



% --- Executes on button press in radiobutton_Nomenclature2.
function radiobutton_Nomenclature2_Callback(hObject, eventdata, handles)
handles = readParams(handles);
handles.oldNomenclature = 2; 

set(handles.text_A1,'String','C1,2a (nm):');
set(handles.text_Phi22,'String','C1,2b (nm):');
set(handles.text_C1,'String','C1 (nm):');

set(handles.text_A2,'String','C2,3a (nm):');
set(handles.text_Phi33,'String','C2,3b (nm):');
set(handles.text_B2,'String','C2,1a (nm):');
set(handles.text_Phi31,'String','C2,1b (nm):');

set(handles.text_A3,'String','C3,4a (um):');
set(handles.text_Phi44,'String','C3,4b (um):');
set(handles.text_S3,'String','C3,2a (um):');
set(handles.text_Phi42,'String','C3,2b (um):');
set(handles.text_C3,'String','C3 (um):');

set(handles.text_A4,'String','C4,5a (um):');
set(handles.text_Phi55,'String','C4,5b (um):');
set(handles.text_D4,'String','C4,3a (um):');
set(handles.text_Phi53,'String','C4,3b (um):');
set(handles.text_B4,'String','C4,1a (um):');
set(handles.text_Phi51,'String','C4,1b (um):');

set(handles.text_A5,'String','C5,6a (mm):');
set(handles.text_Phi66,'String','C5,6b (mm):');
set(handles.text_D5,'String','C5,4a (mm):');
set(handles.text_Phi64,'String','C5,4b (mm):');
set(handles.text_B5,'String','C5,2a (mm):');
set(handles.text_Phi62,'String','C5,2b (mm):');
set(handles.text_C5,'String','C_5 (mm):');

set(handles.radiobutton_deg,'Enable','off');
set(handles.radiobutton_rad,'Enable','off');
writeParams(handles);
guidata(hObject, handles);



% --- Executes on button press in pushbutton_Test.
function pushbutton_Test_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'createNewFigure')
    handles.createNewFigure = 0;
end


Nx = handles.params(1);    % number of pixels in x-direction   
Ny = handles.params(2);    % number of pixels in y-direction   
dx = handles.params(3);    % sampling in x-direction (unit: A) 
dy = handles.params(4);    % sampling in x-direction (unit: A) 
dkx = 1/(Nx*dx);           % rec. space sampling in x-direction
dky = 1/(Ny*dy);           % rec. space sampling in x-direction
handles.params(5) = str2num(get(handles.edit_HighVoltage,'String'));
lambda = wavelength(handles.params(5));
Delta = 10*str2double(get(handles.edit_Delta,'String'));
alpha = 1e-3*str2num(get(handles.edit_Alpha,'String'));
qmax  = sin(alpha)/lambda; 

% Beam tilt: 
tiltKx = 0;
tiltKy = 0;
[kx,ky] = meshgrid(dkx*(-Nx/2+[0:Nx-1])+tiltKx,dky*(-Ny/2+[0:Ny-1])+tiltKy);
% kx2 = kx.^2;
% ky2 = ky.^2;
% k2 = kx2+ky2;
k2 = (kx.^2+ky.^2);
theta2 = k2*lambda*lambda;
ktheta = sqrt(theta2);
kphi = atan2(ky,kx);
ktm = asin(qmax*lambda);
% probe = handles.waveRec;
% probe(find(ktheta > handles.parentHandles.kObjApert)) = 0;

handles = readParams(handles);
guidata(hObject, handles);
a = handles.a;          % Array of amplitudes of aberrations_TEM (according to nomenclature 3)  
phi = handles.phi;      % Array of angles of aberrations_TEM   (according to nomenclature 3)    
c = handles.c;          % List of symmetric aberrations_TEM  (c(2) = def, c(4) = Cs, c(6) = C5) 




chi   = 2*pi/lambda*(1/2*(a(2,2).*cos(2*(kphi-phi(2,2)))+c(2)).*theta2+...
                     1/3*(a(3,3).*cos(3*(kphi-phi(3,3)))+a(3,1).*cos(1*(kphi-phi(3,1)))).*ktheta.^3+...
                     1/4*(a(4,4).*cos(4*(kphi-phi(4,4)))+a(4,2).*cos(2*(kphi-phi(4,2)))+c(4)).*theta2.^2+...
                     1/5*(a(5,5).*cos(5*(kphi-phi(5,5)))+a(5,3).*cos(3*(kphi-phi(5,3)))+a(5,1).*cos(1*(kphi-phi(5,1)))).*ktheta.^5+...
                     1/6*(a(6,6).*cos(6*(kphi-phi(6,6)))+a(6,4).*cos(4*(kphi-phi(6,4)))+a(6,2).*cos(2*(kphi-phi(6,2)))+c(6)).*theta2.^3);
% compute the probe
ctf     = exp(-i*chi);
% figure; subplot(1,2,1); imagesc(theta2); subplot(1,2,2); imagesc(chi); colorbar

if (handles.parentHandles.ObjApertureTheta > 0)
    handles.parentHandles.ObjApertureEdge = 1e-3*str2num(get(handles.parentHandles.edit_ObjApertureEdge,'String'));
    ctf(find(ktheta >= (handles.parentHandles.ObjApertureTheta+handles.parentHandles.ObjApertureEdge))) = 0;
    ind = find((ktheta < (handles.parentHandles.ObjApertureTheta+handles.parentHandles.ObjApertureEdge)) & (ktheta > handles.parentHandles.ObjApertureTheta));
    ctf(ind) = 0.5*ctf(ind).*(1+cos(pi*(ktheta(ind)-handles.parentHandles.ObjApertureTheta)/handles.parentHandles.ObjApertureEdge));
    clear ind
end

% Compute the damping envelopes:
% This is what I used in simImageFromWave - check this!!!:
% Etemp = exp(-2.0*(1.0/handles.lambda*handles.Delta*theta2).^2);  

pl      = pi*lambda;
% fa      = a(2,2).*cos(2*(kphi-phi(2,2)))+c(2);
%Espat = exp(-sign(alpha)*(pi*alpha*(fa+lambda^2*c(4).*k2)).^2.*k2);
% Espat = exp(-(alpha/2lambda)^2*[dchi(q)/dq]^2)
% I will use: [dchi(q)/dq]^2 = [dchi/d|q|]^2+[dchi/|q|dphi]^2
% first let's compute the q-derivative (the division by lambda is removed in both cases):
dchi_dq = ((a(2,2).*cos(2*(kphi-phi(2,2)))+c(2)).*ktheta+...
           (a(3,3).*cos(3*(kphi-phi(3,3)))+a(3,1).*cos(1*(kphi-phi(3,1)))).*theta2+...
           (a(4,4).*cos(4*(kphi-phi(4,4)))+a(4,2).*cos(2*(kphi-phi(4,2)))+c(4)).*ktheta.^3+...
           (a(5,5).*cos(5*(kphi-phi(5,5)))+a(5,3).*cos(3*(kphi-phi(5,3)))+a(5,1).*cos(1*(kphi-phi(5,1)))).*theta2.^2+...
           (a(6,6).*cos(6*(kphi-phi(6,6)))+a(6,4).*cos(4*(kphi-phi(6,4)))+a(6,2).*cos(2*(kphi-phi(6,2)))+c(6)).*ktheta.^5)*2*pi;
dchi_dphi=(1/2*(2*a(2,2).*sin(2*(kphi-phi(2,2)))).*ktheta+...
           1/3*(3*a(3,3).*sin(3*(kphi-phi(3,3)))+1*a(3,1).*sin(1*(kphi-phi(3,1)))).*theta2+...
           1/4*(4*a(4,4).*sin(4*(kphi-phi(4,4)))+2*a(4,2).*sin(2*(kphi-phi(4,2)))).*ktheta.^3+...
           1/5*(5*a(5,5).*sin(5*(kphi-phi(5,5)))+3*a(5,3).*sin(3*(kphi-phi(5,3)))+1*a(5,1).*sin(1*(kphi-phi(5,1)))).*theta2.^2+...
           1/6*(6*a(6,6).*sin(6*(kphi-phi(6,6)))+4*a(6,4).*sin(4*(kphi-phi(6,4)))+2*a(6,2).*sin(2*(kphi-phi(6,2)))).*ktheta.^5)*(-2*pi);
Espat = exp(-sign(alpha)*(alpha/(2*lambda)).^2.*(dchi_dq.^2+dchi_dphi.^2));       
Etemp = exp(-sign(Delta)*(0.5*pi/lambda*Delta*theta2).^2);
Etemp(find(Etemp > 5)) = 5;
Espat(find(Espat > 5)) = 5;

% figure; plot(handles.dkx*10*[-Nx/2:Nx/2-1],[fftshift(Espat(:,1)) fftshift(Etemp(:,1))]); pause(0.1);
% Apply the damping envelopes as well:
ctf = Etemp.*Espat.*ctf;

if (handles.parentHandles.TiltX ~= 0) || (handles.parentHandles.TiltY ~= 0)
    tx = round(handles.parentHandles.TiltX/(handles.parentHandles.lambda*handles.parentHandles.dkx));
    ty = round(handles.parentHandles.TiltY/(handles.parentHandles.lambda*handles.parentHandles.dky));
    probe = abs(ifft2(circshift(handles.waveRec,[ty tx]).*ifftshift(ctf))).^2;
else
    tx = 0;
    ty = 0;
    probe = abs(ifft2(handles.waveRec.*ifftshift(ctf))).^2;
end

if get(handles.parentHandles.checkbox_Vibration,'Value')
    % convert vibration from pm to A:
    vibration_Axis1 = 0.1*str2double(get(handles.parentHandles.edit_Vibration_Axis1,'String'));
    vibration_Axis2 = 0.1*str2double(get(handles.parentHandles.edit_Vibration_Axis2,'String'));
    vibration_Angle = pi/180*str2double(get(handles.parentHandles.edit_Vibration_Angle,'String'));
    vib1 = (pi*vibration_Axis1/(-log(0.5))).^2;
    vib2 = (pi*vibration_Axis2/(-log(0.5))).^2;
    probe = real(ifft2(fft2(probe).*ifftshift(exp(-vib1*k2.*abs(cos(kphi-vibration_Angle))-vib2*k2.*abs(sin(kphi-vibration_Angle))))));
end


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
    if (handles.createNewFigure == 0)
        axes(handles.axes_Probe);
        if get(handles.checkbox_Surf,'Value')
            surf(dx*0.1*([0:Nx-1]),dy*0.1*([0:Ny-1]),probe);
            shading interp;
        else
            imagesc(dx*0.1*([0:Nx-1]),dy*0.1*([0:Ny-1]),probe);
            set(gca,'YDir','normal');
            axis equal; axis tight;
        end
        if get(handles.checkbox_Color,'Value')
            colormap('default');
        else
            colormap('gray');
        end
        xlabel('x in nm');
        ylabel('y in nm');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display CTF and or phase map    
    if (handles.createNewFigure == 0)
        axes(handles.axes_Phasemap);
    else
        figure;
    end
    limLow = str2num(get(handles.edit_CTF_low,'String'));
    limHigh = str2num(get(handles.edit_CTF_high,'String'));
    if get(handles.radiobutton_realPart,'Value')
        imagesc(dkx*10*(-Nx/2+[0:Nx-1]),dky*10*(-Ny/2+[0:Ny-1]),real(ctf));
        title('Contrast Transfer Function (real part)')
    end
    if get(handles.radiobutton_imagPart,'Value')
        imagesc(dkx*10*(-Nx/2+[0:Nx-1]),dky*10*(-Ny/2+[0:Ny-1]),imag(ctf));
        title('Contrast Transfer Function (imaginary part)')
    end
    if get(handles.radiobutton_phasePart,'Value')
        imagesc(dkx*10*(-Nx/2+[0:Nx-1]),dky*10*(-Ny/2+[0:Ny-1]),angle(ctf));
        title('Contrast Transfer Function (phase)')
    end
    if get(handles.radiobutton_amplitudePart,'Value')
        imagesc(dkx*10*(-Nx/2+[0:Nx-1]),dky*10*(-Ny/2+[0:Ny-1]),abs(ctf));
        title('Contrast Transfer Function (amplitude)')
    end
    if get(handles.radiobutton_chi,'Value')
        imagesc(dkx*10*(-Nx/2+[0:Nx-1]),dky*10*(-Ny/2+[0:Ny-1]),chi);
        title('chi(k)')
        caxis([limLow limHigh])
    end
    set(gca,'YDir','normal');
    xlabel('kx in 1/nm');
    ylabel('ky in 1/nm');    

    if get(handles.checkbox_Color,'Value')
        colormap('default');
    else
        colormap('gray');
    end
    axis equal; axis tight;
    handles.createNewFigure = 0;
    guidata(hObject, handles);
end


function handles = readParams(handles)
a = zeros(6,6);    
phi = zeros(6,6);
c   = zeros(6,1);
% Spherical aberration:
% a(2,0) = 10*str2num(get(handles.edit_C1,'String'));
c(2)     = 10*str2num(get(handles.edit_C1,'String'));
% Astigmatism
a(2,2)   = 10*str2num(get(handles.edit_A1,'String'));
phi(2,2) = str2num(get(handles.edit_Phi22,'String'));
% Coma
a(3,1)   = 10*str2num(get(handles.edit_B2,'String'));
phi(3,1) = str2num(get(handles.edit_Phi31,'String'));
% 3-fold astigmatism:
a(3,3)   = 10*str2num(get(handles.edit_A2,'String'));
phi(3,3) = str2num(get(handles.edit_Phi33,'String'));     
% 4-fold astigmatism
a(4,4)   = 1e4*str2num(get(handles.edit_A3,'String')); 
phi(4,4) = str2num(get(handles.edit_Phi44,'String'));
% star aberration
a(4,2)   = 1e4*str2num(get(handles.edit_S3,'String')); 
phi(4,2) = str2num(get(handles.edit_Phi42,'String')); 
% Spherical aberration:
% a(4,0)   = 1e4*str2num(get(handles.edit_C3,'String')));
c(4)     = 1e4*str2num(get(handles.edit_C3,'String'));

% 5-fold astigmatism:
a(5,5)   = 1e4*str2num(get(handles.edit_A4,'String'));
phi(5,5) = str2num(get(handles.edit_Phi55,'String'));     
%
a(5,3)   = 1e4*str2num(get(handles.edit_D4,'String'));
phi(5,3) = str2num(get(handles.edit_Phi53,'String'));     
%
a(5,1)   = 1e4*str2num(get(handles.edit_B4,'String'));
phi(5,1) = str2num(get(handles.edit_Phi51,'String'));

% 6-fold astigmatism
a(6,6)   = 1e7*str2num(get(handles.edit_A5,'String'));
phi(6,6) = str2num(get(handles.edit_Phi66,'String'));     
%
a(6,4)   = 1e7*str2num(get(handles.edit_D5,'String'));
phi(6,4) = str2num(get(handles.edit_Phi64,'String'));     
%
a(6,2)   = 1e7*str2num(get(handles.edit_B5,'String'));
phi(6,2) = str2num(get(handles.edit_Phi62,'String'));
% C5:
c(6)     = 1e7*str2num(get(handles.edit_C5,'String'));

% HighVoltage:
handles.params(5) = str2num(get(handles.edit_HighVoltage,'String'));

% if get(handles.radiobutton_Nomenclature2,'Value')  % do a converion for Krivanek's nomenclature:
if handles.oldNomenclature == 2; 
    phi(2:3,:) = phi(2:3,:)*10;
    phi(4:5,:) = phi(4:5,:)*1e4;
    phi(6,:) = phi(6,:)*1e7;
    v = a+i*phi;
    a = abs(v);
    phi = angle(v);
else
   % adjust scale of angles
   if handles.oldDeg
       phi = phi*pi/180;
   end
end


handles.a   = a;
handles.phi = phi;
handles.c   = c;


% --- Executes on button press in pushbutton_Scherzer.
function pushbutton_Scherzer_Callback(hObject, eventdata, handles)
handles = readParams(handles);
C3  = handles.c(4);  % Cs in A
HighVoltage = handles.params(5);
def = -sign(C3)*sqrt(1.5*abs(C3)*wavelength(HighVoltage)); % defocus in A
handles.c(2) = def;

writeParams(handles);


function writeParams(handles)
a   = handles.a;    
phi = handles.phi;
c   = handles.c;

if get(handles.radiobutton_Nomenclature2,'Value')  % do a converion for Krivanek's nomenclature:
    v = a.*exp(i*phi);
    a = real(v);
    phi = imag(v);
    phi(2:3,:) = phi(2:3,:)*0.1;
    phi(4:5,:) = phi(4:5,:)*1e-4;
    phi(6,:) = phi(6,:)*1e-7;
else
   % adjust scale of angles
   if get(handles.radiobutton_deg,'Value')
       phi = phi*180/pi;
   end
end


% Spherical aberration:
% a(2,0) = 10*str2num(get(handles.edit_C1,'String'));
set(handles.edit_C1,'String',sprintf('%.3f',0.1*c(2)));   % defocus in nm
% Astigmatism
set(handles.edit_A1,'String',sprintf('%.3f',0.1*a(2,2)));  % astigmatism in nm
set(handles.edit_Phi22,'String',sprintf('%.3f',phi(2,2)));
% Coma
set(handles.edit_B2,'String',sprintf('%.3f',0.1*a(3,1)));
set(handles.edit_Phi31,'String',sprintf('%.3f',phi(3,1)));
% 3-fold astigmatism:
set(handles.edit_A2,'String',sprintf('%.3f',0.1*a(3,3)));
set(handles.edit_Phi33,'String',sprintf('%.3f',phi(3,3)));     
% 4-fold astigmatism
set(handles.edit_A3,'String',sprintf('%.3f',1e-4*a(4,4))); 
set(handles.edit_Phi44,'String',sprintf('%.3f',phi(4,4)));
% star aberration
set(handles.edit_S3,'String',sprintf('%.3f',1e-4*a(4,2))); 
set(handles.edit_Phi42,'String',sprintf('%.3f',phi(4,2)));
% Spherical aberration:
% a(4,0)   = 1e4*str2num(get(handles.edit_C3,'String')));
set(handles.edit_C3,'String',sprintf('%.3f',1e-4*c(4)));

% 5-fold astigmatism:
set(handles.edit_A4,'String',sprintf('%.3f',1e-4*a(5,5)));
set(handles.edit_Phi55,'String',sprintf('%.3f',phi(5,5)));     
%
set(handles.edit_D4,'String',sprintf('%.3f',1e-4*a(5,3)));
set(handles.edit_Phi53,'String',sprintf('%.3f',phi(5,3)));     
%
set(handles.edit_B4,'String',sprintf('%.3f',1e-4*a(5,1)));
set(handles.edit_Phi51,'String',sprintf('%.3f',phi(5,1)));

% 6-fold astigmatism
set(handles.edit_A5,'String',sprintf('%.3f',1e-7*a(6,6)));
set(handles.edit_Phi66,'String',sprintf('%.3f',phi(6,6)));     
%
set(handles.edit_D5,'String',sprintf('%.3f',1e-7*a(6,4)));
set(handles.edit_Phi64,'String',sprintf('%.3f',phi(6,4)));     
%
set(handles.edit_B5,'String',sprintf('%.3f',1e-7*a(6,2)));
set(handles.edit_Phi62,'String',sprintf('%.3f',phi(6,2)));
% C5:
set(handles.edit_C5,'String',sprintf('%.3f',1e-7*c(6)));

% HighVoltage:
set(handles.edit_HighVoltage,'String',handles.params(5));




function edit_Alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_Alpha as a double


% --- Executes during object creation, after setting all properties.
function edit_Alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








% --- Executes on button press in radiobutton_deg.
function radiobutton_deg_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_deg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_deg
handles = readParams(handles);
handles.oldDeg = 1;
writeParams(handles);
guidata(hObject, handles);




% --- Executes on button press in radiobutton_rad.
function radiobutton_rad_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_rad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_rad
handles = readParams(handles);
handles.oldDeg = 0;
writeParams(handles);
guidata(hObject, handles);




% --- Executes on button press in pushbutton_ZoomIn.
function pushbutton_ZoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Probe);
%    imagesc(dx*(-Nx/2+[0:Nx-1]),dy*(-Ny/2+[0:Ny-1]),fftshift(abs(ifft2(ifftshift(probe))).^2));
%    colormap('default');
x = xlim;
xlim(0.5*x);
y = ylim;
ylim(0.5*y);
axes(handles.axes_Phasemap);
x = xlim;
xlim(0.5*x);
y = ylim;
ylim(0.5*y);
%    imagesc(dkx*(-Nx/2+[0:Nx-1]),dky*(-Ny/2+[0:Ny-1]),angle(probe));    


% --- Executes on button press in pushbutton_ZoomOut.
function pushbutton_ZoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_Probe);
x = xlim;
xlim(2*x);
y = ylim;
ylim(2*y);
axes(handles.axes_Phasemap);
x = xlim;
xlim(2*x);
y = ylim;
ylim(2*y);




% --- Executes on button press in checkbox_Color.
function checkbox_Color_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Color








function edit_HighVoltage_Callback(hObject, eventdata, handles)
% hObject    handle to edit_HighVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_HighVoltage as text
%        str2double(get(hObject,'String')) returns contents of edit_HighVoltage as a double


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




% --- Executes on button press in checkbox_Surf.
function checkbox_Surf_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Surf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Surf





function edit_Delta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Delta as text
%        str2double(get(hObject,'String')) returns contents of edit_Delta as a double


% --- Executes during object creation, after setting all properties.
function edit_Delta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Delta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_CTF_low_Callback(hObject, eventdata, handles)
pushbutton_Test_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit_CTF_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CTF_low (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_CTF_high_Callback(hObject, eventdata, handles)
pushbutton_Test_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function edit_CTF_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CTF_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton_CTF_NewFigure.
function pushbutton_CTF_NewFigure_Callback(hObject, eventdata, handles)
handles.createNewFigure = 1;
pushbutton_Test_Callback(hObject, eventdata, handles);



% --- Executes on button press in radiobutton_realPart.
function radiobutton_realPart_Callback(hObject, eventdata, handles)
pushbutton_Test_Callback(hObject, eventdata, handles);



% --- Executes on button press in radiobutton_imagPart.
function radiobutton_imagPart_Callback(hObject, eventdata, handles)
pushbutton_Test_Callback(hObject, eventdata, handles);



% --- Executes on button press in radiobutton_phasePart.
function radiobutton_phasePart_Callback(hObject, eventdata, handles)
pushbutton_Test_Callback(hObject, eventdata, handles);



% --- Executes on button press in radiobutton_amplitudePart.
function radiobutton_amplitudePart_Callback(hObject, eventdata, handles)
pushbutton_Test_Callback(hObject, eventdata, handles);



% --- Executes on button press in radiobutton_chi.
function radiobutton_chi_Callback(hObject, eventdata, handles)
pushbutton_Test_Callback(hObject, eventdata, handles);
