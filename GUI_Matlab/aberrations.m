function varargout = aberrations(varargin)
% ABERRATIONS M-file for aberrations.fig
%      ABERRATIONS by itself, creates a new ABERRATIONS or raises the
%      existing singleton*.
%
%      H = ABERRATIONS returns the handle to a new ABERRATIONS or the handle to
%      the existing singleton*.
%
%      ABERRATIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABERRATIONS.M with the given input arguments.
%
%      ABERRATIONS('Property','Value',...) creates a new ABERRATIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before aberrations_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to aberrations_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help aberrations

% Last Modified by GUIDE v2.5 23-Nov-2011 10:44:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @aberrations_OpeningFcn, ...
                   'gui_OutputFcn',  @aberrations_OutputFcn, ...
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

% --- Executes just before aberrations is made visible.
function aberrations_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to aberrations (see VARARGIN)

% Choose default command line output for aberrations
handles.oldNomenclature = 1;
handles.oldDeg = 1;
set(handles.radiobutton_Nomenclature1,'Value',1);
set(handles.radiobutton_deg,'Value',1);
handles.a      = varargin{1};
handles.phi    = varargin{2};
handles.c      = varargin{3};
handles.params = varargin{4};
set(handles.edit_Alpha,'String',handles.params(6));
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

% UIWAIT makes aberrations wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = aberrations_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles = readParams(handles);
handles.params(6) = str2num(get(handles.edit_Alpha,'String'));

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

Nx = handles.params(1);    % number of pixels in x-direction   
Ny = handles.params(2);    % number of pixels in y-direction   
dx = handles.params(3);    % sampling in x-direction (unit: A) 
dy = handles.params(4);    % sampling in x-direction (unit: A) 
dkx = 1/(Nx*dx);           % rec. space sampling in x-direction
dky = 1/(Ny*dy);           % rec. space sampling in x-direction
lambda = wavelength(handles.params(5));
alpha = 1e-3*str2num(get(handles.edit_Alpha,'String'));
qmax  = sin(alpha)/lambda; 

% Beam tilt: 
tiltKx = 0;
tiltKy = 0;

[kx,ky] = meshgrid(dkx*(-Nx/2+[0:Nx-1])+tiltKx,dky*(-Ny/2+[0:Ny-1])+tiltKy);
% kx2 = kx.^2;
% ky2 = ky.^2;
% k2 = kx2+ky2;
ktheta = asin(sqrt(kx.^2+ky.^2)*lambda);
kphi = atan2(ky,kx);
ktm = asin(qmax*lambda);
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




chi   = 2*pi/lambda*(1/2*(a(2,2).*cos(2*(kphi-phi(2,2)))+c(2)).*ktheta.^2+...
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
        zlabel('intensity (arb. units)');
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
    xlabel('x in A');
    ylabel('y in A');
    axes(handles.axes_Phasemap);
    imagesc(dkx*(-Nx/2+[0:Nx-1]),dky*(-Ny/2+[0:Ny-1]),angle(probe));    
    set(gca,'YDir','normal');
    if get(handles.checkbox_Color,'Value')
        colormap('default');
    else
        colormap('gray');
    end
    xlabel('k_x in 1/A');
    ylabel('k_y in 1/A');

end
handles.probe = probe;
guidata(hObject, handles);


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




% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
[filename, pathname] = uiputfile('*.img', 'Save spot image as');
if isequal(filename,0) return; end

Nx = handles.params(1);    % number of pixels in x-direction   
Ny = handles.params(2);    % number of pixels in y-direction   
dx = handles.params(3);    % sampling in x-direction (unit: A) 
dy = handles.params(4);    % sampling in x-direction (unit: A) 

img = fftshift(abs(ifft2(ifftshift(handles.probe))).^2);
binwrite2D(img,fullfile(pathname, filename),dx,dy,1,0,0);


