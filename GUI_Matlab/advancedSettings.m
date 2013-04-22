function varargout = advancedSettings(varargin)
% ADVANCEDSETTINGS M-file for advancedSettings.fig
%      ADVANCEDSETTINGS by itself, creates a new ADVANCEDSETTINGS or raises the
%      existing singleton*.
%
%      H = ADVANCEDSETTINGS returns the handle to a new ADVANCEDSETTINGS or the handle to
%      the existing singleton*.
%
%      ADVANCEDSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCEDSETTINGS.M with the given input arguments.
%
%      ADVANCEDSETTINGS('Property','Value',...) creates a new ADVANCEDSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before advancedSettings_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to advancedSettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help advancedSettings

% Last Modified by GUIDE v2.5 03-Jan-2008 11:38:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @advancedSettings_OpeningFcn, ...
                   'gui_OutputFcn',  @advancedSettings_OutputFcn, ...
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

% --- Executes just before advancedSettings is made visible.
function advancedSettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to advancedSettings (see VARARGIN)

% Choose default command line output for advancedSettings
handles.output = 'Yes';

% Update handles structure
guidata(hObject, handles);

% Insert custom Title and Text if specified by the user
% Hint: when choosing keywords, be sure they are not easily confused 
% with existing figure properties.  See the output of set(figure) for
% a list of figure properties.
handles.data = varargin{1};

set(handles.edit_SaveLevel,'String',sprintf('%d',handles.data.saveLevel));
set(handles.edit_PrintLevel,'String',sprintf('%d',handles.data.printLevel));
set(handles.edit_AtomRadius,'String',sprintf('%.1f',handles.data.atomRadius));
set(handles.checkbox_SavePotential,'Value',handles.data.savePotential);
set(handles.checkbox_SaveTotalPotential,'Value',handles.data.saveTotalPotential);
set(handles.edit_propProgInterval,'String',sprintf('%d',handles.data.propProgInterval));
set(handles.edit_potProgInterval,'String',sprintf('%d',handles.data.potProgInterval));

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


% Make the GUI modal
set(handles.AdvancedSettings,'WindowStyle','modal');

% UIWAIT makes advancedSettings wait for user response (see UIRESUME)
uiwait(handles.AdvancedSettings);

% --- Outputs from this function are returned to the command line.
function varargout = advancedSettings_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.data.saveLevel = floor(str2num(get(handles.edit_SaveLevel,'String')));
handles.data.printLevel = floor(str2num(get(handles.edit_PrintLevel,'String')));
handles.data.atomRadius = abs(str2num(get(handles.edit_AtomRadius,'String')));
handles.data.savePotential = get(handles.checkbox_SavePotential,'Value');
handles.data.saveTotalPotential = get(handles.checkbox_SaveTotalPotential,'Value');
handles.data.propProgInterval = floor(str2num(get(handles.edit_propProgInterval,'String')));
handles.data.potProgInterval = floor(str2num(get(handles.edit_potProgInterval,'String')));

varargout{1} = handles.data;

% The figure can be deleted now
delete(handles.AdvancedSettings);



% --- Executes when user attempts to close AdvancedSettings.
function AdvancedSettings_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to AdvancedSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.AdvancedSettings, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.AdvancedSettings);
else
    % The GUI is no longer waiting, just close it
    delete(handles.AdvancedSettings);
end


% --- Executes on key press over AdvancedSettings with no controls selected.
function AdvancedSettings_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to AdvancedSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for "enter" or "escape"
if isequal(get(hObject,'CurrentKey'),'escape')
    % User said no by hitting escape
    handles.output = 'No';
    
    % Update handles structure
    guidata(hObject, handles);
    
    uiresume(handles.AdvancedSettings);
end    
    
if isequal(get(hObject,'CurrentKey'),'return')
    uiresume(handles.AdvancedSettings);
end    


% --- Executes on button press in pushbutton_Ok.
function pushbutton_Ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.AdvancedSettings);



% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.output = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.AdvancedSettings);
% delete(handles.AdvancedSettings);


function edit_PrintLevel_Callback(hObject, eventdata, handles)
value = floor(str2double(get(hObject,'String')));
if value < 0, value = 0; end
if value > 5, value = 5; end
set(hObject,'String',sprintf('%d',value));


% --- Executes during object creation, after setting all properties.
function edit_PrintLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PrintLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SaveLevel_Callback(hObject, eventdata, handles)
value = floor(str2double(get(hObject,'String')));
if value < 0, value = 0; end
if value > 3, value = 3; end
set(hObject,'String',sprintf('%d',value));


% --- Executes during object creation, after setting all properties.
function edit_SaveLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SaveLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_SavePotential.
function checkbox_SavePotential_Callback(hObject, eventdata, handles)



function edit_AtomRadius_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AtomRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AtomRadius as text
%        str2double(get(hObject,'String')) returns contents of edit_AtomRadius as a double
value = str2double(get(hObject,'String'));
if value < 1, value = 1; end
if value > 15, value = 15; end
set(hObject,'String',sprintf('%.1f',value));

% --- Executes during object creation, after setting all properties.
function edit_AtomRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AtomRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_propProgInterval_Callback(hObject, eventdata, handles)
value = floor(str2double(get(hObject,'String')));
if value < 1, value = 1; end
if value > 10000, value = 10000; end
set(hObject,'String',sprintf('%d',value));


% --- Executes during object creation, after setting all properties.
function edit_propProgInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_propProgInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_potProgInterval_Callback(hObject, eventdata, handles)
value = floor(str2double(get(hObject,'String')));
if value < 1, value = 1; end
if value > 100000, value = 100000; end
set(hObject,'String',sprintf('%d',value));


% --- Executes during object creation, after setting all properties.
function edit_potProgInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_potProgInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_SaveTotalPotential.
function checkbox_SaveTotalPotential_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_SaveTotalPotential (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_SaveTotalPotential


