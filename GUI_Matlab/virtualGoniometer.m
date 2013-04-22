function varargout = virtualGoniometer(varargin)
% VIRTUALGONIOMETER M-file for virtualGoniometer.fig
%      VIRTUALGONIOMETER, by itself, creates a new VIRTUALGONIOMETER or raises the existing
%      singleton*.
%
%      H = VIRTUALGONIOMETER returns the handle to a new VIRTUALGONIOMETER or the handle to
%      the existing singleton*.
%
%      VIRTUALGONIOMETER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIRTUALGONIOMETER.M with the given input arguments.
%
%      VIRTUALGONIOMETER('Property','Value',...) creates a new VIRTUALGONIOMETER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before virtualGoniometer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to virtualGoniometer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help virtualGoniometer

% Last Modified by GUIDE v2.5 28-Apr-2011 10:50:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @virtualGoniometer_OpeningFcn, ...
                   'gui_OutputFcn',  @virtualGoniometer_OutputFcn, ...
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


% --- Executes just before virtualGoniometer is made visible.
function virtualGoniometer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to virtualGoniometer (see VARARGIN)

% Choose default command line output for virtualGoniometer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes virtualGoniometer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = virtualGoniometer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_a_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_a as text
%        str2double(get(hObject,'String')) returns contents of edit_a as a double


% --- Executes during object creation, after setting all properties.
function edit_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_b_Callback(hObject, eventdata, handles)
% hObject    handle to edit_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_b as text
%        str2double(get(hObject,'String')) returns contents of edit_b as a double


% --- Executes during object creation, after setting all properties.
function edit_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_c_Callback(hObject, eventdata, handles)
% hObject    handle to edit_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_c as text
%        str2double(get(hObject,'String')) returns contents of edit_c as a double


% --- Executes during object creation, after setting all properties.
function edit_c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_c (see GCBO)
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



function edit_beta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beta as text
%        str2double(get(hObject,'String')) returns contents of edit_beta as a double


% --- Executes during object creation, after setting all properties.
function edit_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gamma as text
%        str2double(get(hObject,'String')) returns contents of edit_gamma as a double


% --- Executes during object creation, after setting all properties.
function edit_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function chi2 = refineMatrix(H3,H,alpha,beta,gamma)
H(3,:) = H3.';
alphaP = acos(sum(H(2,:).*H(3,:))/sqrt(sum(H(2,:).^2)*sum(H(3,:).^2)))*180/pi;
betaP =  acos(sum(H(1,:).*H(3,:))/sqrt(sum(H(1,:).^2)*sum(H(3,:).^2)))*180/pi;
% fprintf('Angles of matrix: alpha = %.2f°, beta = %.2f°, gamma = %.2f°\n',alphaP,betaP,gammaP);
chi2 = 100000*(alphaP-alpha)^2+(betaP-beta)^2;



% --- Executes on button press in pushbutton_Matrix.
function pushbutton_Matrix_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = str2double(get(handles.edit_a,'String'));
b = str2double(get(handles.edit_b,'String'));
c = str2double(get(handles.edit_c,'String'));

alpha = str2double(get(handles.edit_alpha,'String'))*pi/180;
beta  = str2double(get(handles.edit_beta,'String'))*pi/180;
gamma = str2double(get(handles.edit_gamma,'String'))*pi/180;

if ((abs(alpha) < 1e-4) || (abs(beta) < 1e-4) || (abs(gamma) < 1e-4))
    errordlg('At least one of the angles is too small!')
    return
end

H00 = a;
H01 = 0;
H02 = 0;
H10 = b*cos(gamma);
H11 = b*sin(gamma);
H12 = 0;
H20 = 0;
H21 = 0;
H22 = c;
H = [H00 H01 H02;H10 H11 H12; H20 H21 H22];

% start a refinement to search for the correct matrix entries in the 3rd
% row:
if (abs(alpha-pi/2) > 1e-4) ||  (abs(beta-pi/2) > 1e-4)
    H3 = H(3,:).';
    [H3,chi2] = fminsearch(@(H3) refineMatrix(H3,H,alpha*180/pi,beta*180/pi,gamma*180/pi), H3);
    fprintf('found structure matrix with chi2 = %f\n',chi2);
    H(3,:) = H3.'*c/sqrt(sum(H3.^2));
    H20 = H(3,1);
    H21 = H(3,2);
    H22 = abs(H(3,3));
end

set(handles.edit_H00,'String',sprintf('%.3f',H00));
set(handles.edit_H01,'String',sprintf('%.3f',H01));
set(handles.edit_H02,'String',sprintf('%.3f',H02));

set(handles.edit_H10,'String',sprintf('%.3f',H10));
set(handles.edit_H11,'String',sprintf('%.3f',H11));
set(handles.edit_H12,'String',sprintf('%.3f',H12));

set(handles.edit_H20,'String',sprintf('%.3f',H20));
set(handles.edit_H21,'String',sprintf('%.3f',H21));
set(handles.edit_H22,'String',sprintf('%.3f',H22));

alpha = acos(sum(H(2,:).*H(3,:))/sqrt(sum(H(2,:).^2)*sum(H(3,:).^2)))*180/pi;
beta =  acos(sum(H(1,:).*H(3,:))/sqrt(sum(H(1,:).^2)*sum(H(3,:).^2)))*180/pi;
gamma = acos(sum(H(1,:).*H(2,:))/sqrt(sum(H(1,:).^2)*sum(H(2,:).^2)))*180/pi;
fprintf('Angles of matrix: alpha = %.2f°, beta = %.2f°, gamma = %.2f°\n',alpha,beta,gamma);
a = sqrt(sum(H(1,:).^2));
b = sqrt(sum(H(2,:).^2));
c = sqrt(sum(H(3,:).^2));
fprintf('Vector lengths: a= %.2f, b = %.2f, c = %.2f\n',a,b,c);

% --- Executes on button press in pushbutton_FindTilts.
function pushbutton_FindTilts_Callback(hObject, eventdata, handles)

format compact
H = zeros(3);
H(1,1) = str2double(get(handles.edit_H00,'String'));
H(1,2) = str2double(get(handles.edit_H01,'String'));
H(1,3) = str2double(get(handles.edit_H02,'String'));

H(2,1) = str2double(get(handles.edit_H10,'String'));
H(2,2) = str2double(get(handles.edit_H11,'String'));
H(2,3) = str2double(get(handles.edit_H12,'String'));

H(3,1) = str2double(get(handles.edit_H20,'String'));
H(3,2) = str2double(get(handles.edit_H21,'String'));
H(3,3) = str2double(get(handles.edit_H22,'String'));

refZone = zeros(3,1);
zone = zeros(3,1);
refZone(1) =  str2double(get(handles.edit_href,'String'));
refZone(2) =  str2double(get(handles.edit_kref,'String'));
refZone(3) =  str2double(get(handles.edit_lref,'String'));

zone(1) =  str2double(get(handles.edit_h,'String'));
zone(2) =  str2double(get(handles.edit_k,'String'));
zone(3) =  str2double(get(handles.edit_l,'String'));

[rotAngles,Mrot] = rotateZoneAxis(zone,H.',refZone);


set(handles.edit_tiltx,'String',sprintf('%.3f',rotAngles(1)));
set(handles.edit_tilty,'String',sprintf('%.3f',rotAngles(2)));
set(handles.edit_tiltz,'String',sprintf('%.3f',rotAngles(3)));

Mrot

function edit_H00_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H00 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H00 as text
%        str2double(get(hObject,'String')) returns contents of edit_H00 as a double


% --- Executes during object creation, after setting all properties.
function edit_H00_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H00 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_H01_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H01 as text
%        str2double(get(hObject,'String')) returns contents of edit_H01 as a double


% --- Executes during object creation, after setting all properties.
function edit_H01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_H02_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H02 as text
%        str2double(get(hObject,'String')) returns contents of edit_H02 as a double


% --- Executes during object creation, after setting all properties.
function edit_H02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_H10_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H10 as text
%        str2double(get(hObject,'String')) returns contents of edit_H10 as a double


% --- Executes during object creation, after setting all properties.
function edit_H10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_H11_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H11 as text
%        str2double(get(hObject,'String')) returns contents of edit_H11 as a double


% --- Executes during object creation, after setting all properties.
function edit_H11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_H12_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H12 as text
%        str2double(get(hObject,'String')) returns contents of edit_H12 as a double


% --- Executes during object creation, after setting all properties.
function edit_H12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_H20_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H20 as text
%        str2double(get(hObject,'String')) returns contents of edit_H20 as a double


% --- Executes during object creation, after setting all properties.
function edit_H20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_H21_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H21 as text
%        str2double(get(hObject,'String')) returns contents of edit_H21 as a double


% --- Executes during object creation, after setting all properties.
function edit_H21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_H22_Callback(hObject, eventdata, handles)
% hObject    handle to edit_H22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_H22 as text
%        str2double(get(hObject,'String')) returns contents of edit_H22 as a double


% --- Executes during object creation, after setting all properties.
function edit_H22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_H22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_h_Callback(hObject, eventdata, handles)
% hObject    handle to edit_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_h as text
%        str2double(get(hObject,'String')) returns contents of edit_h as a double


% --- Executes during object creation, after setting all properties.
function edit_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_k_Callback(hObject, eventdata, handles)
% hObject    handle to edit_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_k as text
%        str2double(get(hObject,'String')) returns contents of edit_k as a double


% --- Executes during object creation, after setting all properties.
function edit_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lref_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lref as text
%        str2double(get(hObject,'String')) returns contents of edit_lref as a double


% --- Executes during object creation, after setting all properties.
function edit_lref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
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



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tilty_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tilty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tilty as text
%        str2double(get(hObject,'String')) returns contents of edit_tilty as a double


% --- Executes during object creation, after setting all properties.
function edit_tilty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tilty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function edit_href_Callback(hObject, eventdata, handles)
% hObject    handle to edit_href (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_href as text
%        str2double(get(hObject,'String')) returns contents of edit_href as a double


% --- Executes during object creation, after setting all properties.
function edit_href_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_href (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_l_Callback(hObject, eventdata, handles)
% hObject    handle to edit_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_l as text
%        str2double(get(hObject,'String')) returns contents of edit_l as a double


% --- Executes during object creation, after setting all properties.
function edit_l_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_l (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kref_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kref as text
%        str2double(get(hObject,'String')) returns contents of edit_kref as a double


% --- Executes during object creation, after setting all properties.
function edit_kref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kref (see GCBO)
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



function edit_tiltz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tiltz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tiltz as text
%        str2double(get(hObject,'String')) returns contents of edit_tiltz as a double


% --- Executes during object creation, after setting all properties.
function edit_tiltz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tiltz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


