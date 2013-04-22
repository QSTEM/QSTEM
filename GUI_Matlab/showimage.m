function varargout = showimage(varargin)
% %compile using 'mcc -m showimage'
% SHOWIMAGE M-file for showimage.fig
%      SHOWIMAGE, by itself, creates a new SHOWIMAGE or raises the existing
%      singleton*.
%
%      H = SHOWIMAGE returns the handle to a new SHOWIMAGE or the handle to
%      the existing singleton*.
%
%      SHOWIMAGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWIMAGE.M with the given input arguments.
%
%      SHOWIMAGE('Property','Value',...) creates a new SHOWIMAGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before showimage_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to showimage_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help showimage

% Last Modified by GUIDE v2.5 29-Jun-2010 09:59:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showimage_OpeningFcn, ...
                   'gui_OutputFcn',  @showimage_OutputFcn, ...
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


% --- Executes just before showimage is made visible.
function showimage_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showimage (see VARARGIN)

% Choose default command line output for showimage
handles.output = hObject;
handles.directoryname = pwd();
set(handles.edit_directoryname,'String',handles.directoryname);

if (length(varargin) > 0)
    inputDir = varargin{1};
    filename = '';
    if isempty(findstr(inputDir,':\'))
        inputDir = fullfile(pwd(),inputDir);
    end
    if ~isdir(inputDir)
        % Check to see whether file exists:
        fid = fopen(inputDir);
        if fid < 0
            inputDir = pwd;
        else
            fclose(fid)
            [inputDir,name,ext] = fileparts(inputDir);
            filename = [name,ext];
        end
    end
    handles.directoryname = inputDir;
    set(handles.edit_directoryname,'String',handles.directoryname);
    handles = pushbutton_SelectDir_Callback(hObject, eventdata, handles,0,filename);
    % if file:
    % handles = listbox_image_Callback(hObject, eventdata, handles,1);

end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes showimage wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = showimage_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_directoryname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_directoryname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_directoryname as text
%        str2double(get(hObject,'String')) returns contents of edit_directoryname as a double


% --- Executes during object creation, after setting all properties.
function edit_directoryname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_directoryname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_SelectDir.
function handles = pushbutton_SelectDir_Callback(hObject, eventdata, handles,askUser,filename)
if (nargin < 4)
    askUser = 1;
end
if (nargin < 5)
    filename = '';
end

handles.directoryname = get(handles.edit_directoryname,'String');
if askUser
    handles.directoryname = uigetdir(handles.directoryname);
    if ischar(handles.directoryname)
        set(handles.edit_directoryname,'String',handles.directoryname);
    end
end
if ~ischar(handles.directoryname)
    return
end

% handles.directoryname
handles.files = dir(sprintf('%s\\*.img',handles.directoryname));
% handles.files.name
set(handles.listbox_image,'String',{handles.files.name});
set(handles.listbox_image,'Value',1);
handles = listbox_image_Callback(hObject, eventdata, handles,0,filename);
guidata(hObject, handles);

% --- Executes on button press in pushbutton_NewWindow.
function pushbutton_NewWindow_Callback(hObject, eventdata, handles)
handles = listbox_image_Callback(hObject, eventdata, handles,1);
guidata(hObject, handles);

% --- Executes on selection change in listbox_image.
function handles = listbox_image_Callback(hObject, eventdata, handles,showInNewWindow,fname)
if nargin < 4
    showInNewWindow = 0;
end
if (nargin < 5)
    fname = '';
end

if isempty(fname)
    filelist = get(handles.listbox_image,'String');
    % filename = sprintf('%s\\%s',handles.directoryname,filelist{get(handles.listbox_image,'Value')});
    if isempty(filelist)
        filename = '';
    else
        filename = fullfile(handles.directoryname,filelist{get(handles.listbox_image,'Value')});
    end
else
    filelist = get(handles.listbox_image,'String');
    for j=1:length(filelist)
        if strcmp(filelist{j},fname)
            set(handles.listbox_image,'Value',j);
            break
        end
    end
    filename = fullfile(handles.directoryname,fname);
end
if isempty(filename)
    return
end
[img,t,dx,dy,params,comment] = binread2D(filename);
[Ny,Nx] = size(img);
realFlag = isreal(img);
ss = 0;
STEMimage = 0;

comment = char(comment);
if length(comment) == 10
    if sum(abs(comment-'STEM image')) == 0
        % fprintf('STEM image\n')
        ss = abs(str2num(get(handles.edit_SourceSize,'String')));  
        sourceSize = ss;
        % scale the source size to produce a Gaussian of FWHM = ss:
        ss = ss/(sqrt(-log(0.5))*2);  % = ss/1.66510922231540;
        STEMimage = 1;            
    end
end

if STEMimage
    set(handles.text_ImageProperties,'String',sprintf('Thickness: %.1fA\nSampling: %.1f x %.1fA, %d runs averaged over\n%s',t,dx,dy,params(1),comment));
    set(handles.edit_SourceSize,'Enable','on');
    set(handles.pushbutton_Quantify,'Enable','on');

    if (0)
        set(handles.edit_ReplicateX,'Enable','on');
        set(handles.edit_ReplicateY,'Enable','on');
    end
else
    set(handles.text_ImageProperties,'String',sprintf('Thickness: %.1fA\nSampling: %.1f x %.1fA\n%s',t,dx,dy,comment));
    set(handles.edit_SourceSize,'Enable','off');
    set(handles.pushbutton_Quantify,'Enable','off');
    if (0)
        set(handles.edit_ReplicateX,'Enable','off');
        set(handles.edit_ReplicateY,'Enable','off');
    end
end


handles.Oversampling = abs(str2num(get(handles.edit_Oversampling,'String')));
if handles.Oversampling*max(Nx,Ny) > 1024
    handles.Oversampling = floor(1024/max(Nx,Ny));
end
if handles.Oversampling < 1
    handles.Oversampling = 1;
end
set(handles.edit_Oversampling,'String',handles.Oversampling);

% oversample image, if desired:
scale = 1;
if handles.Oversampling > 1
    NyO = round(handles.Oversampling*Ny);
    NxO = round(handles.Oversampling*Nx);
    imgft = zeros(NyO,NxO);
    NxOMid = floor(handles.Oversampling*Nx/2)+1;
    NyOMid = floor(handles.Oversampling*Ny/2)+1;
    NxMid = floor(Nx/2)+1;
    NyMid = floor(Ny/2)+1;
    imgft(NyOMid-NyMid+[1:Ny],NxOMid-NxMid+[1:Nx]) = fftshift(fft2(img));
    % Apply effective source size for STEM images:
    if ss > 0
        [qx,qy] = meshgrid((-NxOMid+[1:NxO])/(Nx*dx),(-NyOMid+[1:NyO])/(Ny*dy));
        img = ifft2(ifftshift(imgft.*exp(-((pi*ss)^2*(qx.^2+qy.^2)))));        
    else
        img = ifft2(ifftshift(imgft));
    end
    if realFlag
       img = real(img); 
    end
    scale = 1/(Nx*Ny);
    Nx = round(handles.Oversampling*Nx);
    Ny = round(handles.Oversampling*Ny);
    dx = dx/handles.Oversampling;
    dy = dy/handles.Oversampling;
    scale = (scale*Nx*Ny);
else
    if ss > 0
        NxMid = floor(Nx/2)+1;
        NyMid = floor(Ny/2)+1;
        [qx,qy] = meshgrid((-NxMid+[1:Nx])/(Nx*dx),(-NyMid+[1:Ny])/(Ny*dy));
        img = real(ifft2(fft2(img).*ifftshift(exp(-((pi*ss)^2*(qx.^2+qy.^2))))));        
    end    
end

% decide which part of a complex image to use:
if ~isreal(img)
    set(handles.listbox_Complex,'Enable','on')
    % fprintf('Complex choice: %d\n',get(handles.listbox_Complex,'Value'));
    switch get(handles.listbox_Complex,'Value')
        case 1  % Modulus
            img = abs(img);
        case 2  % Phase
            img = angle(img);
        case 3  % real
            img = real(img);
        case 4  % imag
            img = imag(img);
        case 5  % imag
            img = abs(img).^2;
    end        
else
    set(handles.listbox_Complex,'Enable','off')
    if (0)
    switch get(handles.listbox_Complex,'Value')
        case 1  % Modulus
            img = abs(img);
        case 2  % Phase
            img = angle(img);
        case 3  % real
            img = real(img);
        case 4  % imag
            img = imag(img);
        case 5  % imag
            img = abs(img).^2;
    end            
    end
end



if get(handles.checkbox_Log,'Value')
    img = log(1+img-min(min(img)));
end

replicateX = str2num(get(handles.edit_ReplicateX,'String'));
replicateY = str2num(get(handles.edit_ReplicateY,'String'));
if ((replicateX > 1) || (replicateY > 1))
   img = repmat(img,replicateY,replicateX); 
end

img = img*scale;

if showInNewWindow
    figure
    if get(handles.checkbox_Surfplot,'Value')
        surf(dx*[0:replicateX*Nx-1],dy*[0:replicateY*Ny-1],img);
        colormap('default');
        shading interp
    else
        imagesc(dx*[0:replicateX*Nx-1],dy*[0:replicateY*Ny-1],img);
        set(gca,'YDir','normal');
        colormap('gray');
        axis equal; axis tight;
    end
else
    if get(handles.checkbox_Surfplot,'Value')
        surf(dx*[0:replicateX*Nx-1],dy*[0:replicateY*Ny-1],img,'Parent',handles.axes_Image);
        colormap('default');
        shading interp
    else
        imagesc(dx*[0:replicateX*Nx-1],dy*[0:replicateY*Ny-1],img,'Parent',handles.axes_Image);
        set(gca,'YDir','normal');
        colormap('gray');
        axis equal; axis tight;
    end
end
colorbar
handles.imgFinal = img;
handles.dx = dx;
handles.dy = dy;
handles.t = t;
handles.sourceSize = abs(str2num(get(handles.edit_SourceSize,'String'))); 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox_image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Oversampling_Callback(hObject, eventdata, handles)
handles = listbox_image_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_Oversampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Oversampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_Log.
function checkbox_Log_Callback(hObject, eventdata, handles)
handles = listbox_image_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in checkbox_Surfplot.
function checkbox_Surfplot_Callback(hObject, eventdata, handles)
handles = listbox_image_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on selection change in listbox_Complex.
function listbox_Complex_Callback(hObject, eventdata, handles)
handles = listbox_image_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox_Complex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_Complex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_SourceSize_Callback(hObject, eventdata, handles)
handles = listbox_image_Callback(hObject, eventdata, handles);
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
Nx = header(5);
Ny = header(4);
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
    img = fread(fid,[Nx,Ny],'float32');
case 4
    complexFlag = 0;
    doubleFlag = 1;
    if printFlag
    fprintf('64-bit real data, %.3fMB)\n',Nx*Ny*8/1048576);
    end
    img = fread(fid,[Nx,Ny],'float64');
case 2
    complexFlag = 1;
    doubleFlag = 0;
    if printFlag
    fprintf('32-bit complex data, %.3fMB)\n',Nx*Ny*8/1048576);
    end
    % img = fread(fid,[Nx,2*Ny],'float32');
    % img = img(1:Nx,1:2:2*Ny-1)+i*img(1:Nx,2:2:2*Ny);
    img = fread(fid,[2*Nx,Ny],'float32');
    img = img(1:2:2*Nx-1,1:Ny)+i*img(2:2:2*Nx,1:Ny);
case 6
    complexFlag = 1;
    doubleFlag = 1;
    if printFlag
    fprintf('64-bit complex data, %.3fMB)\n',Nx*Ny*16/1048576);
    end
    img = fread(fid,[2*Nx,Ny],'float64');
    img = img(1:2:2*Nx-1,1:Ny)+i*img(2:2:2*Nx,1:Ny);
case 1
    complexFlag = 0;
    doubleFlag = 0;
    if printFlag
        fprintf('16-bit integer data, %.3fMB)\n',Nx*Ny*8/1048576);
    end
    img = fread(fid,[Nx,Ny],'int16');    
case 5
    complexFlag = 0;
    doubleFlag = 0;
    if printFlag
    fprintf('32-bit integer data)\n');
    end
    img = fread(fid,[Nx,Ny],'int32');    
end
% imagesc(img); colormap('gray'); axis equal; axis tight;

fclose(fid);




% function f=Gaussian(sigma,x0,x,[y])
% This computes a normal distribution with a Given variance sigma^2

function f=Gaussian(sigma,x0,x,y)
% decide whether we do this in 1, or 2 dimensions:
[Ny,Nx] = size(x);
if (Nx==1) || (Ny==1)
    dx = x(2)-x(1);
    f = dx/(sigma*sqrt(2*pi))*exp(-(x-x0).^2./(2*sigma^2));
else
    dx = max(x(2,1)-x(1,1),x(1,2)-x(1,1));
    dy = max(y(2,1)-y(1,1),y(1,2)-y(1,1));
    f = dx*dy/(sigma^2*2*pi)*exp(-((x-x0(1)).^2+(y-x0(2)).^2)./(2*sigma^2));
end





function edit_ReplicateX_Callback(hObject, eventdata, handles)
handles = listbox_image_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_ReplicateX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ReplicateX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ReplicateY_Callback(hObject, eventdata, handles)
handles = listbox_image_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_ReplicateY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ReplicateY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end










% --- Executes on button press in pushbutton_Quantify.
function pushbutton_Quantify_Callback(hObject, eventdata, handles)
% handles = scaleSTEMImage(handles);
% guidata(hObject, handles);
if isfield(handles,'sourceSize')
    scaleSTEMImage(handles);
end

