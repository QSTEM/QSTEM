function varargout = DisplayModelProperties(varargin)
% DISPLAYMODELPROPERTIES M-file for DisplayModelProperties.fig
%      DISPLAYMODELPROPERTIES, by itself, creates a new DISPLAYMODELPROPERTIES or raises the existing
%      singleton*.
%
%      H = DISPLAYMODELPROPERTIES returns the handle to a new DISPLAYMODELPROPERTIES or the handle to
%      the existing singleton*.
%
%      DISPLAYMODELPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAYMODELPROPERTIES.M with the given
%      input arguments.
%
%      DISPLAYMODELPROPERTIES('Property','Value',...) creates a new DISPLAYMODELPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DisplayModelProperties_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DisplayModelProperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DisplayModelProperties

% Last Modified by GUIDE v2.5 05-Apr-2010 21:05:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DisplayModelProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @DisplayModelProperties_OutputFcn, ...
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


% --- Executes just before DisplayModelProperties is made visible.
function DisplayModelProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DisplayModelProperties (see VARARGIN)

% Choose default command line output for DisplayModelProperties
handles.output = hObject;
% The input will consist of a row vector of arguments
if size(varargin,2) > 0
    if isstruct(varargin{1})
        if isfield(varargin{1},'coords')
            handles.a = varargin{1}.coords;  % [coords aType DW occ charge];
            handles.Mm  = diag([max(varargin{1}.boundingBox(:,1)) - min(varargin{1}.boundingBox(:,1)), 
                max(varargin{1}.boundingBox(:,2)) - min(varargin{1}.boundingBox(:,2)),
                max(varargin{1}.boundingBox(:,3)) - min(varargin{1}.boundingBox(:,3))]);
        else
            handles.posFilePath = varargin{1}.configPath;  %
            handles.posFileName = varargin{1}.posFileName;  %
            filename = fullfile(handles.posFilePath,handles.posFileName);
            fid = fopen(filename);
            if fid ~= -1
                handles.fileName = filename;
                fclose(fid);
            else
                handles.fileName = handles.posFileName;
            end
            [coords,aType,Mm,DW,occ,charge] =  readCFG(handles.fileName);
            if isempty(coords)
                msgbox(sprintf('Could not open .cfg file %s!',handles.fileName));
            else
                handles.a = [aType coords DW occ charge];
                handles.Mm = Mm;
            end
            clear coords aType Mm DW occ charge
        end
    end
end
if (0)
    load fel0.mat
    handles.V0_neutral = V0_neutral;
    handles.V0_charged = V0_charged;
else
    handles.V0_neutral = [0.5286    0.4173    3.2556    3.0383    2.7850    2.4709    2.2034    1.9839 ...
    1.8017    1.6494    4.7387    5.1781    5.8665    5.7535    5.4738    5.1640 ...
    4.8566    4.5712    8.8966    9.8371    9.2502    8.7270    8.2633    6.9550 ...
    7.4772    7.1403    6.8305    6.5466    6.2851    6.0434    7.1413    7.3741 ...
    7.3686    7.2661    7.0805    6.8811   11.6690   13.0052   12.5953   12.1094 ...
   10.6991   10.2830   10.8336    9.5521    9.2244    7.5645    8.6414    9.2024 ...
   10.5928   11.0376   11.1772   11.1747   10.9811   10.7652   16.3483   18.1123 ...
   17.7054   17.2940   16.8631   16.4897   16.1899   15.7949   15.4700   15.2591 ...
   14.8526   14.5577   14.2730   13.9978   13.7312   13.4728   13.4834   13.2202 ...
   12.9316   12.6434   12.3627   12.0917   11.7793   10.8225   10.5466   10.9356 ...
   12.8045   13.5344   13.8911   14.0629   13.7529   13.4546   18.5513   20.3926 ...
   20.4678   20.2022   19.6404   19.2909];
   handles.V0_charged = [0  0   -3.0987   -1.4783    0     0         0   -1.0576 ...
   -0.8142         0   -3.6102   -2.1740   -1.7424         0         0         0 ...
   -1.5292         0   -5.4666   -3.5648   -2.3475   -1.7210   -1.3395   -1.2596 ...
   -2.3109   -2.1647   -2.0355   -1.9204   -1.8172   -1.7239         0         0 ...
         0         0   -1.7931         0   -6.1346   -4.1850   -2.8717   -2.1581 ...
   -1.5244   -1.2554         0         0         0   -0.6757   -1.7102   -2.0350 ...
         0         0         0         0   -2.2140         0   -7.3329   -5.1562 ...
   -3.6071   -2.7836   -3.3908   -3.3001         0   -3.1356   -3.0601   -3.0212 ...
   -2.9177   -2.8512   -2.7873   -2.7258   -2.6664   -2.6090   -2.6405         0 ...
         0         0         0         0         0         0         0         0 ...
         0         0         0         0         0         0         0         0 ...
         0         0         0         0];    
end

% Update handles structure
guidata(hObject, handles);
pushbutton_UpdateView_Callback(hObject, eventdata, handles)

% UIWAIT makes DisplayModelProperties wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DisplayModelProperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_SelectView.
function listbox_SelectView_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_SelectView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_SelectView contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_SelectView

function edit_PSFwidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PSFwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PSFwidth as text
%        str2double(get(hObject,'String')) returns contents of edit_PSFwidth as a double


% --- Executes during object creation, after setting all properties.
function edit_PSFwidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PSFwidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Nbins_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Nbins as text
%        str2double(get(hObject,'String')) returns contents of edit_Nbins as a double


% --- Executes during object creation, after setting all properties.
function edit_Nbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Nbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_UnitArea_Callback(hObject, eventdata, handles)
% hObject    handle to edit_UnitArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_UnitArea as text
%        str2double(get(hObject,'String')) returns contents of edit_UnitArea as a double


% --- Executes during object creation, after setting all properties.
function edit_UnitArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_UnitArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function listbox_SelectView_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_SelectView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_LoadModel.
function pushbutton_LoadModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LoadModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if isfield(handles,'a'), rmfield(handles,'a'); end
[filename, pathname] = uigetfile({'*.cfg', 'CFG file (*.cfg)';...
    '*.*', 'All files (*.*)'},...
    'Select a file');
if isequal(filename,0) || isequal(pathname,0)
    % disp('User pressed cancel')
else
    handles.fileName = fullfile(pathname,filename);
    [coords,aType,Mm,DW,occ,charge] =  readCFG(handles.fileName);
    if isempty(coords)
        msgbox(sprintf('Could not open .cfg file %s!',handles.fileName));
    else
        handles.a = [aType coords DW occ charge];
        handles.Mm = Mm;
    end
    clear coords aType Mm DW occ charge
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_UpdateView.
function pushbutton_UpdateView_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_UpdateView (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'a')
   % msgbox('Please load a model first') 
   return;
end

% Determine what the user wants to see:
handles.displayMode = get(handles.listbox_SelectView,'Value');
handles.atomSize    = 30;
handles.Nedge       = 100;
handles.PSFwidth    = 10*str2double(get(handles.edit_PSFwidth,'String'));  % convert from nm to A
% scale the PSFwidth to produce a Gaussian of FWHM = PSFwidth:
handles.PSFwidth  = handles.PSFwidth /(sqrt(-log(0.5))*2);
handles.PSFwidth  = (pi*handles.PSFwidth)^2;
handles.Nbins       = round(str2double(get(handles.edit_Nbins,'String')));
handles.unitArea    = str2double(get(handles.edit_UnitArea,'String'));
handles.thickness   = str2double(get(handles.edit_thickness,'String'));
E0                  = str2double(get(handles.edit_HighVoltage,'String'));
beta                = str2double(get(handles.edit_beta,'String'));
simpleAverage       = get(handles.checkbox_SimpleAverage,'Value');



coords = handles.a(:,2:4);
aType  = handles.a(:,1);
DW     = handles.a(:,5);
occ    = handles.a(:,6);
charge = handles.a(:,7);
Mm     = handles.Mm;
Natom  = size(handles.a,1); 
% use nomenclature from Malis et al.:
F      = (1+(E0/1022))/((1+(E0/512))^2);

% handles.a(1:10,:)
% Setup arrays for 2D plot generation
Nbins    = handles.Nbins;
qArray2D = zeros(Nbins);
Zcount   = zeros(Nbins);
dx       = Mm(1,1)/(Nbins-1);
dy       = Mm(2,2)/(Nbins-1);

pos = get(handles.axes_Model,'Position');

axes(handles.axes_Model);
cla(handles.axes_Model,'reset')
% Fill the window, depending on the desired display mode:
switch handles.displayMode
    case 1
        scatter(coords(:,1),coords(:,2),handles.atomSize,aType,'filled');
        axis equal; axis tight;
    case 2
        % collect all the charges:
        for j=1:Natom
            ix = round(coords(j,2)/dy)+1;
            iy = round(coords(j,1)/dx)+1;
            if ix > Nbins, ix = Nbins; end
            if iy > Nbins, iy = Nbins; end
            if (ix > 0) && (iy > 0)    
                qArray2D(ix,iy) = qArray2D(ix,iy)+charge(j)*occ(j);
            end
        end
        fprintf('Total charge: %g\n',sum(sum(qArray2D)));

    case 3
        % collect all the Znumbers:
        for j=1:Natom
            ix = round(coords(j,2)/dy)+1;
            iy = round(coords(j,1)/dx)+1;
            if ix > Nbins, ix = Nbins; end
            if iy > Nbins, iy = Nbins; end
            if (ix > 0) && (iy > 0)
                qArray2D(ix,iy) = qArray2D(ix,iy)+aType(j)*occ(j);
            end
        end
    case 4
        % collect all the Znumbers^2:
        for j=1:Natom
            ix = round(coords(j,2)/dy)+1;
            iy = round(coords(j,1)/dx)+1;
            if ix > Nbins, ix = Nbins; end
            if iy > Nbins, iy = Nbins; end
            if (ix > 0) && (iy > 0)
                qArray2D(ix,iy) = qArray2D(ix,iy)+(aType(j)*occ(j))^2;
            end
        end
    case 5
        % collect all the (DW*Znumbers)^2:
        for j=1:Natom
            ix = round(coords(j,2)/dy)+1;
            iy = round(coords(j,1)/dx)+1;
            if ix > Nbins, ix = Nbins; end
            if iy > Nbins, iy = Nbins; end
            if (ix > 0) && (iy > 0)
                qArray2D(ix,iy) = qArray2D(ix,iy)+(aType(j)*DW(j)*occ(j))^2;
            end
        end
    case 6
        % Potential (neutral atoms)
        for j=1:Natom
            ix = round(coords(j,2)/dy)+1;
            iy = round(coords(j,1)/dx)+1;
            if ix > Nbins, ix = Nbins; end
            if iy > Nbins, iy = Nbins; end
            if (ix > 0) && (iy > 0)
                qArray2D(ix,iy) = qArray2D(ix,iy)+handles.V0_neutral(aType(j))*occ(j);
            end
        end
        % normalize by thickness: 
        qArray2D = qArray2D*0.479/Mm(3,3);
    case 7
        % Potential (charged atoms)
        for j=1:Natom
            ix = round(coords(j,2)/dy)+1;
            iy = round(coords(j,1)/dx)+1;
            if ix > Nbins, ix = Nbins; end
            if iy > Nbins, iy = Nbins; end
            if (ix > 0) && (iy > 0)
                qArray2D(ix,iy) = qArray2D(ix,iy)+(handles.V0_neutral(aType(j))+handles.V0_charged(aType(j))*charge(j))*occ(j);
            end
        end
        qArray2D = qArray2D*0.479/Mm(3,3);
	case 8
		% Displacement map
		DisplacementParams(handles);
		return;
    case 9
        % Inelastic mean free path
        for j=1:Natom
            ix = round(coords(j,2)/dy)+1;
            iy = round(coords(j,1)/dx)+1;
            if ix > Nbins, ix = Nbins; end
            if iy > Nbins, iy = Nbins; end
            if (ix > 0) && (iy > 0)
                if simpleAverage
                    % First, compute Z_eff = sum(Z)/length(Z)
                    qArray2D(ix,iy) = qArray2D(ix,iy)+aType(j)*occ(j);
                    Zcount(ix,iy)   = Zcount(ix,iy)+occ(j);
                else
                    % First, compute Z_eff = sum(Z^(1+0.3))/sum(Z^0.3)
                    qArray2D(ix,iy) = qArray2D(ix,iy)+(aType(j)^1.3)*occ(j);
                    Zcount(ix,iy)   = Zcount(ix,iy)+(aType(j)^0.3)*occ(j);
                end
            end
        end
        % Compute the effective Z number in thos eplaces, where it is defined:
        ind = find(Zcount > 0);
        qArray2D(ind) = qArray2D(ind)./Zcount(ind);
        Zmean = sum(qArray2D(ind))/length(ind);
        
        % qArray2D = griddata(ix,iy,qArray2D(ind),ix2d,iy2d);
        % qArray2D = interp2(ix,iy,qArray2D(ind),ix2d,iy2d,'linear',Zmean);
        qArray2D(1:2,1:Nbins) = 0;
        qArray2D(Nbins-1:Nbins,1:Nbins) = 0;
        qArray2D(1:Nbins,Nbins-1:Nbins) = 0;
        qArray2D(1:Nbins,1:2) = 0;
        
        qArray2D(1,1:Nbins) = Zmean;
        qArray2D(Nbins,1:Nbins) = Zmean;
        qArray2D(1:Nbins,Nbins) = Zmean;
        qArray2D(1:Nbins,1) = Zmean;
        figure; imagesc(qArray2D); pause(0.5);
        axes(handles.axes_Model);
        % We now need to interpolate over those pixels which have not been
        % calculated:
        ind = find(qArray2D ~= 0);
        [ix,iy] = ind2sub(size(qArray2D),ind);
        [ix2d,iy2d] = meshgrid([1:Nbins],[1:Nbins]);
        %if (ix > 0) && (iy > 0)
            qArray2D = griddata(ix,iy,qArray2D(ind),ix2d,iy2d,'nearest');
        %end
        if (0)
        % Next, compute Em = 7.6eV*Z_eff^0.36:
        qArray2D = 7.6*qArray2D.^0.36;
        
        % lambda = 106*F*E0/(Em*log(2*beta*E0/Em))
        qArray2D = 106*F*E0./(qArray2D.*log(2*beta*E0./qArray2D));        
        end
    case 9
        % t/lambda map
        for j=1:Natom
            ix = round(coords(j,2)/dy)+1;
            iy = round(coords(j,1)/dx)+1;
            if ix > Nbins, ix = Nbins; end
            if iy > Nbins, iy = Nbins; end
            if (ix > 0) && (iy > 0)    
                qArray2D(ix,iy) = qArray2D(ix,iy)+occ(j);
            end
        end
                
end

% The common display code that converts to 1D plots from 2D arrays:
if handles.displayMode > 1
    % The real space x-axis:
    x = dx*([0:Nbins-1].');
    y = dy*([0:Nbins-1].');

    % Setup arrays for 2D plot generation
    dkx = 1/Mm(1,1);
    dky = 1/Mm(2,2);
    dx  = Mm(1,1)/Nbins;
    dy  = Mm(2,2)/Nbins;
    
    [kx,ky]  = meshgrid(dkx*[-Nbins/2:Nbins/2-1],dky*[-Nbins/2:Nbins/2-1]);
    k2       = fftshift(kx.^2+ky.^2);
    % k2(1,1)

    % Smoothen the 2D array:
    qArray2D = real(ifft2(fft2(qArray2D).*exp(-handles.PSFwidth*k2)));
    % The scaling should also relate the area of a single unit cell with
    % that of the whole array

    % the scale defines how many pixels are in 1 nm^2:
    scale    = 100.0/(dx*dy);  % in 1/nm^2
    scale    = scale*handles.unitArea/handles.thickness;
    qArray2D = scale*qArray2D;
    
    
    
    % Display:
    % figure
    imagesc(0.1*x,0.1*y,qArray2D);
    % set(gca,'YDir','normal','FontSize',14);
    set(gca,'YDir','normal');
    xlabel('x in nm');
    ylabel('y in nm');
    colorbar;
    axis equal; axis tight;  
    % xlim([0.5 0.092*max(x)]);
    % ylim([0.5 0.092*max(y)]);
    switch handles.displayMode
        case 2
            title(sprintf('Average number of electrons per %g nm^2',handles.unitArea));
        case 3
            title(sprintf('Average Z number per %g nm^2',handles.unitArea));
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show the 1D integrated plot:    
    qArray = sum(qArray2D(handles.Nedge+1:Nbins-handles.Nedge,:),1)/(Nbins-2*handles.Nedge);
    axes(handles.axes_Plot);
    plot(x,qArray);
    xlim([x(1) x(end)]);
    
    handles.qArray2D = qArray2D;
end
set(handles.axes_Model,'Position',pos);

posPlot = get(handles.axes_Plot,'TightInset');
posImg = get(handles.axes_Model,'TightInset');

guidata(hObject, handles);






% --- Executes on button press in pushbutton_SaveData.
function pushbutton_SaveData_Callback(hObject, eventdata, handles)

pushbutton_UpdateView_Callback(hObject, eventdata, handles)
% save handles.qArray2D
if (handles.displayMode > 1)
    [filename, pathname] = uiputfile('*.img', 'Save wave function as');
    if isequal(filename,0) return; end

    Nbins    = handles.Nbins;
    Mm     = handles.Mm;
    dx       = Mm(1,1)/(Nbins-1);
    dy       = Mm(2,2)/(Nbins-1);
    
    binwrite2D(handles.qArray2D,fullfile(pathname, filename),dx,dy,handles.thickness,0,0); 
end




function edit_thickness_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thickness as text
%        str2double(get(hObject,'String')) returns contents of edit_thickness as a double


% --- Executes during object creation, after setting all properties.
function edit_thickness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thickness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





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




% --- Executes on button press in checkbox_SimpleAverage.
function checkbox_SimpleAverage_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_SimpleAverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_SimpleAverage





