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
function [img,t,dx,dy] = binread2D(fileName,printFlag,flag)

headerlength = 4;
if nargin < 2
    printFlag = 1;
end
% open the file and define the file ID (fid):
fid=fopen(fileName,'rb','ieee-le');

% header = [headersize(bytes) paramSize commentSize Nx Ny complFlag doubleFlag dataSize version]

header = fread(fid,8,'int32');
Nx = header(5);
Ny = header(4);
t = fread(fid,1,'float64');
dx = fread(fid,1,'float64');
dy = fread(fid,1,'float64');
% header size: 8*4+3*8=

% read additional parameters from file, if any exist:
paramSize = header(2);
if (paramSize > 0)
    params = fread(fid,paramSize,'float64');
    params
end

% read comments from file, if any exist:
commentSize = header(3);
if (commentSize > 0)
    comment = fread(fid,commentSize,'char');
    % fprintf('Comment: %s\n',comment);
end

    % flag = 0;

integerFlag = 0;
complexFlag = header(6);
doubleFlag = (header(7) == 8*(complexFlag+1));
if (nargin <3)
    flag = integerFlag+doubleFlag*4+complexFlag*2;    
end

% fprintf('Header Size: %d bytes, t: %g, dx: %g, dy: %g\n',header(1),t,dx,dy);


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

% img = img.';

