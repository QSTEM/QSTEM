% function binwrite2D(data,fileName,dx,dy,t,doubleFlag,intFlag)
% this function will write 2-dimensional binary data to the file specified
% by fileName
% 
function binwrite2D(data,fileName,dx,dy,t,doubleFlag,intFlag)

if nargin < 5
    t = 0;
end

if nargin < 7
    intFlag = 0;
end
if nargin < 6
   doubleFlag = 0;
end

if (0) % old version
    headerlength = 4;
    % data = data.';

    % fid=fopen('c:\my documents\physics\thesis\stem\zoom10_mag64k_ROI2_box6.bin','r');
    fid=fopen(fileName,'w');
    [Nx,Ny] = size(data);
    
    cFlag = 2*(1-isreal(data));
    flag = intFlag+cFlag+4*doubleFlag;
    
    
    fwrite(fid,[Ny Nx],'int32');
    fwrite(fid,t,'float32');
    fwrite(fid,flag,'int32');
    fprintf('binwrite2D %s: %d x %d pixels (',fileName,Nx,Ny);
else
    format compact
    % header = [headersize(bytes) paramSize commentSize Nx Ny complFlag doubleFlag dataSize version]
    if nargin < 3
        dx = 1;
        dy = 1;
    end
    [Nx,Ny] = size(data);
    complexFlag = 1-isreal(data);
    dataSize = (doubleFlag+1)*4*(complexFlag+1);

    fid=fopen(fileName,'w+');
    if (fid < 0)
       error('could not write to file %s\n',fileName);
       return
    end
        
        
    header = [4*8+3*8 0 0 Ny Nx complexFlag dataSize 1];
    fwrite(fid,header,'int32');
    fwrite(fid,t,'float64');
    fwrite(fid,dx,'float64');
    fwrite(fid,dy,'float64');
    flag = intFlag+2*complexFlag+4*doubleFlag;
end
switch bitand(flag,7)
case 0
    % complexFlag = 0;
    % doubleFlag = 0;
    fprintf('(32-bit real data)\n');
    fwrite(fid,data,'float32');
case 4
    % complexFlag = 0;
    % doubleFlag = 1;
    fprintf('64-bit real data, %.3fMB)\n',Nx*Ny*8/1048576);
    fwrite(fid,data,'float64');
case 2
    % complexFlag = 1;
    % doubleFlag = 0;
    fprintf('32-bit complex data)\n');
    % img(1:Nx,1:2:2*Ny-1) = real(data);
    % img(1:Nx,2:2:2*Ny) = imag(data);
    img(1:2:2*Nx-1,1:Ny) = real(data);
    img(2:2:2*Nx,1:Ny)   = imag(data);
    % img = fread(fid,[2*Nx,Ny],'float32');
    % img = img(1:2:2*Nx-1,1:Ny)+i*img(2:2:2*Nx,1:Ny);

    fwrite(fid,img,'float32');
case 6
    % complexFlag = 1;
    % doubleFlag = 1;
    fprintf('64-bit complex data)\n');
    img(1:Nx,1:2:2*Ny-1) = real(data);
    img(1:Nx,2:2:2*Ny) = imag(data);
    fwrite(fid,img,'float64');
case 1
    % complexFlag = 0;
    % doubleFlag = 0;
    fprintf('16-bit integer data, %.3fMB)\n',Nx*Ny*2/1048576);
    fwrite(fid,data,'int16');
case 5
    % complexFlag = 0;
    % doubleFlag = 1;
    fprintf('32-bit integer data)\n');
    fwrite(fid,data,'int32');    
end
% imagesc(img); colormap('gray'); axis equal; axis tight;

fclose(fid);


