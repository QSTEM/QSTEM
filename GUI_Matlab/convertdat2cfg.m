function [data2,data1] = convertdat2cfg(fileName,shift,writeFlag);

if nargin < 3
    writeFlag = 1;
end

if nargin <1
    fileName = '/home/koch/Si3N4/Si3N4_12CaO-SiO2_1nm.dat';
end
if nargin < 2
    shift = [0 0 0.62];
end


cfgName = fileName;
ln = length(fileName);
cfgName(ln-2:ln) = 'cfg';

fid = fopen(fileName,'r','ieee-be.l64');
fseek(fid,0,-1);
[s,count1] = fread(fid,2,'int32');
N = s(2);
[data1,count1] = fread(fid,N,'int32');
% fseek(fid,hex2dec('0000f810')+8,-1);
% [data2,count2] = fread(fid,(hex2dec('0006c890')-hex2dec('0000f818'))/8,'float64');
[data2,count2] = fread(fid,[3 N],'float64');
[cell,count2] = fread(fid,3,'float64');
fprintf('N=%d, first and last integer: %d .. %d\n',N,s(1),fread(fid,1,'int32'));
fclose(fid);

% assume that the units are in cm!!!
cell = cell*1e8;
data2 = data2*1e8;

fprintf('positions: %g..%g, %g..%g, %g..%g\n',...
    min(data2(1,:))/cell(1),max(data2(1,:))/cell(1),...
    min(data2(2,:))/cell(2),max(data2(2,:))/cell(2),...
    min(data2(3,:))/cell(3),max(data2(3,:))/cell(3));
    
fprintf('cell dimensions: %g %g %g\n',cell(1),cell(2),cell(3));
% N2 = 400; scatter3(data2(1,1:N2),data2(2,1:N2),data2(3,1:N2),repmat(30,1,N2),data1(1:N2),'filled');
Mm = [cell(1)       0         0
           0  cell(2)         0
           0        0   cell(3)];
       
if writeFlag
    writeCFG(cfgName,Mm,N,data1,data2.',shift);
    fprintf('wrote cfg file %s\n',cfgName);
end
% cut away the amorphous layer in ae using command 'B' and then 0 0 -1 0 0 0.8 !!!
