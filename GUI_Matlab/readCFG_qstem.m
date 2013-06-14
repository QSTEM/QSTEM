% function [coords, aType, Mm, DW, charge] = readCFG(fileName)
%
% fileName: name of input file
% shift:    offset to atomic positions
% mode:     1=cartesian coordinates, 0 = fractional coordinates
function [coords, aType, Mm, DW, occ, charge] = readCFG(fileName,shift,mode)

if nargin < 2
    shift = [0 0 0];
end
if nargin < 3
    mode = 1;
end

coords = [];
aType = [];
DW = [];
occ = [];
charge = [];
Mm = eye(3);

if (1)
    clear readCFG_mex
    [pathname,filen,ext] = fileparts(fileName);
    filename2 = [filen,ext];
    oldPath = pwd();
    cd(pathname);
    [a,Mm] = readCFG_mex(filename2);
    if sum(sum(isnan(a))) > 0 || sum(sum(isnan(Mm))) > 0
        clear readCFG_mex
        fprintf('reading file a second time!\n');
        [a,Mm] = readCFG_mex(filename2);
    end
    
    Mm = Mm.';
    
    cd(oldPath);
    if prod(size(a)) == 1
        msgbox(sprintf('Could not open CFG file %s\n',fileName));  
        return;
    end
    
    
    
    box = [0 0 0;1 0 0; 1 1 0; 0 1 0;0 0 1;1 0 1; 1 1 1; 0 1 1];
    box = (Mm*(box.')).';
    NcellPotOffset = min(box,[],1); % +[handles.PotentialOffsetX handles.PotentialOffsetY 0];

    
    coords = double(a(:,1:3))+repmat(NcellPotOffset,size(a,1),1);    
    aType  = double(a(:,4));
    DW     = double(a(:,5));
    occ    = double(a(:,6));
    charge = double(a(:,7));
% fprintf('Minimum in z-direction: %f\n',min(coords(:,3)));
    if mode ~= 1
        % coords = coords*(inv(Mm).');
        coords = (inv(Mm)*(coords.')).';
    end
    
    clear a   
%    fprintf('Minimum in z-direction: %f\n',min(coords(:,3)));
else    
    fprintf('reading the slow way ...\n');
    atomType = 0;
    name = ['Si','O ','Ca','N '];
    mass = [28,16,40,14];
    
    fid = fopen(fileName,'r');
    if fid == -1
        msgbox(sprintf('Could not open CFG file %s\n',fileName));
        return
    end
    buf = fgets(fid);
    N = str2num(buf(findstr(buf,'=')+1:end));
    fgets(fid);  % size unit
    
    for ix=1:3
        for iy=1:3
            buf = fgets(fid);
            buf = buf(findstr(buf,'=')+2:end);
            i2 = findstr(buf,' ');
            if isempty(i2)
                Mm(ix,iy) = str2num(buf);
            else
                Mm(ix,iy) = str2num(buf(1:i2(1)));
            end
        end
    end
    fgets(fid);
    buf = fgets(fid);
    entryCount = str2num(buf(findstr(buf,'=')+1:end));
    coords     = zeros(N,3);
    aType      = zeros(N,1);
    DW         = zeros(N,1);
    occ        = zeros(N,1);
    charge     = zeros(N,1);
    
    
    for ix=1:N
        buf = fgets(fid);
        num1 = str2num(buf);
        if length(num1) < entryCount
            atom = round(num1(1)/2);
            fgets(fid);
            buf = fgets(fid);
            num1 = str2num(buf);
        end
        aType(ix) = atom;
        % coords(ix,1:3) = (Mm*num1(1:3)')';
        % coords(ix,1:3) = (num1(1:3)*(Mm'));
        coords(ix,1:3) = mod(num1(1:3),1);
        if length(num1) > 3
            DW(ix)         = num1(4);
        end
    end
    fclose(fid);
    Mm = Mm.';
    if mode == 1
        coords = coords*(Mm);
    end
    
end
    