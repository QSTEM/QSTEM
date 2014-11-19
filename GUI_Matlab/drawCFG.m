function drawcfg(fileName,flag)
% Author: Christoph T. Koch
% This function will read and display the content of a .cfg file (atom
% positions).


% If no fileName is supplied by the user, open a file selection box:
if nargin < 1
   [filename, pathname, filterindex] = uigetfile('*.cfg', 'Select CFG file'); 
   fileName = [pathname filename];
end

% if no second argument is supplied, make only show one unit cell
if nargin < 2
    flag = 1;
end

% Read the CFG file:
[coords, aType, Mm, ~, ~, ~] = readCFG(fileName);

% extract number of atoms:
Natom = size(coords,1);
% fprintf('There are %d atoms in this structure\n',Natom);


if flag
    % add atoms at boundary in y-direction:
    Xtal2 = [aType coords];
    clear coords aType
    
    % Xtal2 = [Xtal2; Xtal2+repmat([0 (Mm*[0;1;0]).'],Natom,1); Xtal2+repmat([0 (Mm*[0;-1;0]).'],Natom,1)];
    Xtal2 = [Xtal2; Xtal2+repmat([0 0 1 0],Natom,1); Xtal2+repmat([0 0 -1 0],Natom,1)];
    ind = find(sum(Xtal2(:,2:4) > 1.01,2) == 0);
    Xtal2 = Xtal2(ind,:);
    ind = find(sum(Xtal2(:,2:4) < 0,2) == 0);
    Xtal2 = Xtal2(ind,:);
    Natom = size(Xtal2,1);
    
    % add atoms at boundary in x-direction:
    % Xtal2 = [Xtal2; Xtal2+repmat([0 (Mm*[1;0;0]).'],Natom,1); Xtal2+repmat([0 (Mm*[-1;0;0]).'],Natom,1)];
    Xtal2 = [Xtal2; Xtal2+repmat([0 1 0 0],Natom,1); Xtal2+repmat([0 -1 0 0],Natom,1)];
    ind = find(sum(Xtal2(:,2:4) > 1.01,2) == 0);
    Xtal2 = Xtal2(ind,:);
    ind = find(sum(Xtal2(:,2:4) < 0,2) == 0);
    Xtal2 = Xtal2(ind,:);
    Natom = size(Xtal2,1);

    
    % add atoms at boundary in z-direction:
    %Xtal2 = [Xtal2; Xtal2+repmat([0 (Mm*[0;0;1]).'],Natom,1); Xtal2+repmat([0 (Mm*[0;0;-1]).'],Natom,1)];
    Xtal2 = [Xtal2; Xtal2+repmat([0 0 0 1],Natom,1); Xtal2+repmat([0 0 0 -1],Natom,1)];
    ind = find(sum(Xtal2(:,2:4) > 1.01,2) == 0);
    Xtal2 = Xtal2(ind,:);
    ind = find(sum(Xtal2(:,2:4) < 0,2) == 0);
    Xtal2 = Xtal2(ind,:);
    Natom = size(Xtal2,1);
    
    % show the structure:
    Xtal2(:,2:4) = Xtal2(:,2:4)*(Mm.');
    drawXtal(Xtal2,Mm);
else
    drawXtal([aType coords],Mm);    
end
