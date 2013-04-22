% function writeCFG(fileName,Mm,aType,coords,shift,mode,DW,charge)
%
% fileName:   name of cfg file
% Mm:         metric matrix (for orthogonal unit cells or super structures
%             Mm = diag([ax by cz]), where ax, by, cz are the size of your
%             unit cell / super structure in A.
% aType:      column vector of Z numbers
% coords:     N x 3 matrix of atomic positions [x1 y1 z1; x2 y2 z2; ...]
% shift:      can be used to define an offset to atomic positions (default
%             is [0 0 0]
% mode:       0 = atom positions in coords are fractional coordinates
%             1 = atom positions in coords are cartesian coordinates
% DW:         Debye-Waller factor: this can be a single number or a vector
%             with a separate value for each atom
% charge:     vector with value of charge per atom
function writeCFG(fileName,Mm,aType,coords,shift,mode,DW,charge)

if nargin < 5
    shift = [0 0 0];
end
if nargin < 6
    mode = 1;
end
N = size(coords,1);
Mminv = inv(Mm);
atomType = 0;
name = ['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl',...
    'Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
    'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
    'I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
    'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
    'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'];
mass = 2*[1:length(name)];

fid = fopen(fileName,'w');

fprintf(fid,'Number of particles = %d\n',N);
fprintf(fid,'A = 1.0 Angstrom (basic length-scale)\n');
for ix=1:3
    for iy=1:3
        fprintf(fid,'H0(%d,%d) = %g A\n',ix,iy,Mm(ix,iy));
    end
end
fprintf(fid,'.NO_VELOCITY.\n');
if nargin < 8
    fprintf(fid,'entry_count = 3\n');
end
if nargin == 8
    fprintf(fid,'entry_count = 5\n');
end
if nargin > 8
    fprintf(fid,'entry_count = 6\n');
end

for ix=1:N
    at = aType(ix);
    if at ~= atomType
        fprintf(fid,'%g\n%s\n',mass(at),name(2*at-1:2*at));  % specify type and mass of first atom (all others are the same)
        atomType = at;
    end
    if mode == 1
        fprintf(fid,'%.5f %.5f %.5f',mod(coords(ix,:)*Mminv+shift,1));
    else
        fprintf(fid,'%.5f %.5f %.5f',coords(ix,:));        
    end    
    if nargin > 7
        if length(DW) > 1
            fprintf(fid,' %.2f 1',DW(ix));
        else
            fprintf(fid,' %.2f 1',DW);            
        end
    end
    if nargin > 8
        fprintf(fid,' %.2f',charge(ix));                
    end
    fprintf(fid,'\n');
end
fclose(fid);