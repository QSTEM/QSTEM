% function writeCFG(fileName,Mm,aType,coords,shift,mode,DW,occ,charge)
%
% fileName:   name of cfg file
% Mm:         metric matrix (for orthogonal unit cells or super structures
%             Mm = diag([ax by cz]), where ax, by, cz are the size of your
%             unit cell / super structure in A.
% aType:      column vector of Z numbers
% coords:     N x 3 matrix of atomic positions [x1 y1 z1; x2 y2 z2; ...]
% shift:      can be used to define an offset to atomic positions (default
%             is [0 0 0].  Will only be used of the atomic coordinates are
%             cartesian.
% mode:       0 = atom positions in coords are fractional coordinates
%             1 = atom positions in coords are cartesian coordinates
% DW:         Debye-Waller factor: this can be a single number or a vector
%             with a separate value for each atom
% occ:        vector with value of occupancy for each
% charge:     vector with value of charge per atom
function success = writeCFG(fileName,Mm,aType,coords,shift,mode,DW,occ,charge)

if nargin < 5
    shift = [0 0 0];
end
if nargin < 6
    mode = 1;
end
success = 1;
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
if (fid < 0)
	success = 0;
	return;
end

fprintf(fid,'Number of particles = %d\n',N);
fprintf(fid,'A = 1.0 Angstrom (basic length-scale)\n');
for ix=1:3
    for iy=1:3
        fprintf(fid,'H0(%d,%d) = %10.3f A\n',ix,iy,Mm(iy,ix));
    end
end
fprintf(fid,'.NO_VELOCITY.\n');
if nargin < 7               % 6 or fewer arguments means no DW, occ, or charge
    fprintf(fid,'entry_count = 3\n');
end
if nargin == 7              % 7 args means just DW included
    fprintf(fid,'entry_count = 4\n');
end
if nargin == 8              % 8 args means DW and occ included
    fprintf(fid,'entry_count = 5\n');
end
if nargin > 8               % 9 args means DW, occ, and charge all included
    fprintf(fid,'entry_count = 6\n');
end

for ix=1:N
    at = aType(ix);
    if at ~= atomType
        fprintf(fid,'%g\n%s\n',mass(at),name(2*at-1:2*at));  % specify type and mass of first atom (all others are the same)
        atomType = at;
    end
    if mode == 1
        fprintf(fid,'%.6f %.6f %.6f',mod((coords(ix,:)+shift)*Mminv,1));
    else
        fprintf(fid,'%.6f %.6f %.6f',coords(ix,:));        
    end    
    if nargin > 6           % nargs is at least 7, so write DW
        if length(DW) > 1
            fprintf(fid,' %.4f',DW(ix));
        else
            fprintf(fid,' %.4f',DW);            
        end
    end
    if nargin > 7           % nargs is at least 8, so write occ, as well
        fprintf(fid,' %.4f',occ(ix));
    end
    if nargin > 8           % nargs is at least 9, so write charge, too
        fprintf(fid,' %.4f',charge(ix));                
    end
    fprintf(fid,'\n');
end
fclose(fid);