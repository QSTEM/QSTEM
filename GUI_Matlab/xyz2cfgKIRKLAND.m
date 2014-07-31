function xyz2cfgKIRKLAND(filename)               % filename without extension

% written by Yaron Kauffmann - 01.11.08
% reads xyz files containing atomic positions in cartesian coordinates and converts to cfg format 
% readable by QSTEM (http://www.mf.mpg.de/en/organisation/hsm/koch/stem/index.html).
% the xyz file format should be as follows:
% line 1 : the number of atoms (natom)
% line 2 : empty or one integer number (no text)
% line 3 ... : atom type and x y z positions in cartesian coordinates
% SORTCELL function was downloaded from http://www.mathworks.de/matlabcentral/fileexchange/13770.
% writeCFG function was provided by Christoph Koch (developer of QSTEM).

% Robert A. McLeod comments:
% This doesn't seem to work for Kirkland style data, which has the
% style:
% 
%   Title
%   373.4902 373.4902 373.4902 -- unit cell size
%   31 225.37 197.95 183.38 1 80 - Z-number X Y Z occupancy unknown D-W factor
%   etc. 
%   -1 (termination of file)
%
% So I need to re-write this to work for these wonky .xyz files


    % Search for .xyz file extension.  If not present, append it.
    location = regexp( filename, '\.xyz', 'ONCE' );
    if( isempty(location) )
        filename_ext = [filename, '.xyz' ];
    end
    disp( ['Attempting to open : ', filename_ext ] )
    fid = fopen( filename_ext, 'r' );
  
    %% Start of RAM modification
    % TO DO RAM: modify code to automatically tell the difference between
    % normal and Kirkland .XYZ files
    % Read title
    xyztitle = textscan( fid, '%s', 1 ); 
    disp( [ 'Converting Kirkland-style ', filename, '.XYZ to QSTEM .CFG file for: ', xyztitle ] )
    % Read in dimensions of the Kirland cell and convert from cell to
    % matrix
    celldim = cell2mat( textscan( fid, '%f %f %f', 1 ) );
    % Read in all remaining lines with the proper format
    atominfo = textscan( fid, '%d %f %f %f %f %f' ); % Z-number X Y Z unknown(atom count?) unknown
  
    fclose(fid);
    
    % Convert from cell arrays to normal arrays
    aType = atominfo{1}; 
    aType = aType / 2; % This is atomic number * 2 for some reason in Kirkland
    aType(end) = []; % There is a trailing -1, delete it.
    natom = numel(aType); % Find number of atoms
    
    % Find atomic coordinates
    coords = [atominfo{2} atominfo{3} atominfo{4}];
    coords(end,:) = []; % Trailing NaNs associated with -1 terminator
    
    
    % Apparently we want to sort by atom type
    [aType, sortindex] = sort( aType );
    coords = coords(sortindex,:);  
  
    Mm = diag( celldim );
	% Mm=diag(max(coords)-min(coords)) ;
    % coords=coords-repmat(min(coords),natom,1);

    % aType is already in the form of Z-number so no need to modify for Kirkland files 
    % Call C. Koch's function
    writeCFG(filename_ext(1:end-4),Mm,natom,aType,coords, 0)
  
    
end

  
% function writeCFG(fileName,Mm,N,aType,coords,shift,mode)
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
function writeCFG(filename,Mm,N,aType,coords,shift,mode)

    mode=0;
    if nargin < 6
        shift = [0 0 0];
    end
    if nargin < 7
        mode = 1;
    end

    Mminv = inv(Mm);
    atomType = 0;
    % RAM: is this 'name' supposed to be a cell array?
    name = {'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl',...
        'Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
        'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
        'I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
        'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
        'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'};
    mass = 2.*(1:length(name));

    % RAM: Output filename to command-line
    disp( ['Writing config file to : ', num2str( filename ), '.cfg'] )
    fid = fopen([filename,'.cfg'],'w');

    fprintf(fid,'Number of particles = %d\n',N);
    fprintf(fid,'A = 1.0 Angstrom (basic length-scale)\n');
    for ix=1:3
        for iy=1:3
            fprintf(fid,'H0(%d,%d) = %g A\n',ix,iy,Mm(ix,iy));
        end
    end
    fprintf(fid,'.NO_VELOCITY.\n');
    fprintf(fid,'entry_count = 3\n');

    for ix=1:N
        at = aType(ix);
        if at ~= atomType
            % RAM: This fprintf didn't give the correct string...
            fprintf(fid,'%g\n%s\n',mass(at),name{at});  % specify type and mass of first atom (all others are the same)
            atomType = at;
        end
        if mode == 1
            % RAM: someone forced the mode = 0 above, so shift does nothing
    %        fprintf(fid,'%g %g %g\n',mod(coords(ix,:)*Mminv+shift,1));
            fprintf(fid,'%g %g %g\n',coords(ix,:)*Mminv+shift);
        else
            fprintf(fid,'%g %g %g\n',coords(ix,:));        
        end    
    end
    fclose(fid);
    disp( 'Finished' );
end