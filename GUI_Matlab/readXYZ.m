function [coord,aType,Mm] = readXYZ( filename )

    % Robert A. McLeod
    % 14 April 2014
    % for reading in proteins from Gromacs

    %% Open file
    location = regexp( filename, '\.xyz', 'ONCE' );
    if( isempty(location) )
        filename_ext = [filename, '.xyz' ];
    else
        filename_ext = filename;
    end
    disp( ['Attempting to open : ', filename_ext ] )
    fHand = fopen( filename_ext, 'r' );
    
    if( fHand == -1 )
        error( [ 'Error, ', filename_ext', ' not found.' ] )
    end
    
    %% Parse .XYZ file
    elementName = {'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl',...
        'Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
        'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',...
        'I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
        'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',...
        'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'};
    
    % Read atom count
    nAtom = str2double( fgetl( fHand ) );
    % Read title
    xyztitle = fgetl( fHand );
    
    data = textscan(fHand,'%s %f32 %f32 %f32', nAtom);
    fclose(fHand); 

    % Pull off atom names
    atom_names = data{1};
    % Convert to z-number
    aType = zeros( [nAtom 1] );
    for J = 1:nAtom
        aType(J) = find(strcmpi( atom_names{J}, elementName) );
    end

    % Find coords
    coord = [ data{2}, data{3}, data{4} ];

    % Sort according to atom type
    [aType, aIndex] = sort( aType );
    coord = coord(aIndex,:);

    % Apparently we want to sort by atom type
    [aType, sortindex] = sort( aType );
    coord = coord(sortindex,:);  
  
    celldim = max(coord, [], 1 ) - min( coord, [], 1 );
    Mm = diag( celldim );
end