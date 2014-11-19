function [coord, aType, Mm, DW, occ, charge] = readCFG( filename )

    % Robert A. McLeod
    % 14 April 2014
    % Last updated: 19 November 2014 by Joshua A. Taillon
    % This is designed to replace readCFG_qstem, which likes to break.

    %% Open file
    location = regexp( filename, '\.cfg', 'ONCE' );
    if( isempty(location) )
        filename_ext = [filename, '.cfg' ];
    else
        filename_ext = filename;
    end
    disp( ['Attempting to open : ', filename_ext ] )
    fHand= fopen( filename_ext, 'r' );

    if( fHand == -1 )
        error( [ 'Error, ', filename_ext', ' not found.' ] )
    end

    %% Parse .CFG header
    % first line, number of atoms
    buffer = fgetl( fHand );
    nAtom = str2double( buffer( regexp( buffer, '\d' ) ) );

    % second line, length scale (ignored)
    fgetl( fHand );

    % 3rd-11th lines, the Mm matrix
    Mm = zeros( [3 3] );
    for J = 1:9
        buffer = fgetl( fHand );

        [startM, endM] = regexp( buffer, '[-+]?\d*\.?\d*', 'start', 'end' );
        % We are only interested in the last match, the actual number (and not
        % the matrix indices
        Mm(J) = str2double( buffer( startM(end):endM(end) )  );
    end

    fgetl( fHand ); % No velocity
    buffer = fgetl( fHand ); % entry_count (tells us what factors to use)
    entries = str2double( buffer( regexp( buffer, '\d' ) ) );
    
    %% End of .CFG header, start of body
    % Preallocate arrays based on number of entries
    coord = zeros( [nAtom 3] );
    aType = zeros( [nAtom 1] );
    charge = [];
    occ = [];
    DW = [];
    if ( entries > 5)
        charge = zeros( [nAtom 1] );
        occ = zeros( [nAtom 1] );
        DW = zeros( [nAtom 1] );
    elseif (entries > 4)
        occ = zeros( [nAtom 1] );
        DW = zeros( [nAtom 1] );
    elseif (entries > 3)
        DW = zeros( [nAtom 1] );
    end

    % Try to scan a single number (Z*2) followed by a end-of-line, followed
    % by a character
    atomZ = textscan( fHand, '%d' );
    atomZ = atomZ{1}/2;
    % Check for break
    %     if( isempty( newAtom ) || isnan(newAtom) )
    %         break;
    %     end
    fgetl( fHand ); % next line, string, discard
    atomIndex = 1;
    while( atomIndex <= nAtom )
        coordScan = textscan( fHand, '%f%f%f%f%f%f', 'CollectOutput', 1 ); % Scan until it hits the next non-match (i.e. a new atom)
        coordScan = coordScan{1};
        coordScan = coordScan(:, 1:entries);
        numNewAtoms = size( coordScan,1 );

        if( numNewAtoms <= 0 )
            break;
        end

        if( isnan( coordScan(end,2) ) )
            % Hit a new atom type
            numNewAtoms = numNewAtoms - 1;
%             disp( [ '1: Found ', num2str(numNewAtoms), ' of atomic number Z = ', num2str(atomZ) ] )

            nextAtomZ = coordScan( end, 1 )/2;
            coordScan = coordScan(1:end-1,:);
            

            fgetl( fHand ); % next line, string, discard
        else
%             disp( [ '2: Found ', num2str(numNewAtoms), ' atoms of atomic number Z = ', num2str(atomZ) ] )
        end

        % Add to list
        aType(atomIndex:atomIndex+numNewAtoms-1) = atomZ;
        coord(atomIndex:atomIndex+numNewAtoms-1,:) = coordScan(:,1:3);
        if (entries > 3)
            DW(atomIndex:atomIndex+numNewAtoms-1) = coordScan(:,4);
        end
        if (entries > 4)
            occ(atomIndex:atomIndex+numNewAtoms-1) = coordScan(:,5);
        end
        if (entries > 5)
            charge(atomIndex:atomIndex+numNewAtoms-1) = coordScan(:,6);
        end
        atomIndex = atomIndex + numNewAtoms;

        atomZ = nextAtomZ;
    end
    disp(['Read in ' num2str(nAtom) ' atoms from ' filename_ext '.'])
    fclose( fHand );
end