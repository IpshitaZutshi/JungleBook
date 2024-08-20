function filesizes = get_filesizes( filepaths )
% get file sizes (in bytes) of files 
% filepaths - paths to files of interest (char or cell array)
%
% a single path may be passed as a string, multiple paths must be passed as a cell array 

if ischar( filepaths )
    filepaths = {filepaths};
end

filesizes = nan( length(filepaths), 1 );
for kp = 1:length( filepaths )
    s = dir( filepaths{ kp } );
    filesizes( kp ) = s.bytes;
end

end