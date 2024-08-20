function totaldirsize = compute_directorySizes( targetpath,  varargin )
% compute the size of each directory in targetpath
% 
% recurses through entire subdirectory tree in target path...
% outputs are saved in a text file in the targetpath- each line holds the
% directory name and the size it occupies on disk
% 
% input:
%   targetpath       - path with directories to assess
% output: 
%   totaldirsize     - total directory size (in units specified in the input)
%
% if directories are large- as expected on NYU share- this will take a to while

p = inputParser;
addParameter( p,'outputfile','file_sizes.txt',@ischar );  % name of file where to store data sizes
addParameter( p,'filesize_unit','GB',@ischar );           % file size unit- GB or MB

parse( p,varargin{:} );

outputfile = p.Results.outputfile;
filesize_unit = p.Results.filesize_unit;
if ~ischar( targetpath )
    error('''targetpath'' must be a string, and valid path')
end

clc
if strcmp( filesize_unit, 'GB' )
    denom = 1e9;
elseif strcmp( filesize_unit, 'MB' )
    denom = 1e6;
end

%% get all folders in desired path

cd( targetpath )
content = dir;
content_filenames = { content.name };

% directories
dirs_logical = [ content.isdir ];
dirs_names = content_filenames( dirs_logical );
dirs_names( 1:2 ) = [];   % excludes '.' and '..', the current and parent directory 


ndirs = length( dirs_names );

diary( outputfile )
% loop through all directories
for d = 1:ndirs
    
    % progress so you know where things stand... 
    fprintf('directory %d/%d...\n',d, length( dirs_names ) )
    cd( dirs_names{d} )
    % call to the powershell- this does the heavy lifting...
    system('powershell -command "Get-ChildItem -Recurse | Measure-Object -Sum Length"');
    cd ..

end

diary off
fprintf('done ! \n')


%%

dirsizes = nan( ndirs, 1 );

counter = 1;
fid = fopen( outputfile );
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    % line with file size on it
    if length(tline) >=3 && strcmp( tline(1:3), 'Sum' )
        % parse string that holds the filesize
        cmp = strsplit(tline, ' ');
        dirsizes( counter ) = str2double( cmp{ 3 } ) ./ denom;
        counter = counter+1;
    end
end
fclose(fid);
% write contents to disk
fid = fopen( outputfile,'w' );
for ln = 1:ndirs
    
    line = sprintf('%s\t\t%.2f %s\n', dirs_names{ln}, dirsizes(ln), filesize_unit);
    fwrite( fid,line );
end
totaldirsize = sum( dirsizes );
line = sprintf('Title size on disk: %.2f %s', totaldirsize, filesize_unit );
fwrite( fid,line );
fclose(fid);

end