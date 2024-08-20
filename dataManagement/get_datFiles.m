
function [ datfiles, datfile_sizes, is_symlink ] = get_datFiles( basepath, varargin )
% given basepath, return paths to all .dat files
%
% input:
%       basepath
%       varargin
%           subset :  'cat' - only concatenated ones in basepath folder
%                     'subsession' - only the ones in subsession fodlers (component .dat files)
%                     'all'
% output:
%       datfiles        - paths to .dat files
%       datfile_sizes   - file sizes
%       is_symlink      - logical, true if it's a symlink (data is in cold storage)
%                        CAVEAT: symlinks are identified based on file size, less than 1MB ; this might cause weirdness if there
%                                were a .dat file brief enough to be this small, but I have a hard time imagining this case in practice
%
% 
% NOTE: subsessions are expected to have the basename in their folder name
%       e.g., if basename is Cicero, a valid subsession folder might be
%       Cicero_10302024

p = inputParser;
addParameter( p,'subset','all',@ischar );  % 'all', 'cat', 'subsess'        

parse( p,varargin{:} );

subset = p.Results.subset;
basename = bz_BasenameFromBasepath( basepath );

currdir = pwd; % save current directoryy
cd( basepath )
% get .dat files inside basepath
dat_sessdir = dir( '*.dat' );
dat_sessdir_path = cellfun(@(x) fullfile( basepath, x ), { dat_sessdir.name } , 'UniformOutput', false)';

if strcmp(subset, 'cat')
    datfiles = dat_sessdir_path ;
    datfile_sizes = get_filesizes( datfiles );
    cd( currdir ) % switch to original path
    is_symlink = datfile_sizes < 1000;
    return
end

content = dir( basepath );
dir_logical = [ content.isdir ];
dir_name = {content.name}; dir_name( ~dir_logical ) = [];
dir_name( 1:2 ) = [];
% look for folders with basename included in name
subsess_indicator = cellfun(@(x) any(regexp(x, basename)), dir_name );
subsess_dir_name = dir_name( subsess_indicator );

% store all .dat files within a basename
dat_subssDir_path = [];
for kp = 1:length( subsess_dir_name )

    cd( subsess_dir_name{kp} )
    dat_subssDir = dir( '*.dat' );
    dat_subssDir_path = [ dat_subssDir_path ; ...
            cellfun(@(x) fullfile( basepath, subsess_dir_name{kp}, x ), { dat_subssDir.name } , 'UniformOutput', false)' ];
    cd ..
end

if strcmp(subset, 'subsession')
    datfiles = dat_subssDir_path ;
elseif strcmp(subset, 'all')
    datfiles = [ dat_sessdir_path ; dat_subssDir_path ];
end
datfile_sizes = get_filesizes( datfiles );
is_symlink = datfile_sizes < 1000;
cd( currdir ) % switch to original path

end


