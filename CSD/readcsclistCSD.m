function cscData = readcsclistCSD(cscList,cscDirs)
% READCSCLIST   Read list of CSC files and store raw data as a structure array
%
% CSCDATA = READCSCLIST(CSCLIST,CSCDIRS) returns an MxN structure array of
% raw data, where N is the number of listed CSC files, M is the number of
% directories in which these CSC files are stores, and CSCDATA(I,J) is the
% raw data corresponding to the Jth listed tetrode in the Ith directory.
%           CSCLIST can be a cell array of CSC file names or a double array 
%           of tetrode numbers
%               e.g. 
%                   CSCLIST = {'CSC5.ncs','CSC7.ncs'}
%                       and
%                   CSCLIST = [5 7] 
%               will return precisely the same output.
%
%           CSCDIRS can be either a cell array of directories where the CSC
%           files can be found or a character string specifying a single
%           directory.
%
% Each array element of CSCDATA contains the following data fields:
%
%       - TT
%         Tetrode number (e.g. if the file name is 'CSC12.ncs', then TT is 12)
%
%       - SAMPLE
%         A 512xR array containg the sample trace, where R is the
%         number of records and 512 is the number of samples per record.
%         SAMPLE is returned in units of bits.
%
%       - TIMESTAMPS
%         An R vector containing the timestamp associated with the first
%         sample in each of the R records.
%
%       - FS
%         An R vector of sampling frequencies for each of the R records
%
%       - NVALIDSAMP
%         An R vector containing the number of valid samples for that
%         record. NOTE: THE NUMBER OF SAMPLES IS STILL 512
%         
%       - CH
%         An R vector containing the channel number on which each record
%         was recorded
%
%       - BITS2VOLTS
%         The conversion factor from bits to volts, extracted from the CSC
%         file header
%
%       - LOCATION
%         The path to the location of the CSC file (derived from
%         CSCFILEPATH)
%
%       - FILE
%         The name of the CSC file (derived from CSCFILEPATH)
%         
% See also: GETCSCDATA, FORMATCSC, READCSCFILE, NLX2MATCSC_V3

if isnumeric(cscList)
    cscList = arrayfun(@(x)['CSC' num2str(x) '.ncs'],cscList,'uniformoutput',0);
elseif ischar(cscList)
    cscList = {cscList};
end

nFiles = length(cscList);
nDirs = length(cscDirs);
for iD = 1:nDirs
    for iF = 1:nFiles
        cscFilePath = [cscDirs{iD} filesep cscList{iF}];
        cscData(iD,iF) = readcscfile(cscFilePath);
    end
end