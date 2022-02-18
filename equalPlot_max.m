

function rateMaps = equalPlot_max(positions,tSp)

sLength = 100; % Side length in cm
bins = 50; % Number of bins
h = 2; % Smoothing factor when calculating the ratemap
binWidth = sLength/bins;
segment = [0,inf];
suppressRedraw = 0;


mapAxis = (-sLength/2+binWidth/2):binWidth:(sLength/2-binWidth/2);

t = positions(:,1);
x = positions(:,2);
if size(positions,2) > 2
    y = positions(:,3);
end

% Do median filtering to supress off-values
%[x, y, t] = medianFilter(x, y, t);

% Smoothing the position samples with a simple mean filter
for cc=8:length(x)-7
    x(cc) = nanmean(x(cc-7:cc+7));
    y(cc) = nanmean(y(cc-7:cc+7));
end
    
% With many off-values in a row the median filter doesn't remove all, must
% apply an additional off-value filter.
%[x, y, t] = offValueFilter(x, y, t);
    
[x,y] = centreBox(x,y);
    
% Calculate the time map for this session
timeMaps = findTimeMap(x,y,t,bins,sLength,binWidth);
    
% Store path
pathCoord{1,1} = x;
pathCoord{1,2} = y;
pathCoord{1,3} = t;
% END COMPUTING POSITION DATA
    
% Get position to spikes
[spkx,spky] = spikePos(tSp,x,y,t);
% Calculate rate map
map = ratemap(tSp,spkx,spky,x,y,t,h,mapAxis);

rateMaps = map;
peak = max(max(map));
timeStamps = tSp;


fieldAxis = linspace(-sLength/2,sLength/2,bins);

% Remove unvisited parts of the box from the map
h1 = figure; 
if suppressRedraw
    set(h1,'Visible','off');
 %   set(h2,'Visible','off');
end

drawfield(rateMaps,mapAxis,'jet',peak,peak);
set(gcf,'Color', [1 1 1]);
drawnow;

h2 = figure;
% Plot path of rat with spikes
%set(h2,'CurrentFigure',h2)
if suppressRedraw
    set(h2,'Visible','off');
end
h2 = plot(pathCoord{1,1},pathCoord{1,2}); % Path
set(h2,'color',[0.5 0.5 0.5]); % Set color of the path to gray
%title('path with spikes');
hold on;
% spkInd = spkIndex(pathCoord{1,3},timeStamps);
% spkX = pathCoord{1,1}(spkInd);
% spkY = pathCoord{1,2}(spkInd);
plot(spkx,spky,'r.'); % Spikes
set(gcf,'Color', [1 1 1]);
hold off;
axis image;
axis([-sLength/2,sLength/2,-sLength/2,sLength/2]);
axis off;
drawnow;





%__________________________________________________________________________
%
% Field functions
%__________________________________________________________________________

% timeMap will find the amount of time the rat spends in each bin of the
% box. The box is divided into the same number of bins as in the rate map.
% posx, posy, post: Path coordinates and time stamps.
% bins: Number of bins the box is diveded int (bins x bins)
% extremal: minimum and maximum positions.
function timeMap = findTimeMap(posx,posy,post,bins,sLength,binWidth)

% Duration of trial
duration = post(end)-post(1);
% Average duration of each position sample
sampDur = duration/length(posx);

% X-coord for current bin position
pcx = -sLength/2-binWidth/2;
timeMap = zeros(bins);
% Find number of position samples in each bin
for ii = 1:bins
    % Increment the x coordinate
    pcx = pcx + binWidth;
    I = find(posx >= pcx & posx < pcx+binWidth);
    % Y-coord for current bin position
    pcy = -sLength/2-binWidth/2;;
    for jj=1:bins
        % Increment the y coordinate
        pcy = pcy + binWidth;
        J = find(posy(I) >= pcy & posy(I) < pcy+binWidth);
        % Number of position samples in the current bin
        timeMap(jj,ii) = length(J);
    end
end
% Convert to time spent in bin
timeMap = timeMap * sampDur;


% Calculates the rate map.
function map = ratemap(ts,spkx,spky,posx,posy,post,h,mapAxis)
invh = 1/h;
map = zeros(length(mapAxis),length(mapAxis));
yy = 0;
for y = mapAxis
    yy = yy + 1;
    xx = 0;
    for x = mapAxis
        xx = xx + 1;
        map(yy,xx) = rate_estimator(ts,spkx,spky,x,y,invh,posx,posy,post);
    end
end

% Calculate the rate for one position value
function r = rate_estimator(ts,spkx,spky,x,y,invh,posx,posy,post)
% edge-corrected kernel density estimator

conv_sum = sum(gaussian_kernel(((spkx-x)*invh),((spky-y)*invh)));
edge_corrector =  trapz(post,gaussian_kernel(((posx-x)*invh),((posy-y)*invh)));
edge_corrector(edge_corrector<0.15) = NaN;
r = (conv_sum / (edge_corrector + 0.1)) + 0.1; % regularised firing rate for "wellbehavedness"
                                                       % i.e. no division by zero or log of zero
                                                       
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,y)
% k(u) = ((2*pi)^(-length(u)/2)) * exp(u'*u)
r = 0.15915494309190 * exp(-0.5*(x.*x + y.*y));


% Finds the position to the spikes
function [spkx,spky] = spikePos(ts,posx,posy,post)
N = length(ts);
spkx = zeros(N,1);
spky = zeros(N,1);
for ii = 1:N
    tdiff = (post-ts(ii)).^2;
    [m,ind] = min(tdiff);
    spkx(ii) = posx(ind(1));
    spky(ii) = posy(ind(1));
end

%__________________________________________________________________________
%
% Function for modifing the path
%__________________________________________________________________________

% Median filter for positions
function [posx,posy,post] = medianFilter(posx,posy,post)
N = length(posx);
x = [NaN*ones(15,1); posx'; NaN*ones(15,1)];
y = [NaN*ones(15,1); posy'; NaN*ones(15,1)];
X = ones(1,N);
Y = ones(1,N);
for cc=16:N+15
    lower = cc-15;
    upper = cc+15;
    X(cc-15) = nanmedian(x(lower:upper));
    Y(cc-15) = nanmedian(y(lower:upper));
end
index = find(isfinite(X));
posx=X(index);
posy=Y(index);
post=post(index);

% Removes position values that are a result of bad tracking
function [posx,posy,post] = offValueFilter(posx,posy,post)

x1 = posx(1:end-1);
x2 = posx(2:end);
y = posy(2:end);
t = post(2:end);

N = 1;
while N > 0
    len = length(x1);
    dist = abs(x2 - x1);
    % Finds were the distance between two neighbour samples are more than 2
    % cm.
    ind = find(dist > 2);
    N = length(ind);
    if N > 0
        if ind(1) ~= len
            x2(:,ind(1)) = [];
            x1(:,ind(1)+1) = [];
            y(:,ind(1)) = [];
            t(:,ind(1)) = [];
        else
            x2(:,ind(1)) = [];
            x1(:,ind(1)) = [];
            y(:,ind(1)) = [];
            t(:,ind(1)) = [];
        end
    end
end
x = x2(2:end);
t = t(2:end);
y1 = y(1:end-1);
y2 = y(2:end);

N = 1;
while N > 0
    len2 = length(y1);
    dist = abs(y2 - y1);
    ind = find(dist > 3);
    N = length(ind);
    if N > 0
        if ind(1) ~= len2
            y2(:,ind(1)) = [];
            y1(:,ind(1)+1) = [];
            x(:,ind(1)) = [];
            t(:,ind(1)) = [];
        else
            y2(:,ind(1)) = [];
            y1(:,ind(1)) = [];
            x(:,ind(1)) = [];
            t(:,ind(1)) = [];
        end
    end
end
posx = x;
posy = y2;
post = t;

% Centre the path/box for the two sessions. Both boxes are set according to
% the path/coordinates to the box for the first session.
function [posx,posy] = centreBox(posx,posy)

% Find border values for box for session 1
maxX = max(posx);
minX = min(posx);
maxY = max(posy);
minY = min(posy);

% Set the corners of the reference box
NE = [maxX, maxY];
NW = [minX, maxY];
SW = [minX, minY];
SE = [maxX, minY];

% Get the centre coordinates of the box
centre = findCentre(NE,NW,SW,SE);
% Centre both boxes according to the coordinates to the first box
posx = posx - centre(1);
posy = posy - centre(2);


% Calculates the centre of the box from the corner coordinates
function centre = findCentre(NE,NW,SW,SE);

% The centre will be at the point of interception by the corner diagonals
a = (NE(2)-SW(2))/(NE(1)-SW(1)); % Slope for the NE-SW diagonal
b = (SE(2)-NW(2))/(SE(1)-NW(1)); % Slope for the SE-NW diagonal
c = SW(2);
d = NW(2);
x = (d-c+a*SW(1)-b*NW(1))/(a-b); % X-coord of centre
y = a*(x-SW(1))+c; % Y-coord of centre
centre = [x,y];


%__________________________________________________________________________
%
% Additional graphics function
%__________________________________________________________________________

function drawfield(map,faxis,cmap,maxrate_actual,maxrate_draw)
   
   % This function will calculate an RGB image from the rate
   % map. We do not just call image(map) and caxis([0 maxrate]),
   % as it would plot unvisted parts with the same colour code
   % as 0 Hz firing rate. Instead we give unvisited bins
   % their own colour (e.g. gray or white).

   maxrate_actual = ceil(maxrate_actual);
   maxrate_draw = ceil(maxrate_draw);
   if maxrate_draw < 1
      maxrate_draw = 1;
   end
   n = size(map,1);
   plotmap = ones(n,n,3);
   for jj = 1:n
      for ii = 1:n
         if isnan(map(jj,ii))
            plotmap(jj,ii,1) = 1; % give the unvisited bins a gray colour
            plotmap(jj,ii,2) = 1;
            plotmap(jj,ii,3) = 1;
         else
             if (map(jj,ii) > maxrate_draw)
                 plotmap(jj,ii,1) = 1; % give the unvisited bins a gray colour
                 plotmap(jj,ii,2) = 0;
                 plotmap(jj,ii,3) = 0;
             else    
                 rgb = pixelcolour(map(jj,ii),maxrate_draw,cmap);
                 plotmap(jj,ii,1) = rgb(1);
                 plotmap(jj,ii,2) = rgb(2);
                 plotmap(jj,ii,3) = rgb(3);
             end
         end
      end
   end
   image(faxis,faxis,plotmap);
   set(gca,'YDir','Normal');
   axis image;
   axis off


function rgb = pixelcolour(map,maxrate,cmap)

   % This function calculates a colour for each bin
   % in the rate map.

   cmap1 = ...
    [    0         0    0.5625; ...
         0         0    0.6875; ...
         0         0    0.8125; ...
         0         0    0.9375; ...
         0    0.0625    1.0000; ...
         0    0.1875    1.0000; ...
         0    0.3125    1.0000; ...
         0    0.4375    1.0000; ...
         0    0.5625    1.0000; ...
         0    0.6875    1.0000; ...
         0    0.8125    1.0000; ...
         0    0.9375    1.0000; ...
    0.0625    1.0000    1.0000; ...
    0.1875    1.0000    0.8750; ...
    0.3125    1.0000    0.7500; ...
    0.4375    1.0000    0.6250; ...
    0.5625    1.0000    0.5000; ...
    0.6875    1.0000    0.3750; ...
    0.8125    1.0000    0.2500; ...
    0.9375    1.0000    0.1250; ...
    1.0000    1.0000         0; ...
    1.0000    0.8750         0; ...
    1.0000    0.7500         0; ...
    1.0000    0.6250         0; ...
    1.0000    0.5000         0; ...
    1.0000    0.3750         0; ...
    1.0000    0.2500         0; ...
    1.0000    0.1250         0; ...
    1.0000         0         0; ...
    0.8750         0         0; ...
    0.7500         0         0; ...
    0.6250         0         0 ];
   
   cmap2 = ...
   [0.0417         0         0; ...
    0.1250         0         0; ...
    0.2083         0         0; ...
    0.2917         0         0; ...
    0.3750         0         0; ...
    0.4583         0         0; ...
    0.5417         0         0; ...
    0.6250         0         0; ...
    0.7083         0         0; ...
    0.7917         0         0; ...
    0.8750         0         0; ...
    0.9583         0         0; ...
    1.0000    0.0417         0; ...
    1.0000    0.1250         0; ...
    1.0000    0.2083         0; ...
    1.0000    0.2917         0; ...
    1.0000    0.3750         0; ...
    1.0000    0.4583         0; ...
    1.0000    0.5417         0; ...
    1.0000    0.6250         0; ...
    1.0000    0.7083         0; ...
    1.0000    0.7917         0; ...
    1.0000    0.8750         0; ...
    1.0000    0.9583         0; ...
    1.0000    1.0000    0.0625; ...
    1.0000    1.0000    0.1875; ...
    1.0000    1.0000    0.3125; ...
    1.0000    1.0000    0.4375; ...
    1.0000    1.0000    0.5625; ...
    1.0000    1.0000    0.6875; ...
    1.0000    1.0000    0.8125; ...
    1.0000    1.0000    0.9375];
   
   if strcmp(cmap,'jet')
      steps = (31*(map/maxrate))+1;
      steps = round(steps);
      rgb = cmap1(steps,:);
   else
      steps = (31*(map/maxrate))+1;
      steps = round(steps);
      rgb = cmap2(steps,:);
   end
   
% %__________________________________________________________________________
% %
% %      Functions for reading Mclust data
% %__________________________________________________________________________
%    
% function S = LoadSpikes(tfilelist, path, NaFile)
% % tfilelist:    List of t-files. Each file contains a cluster of spikes
% %               from a cell.
% % path:         Path to the directory were the t-files are stored
% % NaFile:       List of file names in tfilelist that don't exist
% %               in the current directory
% %
% % inp: tfilelist is a cellarray of strings, each of which is a
% % tfile to open.  Note: this is incompatible with version unix3.1.
% % out: Returns a cell array such that each cell contains a ts 
% % object (timestamps which correspond to times at which the cell fired)
% %
% % Edited by: Raymond Skjerpeng
% 
% 
% %-------------------
% % Check input type
% %-------------------
% 
% if ~isa(tfilelist, 'cell')
%    error('LoadSpikes: tfilelist should be a cell array.');
% end
% 
% 
% % Number of file names in tfilelist
% nFiles = length(tfilelist);
% % Actual number of files to be loaded
% anFiles = nFiles - length(NaFile);
% 
% 
% %--------------------
% % Read files
% %--------------------
% fprintf(2, 'Reading %d files.', anFiles);
% 
% % for each tfile
% % first read the header, then read a tfile 
% % note: uses the bigendian modifier to ensure correct read format.
% 
% 
% S = cell(nFiles, 1);
% for iF = 1:nFiles
%     %DisplayProgress(iF, nFiles, 'Title', 'LoadSpikes');
%     tfn = tfilelist{iF};
%     % Check if file exist
%     if length(strmatch(char(tfn),NaFile,'exact'))>0
%         S{iF} = -1; % Set this as default for each file that doesn't exist
%     else
%         tfn = strcat(strcat(path,'\'),tfn); % Path to file + file name
%         if ~isempty(tfn)
%             tfp = fopen(tfn, 'rb','b');
%             if (tfp == -1)
%                 warning([ 'Could not open tfile ' tfn]);
%             end
% 
%             ReadHeader(tfp);    
%             S{iF} = fread(tfp,inf,'uint32');	%read as 32 bit ints
%             S{iF} = ts(S{iF});
%             fclose(tfp);
%         end 	% if tfn valid
%     end
% end		% for all files
% fprintf(2,'\n');


% function F = ReadFileList(fn)
% 
% % F = ReadFileList(fn)
% %
% % INPUTS: 
% %   fn -- an ascii file of filenames, 1 filename per line
% % 
% % OUTPUTS:
% %   F -- a cell array of filenames suitable for use in programs
% %        such as LoadSpikes
% %
% % Now can handle files with headers
% 
% % ADR 1998
% % version L4.1
% % status: PROMOTED
% 
% % v4.1 added ReadHeader
% 
% [fp,errmsg] = fopen(fn, 'rt');
% if (fp == -1)
%    error(['Could not open "', fn, '". ', errmsg]);
% end
% 
% ReadHeader(fp);
% ifp = 1;
% while (~feof(fp))
%    F{ifp} = fgetl(fp);
%    ifp = ifp+1;
% end
% fclose(fp);
% 
% F = F';
% 
% 
% function H = ReadHeader(fp)
% % H = ReadHeader(fp)
% %  Reads NSMA header, leaves file-read-location at end of header
% %  INPUT: 
% 
% %      fid -- file-pointer (i.e. not filename)
% %  OUTPUT: 
% %      H -- cell array.  Each entry is one line from the NSMA header
% % Now works for files with no header.
% % ADR 1997
% % version L4.1
% % status: PROMOTED
% % v4.1 17 nov 98 now works for files sans header
% %
% % May 2010, modified by Emily to work with new Matlab
% %---------------
% 
% % Get keys
% beginheader = '%%BEGINHEADER';
% endheader = '%%ENDHEADER';
% 
% iH = 1; H = {};
% curfpos = ftell(fp);
% 
% % look for beginheader
% headerLine = my_fgetl(fp);
% if strcmp(headerLine, beginheader)
%    H{1} = headerLine;
%    while ~feof(fp) & ~strcmp(headerLine, endheader)     
%       headerLine = my_fgetl(fp);
%       iH = iH+1;
%       H{iH} = headerLine;
%    end
% else % no header
%    fseek(fp, curfpos, 'bof');
% end
% 
% function tline = my_fgetl(fid)
% 
% try
%     [tline,lt] = fgets(fid);
%     tline = tline(1:end-length(lt));
%     fseek(fid, -(length(lt)-1), 'cof');
%     
%     if isempty(tline)
%         tline = '';
%     end
% 
% catch exception
%     if nargin ~= 1
%         error (nargchk(1,1,nargin,'struct'))
%     end
%     throw(exception);
% end
% 
% function spkInd = spkIndex(post,ts)
% M = length(ts);
% spkInd = zeros(M,1);
% for ii=1:M
%     tdiff = (post-ts(ii)).^2;
%     [m,ind] = min(tdiff);
%     spkInd(ii) = ind;
% end