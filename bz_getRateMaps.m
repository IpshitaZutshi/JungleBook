function [firingMaps] = bz_getRateMaps(positions,spikes,varargin)

% USAGE
% [firingMaps] = bz_getRateMaps(positions,spikes,varargin)
% Calculates averaged firing map for a set of linear postions 
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%   <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'lin'             % whether to make a linear or 2D map
%     'speedThresh'		speed threshold to compute firing rate
%     'binSize'			% Spatial bin size (cm)
%     'xRange'			% Range of x coordinate values
%     'yRange'          % Range of y coordinate value
%     'minOccupancy' 	minimum time a bin needs to be occupied
%     'posFs'           sampling rate
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'saveMat'   		- logical (default: false) that saves firingMaps file

%% parse inputs
p=inputParser;
addParameter(p,'speedThresh',0.1,@isnumeric);
addParameter(p,'lin',true,@islogical);
addParameter(p,'binSize',2.5,@isnumeric);
addParameter(p,'xRange',[0 6],@isnumeric);
addParameter(p,'yRange',[0 125],@isnumeric);
addParameter(p,'minOccupancy',0.01,@isnumeric);
addParameter(p,'posFs',30,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'saveMat',false,@islogical);

parse(p,varargin{:});
speedThresh = p.Results.speedThresh;
lin = p.Results.lin;
binSize = p.Results.binSize;
xRange = p.Results.xRange;
yRange = p.Results.yRange;
minOccupancy = p.Results.minOccupancy;
posFs = p.Results.posFs;
order = p.Results.orderKalmanVel;
saveMat = p.Results.saveMat;

% number of conditions
if iscell(positions)
    conditions = length(positions); 
elseif isvector(positions)
    conditions = 1;
end

%% Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    if ~isempty(positions{iCond})
        post = positions{iCond}(:,1);
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        % If all nans just skip this 
        if ~isnan(posx)
            [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);

            % Absolute speed
            if sum(~isnan(posx))>0
                v = sqrt(vx.^2+vy.^2);
                % Compute timestamps where speed is under threshold
                positions{iCond}(:,2) = posx;
                positions{iCond}(:,3) = posy;
                positions{iCond}(v<speedThresh,:) = [];
            end
        end
    end

end

if lin
    mapSize(1) = 1;
    mapSize(2) = round((yRange(2)-yRange(1))./binSize);
else
    mapSize(1) = (xRange(2)-xRange(1))./binSize;
    mapSize(2) = (yRange(2)-yRange(1))./binSize;
end

binSizeX = diff(xRange)/mapSize(1);
binSizeY = diff(yRange)/mapSize(2);
    
mapAxisX = (xRange(1)+binSizeX/2):binSizeX:(xRange(2)-binSizeX/2);
mapAxisY = (yRange(1)+binSizeY/2):binSizeY:(yRange(2)-binSizeY/2);

%% Build rate maps
for unit = 1:length(spikes.times)
    tSp = spikes.times{unit};
    for c = 1:conditions
        t = positions{c}(:,1);
        x = positions{c}(:,2);
        y = positions{c}(:,3);        
        countMap = getTimeMap(x,y,xRange,yRange,mapSize); % get occupancy map;
        timeMap = countMap/posFs;        
        [xSp,ySp,~,iSp] = bz_spikePos2(tSp,x,y,t);
        spikeMap = getTimeMap(xSp,ySp,xRange,yRange,mapSize); 
        rateMap = bz_calRateMap(spikeMap,timeMap,minOccupancy);
        visited = visitedBins(x,y,mapAxisX, mapAxisY);
        visited = visited | timeMap>0;
        rateMap(visited==0) = NaN;
        map{unit}{c}.z = rateMap;
        map{unit}{c}.count = spikeMap;
        map{unit}{c}.time = timeMap;
        posPDF = timeMap/sum(sum(timeMap));        
        [map{unit}{c}.information,map{unit}{c}.sparsity,map{unit}{c}.selectivity] = bz_mapstat(rateMap,posPDF);
    end
end

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
%firingMaps.sessionName = spikes.sessionName;
try
firingMaps.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end
firingMaps.params.speedThresh = speedThresh;
firingMaps.params.binSize = binSize;
firingMaps.params.minOccupancy = minOccupancy;

for unit = 1:length(spikes.times)
    for c = 1:conditions
        if ~isempty(map{unit}{c})
            firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
            firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
            firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
            firingMaps.information{unit,1}{c} = map{unit}{c}.information;
            firingMaps.selectivity{unit,1}{c} = map{unit}{c}.selectivity;
            firingMaps.sparsity{unit,1}{c} = map{unit}{c}.sparsity;            
            firingMaps.peakRate{unit,1}{c} = nanmax(map{unit}{c}.z); 
            firingMaps.avgRate{unit,1}{c} = nanmean(map{unit}{c}.z); 
        else
            firingMaps.rateMaps{unit,1}{c} = [];
            firingMaps.countMaps{unit,1}{c} = [];
            firingMaps.occupancy{unit,1}{c} = [];
            firingMaps.information{unit,1}{c} = nan;
            firingMaps.selectivity{unit,1}{c} = nan;
            firingMaps.sparsity{unit,1}{c} = nan;  
            firingMaps.peakRate{unit,1}{c} = nan; 
            firingMaps.avgRate{unit,1}{c} = nan;
        end
    end
end

if saveMat
   save([spikes.basename '.rateMaps.cellinfo.mat'],'firingMaps');     
end

end

function countMap = getTimeMap(x,y,xRange,yRange,mapSize)
    % suppressedges - set position samples outside the defined box equal to the
    % border
    ind = x < xRange(1);
    x(ind) = xRange(1);
    ind = x > xRange(2);
    x(ind) = xRange(2);

    ind = y < yRange(1);
    y(ind) = yRange(1);
    ind = y > yRange(2);
    y(ind) = yRange(2);

    binSizeX = diff(xRange)/mapSize(1);
    binSizeY = diff(yRange)/mapSize(2);
    xEdges = xRange(1):binSizeX:xRange(2);
    yEdges = yRange(1):binSizeY:yRange(2);

    [~,xBin] = histc(x,xEdges); xBin(xBin==0) = 1;
    [~,yBin] = histc(y,yEdges); yBin(yBin==0) = 1;
    xBin(xBin>mapSize(1)) = mapSize(1);
    yBin(yBin>mapSize(2)) = mapSize(2);
    nSamples = length(x);
    countMap = zeros(mapSize);

    for s = 1:nSamples
        countMap(xBin(s),yBin(s)) = countMap(xBin(s),yBin(s))+1;
    end
end

function visited = visitedBins(posx,posy,mapAxisX, mapAxisY)
    % Number of bins in each direction of the map
    Nx = length(mapAxisX);
    Ny = length(mapAxisY);
    binWidth = mapAxisY(2) - mapAxisY(1);

    visited = zeros(Nx,Ny);
    if ~isnan(posy)
        iRegion = convhull(posx,posy);

        for ii = 1:Nx
            for jj = 1:Ny
                px = mapAxisX(ii);
                py = mapAxisY(jj);
                if inpolygon(px,py,posx(iRegion),posy(iRegion))
                    distance = sqrt( (px-posx).^2 + (py-posy).^2 );            
                    if min(distance) <= binWidth
                        visited(ii,jj) = 1;
                    end
                end
            end
        end
    end
end
