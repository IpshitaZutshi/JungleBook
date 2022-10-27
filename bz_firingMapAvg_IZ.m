function [firingMaps] = bz_firingMapAvg_IZ(positions,spikes,varargin)

% USAGE
% [firingMaps] = bz_firingMapAvg(positions,spikes,varargin)
% Calculates averaged firing map for a set of linear postions 
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%      or
%   behavior  - buzcode format behavior struct - NOT YET IMPLEMENTED
%   <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'			smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'			number of bins (default = 50)
%     'speedThresh'		speed threshold to compute firing rate
%     'minTime'			minimum time spent in each bin (in s, default = 0)
%     'mode'			interpolate' to interpolate missing points (< minTime),
%                   	or 'discard' to discard them (default)
%     'maxDistance'		maximal distance for interpolation (default = 5)
%     'maxGap'			z values recorded during time gaps between successive (x,y)
%                   	samples exceeding this threshold (e.g. undetects) will not
%                	    be interpolated; also, such long gaps in (x,y) sampling
%                 	    will be clipped to 'maxGap' to compute the occupancy map
%                 	    (default = 0.100 s)
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'saveMat'   		- logical (default: false) that saves firingMaps file
%     'CellInspector'  	- logical (default: false) that creates an otuput
%                   	compatible with CellInspector

%
%
% OUTPUT
%
%   firingMaps - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%
% Antonio FR, 10/2019

%% parse inputs
p=inputParser;
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'speedThresh',0.1,@isnumeric);
addParameter(p,'nBins',55,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0.1,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'CellInspector',false,@islogical);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'downsample',false,@islogical);
addParameter(p,'plotFig',true,@islogical);

parse(p,varargin{:});
smooth = p.Results.smooth;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
saveMat = p.Results.saveMat;
CellInspector = p.Results.CellInspector;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
order = p.Results.orderKalmanVel;
plotFig = p.Results.plotFig;
downsample = p.Results.downsample;

% number of conditions
  if iscell(positions)
     conditions = length(positions); 
  elseif isvector(positions)
     conditions = 1;
  end
  %%% TODO: conditions label
  
%% Calculate
%Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    if ~isempty(positions{iCond})
        post = positions{iCond}(:,1);
        % - 1D 
        if size(positions{iCond},2)==2
            posx = positions{iCond}(:,2);
    %         for cc=8:length(posx)-7
    %             posx(cc) = nanmedian(posx(cc-7:cc+7));
    %         end
        if sum(~isnan(posx))>0
            [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
        end
        elseif size(positions{iCond},2)==3
            posx = positions{iCond}(:,2);
            posy = positions{iCond}(:,3);
    %         for cc=8:length(posx)-7
    %             posx(cc) = nanmedian(posx(cc-7:cc+7));
    %             posy(cc) = nanmedian(posy(cc-7:cc+7));
    %         end
            [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
        else
            warning('This is not a linear nor a 2D space!');
        end
        % Absolute speed
        if sum(~isnan(posx))>0
            v = sqrt(vx.^2+vy.^2);

            % Compute timestamps where speed is under threshold
            positions{iCond}(:,2) = posx;
            if size(positions{iCond},2) > 2
                positions{iCond}(:,3) = posy;
            end
            positions{iCond}(v<speedThresh,:) = [];
        end
    end

end


% get firing rate maps
% if plotFig
%     f1 = figure;
%     f2 = figure;
% end
for unit = 1:length(spikes.times)
    if plotFig     
        figure;
        set(gcf,'Position',[100 100 1800 400])
    end
    for c = 1:conditions
        if ~isempty(positions{c})
            map{unit}{c} = Map(positions{c},spikes.times{unit},'smooth',smooth,'minTime',minTime,...
               'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
            %map{unit}{c}.z(map{unit}{c}.z==0) = nan;

            if plotFig     
                %Plot the trajectory
                if nanmean(nanmean(map{unit}{c}.z)) < 10
%                     if mod(c,2) == 0
%                         loc = 2;
%                     else
%                         loc = 1;
%                     end
%                     subplot(3,conditions,[loc loc+conditions])
                    subplot(3,conditions,(2*conditions)+c)
                    if size(positions{c},2) == 2 %Its a linear map
                        plot(positions{c}(:,2));
                    else
                        a= [positions{c}(:,2)];%,positions{c}(:,3)];
                        b = abs(diff(a));
                        X = positions{c}(2:end,2);
                        X(b>4) = nan;   
                        plot(positions{c}(:,2),positions{c}(:,3),'k')
                        %scatter(positions{c}(:,2),positions{c}(:,3),3,[0.5 0.5 0.5],'.')
                    end
                    hold on
                    spikePos = findspikePos(positions{c},spikes.times{unit});
                    if size(positions{c},2) == 2 %Its a linear map
                        plot(spikePos(:,2),spikePos(:,3),'r.');
                    else
                        plot(spikePos(:,3),spikePos(:,4),'r.');
                    end
                    view([90 -90])
                end
                
                %Plot the rate map
                peak = max(max(map{unit}{c}.z));
                avgRate = nanmean(nanmean(map{unit}{c}.z));
                subplot(3,conditions,(2*conditions)+c)
                if size(positions{c},2) == 2 %Its a linear map
                    imagesc(map{unit}{c}.z)
                    colormap('jet')
                else
                    imagesc(map{unit}{c}.z)
                    drawfield(map{unit}{c}.z,'jet',peak);
                end
                title(strcat('Peakrate:',num2str(peak)));
            end
        else
            map{unit}{c} = [];
        end
    end
    if plotFig     
        if ~isfolder('FiringMap')
            mkdir('FiringMap')
        end
        saveas(gcf,['FiringMap',filesep ,'cell_' num2str(unit) '.png'],'png');
        saveas(gcf,['FiringMap',filesep ,'cell_' num2str(unit) '.fig'],'fig');
        close all;
    end
end
%%% TODO: pass rest of inputs to Map

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
%firingMaps.sessionName = spikes.sessionName;
try
firingMaps.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end

firingMaps.params.smooth = smooth;
firingMaps.params.minTime = minTime;
firingMaps.params.nBins = nBins;
firingMaps.params.maxGap = maxGap;
firingMaps.params.mode = mode;
firingMaps.params.maxDistance = maxDistance;

for unit = 1:length(spikes.times)
    for c = 1:conditions
        if ~isempty(map{unit}{c})
            firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
            firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
            firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
        else
            firingMaps.rateMaps{unit,1}{c} = [];
            firingMaps.countMaps{unit,1}{c} = [];
            firingMaps.occupancy{unit,1}{c} = [];
        end
    end
end

if saveMat
   if ~downsample
       save([spikes.basename '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 
   else
       save([spikes.basename '_DS.firingMapsAvg.cellinfo.mat'],'firingMaps'); 
   end       
end

end

function spikePos = findspikePos(positions,tSp)

for ii = 1:length(tSp)
    [~,idx] = min(abs(positions(:,1)-tSp(ii)));
    spikePos(ii,1) = positions(idx,1);
    spikePos(ii,2) = idx;
    spikePos(ii,3) = positions(idx,2);
    if size(positions,2)>2
        spikePos(ii,4) = positions(idx,3);
    end
end

end
