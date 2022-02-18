% Compile data across all sessions

function PP = SessBehaviorPhasePrecession(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'savePlots',true,@islogical);
addParameter(p,'makePlots',true,@islogical);
addParameter(p,'byCycle',true,@islogical);
addParameter(p,'normalize',true,@islogical);
addParameter(p,'heat',false,@islogical);

parse(p,varargin{:});

expPath = p.Results.expPath;
savePlots = p.Results.savePlots;
makePlots = p.Results.makePlots;
saveMat = p.Results.saveMat;
force = p.Results.force;
byCycle = p.Results.byCycle;
normalize = p.Results.normalize;
heat = p.Results.heat;

%%
if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if exist([expPath '\Summ\phasePrecession.mat'],'file') && ~force 
    disp('Phase Precession already computed! Loading file.');
    load([expPath '\Summ\phasePrecession.mat']);
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:4 %4 fields, left trajectories, right trajectories, stim left and right trajectories
                PP.slope{rr,cc}{zz} = [];
                PP.intercept{rr,cc}{zz} = [];
                PP.r2{rr,cc}{zz} = [];
                PP.p{rr,cc}{zz} = [];                
                PP.boundaries{rr,cc}{zz} = [];
                PP.region{rr,cc}{zz} = [];
                PP.putativeCellType{rr,cc}{zz} = [];                 
            end
        end
    end

    %% Start collecting data
    for ii = 1:size(allSess,1)
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.session.mat'));
        load(file.name);  
        
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            load(strcat(efields{jj},'\',allSess(ii).name,'.phasePrecessionAvg.cellinfo.mat'))
            
            for kk = 1:length(phasePrecession.data)
                for ss = 1:size(phasePrecession.data{kk},1)
                    if makePlots
                        figure
                    end
                    for zz = 1:size(phasePrecession.data{kk},2)
                        if isstruct(phasePrecession.data{kk}{ss,zz})
                            if makePlots
                                plotPrecession(zz,phasePrecession.data{kk}{ss,zz},phasePrecession.stats{kk}{ss,zz},'byCycle',byCycle,'normalize',normalize,'heat',heat);
                            end
                            PP.slope{region, target}{zz} = [PP.slope{region, target}{zz};phasePrecession.stats{kk}{ss,zz}.slope];
                            PP.intercept{region, target}{zz} = [PP.intercept{region, target}{zz};phasePrecession.stats{kk}{ss,zz}.intercept];
                            PP.p{region, target}{zz} = [PP.p{region, target}{zz};phasePrecession.stats{kk}{ss,zz}.p];
                            PP.r2{region, target}{zz} = [PP.r2{region, target}{zz};phasePrecession.stats{kk}{ss,zz}.r2];                            
                            PP.boundaries{region, target}{zz} = [PP.boundaries{region, target}{zz};phasePrecession.stats{kk}{ss,zz}.boundaries];
                            
                            % Assign numerical tag to putative class
                            if strcmp(cell_metrics.putativeCellType{kk},'Narrow Interneuron') == 1
                                PP.putativeCellType{region, target}{zz} = [PP.putativeCellType{region, target}{zz}; 1];
                            elseif strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1
                                PP.putativeCellType{region, target}{zz} = [PP.putativeCellType{region, target}{zz}; 2];
                            elseif strcmp(cell_metrics.putativeCellType{kk},'Wide Interneuron') == 1
                                PP.putativeCellType{region, target}{zz} = [PP.putativeCellType{region, target}{zz}; 3];
                            else 
                                PP.putativeCellType{region, target}{zz} = [PP.putativeCellType{region, target}{zz}; 4];                       
                            end                    

                            % Assign numerical tag to region
                            if strcmp(cell_metrics.brainRegion{kk},'CA1') == 1
                                PP.region{region, target}{zz} = [PP.region{region, target}{zz}; 1];
                            elseif strcmp(cell_metrics.brainRegion{kk},'DG') == 1
                                PP.region{region, target}{zz} = [PP.region{region, target}{zz}; 2];
                            elseif strcmp(cell_metrics.brainRegion{kk},'CA3') == 1
                                PP.region{region, target}{zz} = [PP.region{region, target}{zz}; 3];
                            else 
                                PP.region{region, target}{zz} = [PP.region{region, target}{zz}; 4];                       
                            end                            
                           
                        else
                            continue;
                        end
                    end
                    if savePlots
                        if ~isfolder(strcat(efields{jj},'\PrecessPlots'))
                            mkdir(efields{jj},'PrecessPlots')
                        end
                        saveas(gcf,[efields{jj},filesep,'PrecessPlots',filesep ,'cell_' num2str(kk),'field_',num2str(ss),'.png'],'png');
                        saveas(gcf,[efields{jj},filesep,'PrecessPlots',filesep ,'cell_' num2str(kk),'field_',num2str(ss),'.fig'],'fig');
                        saveas(gcf,[efields{jj},filesep,'PrecessPlots',filesep ,'cell_' num2str(kk),'field_',num2str(ss),'.eps'],'epsc');                
                        close all;
                    end
                end
            end

        end

    end
    if saveMat
        save([expPath '\Summ\phasePrecession.mat'], 'PP');
    end
end

end


function plotPrecession(zz,data,stats,varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'byCycle',true,@islogical);
addParameter(p,'normalize',true,@islogical);
addParameter(p,'heat',false,@islogical);
parse(p,varargin{:});

byCycle = p.Results.byCycle;
normalize = p.Results.normalize;
heat = p.Results.heat;

%%
t = data.position.t;
x = data.position.x;
pp = data.position.phase; 
measure = data.lindat;
phase = data.phidat;
slope = stats.slope;
offset = stats.intercept;

a = subplot(2,4,zz);
nBins = 200;
curve = FiringCurve(data.x,t,'nbins',nBins,'smooth',2,'type','ll');
plot(curve.x,curve.rate,'Color',[85/243 85/243 85/243]);
hold on
idxField = curve.x >=stats.boundaries(1) & curve.x <= stats.boundaries(2);
plot(curve.x(idxField),curve.rate(idxField),'Color',[8/243 133/243 161/243]);
xlim([0 1]);
set(a,'tickdir','out');
ylabel('Firing Rate (Hz)');


if normalize
    measure = measure/max(measure);
    if heat
        xRange = [0 1];
    else
        xRange = [min(measure) max(measure)];
    end
end

if byCycle
    phase = (phase+pi)/(2*pi);
    unit = 'Cycles';
    shift = 1;
    yRange = [0 2];
    offset = (offset+pi)/(2*pi);
else
    slope = slope*2*pi;
    unit = 'Radians';
    shift = 2*pi;
    yRange = [-pi pi];
end

reg = offset+slope*measure;

a = subplot(2,4,zz+4);
if heat
    xRange = [0 1];
    heatPlot([measure measure],[phase phase+shift],xRange,yRange,[12 12],[100 100]);
    lineStyle = 'w-';
    lineWidth = 3;
else
    plot([measure measure],[phase phase+shift],'k.');
    lineStyle = 'r-';
    lineWidth = 0.5;
end
hold on;
plot([measure measure],[reg reg+shift],lineStyle,'LineWidth',lineWidth);
ylabel(['Phase (' unit ')']);
ylim(yRange);
title(strcat('S:',num2str(stats.slope),'p =',num2str(stats.r2)));

end
