function compareHCvsMazeCA3(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'savePlot',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
savePlot = p.Results.savePlot;
saveMat = p.Results.saveMat;
force = p.Results.force;

tag = 'CA3Saline'; 

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ11\Final','IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ28\Final',...
        'IZ29\Final','IZ30\Final','IZ31\Final','IZ32\Final','IZ33\Final','IZ34\Final'};  % To add: IZ23, IZ32, IZ33, IZ34
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'}; 
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'}; % To add: IZ32, IZ33, IZ34
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    %
    reg = {'contramEC','ipsimEC','Both'};
end

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
tar = {'STEM', 'RETURN'};
totCells = 1;

if exist([parentDir 'Compiled\mazeVsHC' tag '.mat'],'file') && ~force 
    disp('Unit data already computed! Loading file.');
    load([parentDir 'Compiled\mazeVsHC' tag '.mat']);
else
    for m = 1:length(mice)

        cd(strcat(parentDir, mice{m}));
        allSess = dir('*_sess*');

        %% Start collecting data
        for ii = 1:size(allSess,1)

           cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
           [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    

           %% Load spikes
           if ~isempty(dir('*Kilosort*')) &&  ~isempty(dir('summ'))
                spikes = bz_LoadPhy('noPrompts',true);
           else
               continue
           end

           %% Load pulses
            if exist([sessionInfo.FileName '.pulses.events.mat'],'file') 
                load([sessionInfo.FileName '.pulses.events.mat']);
            else
                continue
            end

            %% Load cell metrics
            if exist([sessionInfo.FileName '.cell_metrics.cellinfo.mat'],'file') 
                load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
            else
                continue
            end

            file = dir(('*.SessionPulses.Events.mat'));
            load(file.name);
            file = dir(('*.SessionArmChoice.Events.mat'));
            load(file.name);         

            efields = fieldnames(sessionPulses);    
            baseRate = [];
            rateMazeBaseStem = [];
            rateMazeBaseRet = [];
            
            for flds = 1:length(efields)

                region = sessionPulses.(efields{flds}).region; %1 is CA1/CA3, 2 is mec, 3 is both
                target = sessionPulses.(efields{flds}).target; %1 is stem, 2 is return

                rewardTS = sessionArmChoice.(efields{flds}).timestamps;
                startDelay = sessionArmChoice.(efields{flds}).delay.timestamps(1,:)';     
                endDelay = sessionArmChoice.(efields{flds}).delay.timestamps(2,:)';  

                % Get pulse times for maze specific modality in the homecage
                if exist('pulses')
                    if region <= 2
                        pulTr = (pulses.stimComb==region);
                    else
                        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==region);
                    end
                end

                %Only select the pulses that happened in the home cage,
                %i.e., which were 5 seconds long
                homeCagePulse = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
                homeCagePulseidx = homeCagePulse < 5.05 & homeCagePulse > 4.95;
                pulTr = pulTr & homeCagePulseidx;
                st = pulses.intsPeriods(1,pulTr)';           

                for zz = 1:6
                    %Extract relevant intervals for calculating spike rates
                    switch zz
                        case 1  %First, no stim trials, return        
                            startTS = rewardTS(sessionPulses.(efields{flds}).stim(1:(end-1))==0);
                            endTS = startDelay(sessionPulses.(efields{flds}).stim(1:(end-1))==0);
                            events{zz} = [startTS'; endTS'];
                        case 2  %No stim, stem
                            startTS = endDelay(sessionPulses.(efields{flds}).stim(1:(end-1))==0);        
                            endTS = rewardTS(find(sessionPulses.(efields{flds}).stim(1:(end-1))==0)+1);
                            events{zz} = [startTS';endTS'];
                        case 3 %No stim, delay
                            startTS = startDelay(sessionPulses.(efields{flds}).stim(1:(end-1))==0);        
                            endTS = endDelay(sessionPulses.(efields{flds}).stim(1:(end-1))==0); 
                            events{zz} = [startTS';endTS'];  
                        case 4  % Stim, return
                            startTS = rewardTS(sessionPulses.(efields{flds}).stim(1:(end-1))==1);
                            endTS = startDelay(sessionPulses.(efields{flds}).stim(1:(end-1))==1);
                            events{zz} = [startTS';endTS'];                    
                        case 5   % Stim, stem
                            startTS = endDelay(sessionPulses.(efields{flds}).stim(1:(end-1))==1);        
                            endTS = rewardTS(find(sessionPulses.(efields{flds}).stim(1:(end-1))==1)+1);
                            events{zz} = [startTS';endTS'];                      
                        case 6    %stim, delay
                            startTS = startDelay(sessionPulses.(efields{flds}).stim(1:(end-1))==1);        
                            endTS = endDelay(sessionPulses.(efields{flds}).stim(1:(end-1))==1); 
                            events{zz} = [startTS';endTS'];
                    end       
                end

                for jj = 1:size(spikes.UID,2)

                    % Get psth
                    [stccg] = CCG({spikes.times{jj} st},[],'binSize',0.1,'duration',12);
                    
                    if flds == 1
                        compiledData.fr_ratio_hc(totCells) = nanmean(stccg(61:110,2,1))./nanmean(stccg(1:60,2,1));
                        baseRate(jj) = nanmean(stccg(1:60,2,1));
                        compiledData.analogCh(totCells) = 2;
                        compiledData.tarZone(totCells) = target;
                        
                        %add cell metrics
                        compiledData.FR(totCells) = cell_metrics.firingRate(jj); 

                        % Assign numerical tag to putative class
                        if strcmp(cell_metrics.putativeCellType{jj},'Narrow Interneuron') == 1
                            compiledData.putativeClass(totCells) = 1;
                        elseif strcmp(cell_metrics.putativeCellType{jj},'Pyramidal Cell') == 1
                            compiledData.putativeClass(totCells) = 2;
                        elseif strcmp(cell_metrics.putativeCellType{jj},'Wide Interneuron') == 1
                            compiledData.putativeClass(totCells) = 3;
                        else 
                            compiledData.putativeClass(totCells) = 4;                       
                        end                    

                        % Assign numerical tag to region
                        if strcmp(cell_metrics.brainRegion{jj},'CA1') == 1
                            compiledData.region(totCells) = 1;
                        elseif strcmp(cell_metrics.brainRegion{jj},'DG') == 1
                            compiledData.region(totCells) = 2;
                        elseif strcmp(cell_metrics.brainRegion{jj},'CA3') == 1
                            compiledData.region(totCells) = 3;
                        else 
                            compiledData.region(totCells) = 4;                       
                        end

                        %Keep animal id data
                        compiledData.mice(totCells) = m;   

                        % Finally, get rate ratio in the relevant section of the
                        % maze
                        for zz = 1:6
                            [status] = InIntervals(spikes.times{jj},events{zz}');
                            totTime = sum(events{zz}(2,:)-events{zz}(1,:));
                            rateMaze(zz) = sum(status)/totTime;
                        end
                        
                        rateMazeBaseStem(jj) = rateMaze(2);
                        rateMazeBaseRet(jj) = rateMaze(1);
                        
                        if target == 1 % i.e., stem
                            compiledData.fr_ratio_maze(totCells) = rateMaze(5)/rateMaze(2);
                            compiledData.fr_ratio_control(totCells) = rateMaze(4)/rateMaze(1);
                        elseif target == 2
                            compiledData.fr_ratio_maze(totCells) = rateMaze(4)/rateMaze(1);
                            compiledData.fr_ratio_control(totCells) = rateMaze(5)/rateMaze(2);                    
                        end
                        totCells = totCells+1;
                        
                    elseif flds==2
                        for loop = 1:2
                            if loop==1
                                compiledData.fr_ratio_hc(totCells) = nanmean(stccg(1:60,2,1))./baseRate(jj);                        
                                compiledData.analogCh(totCells) = 1;
                            else
                                compiledData.fr_ratio_hc(totCells) = nanmean(stccg(61:110,2,1))./baseRate(jj);                        
                                compiledData.analogCh(totCells) = 3;
                            end
                            compiledData.tarZone(totCells) = target;                        

                            %add cell metrics
                            compiledData.FR(totCells) = cell_metrics.firingRate(jj); 

                            % Assign numerical tag to putative class
                            if strcmp(cell_metrics.putativeCellType{jj},'Narrow Interneuron') == 1
                                compiledData.putativeClass(totCells) = 1;
                            elseif strcmp(cell_metrics.putativeCellType{jj},'Pyramidal Cell') == 1
                                compiledData.putativeClass(totCells) = 2;
                            elseif strcmp(cell_metrics.putativeCellType{jj},'Wide Interneuron') == 1
                                compiledData.putativeClass(totCells) = 3;
                            else 
                                compiledData.putativeClass(totCells) = 4;                       
                            end                    

                            % Assign numerical tag to region
                            if strcmp(cell_metrics.brainRegion{jj},'CA1') == 1
                                compiledData.region(totCells) = 1;
                            elseif strcmp(cell_metrics.brainRegion{jj},'DG') == 1
                                compiledData.region(totCells) = 2;
                            elseif strcmp(cell_metrics.brainRegion{jj},'CA3') == 1
                                compiledData.region(totCells) = 3;
                            else 
                                compiledData.region(totCells) = 4;                       
                            end

                            %Keep animal id data
                            compiledData.mice(totCells) = m;   

                            % Finally, get rate ratio in the relevant section of the
                            % maze
                            for zz = 1:6
                                [status] = InIntervals(spikes.times{jj},events{zz}');
                                totTime = sum(events{zz}(2,:)-events{zz}(1,:));
                                rateMaze(zz) = sum(status)/totTime;
                            end

                            if loop == 1
                                if target == 1 % i.e., stem
                                    compiledData.fr_ratio_maze(totCells) = rateMaze(2)/rateMazeBaseStem(jj);
                                    compiledData.fr_ratio_control(totCells) = rateMaze(1)/rateMazeBaseRet(jj);
                                elseif target == 2
                                    compiledData.fr_ratio_maze(totCells) = rateMaze(1)/rateMazeBaseRet(jj);
                                    compiledData.fr_ratio_control(totCells) = rateMaze(2)/rateMazeBaseStem(jj);                    
                                end
                            elseif loop==2
                                if target == 1 % i.e., stem
                                    compiledData.fr_ratio_maze(totCells) = rateMaze(5)/rateMazeBaseStem(jj);
                                    compiledData.fr_ratio_control(totCells) = rateMaze(4)/rateMazeBaseRet(jj);
                                elseif target == 2
                                    compiledData.fr_ratio_maze(totCells) = rateMaze(4)/rateMazeBaseRet(jj);
                                    compiledData.fr_ratio_control(totCells) = rateMaze(5)/rateMazeBaseStem(jj);                    
                                end                             
                            end
                                
                            totCells = totCells+1;
                        end
                    end
                end

            end
        end
    end
    if saveMat
        save([parentDir 'Compiled\mazeVsHC' tag '.mat'], 'compiledData');
    end
end

col = [0.9856 0.7372 0.2537;...
       0.1986 0.7214 0.6310;...
       0.2305 0.2510 0.6173];
   
figure
set(gcf,'Renderer','painters')
for jj = 1:3 %Across cell types
    statData = [];
    statDataID = [];
    for ii = 1:3 % Across analog channels
        idx = compiledData.analogCh==ii & compiledData.region==1 & compiledData.putativeClass == jj;
        subplot(3,3,3*(jj-1)+1)
        scatter(compiledData.fr_ratio_hc(idx),compiledData.fr_ratio_maze(idx),5,col(ii,:),'filled','MarkerFaceAlpha',0.8)
        hold on   
        refline(1,0)
        xlim([0 2])
        xlabel('Home Cage Ratio')
        ylim([0 2])
        ylabel('Maze Ratio')
        
        subplot(3,3,3*(jj-1)+2)
        scatter(compiledData.fr_ratio_control(idx),compiledData.fr_ratio_maze(idx),5,col(ii,:),'filled','MarkerFaceAlpha',0.8)
        hold on   
        refline(1,0)
        xlim([0 2])
        xlabel('Control Ratio')
        ylim([0 2])
        ylabel('Maze Ratio')
        
        subplot(3,3,3*(jj-1)+3)
        h = cdfplot(compiledData.fr_ratio_maze(idx));
        set(h,'Color',col(ii,:),'LineWidth',2)
        hold on   
        xlim([0 2])
        
        statData = [statData compiledData.fr_ratio_maze(idx)];
        statDataID = [statDataID ones(1,sum(idx))*ii];
    end
    stats{jj,ii} = groupStats(statData,statDataID,'doPlot',false);
    
end

if savePlot
    saveas(gcf,[parentDir 'Compiled\mazeVsHC\mazeVsHC' tag '.png']);
    saveas(gcf,[parentDir 'Compiled\mazeVsHC\mazeVsHC' tag '.eps'],'epsc');
    saveas(gcf,[parentDir 'Compiled\mazeVsHC\mazeVsHC' tag '.fig'],'fig');
    save(strcat(parentDir,'Compiled\mazeVsHC\mazeVsHC',tag,'.mat'),'stats') 
end
end