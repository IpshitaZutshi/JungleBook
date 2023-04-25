% Compile data across all sessions

function [PF, sessCorr] = compileMicePlaceFields(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
addParameter(p,'savePlot',false,@islogical);
addParameter(p,'useZScore',false,@islogical);
addParameter(p,'downsample',false,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;
savePlot = p.Results.savePlot;
useZScore = p.Results.useZScore;
downsample = p.Results.downsample;

tag = 'CA3'; % or mEC
YlGnBu=cbrewer('seq', 'YlGnBu', 11);

if strcmp(tag,'CA1') == 1
    mice = {'IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};%'IZ15\Final'
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ17\Final','IZ18\Final','IZ20\Final'...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
        'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Final'};  % To add: 'IZ15\Final'
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ32\Saline','IZ33\Saline','IZ34\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Both'};
end


if ~isempty(analogEv)
    for ii = 1:numAnalog
        analogCh(ii) = (analogEv-1)+ii;
    end
end

for rr = 1:3
    for cc = 1:2
        for zz = 1:8 %8 fields, left trajectories, right trajectories, stim left and right trajectories (correct & incorrect)
            PF.rateMap{rr,cc}{zz} = [];
            PF.sessCorrMap{rr,cc}{zz} = [];
            PF.peakRate{rr,cc}{zz} = [];
            PF.avgRate{rr,cc}{zz} = [];
            PF.mice{rr,cc}{zz} = [];
            PF.placefield{rr,cc}{zz} = [];
            PF.putativeClass{rr,cc}{zz} = [];
            PF.region{rr,cc}{zz} = [];
        end
    end
end

%% Loop through the mice
for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');

    % Start collecting data
    for ii = 1:size(allSess,1)
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
        load([sessionInfo.FileName '.cell_metrics.cellinfo.mat']);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);
        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return        
            if ~downsample
                load(strcat(efields{jj},'\',allSess(ii).name,'.firingMapsAvg.cellinfo.mat'))
                load(strcat(efields{jj},'\',allSess(ii).name,'.placeFields.cellinfo.mat'))
            else
                load(strcat(efields{jj},'\',allSess(ii).name,'_DS.firingMapsAvg.cellinfo.mat'))
                load(strcat(efields{jj},'\',allSess(ii).name,'_DS.placeFields.cellinfo.mat'))
            end
            
            for kk = 1:length(firingMaps.rateMaps)
                %if ~isnan(firingMaps.rateMaps{kk,1}{1})
                for zz = 1:size(firingMaps.rateMaps{1},2)
                    PF.rateMap{region,target}{zz} = [PF.rateMap{region,target}{zz};firingMaps.rateMaps{kk,1}{1,zz}];   
                    PF.sessCorrMap{region,target}{zz} = [PF.sessCorrMap{region,target}{zz};firingMaps.sessCorrMaps{kk,1}{1,zz}];                                
                    PF.peakRate{region,target}{zz} = [PF.peakRate{region,target}{zz};max(firingMaps.rateMaps{kk,1}{1,zz})];
                    PF.avgRate{region,target}{zz} = [PF.avgRate{region,target}{zz};nanmean(firingMaps.rateMaps{kk,1}{1,zz})];    
                    PF.mice{region,target}{zz} = [PF.mice{region,target}{zz}; m];
                    % Assign numerical tag to putative class
                    if strcmp(cell_metrics.putativeCellType{kk},'Narrow Interneuron') == 1
                        PF.putativeClass{region,target}{zz} = [PF.putativeClass{region,target}{zz}; 1];
                    elseif strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1
                        PF.putativeClass{region,target}{zz} = [PF.putativeClass{region,target}{zz}; 2];
                    elseif strcmp(cell_metrics.putativeCellType{kk},'Wide Interneuron') == 1
                        PF.putativeClass{region,target}{zz} = [PF.putativeClass{region,target}{zz}; 3];
                    else 
                        PF.putativeClass{region,target}{zz} = [PF.putativeClass{region,target}{zz}; 4];                       
                    end                    

                    % Assign numerical tag to region
                    if strcmp(cell_metrics.brainRegion{kk},'CA1') == 1
                        PF.region{region,target}{zz} = [PF.region{region,target}{zz}; 1];
                    elseif strcmp(cell_metrics.brainRegion{kk},'DG') == 1
                        PF.region{region,target}{zz} = [PF.region{region,target}{zz}; 2];
                    elseif strcmp(cell_metrics.brainRegion{kk},'CA3') == 1
                        PF.region{region,target}{zz} = [PF.region{region,target}{zz}; 3];
                    else 
                        PF.region{region,target}{zz} = [PF.region{region,target}{zz}; 4];                       
                    end

                    if zz<5
                        if ~isnan(placeFieldStats.mapStats{kk,1}{1,zz}.x(1))
                            PF.placefield{region,target}{zz} = [PF.placefield{region,target}{zz}; 1];
                        else
                            PF.placefield{region,target}{zz} = [PF.placefield{region,target}{zz}; 0];
                        end
                    else
                        if ~isnan(placeFieldStats.mapStats{kk,1}{1,(zz-4)}.x(1))
                            PF.placefield{region,target}{zz} = [PF.placefield{region,target}{zz}; 1];
                        else
                            PF.placefield{region,target}{zz} = [PF.placefield{region,target}{zz}; 0];
                        end
                    end
                          
                end
            end
        end
    end
end

target = {'STEM', 'RETURN'};
zone = {'LeftOFF','LeftON','RightOFF','RightON'};%,'LeftOFFIN','LeftONIN','RightOFFIN','RightONIN'};

%Plot differently depending on whether plotting mEC only, or manipulations
%with separate analog channels, or CA3
if strcmp(tag,'CA1') == 1 || strcmp(tag,'mECBilateral') == 1 || strcmp(tag,'mEC')==1
    colMat = [85/243 85/243 85/243;...
            224/243 163/243 46/243;...
            8/243 133/243 161/243;...
            56/243 61/243 150/243];  
    for jj = 2%1:length(target)            
        figure(jj)
        set(gcf,'renderer','painters');
        set(gcf,'Position',[50 50 1000 800])
        for ii = 1:length(reg)
            PFbase = [];
            PFstim = [];
            for zz = 1:length(zone)
                % Separate CA1, combine DG & CA3, select principal cells with place
                % fields
                if zz == 1 || zz == 2
                    idxSelect = PF.putativeClass{ii,jj}{1,zz} == 2 & PF.region{ii,jj}{1,zz} == 1 &...
                        (PF.placefield{ii,jj}{1,1} == 1 | PF.placefield{ii,jj}{1,2} == 1);
                    PFtemp = PF.rateMap{ii,jj}{1,zz}(idxSelect,:);
                elseif zz == 3 || zz == 4
                    idxSelect = PF.putativeClass{ii,jj}{1,zz} == 2 & PF.region{ii,jj}{1,zz} == 1 &...
                        (PF.placefield{ii,jj}{1,3} == 1 | PF.placefield{ii,jj}{1,4} == 1);
                    PFtemp = PF.rateMap{ii,jj}{1,zz}(idxSelect,:);
                end
                
                if zz == 1 || zz ==3
                    [maxPF,idx] = max(PFtemp,[],2);
                    [~,sortidx] = sort(idx,'ascend');
                    PFbase = [PFbase;PFtemp(sortidx,:)];
                end
                if zz == 2 || zz == 4
                    PFstim = [PFstim;PFtemp(sortidx,:)];
                end                
            end
 
            subplot(length(reg),4,4*(ii-1)+1)                
            %imagesc(PFtemp(sortidx,:)./maxPF(sortidx))
            if ~useZScore
                imagesc(zscore(PFbase,0,2))
                colormap(YlGnBu)
                colorbar
            else
                temp1 = zscore(PFbase(:,1:63),[],2);
                temp2 = zscore(PFbase(:,64:100),[],2);
                temp = [temp1 temp2];                
                imagesc(temp)
                colormap(YlGnBu)
                colorbar
            end

            caxis([-1 5])
            hold on
            ret = round((size(PFtemp,2)/175)*80);
            stem = round((size(PFtemp,2)/175)*110);
            line([ret ret],[1 size(PFbase,1)],'Color','red','LineWidth',1.5)
            line([stem stem],[1 size(PFbase,1)],'Color','red','LineWidth',1.5)
            title(strcat(target(jj),' Baseline'));

            subplot(length(reg),4,4*(ii-1)+2)                
            %imagesc(PFtemp(sortidx,:)./maxPF(sortidx))
            if ~useZScore
                imagesc(zscore(PFstim,0,2))
                colormap(YlGnBu)
                colorbar
            else
                temp1 = zscore(PFstim(:,1:63),[],2);
                temp2 = zscore(PFstim(:,64:100),[],2);
                temp = [temp1 temp2];                
                imagesc(temp)
                colormap(YlGnBu)
                colorbar
            end

            caxis([-0.5 5])
            hold on
            ret = round((size(PFtemp,2)/175)*80);
            stem = round((size(PFtemp,2)/175)*110);
            line([ret ret],[1 size(PFstim,1)],'Color','red','LineWidth',1.5)
            line([stem stem],[1 size(PFstim,1)],'Color','red','LineWidth',1.5)
            title(strcat(target(jj),' Stim'));
            
            subplot(length(reg),4,4*(ii-1)+3) 
            hold on
%             if ~useZScore
%                 avgRMbase = nanmedian(PFbase,1);
%                 stderrRMbase = 1.2533*(nanstd(PFbase,0,1)/sqrt(size(PFbase,1)));
%             else
%                 temp1 = zscore(PFbase(:,1:63),[],2);
%                 temp2 = zscore(PFbase(:,64:100),[],2);
%                 temp = [temp1 temp2];
%                 avgRMbase = nanmedian(temp,1);
%                 stderrRMbase = 1.2533*(nanstd(temp,0,1)/sqrt(size(temp,1)));                
%             end            
            avgRM = nanmedian(PFbase,1);
            stderrRM = 1.2533*(nanstd(PFbase,0,1)/sqrt(size(PFbase,1)));
            fill([(1:1:length(avgRM))'; (length(avgRM):-1:1)'],[avgRM'-stderrRM';flipud(avgRM'+stderrRM')],colMat(1,:),'linestyle','none','FaceAlpha',0.2);
            hold on
            plot(avgRM,'color',colMat(1,:),'LineWidth',1.5)           

            avgRM = nanmedian(PFstim,1);         
            stderrRM = 1.2533*(nanstd(PFstim,0,1)/sqrt(size(PFstim,1)));
            fill([(1:1:length(avgRM))'; (length(avgRM):-1:1)'],[avgRM'-stderrRM';flipud(avgRM'+stderrRM')],colMat(ii+1,:),'linestyle','none','FaceAlpha',0.2);            hold on
            plot(avgRM,'color',colMat(ii+1,:),'LineWidth',1.5)            
            line([ret ret],[min(avgRM)-0.5 max(avgRM)+0.5],'Color','red','LineWidth',1)
            line([stem stem],[min(avgRM)-0.5 max(avgRM)+0.5],'Color','red','LineWidth',1)
            title(strcat('Firing rate ',target(jj)));
            if ~isempty(avgRM)
                ylim([min(avgRM)-0.5 max(avgRM)+0.5])
                xlim([1 length(avgRM)])
            end
            ylabel('Median pop. FR')
            xlabel('Position along track')      
            
            %Get rate map stability
            Numfields{ii,jj}(1,:) = [sum(PF.placefield{ii,jj}{1,1})/length(PF.placefield{ii,jj}{1,1})...
                sum(PF.placefield{ii,jj}{1,2})/length(PF.placefield{ii,jj}{1,2})... 
                sum(PF.placefield{ii,jj}{1,3})/length(PF.placefield{ii,jj}{1,3})... 
                sum(PF.placefield{ii,jj}{1,4})/length(PF.placefield{ii,jj}{1,4})];

            idxPF = PF.placefield{ii,jj}{1,1} == 1 | PF.placefield{ii,jj}{1,2} == 1 | PF.placefield{ii,jj}{1,3} == 1 | PF.placefield{ii,jj}{1,4} == 1;
            idxSelect = PF.putativeClass{ii,jj}{1,1} == 2 & PF.region{ii,jj}{1,1} == 1 & idxPF == 1;

            % Correlations
%             sessCorr{ii,jj}(:,1) = diag(corr(PF.rateMap{ii,jj}{1,1}(idxSelect,10:ret)',PF.rateMap{ii,jj}{1,2}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
%             sessCorr{ii,jj}(:,2) = diag(corr(PF.rateMap{ii,jj}{1,1}(idxSelect,stem:95)',PF.rateMap{ii,jj}{1,2}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
%             sessCorr{ii,jj}(:,3) = diag(corr(PF.rateMap{ii,jj}{1,3}(idxSelect,10:ret)',PF.rateMap{ii,jj}{1,4}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
%             sessCorr{ii,jj}(:,4) = diag(corr(PF.rateMap{ii,jj}{1,3}(idxSelect,stem:95)',PF.rateMap{ii,jj}{1,4}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));            
            if ~useZScore
                sessCorr{ii,jj}(:,1) = diag(corr(PF.sessCorrMap{ii,jj}{1}(idxSelect,10:ret)',PF.sessCorrMap{ii,jj}{4}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,2) = diag(corr(PF.sessCorrMap{ii,jj}{1}(idxSelect,stem:95)',PF.sessCorrMap{ii,jj}{4}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,3) = diag(corr(PF.sessCorrMap{ii,jj}{5}(idxSelect,10:ret)',PF.sessCorrMap{ii,jj}{8}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,4) = diag(corr(PF.sessCorrMap{ii,jj}{5}(idxSelect,stem:95)',PF.sessCorrMap{ii,jj}{8}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,5) = diag(corr(PF.sessCorrMap{ii,jj}{1}(idxSelect,10:ret)',PF.sessCorrMap{ii,jj}{2}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,6) = diag(corr(PF.sessCorrMap{ii,jj}{1}(idxSelect,stem:95)',PF.sessCorrMap{ii,jj}{2}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,7) = diag(corr(PF.sessCorrMap{ii,jj}{5}(idxSelect,10:ret)',PF.sessCorrMap{ii,jj}{6}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,8) = diag(corr(PF.sessCorrMap{ii,jj}{5}(idxSelect,stem:95)',PF.sessCorrMap{ii,jj}{6}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            else
                sessCorr{ii,jj}(:,1) = diag(corr(zscore(PF.sessCorrMap{ii,jj}{1}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{ii,jj}{4}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,2) = diag(corr(zscore(PF.sessCorrMap{ii,jj}{1}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{ii,jj}{4}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,3) = diag(corr(zscore(PF.sessCorrMap{ii,jj}{5}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{ii,jj}{8}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,4) = diag(corr(zscore(PF.sessCorrMap{ii,jj}{5}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{ii,jj}{8}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,5) = diag(corr(zscore(PF.sessCorrMap{ii,jj}{1}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{ii,jj}{2}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,6) = diag(corr(zscore(PF.sessCorrMap{ii,jj}{1}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{ii,jj}{2}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,7) = diag(corr(zscore(PF.sessCorrMap{ii,jj}{5}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{ii,jj}{6}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
                sessCorr{ii,jj}(:,8) = diag(corr(zscore(PF.sessCorrMap{ii,jj}{5}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{ii,jj}{6}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
                
            end
            if ~useZScore
                for m = unique(PF.mice{ii,jj}{1})'
                    idxMouse  = PF.mice{ii,jj}{1} == m;
                    corrbase(1,:,m) = diag(corr(PF.sessCorrMap{ii,jj}{1}(idxMouse,:),PF.sessCorrMap{ii,jj}{2}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                    corrbase(2,:,m) = diag(corr(PF.sessCorrMap{ii,jj}{5}(idxMouse,:),PF.sessCorrMap{ii,jj}{6}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                    corrstim(1,:,m) = diag(corr(PF.sessCorrMap{ii,jj}{1}(idxMouse,:),PF.sessCorrMap{ii,jj}{4}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                    corrstim(2,:,m) = diag(corr(PF.sessCorrMap{ii,jj}{5}(idxMouse,:),PF.sessCorrMap{ii,jj}{8}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                end
            else 
                clear temp
                for nn = 1:8
                    temp1 = zscore(PF.sessCorrMap{ii,jj}{nn}(idxSelect,1:63),[],2);
                    temp2 = zscore(PF.sessCorrMap{ii,jj}{nn}(idxSelect,64:100),[],2);
                    temp{nn} = [temp1 temp2];
                end
                corrbase(1,:) = diag(corr(temp{1},temp{2},'Rows','pairwise','Type','Spearman'));
                corrbase(2,:) = diag(corr(temp{5},temp{6},'Rows','pairwise','Type','Spearman'));
                corrstim(1,:) = diag(corr(temp{2},temp{4},'Rows','pairwise','Type','Spearman'));
                corrstim(2,:) = diag(corr(temp{6},temp{8},'Rows','pairwise','Type','Spearman'));
            end
%             corrstim(1,:) = diag(corr(PF.rateMap{ii,jj}{1,1},PF.rateMap{ii,jj}{1,2},'Rows','pairwise','Type','Spearman'));
%             corrstim(2,:) = diag(corr(PF.rateMap{ii,jj}{1,3},PF.rateMap{ii,jj}{1,4},'Rows','pairwise','Type','Spearman'));      
            
            corrRM(1,:,:) = nanmean(corrbase,1);         
            corrRM(2,:,:) = nanmean(corrstim,1);
            
            corrPV(1,:,:) = corrRM(2,:,:)-corrRM(1,:,:);
            
            avgCorr = nanmean(corrPV,3);
            stderrCorr = nanstd(corrPV,[],3)./sqrt(size(corrPV,3));
            subplot(length(reg),4,4*(ii-1)+4)    
            col = colMat(ii+1,:);
            fill([(1:1:length(avgCorr))'; (length(avgCorr):-1:1)'],[avgCorr'-stderrCorr';flipud(avgCorr'+stderrCorr')],col,'linestyle','none','FaceAlpha',0.2);                   
            hold on       
            line([ret ret],[-0.4 0.2],'Color','red')
            line([stem stem],[-0.4 0.2],'Color','red')
            plot(avgCorr,'color',col,'LineWidth',1.5);
            ylim([-0.4 0.2])
            xlim([1 length(avgCorr)])        
            xlabel('Position along track')
            ylabel('Correlation')
            
            corrStats(1:size(corrPV,3),2) = nanmean(corrPV(1,65:90,:),2);
            corrStats(1:size(corrPV,3),1) = nanmean(corrPV(1,15:40,:),2);
            for ss = 1:2
                [~,PVstats{ii}.p(ss),~, PVstats{ii}.stats{ss}] = ttest(corrStats(:,ss));
            end
            
            title(strcat('Side:',num2str(PVstats{ii}.p(1)),' Center:',num2str(PVstats{ii}.p(2))));
            clear corrbase corrstim corrPV corrRM
        end
        if savePlot
            if ~useZScore && ~downsample
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'.png'));
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'.eps'),'epsc');
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'.fig'));
            elseif downsample
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'DS.png'));
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'DS.eps'),'epsc');
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'DS.fig'));                
            else
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'ZScore.png'));
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'ZScore.eps'),'epsc');
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField',target{jj},'_',tag,'ZScore.fig'));
            end
        end
    end
    
    sessCorrmat = [];
    sessCorrIDzone = [];
    sessCorrManip = [];
    
    for ii = 1:length(reg)
        for jj = 1%:length(target)

          sessCorrmat = [sessCorrmat;sessCorr{ii,jj}(:,1);sessCorr{ii,jj}(:,3);sessCorr{ii,jj}(:,2);sessCorr{ii,jj}(:,4)];%...
%               sessCorr{ii,jj}(:,5);sessCorr{ii,jj}(:,7);sessCorr{ii,jj}(:,6);sessCorr{ii,jj}(:,8)];
          sessCorrIDzone = [sessCorrIDzone;ones(size(sessCorr{ii,jj}(:,1),1),1)*1;ones(size(sessCorr{ii,jj}(:,3),1),1)*1;...
              ones(size(sessCorr{ii,jj}(:,2),1),1)*2;ones(size(sessCorr{ii,jj}(:,4),1),1)*2];%...
%               ones(size(sessCorr{ii,jj}(:,5),1),1)*3;ones(size(sessCorr{ii,jj}(:,7),1),1)*3;...
%               ones(size(sessCorr{ii,jj}(:,6),1),1)*4;ones(size(sessCorr{ii,jj}(:,8),1),1)*4];
          sessCorrManip = [sessCorrManip; ones(size(sessCorr{ii,jj}(:,1),1),1)*ii;ones(size(sessCorr{ii,jj}(:,3),1),1)*ii;...
              ones(size(sessCorr{ii,jj}(:,2),1),1)*ii;ones(size(sessCorr{ii,jj}(:,4),1),1)*ii];%...
%               ones(size(sessCorr{ii,jj}(:,5),1),1)*ii;ones(size(sessCorr{ii,jj}(:,7),1),1)*ii;...
%               ones(size(sessCorr{ii,jj}(:,6),1),1)*ii;ones(size(sessCorr{ii,jj}(:,8),1),1)*ii];
        end
    end
    
    fig = figure;
    set(gcf,'renderer','painters');
    colMat = [85/243 85/243 85/243;224/243 163/243 46/243;... 
        85/243 85/243 85/243;8/243 133/243 161/243;... 
        85/243 85/243 85/243;56/243 61/243 150/243];
    
    stats = groupStats(sessCorrmat,[sessCorrManip sessCorrIDzone],'inAxis',true,'Color',colMat);
    
    for ii = 1:length(reg)
        sessCorrmat = [];
        sessCorrmat = [sessCorr{ii,1}(:,[1 2]); sessCorr{ii,1}(:,[3 4])];
        [stats.signrank.p(ii),~,stats.signrank.stats{ii}] = signrank(sessCorrmat(:,1),sessCorrmat(:,2));
    end
    stats.PVcorr = PVstats;
    
    if savePlot
        if ~useZScore && ~downsample
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'.png'));
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'.eps'),'epsc');
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'.fig'));
            save(strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'.mat'),'stats');
        elseif downsample
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'DS.png'));
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'DS.eps'),'epsc');
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'DS.fig'));
            save(strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'DS.mat'),'stats');            
        else
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'ZScore.png'));
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'ZScore.eps'),'epsc');
            saveas(fig,strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'ZScore.fig'));
            save(strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'ZScore.mat'),'stats');
        end
    end

elseif strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    
    colMat = [85/243 85/243 85/243;...
            8/243 133/243 161/243;...
            187/243 86/243 149/243;...                   
            56/243 61/243 150/243;...
            90/243 60/243 108/243];     
    figure(1)
    set(gcf,'renderer','painters');
    set(gcf,'Position',[50 50 1200 650])   
    titlename = {'Baseline','mEC','CA3','Both'};
    for jj = 1%:length(target)  
        figure(jj)
        set(gcf,'renderer','painters');
        set(gcf,'Position',[50 50 1000 800])
        for ii = 2:3
            PFbase = [];
            PFstim = []; 
            for zz = 1:4               
                % Separate CA1, select principal cells with place fields
%                 if zz == 1 || zz == 2
%                     idxSelect = PF.putativeClass{2,jj}{1,zz} == 2 & PF.region{2,jj}{1,zz} == 1 &...
%                         (PF.placefield{2,jj}{1,1} == 1 | PF.placefield{2,jj}{1,2} == 1 | PF.placefield{3,jj}{1,1} == 1 | PF.placefield{3,jj}{1,1} == 1);
%                     PFtemp = PF.rateMap{ii,jj}{1,zz}(idxSelect,:);
%                 elseif zz == 3 || zz == 4
%                     idxSelect = PF.putativeClass{2,jj}{1,zz} == 2 & PF.region{2,jj}{1,zz} == 1 &...
%                         (PF.placefield{2,jj}{1,3} == 1 | PF.placefield{2,jj}{1,4} == 1 | PF.placefield{3,jj}{1,3} == 1 | PF.placefield{3,jj}{1,4} == 1);
%                     PFtemp = PF.rateMap{ii,jj}{1,zz}(idxSelect,:);
%                 end
                if zz == 1 || zz == 2
                    idxSelect = PF.putativeClass{2,jj}{1,zz} == 1 & PF.region{2,jj}{1,zz} == 1;
                    PFtemp = PF.rateMap{ii,jj}{1,zz}(idxSelect,:);
                elseif zz == 3 || zz == 4
                    idxSelect = PF.putativeClass{2,jj}{1,zz} == 1 & PF.region{2,jj}{1,zz} == 1;                      
                    PFtemp = PF.rateMap{ii,jj}{1,zz}(idxSelect,:);
                end

                if zz == 1 || zz ==3
                    PFsorting = PF.rateMap{2,jj}{1,zz}(idxSelect,:);
                    [~,idx] = max(PFsorting,[],2);
                    [~,sortidx] = sort(idx,'ascend');
                    PFbase = [PFbase;PFtemp(sortidx,:)];
                end
                if zz == 2 || zz == 4
                    PFstim = [PFstim;PFtemp(sortidx,:)];
                end  
            end
            
            subplot(2,4,2*(ii-1)-1)                
            %imagesc(PFtemp(sortidx,:)./maxPF(sortidx))
            a = ~isnan(PFbase);
            if ~useZScore
                imagesc(zscore(PFbase(a(:,1),1:(end-1)),0,2))
                colormap(YlGnBu)
            else
                temp1 = zscore(PFbase(a(:,1),1:63),[],2);
                temp2 = zscore(PFbase(a(:,1),64:100),[],2);
                temp = [temp1 temp2];                
                imagesc(temp)
                colormap(YlGnBu)
            end  
            caxis([-0.5 5])
            hold on
            ret = round((size(PFtemp,2)/175)*80);
            stem = round((size(PFtemp,2)/175)*110);
            line([ret ret],[1 size(PFbase,1)],'Color','red','LineWidth',1.5)
            line([stem stem],[1 size(PFbase,1)],'Color','red','LineWidth',1.5)
            title(titlename{2*(ii-1)-1});

            subplot(2,4,2*(ii-1))        
            a = ~isnan(PFstim);            
            if ~useZScore
                imagesc(zscore(PFstim(a(:,1),1:(end-1)),0,2))
                colormap(YlGnBu)
            else
                temp1 = zscore(PFstim(a(:,1),1:63),[],2);
                temp2 = zscore(PFstim(a(:,1),64:100),[],2);
                temp = [temp1 temp2];                
                imagesc(temp)
                colormap(YlGnBu)
            end            
            caxis([-0.5 5])
            hold on
            ret = round((size(PFtemp,2)/175)*80);
            stem = round((size(PFtemp,2)/175)*110);
            line([ret ret],[1 size(PFstim,1)],'Color','red','LineWidth',1.5)
            line([stem stem],[1 size(PFstim,1)],'Color','red','LineWidth',1.5)
            title(titlename{2*(ii-1)});
            
            subplot(2,4,5) 
            hold on
            if ~useZScore
                avgRMbase = nanmedian(PFbase,1);
                stderrRMbase = 1.2533*(nanstd(PFbase,0,1)/sqrt(size(PFbase,1)));
            else
                temp1 = zscore(PFbase(:,1:63),[],2);
                temp2 = zscore(PFbase(:,64:100),[],2);
                temp = [temp1 temp2];
                avgRMbase = nanmedian(temp,1);
                stderrRMbase = 1.2533*(nanstd(temp,0,1)/sqrt(size(temp,1)));                    
            end
            avgRM = avgRMbase;%-avgRMbase;
            stderrRM = stderrRMbase;            
            fill([(1:1:length(avgRM))'; (length(avgRM):-1:1)'],[avgRM'-stderrRM';flipud(avgRM'+stderrRM')],colMat(2*(ii-1)-1,:),'linestyle','none','FaceAlpha',0.2);
            hold on
            plot(avgRM,'color',colMat(2*(ii-1)-1,:),'LineWidth',1.5)        
            
            if ~useZScore
                avgRMstim = nanmedian(PFstim,1);
                stderrRMstim = 1.2533*(nanstd(PFstim,0,1)/sqrt(size(PFstim,1)));
            else
                temp1 = zscore(PFstim(:,1:63),[],2);
                temp2 = zscore(PFstim(:,64:100),[],2);
                temp = [temp1 temp2];
                avgRMstim = nanmedian(temp,1);
                stderrRMstim = 1.2533*(nanstd(temp,0,1)/sqrt(size(temp,1)));                
            end
            
            avgRM = avgRMstim;%-avgRMbase;
            stderrRM = stderrRMstim;
            fill([(1:1:length(avgRM))'; (length(avgRM):-1:1)'],[avgRM'-stderrRM';flipud(avgRM'+stderrRM')],colMat(2*(ii-1),:),'linestyle','none','FaceAlpha',0.2);
            hold on
            plot(avgRM,'color',colMat(2*(ii-1),:),'LineWidth',1.5)            
            line([ret ret],[min(avgRM)-0.5 max(avgRM)+0.5],'Color','red','LineWidth',1)
            line([stem stem],[min(avgRM)-0.5 max(avgRM)+0.5],'Color','red','LineWidth',1)
            title(strcat('Firing rate ',target(jj)));
            if ~isempty(avgRM)
                ylim([min(avgRM)-0.5 max(avgRM)+0.5])
                xlim([1 length(avgRM)])
            end
            ylabel('Average pop. FR')
            xlabel('Position along track')       
        end
            
        % Get session correlations
        idxPF = PF.placefield{2,jj}{1,1} == 1 | PF.placefield{2,jj}{1,2} == 1 | PF.placefield{2,jj}{1,3} == 1 | PF.placefield{2,jj}{1,4} == 1;
        idxSelect = PF.putativeClass{2,jj}{1,1} == 1 & PF.region{2,jj}{1,1} == 1;% & idxPF == 1;
        
         % Correlations
         
        %16 column matrix - Baseline_left_return,
        %Baseline_right_return,mec_left_return, mec_right return,
        %ca3_left_return, ca3_right_return, both_left return,
        %both_right_ret both_left return (wrt ca3),
        %both_right_ret(wrt ca3)
        %% First return arm
        if ~useZScore
            sessCorrRet{jj}(:,1) = diag(corr(PF.sessCorrMap{2,jj}{1}(idxSelect,10:ret)',PF.sessCorrMap{2,jj}{2}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,2) = diag(corr(PF.sessCorrMap{2,jj}{5}(idxSelect,10:ret)',PF.sessCorrMap{2,jj}{6}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,3) = diag(corr(PF.sessCorrMap{2,jj}{1}(idxSelect,10:ret)',PF.sessCorrMap{2,jj}{4}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,4) = diag(corr(PF.sessCorrMap{2,jj}{5}(idxSelect,10:ret)',PF.sessCorrMap{2,jj}{8}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,5) = diag(corr(PF.sessCorrMap{2,jj}{2}(idxSelect,10:ret)',PF.sessCorrMap{3,jj}{2}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,6) = diag(corr(PF.sessCorrMap{2,jj}{6}(idxSelect,10:ret)',PF.sessCorrMap{3,jj}{6}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,7) = diag(corr(PF.sessCorrMap{2,jj}{2}(idxSelect,10:ret)',PF.sessCorrMap{3,jj}{4}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,8) = diag(corr(PF.sessCorrMap{2,jj}{6}(idxSelect,10:ret)',PF.sessCorrMap{3,jj}{8}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,9) = diag(corr(PF.sessCorrMap{3,jj}{2}(idxSelect,10:ret)',PF.sessCorrMap{3,jj}{4}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,10) = diag(corr(PF.sessCorrMap{3,jj}{6}(idxSelect,10:ret)',PF.sessCorrMap{3,jj}{8}(idxSelect,10:ret)','Rows','pairwise','Type','Spearman'));

            sessCorr{jj}(:,1) = diag(corr(PF.sessCorrMap{2,jj}{1}(idxSelect,stem:95)',PF.sessCorrMap{2,jj}{2}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,2) = diag(corr(PF.sessCorrMap{2,jj}{5}(idxSelect,stem:95)',PF.sessCorrMap{2,jj}{6}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,3) = diag(corr(PF.sessCorrMap{2,jj}{1}(idxSelect,stem:95)',PF.sessCorrMap{2,jj}{4}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,4) = diag(corr(PF.sessCorrMap{2,jj}{5}(idxSelect,stem:95)',PF.sessCorrMap{2,jj}{8}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,5) = diag(corr(PF.sessCorrMap{2,jj}{2}(idxSelect,stem:95)',PF.sessCorrMap{3,jj}{2}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,6) = diag(corr(PF.sessCorrMap{2,jj}{6}(idxSelect,stem:95)',PF.sessCorrMap{3,jj}{6}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,7) = diag(corr(PF.sessCorrMap{2,jj}{2}(idxSelect,stem:95)',PF.sessCorrMap{3,jj}{4}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,8) = diag(corr(PF.sessCorrMap{2,jj}{6}(idxSelect,stem:95)',PF.sessCorrMap{3,jj}{8}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,9) = diag(corr(PF.sessCorrMap{3,jj}{2}(idxSelect,stem:95)',PF.sessCorrMap{3,jj}{4}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,10) = diag(corr(PF.sessCorrMap{3,jj}{6}(idxSelect,stem:95)',PF.sessCorrMap{3,jj}{8}(idxSelect,stem:95)','Rows','pairwise','Type','Spearman'));
        else
            sessCorrRet{jj}(:,1) = diag(corr(zscore(PF.sessCorrMap{2,jj}{1}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{2,jj}{2}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,2) = diag(corr(zscore(PF.sessCorrMap{2,jj}{5}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{2,jj}{6}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,3) = diag(corr(zscore(PF.sessCorrMap{2,jj}{1}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{2,jj}{4}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,4) = diag(corr(zscore(PF.sessCorrMap{2,jj}{5}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{2,jj}{8}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,5) = diag(corr(zscore(PF.sessCorrMap{2,jj}{2}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{3,jj}{2}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,6) = diag(corr(zscore(PF.sessCorrMap{2,jj}{6}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{3,jj}{6}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,7) = diag(corr(zscore(PF.sessCorrMap{2,jj}{2}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{3,jj}{4}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,8) = diag(corr(zscore(PF.sessCorrMap{2,jj}{6}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{3,jj}{8}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,9) = diag(corr(zscore(PF.sessCorrMap{3,jj}{1}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{3,jj}{4}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorrRet{jj}(:,10) = diag(corr(zscore(PF.sessCorrMap{3,jj}{5}(idxSelect,10:ret),[],2)',zscore(PF.sessCorrMap{3,jj}{8}(idxSelect,10:ret),[],2)','Rows','pairwise','Type','Spearman'));

            sessCorr{jj}(:,1) = diag(corr(zscore(PF.sessCorrMap{2,jj}{1}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{2,jj}{2}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,2) = diag(corr(zscore(PF.sessCorrMap{2,jj}{5}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{2,jj}{6}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,3) = diag(corr(zscore(PF.sessCorrMap{2,jj}{1}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{2,jj}{4}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,4) = diag(corr(zscore(PF.sessCorrMap{2,jj}{5}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{2,jj}{8}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,5) = diag(corr(zscore(PF.sessCorrMap{2,jj}{2}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{3,jj}{2}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,6) = diag(corr(zscore(PF.sessCorrMap{2,jj}{6}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{3,jj}{6}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,7) = diag(corr(zscore(PF.sessCorrMap{2,jj}{2}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{3,jj}{4}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,8) = diag(corr(zscore(PF.sessCorrMap{2,jj}{6}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{3,jj}{8}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,9) = diag(corr(zscore(PF.sessCorrMap{3,jj}{1}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{3,jj}{4}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
            sessCorr{jj}(:,10) = diag(corr(zscore(PF.sessCorrMap{3,jj}{5}(idxSelect,stem:95),[],2)',zscore(PF.sessCorrMap{3,jj}{8}(idxSelect,stem:95),[],2)','Rows','pairwise','Type','Spearman'));
        end
        
        if ~useZScore
            for m = unique(PF.mice{ii,jj}{1})'
                idxMouse  = PF.mice{ii,jj}{1} == m;
                corrbase(1,:,m) = diag(corr(PF.sessCorrMap{2,jj}{1}(idxMouse,:),PF.sessCorrMap{2,jj}{2}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrbase(2,:,m) = diag(corr(PF.sessCorrMap{2,jj}{5}(idxMouse,:),PF.sessCorrMap{2,jj}{6}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrmEC(1,:,m) = diag(corr(PF.sessCorrMap{2,jj}{2}(idxMouse,:),PF.sessCorrMap{2,jj}{4}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrmEC(2,:,m) = diag(corr(PF.sessCorrMap{2,jj}{6}(idxMouse,:),PF.sessCorrMap{2,jj}{8}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrCA3(1,:,m) = diag(corr(PF.sessCorrMap{2,jj}{2}(idxMouse,:),PF.sessCorrMap{3,jj}{2}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrCA3(1,:,m) = diag(corr(PF.sessCorrMap{2,jj}{6}(idxMouse,:),PF.sessCorrMap{3,jj}{6}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrBoth(1,:,m) = diag(corr(PF.sessCorrMap{2,jj}{2}(idxMouse,:),PF.sessCorrMap{3,jj}{4}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrBoth(2,:,m) = diag(corr(PF.sessCorrMap{2,jj}{6}(idxMouse,:),PF.sessCorrMap{3,jj}{8}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrBothCA3(1,:,m) = diag(corr(PF.sessCorrMap{3,jj}{1}(idxMouse,:),PF.sessCorrMap{3,jj}{4}(idxMouse,:),'Rows','pairwise','Type','Spearman'));
                corrBothCA3(2,:,m) = diag(corr(PF.sessCorrMap{3,jj}{5}(idxMouse,:),PF.sessCorrMap{3,jj}{8}(idxMouse,:),'Rows','pairwise','Type','Spearman')); 
            end
        else
            clear temp2 temp3
            for nn = 1:8
                tempA = zscore(PF.sessCorrMap{2,jj}{nn}(:,1:62),[],2);
                tempB = zscore(PF.sessCorrMap{2,jj}{nn}(:,63:100),[],2);
                temp2{nn} = [tempA tempB];
                
                tempA = zscore(PF.sessCorrMap{3,jj}{nn}(:,1:62),[],2);
                tempB = zscore(PF.sessCorrMap{3,jj}{nn}(:,63:100),[],2);
                temp3{nn} = [tempA tempB];                
            end
            corrbase(1,:) = diag(corr(temp2{1},temp2{2},'Rows','pairwise','Type','Spearman'));
            corrbase(2,:) = diag(corr(temp2{5},temp2{6},'Rows','pairwise','Type','Spearman'));
            corrmEC(1,:) = diag(corr(temp2{2},temp2{4},'Rows','pairwise','Type','Spearman'));
            corrmEC(2,:) = diag(corr(temp2{6},temp2{8},'Rows','pairwise','Type','Spearman'));
            corrCA3(1,:) = diag(corr(temp2{2},temp3{2},'Rows','pairwise','Type','Spearman'));
            corrCA3(1,:) = diag(corr(temp2{6},temp3{6},'Rows','pairwise','Type','Spearman'));
            corrBoth(1,:) = diag(corr(temp2{2},temp3{4},'Rows','pairwise','Type','Spearman'));
            corrBoth(2,:) = diag(corr(temp2{6},temp3{8},'Rows','pairwise','Type','Spearman'));
            corrBothCA3(1,:) = diag(corr(temp3{2},temp3{4},'Rows','pairwise','Type','Spearman'));
            corrBothCA3(2,:) = diag(corr(temp3{6},temp3{8},'Rows','pairwise','Type','Spearman')); 
        end

        corrRM(1,:,:) = nanmean(corrbase,1);         
        corrRM(2,:,:) = nanmean(corrmEC,1);
        corrRM(3,:,:) = nanmean(corrCA3,1);
        corrRM(4,:,:) = nanmean(corrBoth,1);
        corrRM(5,:,:) = nanmean(corrBothCA3,1);
        
        corrPV(1,:,:) = corrRM(2,:,:)-corrRM(1,:,:);
        corrPV(2,:,:) = corrRM(3,:,:)-corrRM(1,:,:);
        corrPV(3,:,:) = corrRM(4,:,:)-corrRM(1,:,:);
        corrPV(4,:,:) = corrRM(5,:,:)-corrRM(1,:,:);
  
        subplot(2,4,6) 
        for mm = 1:4%5
            avgCorr = nanmean(corrPV(mm,:,:),3);
            stderrCorr = nanstd(corrPV(mm,:,:),[],3)./sqrt(size(corrPV(mm,:,:),3));
            fill([(1:1:length(avgCorr))'; (length(avgCorr):-1:1)'],[avgCorr'-stderrCorr';flipud(avgCorr'+stderrCorr')],colMat(mm+1,:),'linestyle','none','FaceAlpha',0.2);                   
            hold on                 
            hold on       
            line([ret ret],[-0.45 0.2],'Color','red')
            line([stem stem],[-0.45 0.2],'Color','red')
            plot(avgCorr,'color',colMat(mm+1,:),'LineWidth',1.5);
            
                        
            corrStats(1:size(corrPV(mm,:,:),3),2) = nanmean(corrPV(mm,65:90,:),2);
            corrStats(1:size(corrPV(mm,:,:),3),1) = nanmean(corrPV(mm,15:40,:),2);
            for ss = 1:2
                [~,PVstats{mm}.p(ss),~, PVstats{mm}.stats{ss}] = ttest(corrStats(:,ss));
            end           
            clear corrStats
        end
        ylim([-0.45 0.2])
        xlim([1 length(avgCorr)])        
        xlabel('Position along track')
        ylabel('Correlation')
        title(strcat(num2str(PVstats{1}.p(2)),'.',num2str(PVstats{2}.p(2)),'.',num2str(PVstats{3}.p(2))));        
                  
        %CorrMatReturn
        sessCorrmat = [];
        sessCorrManip = [];
        for ii = 1:10
            sessCorrmat = [sessCorrmat;sessCorrRet{jj}(:,ii)];
            if mod(ii,2) == 0
                sessCorrManip = [sessCorrManip;ones(size(sessCorrRet{jj}(:,ii),1),1)*(ii/2);ones(size(sessCorrRet{jj}(:,ii),1),1)*(ii/2)];
            end
        end
        subplot(2,4,7) 
        stats.return = groupStats(sessCorrmat,sessCorrManip,'inAxis',true,'color',colMat);
        title('Return')       

        %CorrMatStem
        sessCorrmat = [];
        sessCorrManip = [];
        for ii = 1:8%10
            sessCorrmat = [sessCorrmat;sessCorr{jj}(:,ii)];
            if mod(ii,2) == 0
                sessCorrManip = [sessCorrManip;ones(size(sessCorr{jj}(:,ii),1),1)*(ii/2);ones(size(sessCorr{jj}(:,ii),1),1)*(ii/2)];
            end
        end
        subplot(2,4,8) 
        stats.stem = groupStats(sessCorrmat,sessCorrManip,'inAxis',true,'color',colMat);
        title('Stem')
        
        stats.PVcorr = PVstats;
        if savePlot
            if ~useZScore && ~downsample
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'.png'));
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'.eps'),'epsc');
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'.fig'));
                save(strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'.mat'),'stats');
            elseif downsample
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'DS.png'));
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'DS.eps'),'epsc');
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'DS.fig'));
                save(strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'DS.mat'),'stats');
            else
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'ZScore.png'));
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'ZScore.eps'),'epsc');
                saveas(figure(jj),strcat(parentDir,'Compiled\Place Fields\PlaceField_',tag,'ZScore.fig'));
                save(strcat(parentDir,'Compiled\Place Fields\PlaceFieldSessCorr',tag,'ZScore.mat'),'stats');
            end
        end
    end       
end

end
