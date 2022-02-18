% Compile data across all sessions

function PF = SessBehaviorPlaceFields(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
parse(p,varargin{:});

expPath = p.Results.expPath;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;

%%
if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end
allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

if ~isempty(analogEv)
    for ii = 1:numAnalog
        analogCh(ii) = (analogEv-1)+ii;
    end
end

for rr = 1:3
    for cc = 1:2
        for zz = 1:4 %4 fields, left trajectories, right trajectories, stim left and right trajectories
            PF.rateMap{rr,cc}{zz} = [];
            PF.peakRate{rr,cc}{zz} = [];
            PF.avgRate{rr,cc}{zz} = [];
        end
    end
end


%% Start collecting data
for ii = 1:size(allSess,1)

    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);    
    file = dir(('*.SessionPulses.Events.mat'));
    load(file.name);
    
    efields = fieldnames(sessionPulses);    

    for jj = 1:length(efields)
        region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
        target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return
        
        load(strcat(efields{jj},'\',allSess(ii).name,'.firingMapsAvg.cellinfo.mat'))
        for kk = 1:length(firingMaps.rateMaps)
            for zz = 1:4%size(firingMaps.rateMaps{1},2)
                PF.rateMap{region,target}{zz} = [PF.rateMap{region,target}{zz};firingMaps.rateMaps{kk,1}{1,zz}];
                PF.peakRate{region,target}{zz} = [PF.peakRate{region,target}{zz};max(firingMaps.rateMaps{kk,1}{1,zz})];
                PF.avgRate{region,target}{zz} = [PF.avgRate{region,target}{zz};nanmean(firingMaps.rateMaps{kk,1}{1,zz})];
            end
        end
        
    end

end

reg = {'CA1','mEC','Both'};
target = {'STEM', 'RETURN'};
zone = {'LeftOFF','LeftON','RightOFF','RightON'};
for ii = 1:length(reg)
    figure
    set(gcf,'Position',[0 0 1000 1000])
    for jj = 1:length(target)
        for zz = 1:length(zone)
            subplot(3,4,4*(jj-1)+zz)
            if zz <=2
                zzloc = 1;
            else
                zzloc = 3;
            end
            [idxPyr] = PF.avgRate{ii,jj}{1,zzloc}<=10;
            RM = PF.rateMap{ii,jj}{1,zz}(idxPyr,:);
            if zz == 1 || zz ==3
                %exclude INs
                [~,idx] = max(RM,[],2);
                [~,sortidx] = sort(idx,'ascend');
            end
            %imagesc(zscore(PF.rateMap{ii,jj}{1,zz}(sortidx,:),0,2))
            imagesc(zscore(RM(sortidx,:),0,2))
            title(strcat(target(jj),'.',zone(zz)));
            
            subplot(3,4,8+jj)
            if zz ==1 || zz==3
                col = [85/243 85/243 85/243];
            else
                col = [8/243 133/243 161/243];
            end
            %plot(nanmean(PF.rateMap{ii,jj}{1,zz}),'color',col,'LineWidth',1.5)
            plot(nanmean(RM),'color',col,'LineWidth',1.5)
            hold on
            title(strcat('Firing rate ',target(jj)));
            %ylim([0 max(nanmean(PF.rateMap{ii,jj}{1,zz}))+2])
            ylim([0 max(nanmean(RM))+2])
            ylabel('Average pop. FR')
            xlabel('Position along track')           
        end
        
        subplot(3,4,10+jj)        
        [rho, pval] = corr(zscore(PF.rateMap{ii,jj}{1,1},0,2),zscore(PF.rateMap{ii,jj}{1,3},0,2),'Rows','pairwise','Type','Spearman');
        plot(diag(rho),'color',[85/243 85/243 85/243],'LineWidth',1.5)
        hold on
        [rho, pval] = corr(zscore(PF.rateMap{ii,jj}{1,1},0,2),zscore(PF.rateMap{ii,jj}{1,2},0,2),'Rows','pairwise','Type','Spearman');
        plot(diag(rho),'color',[56/243 61/243 150/243],'LineWidth',1.5)
        [rho, pval] = corr(zscore(PF.rateMap{ii,jj}{1,3},0,2),zscore(PF.rateMap{ii,jj}{1,4},0,2),'Rows','pairwise','Type','Spearman');
        plot(diag(rho),'color',[8/243 133/243 161/243],'LineWidth',1.5)
        legend({'BvsB','Left','Right'},'Location','southeast')
        legend('boxoff')
        ylim([0 1.2])
        xlim([1 50])        
        xlabel('Position along track')
        ylabel('Correlation')
        title(strcat('Map Corr ',target(jj)));
    end
    
    saveas(figure(ii),strcat(expPath,'\Summ\PlaceField',reg{ii},'.png'));
    saveas(figure(ii),strcat(expPath,'\Summ\PlaceField',reg{ii},'.eps'),'epsc');
    saveas(figure(ii),strcat(expPath,'\Summ\PlaceField',reg{ii},'.fig'));
end

