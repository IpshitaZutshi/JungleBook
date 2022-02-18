function SessProfileSlopeCorrected(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
parse(p,varargin{:});
frange = [2 300];

expPath = p.Results.expPath;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

analogCh = [64 65];
k = 1;

for ii = 6:14
    fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    file = dir(['*.region.mat']);
    load(file.name);
    lfpCh = region.CA1sp;
    
    [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);

    [pulses] = bz_getAnalogPulses('analogCh',analogCh);
    lfp = bz_GetLFP(lfpCh,'noPrompts', true);
    
     for i =1:3
        if i<=numAnalog
            pulTr = (pulses.stimComb==i);
        else
            pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
        end
        events = pulses.intsPeriods(:,pulTr)';
        events(1,:) =  events(1,:)+1;
        events(2,:) =  events(2,:)-0.5;
        events_pre = events-6;
        
        [specslope_stim,spec_stim] = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',frange,'spectype','fft','ints',events);
        [specslope,spec] = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',frange,'spectype','fft','ints',events_pre);

        subplot(10,3,3*(k-1)+i)
        plot(specslope.freqs,mean(specslope.resid,1),'k')
        hold on
        plot(specslope_stim.freqs,mean(specslope_stim.resid,1),'b')
        ylabel('Corrected power')
        xlabel('Frequency')
        if i ==2
           title(strcat('Channel ',num2str(lfpCh(l))))
        end
        
        SlopeData
     end
    k = k+1;
    
end

reg = {'CA1','mEC','Both'};
zone = {'Stem','Return'};

% for ii = 1:size(BehavData,1)
%     for jj = 1:size(BehavData,2)
%         subplot(3,2,2*(ii-1)+jj)
%         err = std(BehavData{ii,jj},1)/sqrt(size(BehavData,1));
%         errorbar(nanmean(BehavData{ii,jj}),err,'LineWidth',2);
%         xlabel('Trial Blocks')
%         ylabel('Performance')
%         hold on
%         scatter([1:1:4],BehavData{ii,jj}(1,:),'ro','filled');
%         scatter([1:1:4],BehavData{ii,jj}(2,:),'ko','filled');
%         scatter([1:1:4],BehavData{ii,jj}(3,:),'mo','filled');
%         ylim([0.5 1.1])
%         xlim([0 7])
%         title(strcat('Region: ',reg{ii},' Zone:',zone{jj}));
%     end
% end

figure
for ii = 1:size(BehavDataSumm,1)
    for jj = 1:size(BehavDataSumm,2)
        subplot(3,2,2*(ii-1)+jj)
        err = std(BehavDataSumm{ii,jj},1)/sqrt(size(BehavDataSumm,1));
        errorbar(nanmean(BehavDataSumm{ii,jj}),err,'LineWidth',2);
        xlabel('Trial Blocks')
        ylabel('Performance')
        hold on
%         scatter([1:1:2],BehavDataSumm{ii,jj}(1,:),'ro','filled');
%         scatter([1:1:2],BehavDataSumm{ii,jj}(2,:),'ko','filled');
%         scatter([1:1:2],BehavDataSumm{ii,jj}(3,:),'mo','filled');
        ylim([0.5 1.1])
        xlim([0 3])
        title(strcat('Region: ',reg{ii},'      Zone:',zone{jj}));
    end
end

end