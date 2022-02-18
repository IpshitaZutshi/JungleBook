function [BehavData BehavDataSumm] = SessBehavior(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'plotfig',true,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
plotfig = p.Results.plotfig;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');

BehavData{3,2} = [];
BehavDataSumm{3,2} = [];

for ii = 1:size(allSess,1)
    fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    file = dir(['*.SessionPulses.Events.mat']);
    load(file.name);
    efields = fieldnames(sessionPulses);
    for jj = 1:length(efields)
        BehavData{sessionPulses.(efields{jj}).region, sessionPulses.(efields{jj}).target} = ...
            catpad(1,BehavData{sessionPulses.(efields{jj}).region, sessionPulses.(efields{jj}).target}, sessionPulses.(efields{jj}).performanceBlocks(2,:));
        BehavDataSumm{sessionPulses.(efields{jj}).region, sessionPulses.(efields{jj}).target} = ...
            catpad(1,BehavDataSumm{sessionPulses.(efields{jj}).region, sessionPulses.(efields{jj}).target}, sessionPulses.(efields{jj}).performance(2,:));

    end
    
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
if plotfig
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

end