function [sigMat, mutInfo] = compileTuningProportion;

sessName = {'IZ44_sess7_freqlin','IZ39_sess10_freqlin','IZ43_sess4_freqlin','IZ40_sess18_freqlin'...
    'IZ41_sess7_freqlin','IZ41_sess9_freqlin','IZ41_sess11_freqlin','IZ41_sess14_freqlin'};

% sessName = {'IZ41_sess7_freqlin','IZ41_sess14_freqlin',...
%     'IZ43_sess4_freqlin','IZ39_sess10_freqlin'};
sessType = [1 1 1 1 0 0 0 0];
plotcells = 0;
basepath = 'C:\Data\PGAMAnalysis\';

for ss = 1:length(sessName)
    cd(strcat(basepath,sessName{ss}))
    
    load('results_struct.mat');

    % Get neuronID 
    for aa = 1:size(results,2)
        neuronID(aa) = results(aa).neuron;
    end

    % For each neuron, plot the tuning

    mkdir('GAM tuning')
    cd('GAM tuning')

    numCells = max(neuronID);

    sigMat = zeros(numCells,15);
    mutInfo = zeros(numCells,15);
    tuningFreq = [];
    tuningSpace = [];

    for aa = 0:numCells
        idx = find(neuronID == aa);
        if plotcells 
            figure
            set(gcf,'Position',[440 192 10001 605])
        end
        for id = 1:length(idx)       
            dev1 = results(idx(id)).kernel_mCI;
            dev2 = results(idx(id)).kernel_pCI;

            hold on
            if (results(idx(id)).pval < 0.001 && results(idx(id)).pval >0)% && ~isnan(results(idx(id)).mutual_info))
                sigMat(aa+1,id) = 1;
                mutInfo(aa+1,id) = results(idx(id)).mutual_info;
                col = 'b';

%                 if id == 2 % if Yfwd
%                     tuningSpace = [tuningSpace; results(idx(id)).raw_rate_Hz];
%                 elseif id == 6 % if Freq
%                     tuningFreq = [tuningFreq; results(idx(id)).raw_rate_Hz];
%                 end
            else
                col = 'k';
                sigMat(aa+1,id) = 0;
                mutInfo(aa+1,id) = nan;
            end
            if plotcells 
                subplot(6,6,id)
                fill([results(idx(id)).kernel_x flip(results(idx(id)).kernel_x)], [dev1 flip(dev2)],col,'FaceAlpha',.2,'EdgeColor','none')
                hold on
                plot(results(idx(id)).kernel_x, results(idx(id)).kernel,col,'LineWidth',2)
                title([results(idx(id)).variable])

                subplot(6,6,id+18)
                plot(results(idx(id)).x, results(idx(id)).model_rate_Hz,'k','LineWidth',1.5)
                hold on
                plot(results(idx(id)).x, results(idx(id)).raw_rate_Hz,'r','LineWidth',1.5)
                title([results(idx(id)).variable])
                ylabel('Spike/s')
                xlabel('Tuning function axis')
                if id == length(idx)
                    legend({'Model rate','raw rate'},'Location','southwest','box','off')
                end
            end
        end
        if plotcells
            saveas(gcf,['cell_' num2str(aa+1) '.png'],'png');
            saveas(gcf,['cell_' num2str(aa+1) '.fig'],'fig');
            close all;
        end
    end

%     figure
%     set(gcf,'Renderer','painters')
%     subplot(2,3,[1 2 3])
%     propSig = sum(sigMat)./size(sigMat,1);
%     plot(propSig,'+');%bar(propSig)
%     ylim([0 1])
%     xticklabels({'x','yFwd','yRev','yLin','vel','freq','trialEnd','trialRamp','correct','incorrect','lick0','lick1','lick2','lick3','lick4','lick5','lick6','spikehist'});
%     title('Distribution of significant tuning')
% 
%     %For cells tuned to space and frequency, what are they better tuned to
%     sigY = sigMat(:,[2,6]);
%     mutY = mutInfo(:,[2,6]);
%     idxBoth = find(sum(sigY,2)==2);
% 
%     for ii = 1:length(idxBoth)
%         if mutY(idxBoth(ii),1)>mutY(idxBoth(ii),2)
%             space(ii) = 1;
%         else 
%             space(ii) = 0;
%         end
% %     end
% 
%     propDist = [sum(space==1) sum(space==0)]./length(space);
%     subplot(2,3,4)
%     bar(propDist)
%     ylim([0 1])
%     xticklabels({'Space','Tone'})
%     title('Preference for space vs tone')
% 
%     subplot(2,3,5)
%     [~,idx] = max(tuningSpace,[],2);
%     [~,sortidx] = sort(idx);
%     hm = zscore(tuningSpace,[],2);
%     imagesc(results(2).x,1:1:length(sortidx),hm(sortidx,:))
%     title('Tuning to Space')
% 
%     subplot(2,3,6)
%     [~,idx] = max(tuningFreq,[],2);
%     [~,sortidx] = sort(idx);
%     hm = zscore(tuningFreq,[],2);
%     imagesc(results(6).x,1:1:length(sortidx),hm(sortidx,:))
%     title('Tuning to Frequency')
%     
    clear neuronID
    propSigAll(ss,:) = sum(sigMat)./size(sigMat,1);
end

figure
plot(1:2:36,propSigAll(sessType==0,:)','.','MarkerSize',20,'Color',[0.5 0.5,0.5,0.2])
hold on
plot(1:2:36,mean(propSigAll(sessType==0,:)),'+','MarkerSize',20)

plot((1:2:36)+0.7,propSigAll(sessType==1,:)','.','MarkerSize',20,'Color',[216/243 0/143 115/243 0.2])
hold on
plot((1:2:36)+0.7,mean(propSigAll(sessType==1,:)),'+','MarkerSize',20)

for ii = 1:size(propSigAll,2)
    [~,sig(ii)] = ttest(propSigAll(sessType==0,ii),propSigAll(sessType==1,ii));
end

end