function getPeriStimModIdx

[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
load([sessionInfo.FileName '.region.mat']);
pyrCh = region.CA1sp;
[pulses] = bz_getAnalogPulsesSine('analogCh',[64 65]);
phaserange = 2:0.2:20;
amprange = 30:1:300;
numAnalog = 2;

reg = {'CA1/CA3','mEC','Both'};
figure

for i = 1:(numAnalog+1)
    
    if i<=numAnalog
        pulTr = (pulses.stimComb==i);
    else
        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
    end
    events = pulses.intsPeriods(:,pulTr);
    events = events(:,(events(1,:) - 5) > 0);
    lfp = bz_GetLFP(pyrCh,'noPrompts', true);
    %lfp = bz_GetLFP(47,'noPrompts', true);
    
    [comodPost] = bz_ModIndex_IZ(lfp,'intervals',events','flagPlot', false,'phaserange',phaserange,'amprange',amprange);
    [comodPre] = bz_ModIndex_IZ(lfp,'intervals',(events-5)','flagPlot', false,'phaserange',phaserange,'amprange',amprange);
    
    subplot(3,2,2*(i-1)+1)
    imagesc(phaserange(2:end),amprange(2:end),comodPre)
    set(gca,'YDir','normal')
    colormap jet
    colorbar
    xlabel('Frequency phase');
    ylabel('Frequency amplitude');
    title(strcat('Baseline  ',reg(i)));
    
    subplot(3,2,2*(i-1)+2)
    imagesc(phaserange(2:end),amprange(2:end),comodPost)
    set(gca,'YDir','normal')
    colormap jet
    colorbar
    xlabel('Frequency phase');
    ylabel('Frequency amplitude');
    title(strcat('Stimulation  ',reg(i)));
    
end

saveas(gcf,'Summ\Theta_gamma coupling.png');
saveas(gcf,'Summ\Theta_gamma coupling.eps');
saveas(gcf,'Summ\Theta_gamma coupling.fig');