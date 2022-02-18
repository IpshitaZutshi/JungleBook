function plotBilateralEffects

pathDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\IZ26\Final\IZ26_576um_201020_sess4';
analogEv = 64;
numAnalog = 2;
pre = 2;
post = 6;
for ii = 1:numAnalog
    analogCh(ii) = (analogEv-1)+ii;
end

cd(pathDir)
[pulses] = bz_getAnalogPulses('analogCh',analogCh);
[sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
chanIpsi = 14;
chanContra = 126;

lfp = bz_GetLFP([chanIpsi chanContra],'noPrompts', true);
data = lfp.data;
timestamps = lfp.timestamps;

colMat = [85/243 85/243 85/243;
         103/243 189/243 170/243];
 
figure     
for i = 1:3%:(numAnalog+1)   
%     figure(i)
%     set(gcf,'renderer','Painters')    
    %Extract relevant pulses
    if i<=numAnalog
        pulTr = (pulses.stimComb==i);
    else
        pulTr = (pulses.stimPerID'==1 & pulses.stimComb==i);
    end
%     pulseDur = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
%     Pulseidx = pulseDur < 4.95 | pulseDur > 5.05;
%     pulTr = pulTr & Pulseidx;    
%     events = pulses.intsPeriods(:,pulTr);
%     events = round(events*1250);
%     events = events(:,(events(1,:) + (15*1250) <= size(data,1)) & (events(1,:) - (5*1250) > 0));
%     k = 1;
%     for pp = [8 9 13]%1:min(length(events(1,:)),30)    
%         subplot(3,1,k)
%         hold on
%         for kk = 1:2
%             plot(timestamps((events(1,pp)-(1250*pre)):(events(1,pp)+(1250*post))),0.75*(data((events(1,pp)-(1250*pre)):(events(1,pp)+(1250*post)),kk))-(kk-1)*3000,'Color',colMat(kk,:),'LineWidth',1)
%         end
%         xlim([timestamps((events(1,pp)-(1250*pre))) timestamps((events(1,pp)+(1250*post)))])
%         line([timestamps(events(1,pp)) timestamps(events(1,pp))],[2500 (-6000)],'Color','red','LineWidth',2)
%         line([timestamps(events(2,pp)) timestamps(events(2,pp))],[2500 (-6000)],'Color','red','LineWidth',2)
%         k= k+1;
%     end
    pulseDur = pulses.intsPeriods(2,:) - pulses.intsPeriods(1,:);
    Pulseidx = pulseDur < 4.95 | pulseDur > 5.05;
    pulTr = pulTr & Pulseidx;       
    events = pulses.intsPeriods(:,pulTr);
%     events = round(events*1250);
%     events = events(:,(events(1,:) + (15*1250) <= size(data,1)) & (events(1,:) - (5*1250) > 0));  
    
    specslope = bz_PowerSpectrumSlope(lfp,1,0.5,'frange',[1 100],'spectype','fft','ints',events');
    subplot(1,3,i)
    plot(specslope.freqs,mean(10.^specslope.specgram(:,:,1),1),'Color',colMat(1,:),'LineWidth',1.5)
    hold on
    plot(specslope.freqs,mean(10.^specslope.specgram(:,:,2),1),'Color',colMat(2,:),'LineWidth',1.5)
end

freqs = logspace(log10(1),log10(300),200);
freqsIdx = freqs>6&freqs<12;
freqsInt = freqs(freqsIdx);
figure
for ii = 1:3
     for jj = 1
         for kk = [2 5]                
             spec.ipsi = [];
             spec.contra = [];
             for pp = 2:size(PSS.specData{ii,jj}{kk},3) % First array is nans
                spec.ipsi = [spec.ipsi (10.^PSS.specGram{ii,jj}{kk}(:,:,pp))']; 
                spec.contra = [spec.contra (10.^contraPSS.specGram{ii,jj}{kk}(:,:,pp))']; 
             end   
             if kk == 2
                 col = [56/243 61/243 150/243];
             else
                 col = [8/243 133/243 161/243];
             end
             subplot(1,3,ii)
             meanPower = nanmean(spec.ipsi,2);
             stdPower = nanstd(spec.ipsi,[],2)./sqrt(size(spec.ipsi,2));
             dev1 = meanPower - stdPower;
             dev2 = meanPower + stdPower;
             hold on
             fill([freqs flip(freqs)],[dev1' flip(dev2)'],col,'FaceAlpha',.2,'EdgeColor','none');
             plot(freqs,meanPower,'color',col,'LineWidth',1.5);
             [~,maxFreqIdx] = max(meanPower(freqsIdx));
             maxFreq = freqsInt(maxFreqIdx);
             if kk == 5
                line([maxFreq maxFreq],[0 1.5*10^5],'color',col,'LineWidth',1.5) 
             end
             
             hold on
             if kk == 2
                 col = [214/243 126/243 44/243];
             else
                 col = [231/243 199/243 31/243];
             end             
             meanPower = nanmean(spec.contra,2);
             stdPower = nanstd(spec.contra,[],2)./sqrt(size(spec.contra,2));
             dev1 = meanPower - stdPower;
             dev2 = meanPower + stdPower;
             hold on
             fill([freqs flip(freqs)],[dev1' flip(dev2)'],col,'FaceAlpha',.2,'EdgeColor','none');
             plot(freqs,meanPower,'color',col,'LineWidth',1.5);
             ylim([0 1.5*10^5])
             xlim([0 25])
             [~,maxFreqIdx] = max(meanPower(freqsIdx));
             maxFreq = freqsInt(maxFreqIdx);
             if kk == 5
                line([maxFreq maxFreq],[0 1.5*10^5],'color',col,'LineWidth',1.5)                           
             end
         end
     end
end