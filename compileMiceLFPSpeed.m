function compileMiceLFPSpeed

parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

tag = 'CA3';

freqs = logspace(log10(1),log10(150),200);
reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};
speedRange = 0:1:25;
thetaRange = [6.5 11];

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ21\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ11\Final','IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ28\Final',...
        'IZ29\Final','IZ30\Final','IZ31\Final'};  % To add: IZ23, IZ32, IZ33, IZ34
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


compiledLFPSpeed.mice = mice;

for rr = 1:length(reg)
    for cc = 1:length(target)
        for zz = 1:length(zone)
           compiledLFPSpeed.specmean{rr,cc}{zz} = [];
           compiledLFPSpeed.peakThetaFreq{rr,cc}{zz} = [];
           compiledLFPSpeed.avgThetaAmp{rr,cc}{zz} = [];
           
        end
    end
end

if strcmp(tag,'CA3')== 1 || strcmp(tag,'CA3Saline') == 1
    for m = 1:length(mice)
        cd(strcat(parentDir, mice{m},'\Summ'));
        if exist('LFPSpeed.mat','file')
            load(strcat('LFPSpeed.mat'));
        else 
            disp(['LFPSpeed not computed for mouse' mice{m}])
            continue;
        end

        for rr = 2:length(reg)
            for cc = 1%:length(target)
                for zz = 1:length(zone)      
                    for ss = 1:(length(speedRange)-1)
                        idxCurr = LFPSpeed.vel{rr,cc}{1,zz}>=speedRange(ss) & LFPSpeed.vel{rr,cc}{1,zz} < speedRange(ss+1);
                        spectemp(:,ss) = nanmean(10.^(LFPSpeed.specgram{rr,cc}{1,zz}(idxCurr,:)))';

                    end
                    compiledLFPSpeed.specmean{rr,cc}{zz}(:,:,m) = spectemp;%zscore(spectemp,0,'all');
                    freqIdx = find(freqs>thetaRange(1) & freqs <thetaRange(2));
                    freqRange = freqs(freqIdx);
                    [~, idxmax] = max(spectemp(freqIdx,:));
                    spectemp_theta = spectemp(freqIdx,:);
                    compiledLFPSpeed.peakThetaFreq{rr,cc}{zz}(:,:,m) = freqRange(idxmax);
                    compiledLFPSpeed.avgThetaAmp{rr,cc}{zz}(:,:,m) = mean(spectemp_theta);

%                     if zz < 4 
%                         loc = zz;
%                         colidx = m;%2*(m-1)+2;
%                         LStyle = '-';
%                     else
%                         loc = zz-3;
%                         colidx = m;%2*(m-1)+1;
%                         LStyle = '--';
%                     end
%                     subplot(2,6,6*(cc-1)+loc) 
%                     plot(speedRange(1:(end-1)),compiledLFPSpeed.peakThetaFreq{rr,cc}{zz}(:,:,m),'Color',cmap(colidx,:),'LineWidth',1.5,'LineStyle',LStyle);
%                     ylim([4.5 12.5])
%                     xlabel('Running Speed (cm/s)');
%                     ylabel('Peak theta frequency (Hz)');
%                     hold on
%                     subplot(2,6,6*(cc-1)+loc+3)  
%                     plot(speedRange(1:(end-1)),compiledLFPSpeed.avgThetaAmp{rr,cc}{zz}(:,:,m),'Color',cmap(colidx,:),'LineWidth',1.5,'LineStyle',LStyle);
%                     xlabel('Running Speed (cm/s)');
%                     ylabel('Average theta power');
%                     hold on
                end
            end

        end  
    end
    
    figure(1)
    set(gcf,'Position',[50 50 700 700])
    set(gcf,'renderer','painters')   

    figure(2)
    set(gcf,'Position',[50 50 700 700])
    set(gcf,'renderer','painters')
    
    for rr = 2:length(reg)
         for cc = 1%:length(target)
             for zz = 1:length(zone)    
                if zz < 4 
                    loc = zz;
                else
                    loc = zz-3;
                end
                if rr == 2 && zz<4
                    col = [85/243 85/243 85/243];
                elseif rr==2 && zz>=4
                    col = [8/243 133/243 161/243];
                elseif rr == 3 && zz <4
                    col = [224/243 163/243 46/243];
                else
                    col = [56/243 61/243 150/243];
                end
               
                figure(1)
                subplot(2,3,loc) 
                meanlfp = nanmean(compiledLFPSpeed.peakThetaFreq{rr,cc}{zz},3);
                stdlfp = nanstd(compiledLFPSpeed.peakThetaFreq{rr,cc}{zz},[],3)./sqrt(size(compiledLFPSpeed.peakThetaFreq{rr,cc}{zz},3));
                dev1 = meanlfp-stdlfp;
                dev2 = meanlfp+stdlfp;
                fill([speedRange(1:(end-1)) flip(speedRange(1:(end-1)))],[dev1 flip(dev2)],col,'FaceAlpha',.2,'EdgeColor','none')
                hold on
                plot(speedRange(1:(end-1)),meanlfp,'color',col,'LineWidth',1.5);             
                ylim([7 9.5])
                xlim([0 25])
                ylabel('Theta frequency')
                xlabel('Running speed')
                title(zone{loc}(1:(end-1)))


                subplot(2,3,loc+3) 
                meanlfp = nanmean(compiledLFPSpeed.avgThetaAmp{rr,cc}{zz},3);
                stdlfp = nanstd(compiledLFPSpeed.avgThetaAmp{rr,cc}{zz},[],3)./sqrt(size(compiledLFPSpeed.avgThetaAmp{rr,cc}{zz},3));
                dev1 = meanlfp-stdlfp;
                dev2 = meanlfp+stdlfp;
                fill([speedRange(1:(end-1)) flip(speedRange(1:(end-1)))],[dev1 flip(dev2)],col,'FaceAlpha',.2,'EdgeColor','none')
                hold on
                plot(speedRange(1:(end-1)),meanlfp,'color',col,'LineWidth',1.5);             
                ylabel('Theta amplitude')
                xlabel('Running speed')
                ylim([4*10^4 13*10^4])
                xlim([0 25])
                title(zone{loc}(1:(end-1)))
                
                figure(2)
                if zz<4
                    subplot(3,4,4*(loc-1)+2*(rr-2)+1)  
                else
                    subplot(3,4,4*(loc-1)+2*(rr-2)+2)  
                end
                contourf(speedRange(1:(end-1)),freqs,nanmean(compiledLFPSpeed.specmean{rr,cc}{zz},3),40,'LineColor','none');hold on;
                %imagesc(speedRange(1:(end-1)),freqs,nanmean(compiledLFPSpeed.specmean{rr,cc}{zz},3))
                set(gca,'YDir','normal')
                set(gca,'YScale','log')
                colormap jet
                colorbar
                caxis([1.5*10^4 1.6*10^5])
                %caxis([-1 3])
                ylim([5 15])
                xlabel('Running speed (cm/s)');
                ylabel('Frequency (Hz)');
                title(zone{loc}(1:(end-1)))
             end
         end
    end
     saveas(figure(1),strcat(parentDir,'\Compiled\LFPSpeed\compiledLFPSpeed',tag,'.png'));
     saveas(figure(1),strcat(parentDir,'\Compiled\LFPSpeed\compiledLFPSpeed',tag,'.fig'));
     saveas(figure(1),strcat(parentDir,'\Compiled\LFPSpeed\compiledLFPSpeed',tag,'.eps'),'epsc');

     saveas(figure(2),strcat(parentDir,'\Compiled\LFPSpeed\compiledLFPSpeedHeatmap',tag,'.png'));
     saveas(figure(2),strcat(parentDir,'\Compiled\LFPSpeed\compiledLFPSpeedHeatmap',tag,'.fig'));
     saveas(figure(2),strcat(parentDir,'\Compiled\LFPSpeed\compiledLFPSpeedHeatmap',tag,'.eps'),'epsc');
     
    %% Perform statistics 

    range = [5 15];
    figure
    col = [85/243 85/243 85/243;224/243 163/243 46/243;8/243 133/243 161/243;56/243 61/243 150/243];
    for zz = 1:3
        thetaFreqSumm = [];
        thetaAmpSumm = [];
        rangeId = [];
        manipId = [];        

        for kk = 1:size(range,1)             
            for rr = 2:3
                reshapedims = (range(kk,2)-range(kk,1)+1)*size(compiledLFPSpeed.peakThetaFreq{rr,1}{zz},3);
                thetaFreqSumm = [thetaFreqSumm; reshape(compiledLFPSpeed.peakThetaFreq{rr,1}{zz}(1,range(kk,1):range(kk,2),:),[reshapedims 1])];
                thetaAmpSumm = [thetaAmpSumm; reshape(compiledLFPSpeed.avgThetaAmp{rr,1}{zz}(1,range(kk,1):range(kk,2),:),[reshapedims 1])];
                rangeId = [rangeId; ones(reshapedims,1)*kk];
                if rr == 2
                    manipId = [manipId; ones(reshapedims,1)*1];
                else
                    manipId = [manipId; ones(reshapedims,1)*2];
                end
                thetaFreqSumm = [thetaFreqSumm; reshape(compiledLFPSpeed.peakThetaFreq{rr,1}{zz+3}(1,range(kk,1):range(kk,2),:),[reshapedims 1])];
                thetaAmpSumm = [thetaAmpSumm; reshape(compiledLFPSpeed.avgThetaAmp{rr,1}{zz+3}(1,range(kk,1):range(kk,2),:),[reshapedims 1])];
                rangeId = [rangeId; ones(reshapedims,1)*kk];
                if rr == 2
                    manipId = [manipId; ones(reshapedims,1)*3];
                else
                    manipId = [manipId; ones(reshapedims,1)*4];
                end
            end
        end
%         subplot(2,3,zz)
%         statsFreq{zz} = groupStats(thetaFreqSumm,[manipId rangeId],'inAxis',true);
%         subplot(2,3,zz+3)
%         statsAmp{zz} = groupStats(thetaAmpSumm,[manipId rangeId],'inAxis',true);
        
        subplot(2,3,zz)
        statsFreq{zz} = groupStats(thetaFreqSumm,manipId,'inAxis',true,'repeatedMeasures',true,'color',col);
        subplot(2,3,zz+3)
        statsAmp{zz} = groupStats(thetaAmpSumm,manipId,'inAxis',true,'repeatedMeasures',true,'color',col);
    end
    save(strcat(parentDir,'\Compiled\LFPSpeedStats',tag,'.mat'),'statsFreq', 'statsAmp');
end

end