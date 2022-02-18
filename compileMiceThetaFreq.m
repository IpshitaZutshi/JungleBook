function compileMiceThetaFreq(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'passband',[6 10],@isnumeric);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
saveMat = p.Results.saveMat;
force = p.Results.force;
passband = p.Results.passband;

tag = 'mECBilateral'; % or mEC

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ21\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ11\Final','IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final'...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline',...
        'IZ29\Final','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};  
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ26\Final'}; %'IZ24\Final','IZ25\Final',    
    reg = {'contramEC','ipsimEC','Both'};
end

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};
freqs = logspace(log10(1),log10(300),200);

for rr = 1:length(reg)
    for cc = 1:length(target)
        peakFreq{rr,cc} = [];        
        peakAmp{rr,cc} = [];
    end
end

for m = 1:length(mice)

    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('PowerSpectrumSlopePyrChPlus0.mat','file')
        load('PowerSpectrumSlopePyrChPlus0.mat');
    else 
        disp(['Power spectrum not computed for mouse' mice{m}])
        continue;
    end
    
    if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1 || strcmp(tag,'mEC')
        iiRange = [2 3];
    else
        iiRange = [1 2 3];
    end
    for ii = iiRange
        for jj = 1%:length(target)
            freqMat = [];
            powerMat = [];
            for kk = 1:length(zone)               
                freqArr = [];
                powerArr = [];
                if strcmp(mice{m},'IZ29\Final') == 1
                    ppRange = 2:(size(PSS.specGram{ii,jj}{kk},3)-1);
                else
                    ppRange = 2:size(PSS.specGram{ii,jj}{kk},3);
                end
                for pp =  ppRange % First array is nans     
                    specGram = PSS.specGram{ii,jj}{kk}(:,:,pp);%.^10;
                    meanPower = nanmean(specGram,1);
                    freqIdx = find(freqs>passband(1) & freqs <passband(2));
                    freqRange = freqs(freqIdx);
                    [maxTheta, idxmax] = max(meanPower(freqIdx));
                    freqArr = [freqArr; freqRange(idxmax)];  
                    powerArr = [powerArr;maxTheta];
                end
                freqMat(:,kk) = freqArr;
                powerMat(:,kk) = powerArr;
            end
            peakFreq{ii,jj} = [peakFreq{ii,jj};freqMat];
            peakAmp{ii,jj} = [peakAmp{ii,jj};powerMat];
%             peakFreq{ii,jj} = [peakFreq{ii,jj};nanmean(freqMat,1)];
%             peakAmp{ii,jj} = [peakAmp{ii,jj};nanmean(powerMat,1)];            
        end
    end
end

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
        colMat = [85/243 85/243 85/243;...
        8/243 133/243 161/243;...
        224/243 163/243 46/243;...         
        56/243 61/243 150/243];   
    figure
    set(gcf,'Position',[100 100 1300 700])    
    set(gcf,'renderer','painters');  
    
    for jj = 1%:length(target)
        freqMatSide = [peakFreq{2,jj}(:,[1 4]) peakFreq{3,jj}(:,[1 4])];   
        freqMatSideID = reshape(freqMatSide,[size(freqMatSide,1)*size(freqMatSide,2),1]);
        stimIDSide = [ones(size(freqMatSide,1),1);ones(size(freqMatSide,1),1)*2;ones(size(freqMatSide,1),1)*3;ones(size(freqMatSide,1),1)*4];            
        
        freqMatCentral = [peakFreq{2,jj}(:,[2 5]) peakFreq{3,jj}(:,[2 5])];   
        freqMatCentralID = reshape(freqMatCentral,[size(freqMatCentral,1)*size(freqMatCentral,2),1]);
        stimIDCentral = [ones(size(freqMatCentral,1),1);ones(size(freqMatCentral,1),1)*2;ones(size(freqMatCentral,1),1)*3;ones(size(freqMatCentral,1),1)*4];         
        
        AmpMatSide = [peakAmp{2,jj}(:,[1 4]) peakAmp{3,jj}(:,[1 4])];   
        AmpMatSideID = reshape(AmpMatSide,[size(AmpMatSide,1)*size(AmpMatSide,2),1]);

        AmpMatCentral = [peakAmp{2,jj}(:,[2 5]) peakAmp{3,jj}(:,[2 5])];   
        AmpMatCentralID = reshape(AmpMatCentral,[size(AmpMatCentral,1)*size(AmpMatCentral,2),1]);
    end
        
    subplot(2,2,1)
    stats{1} = groupStats(freqMatSideID,stimIDSide,'inAxis',true,'RepeatedMeasures',true,'color',colMat);
    hold on
%     plot(freqMatSide','Color',[0.5 0.5 0.5]);
    ylabel('Peak theta frequency')
    title('Side arm');
    
    subplot(2,2,2)
    stats{2} = groupStats(freqMatCentralID,stimIDCentral,'inAxis',true,'RepeatedMeasures',true,'color',colMat);
    hold on
%     plot(freqMatCentral','Color',[0.5 0.5 0.5]);
    ylabel('Peak theta frequency')
    title('Central arm');
    
    subplot(2,2,3)
    stats{3} = groupStats(AmpMatSideID,stimIDSide,'inAxis',true,'RepeatedMeasures',true,'color',colMat);
    hold on
%     plot(freqMatSide','Color',[0.5 0.5 0.5]);
    ylabel('Peak theta power')
    title('Side arm');
    
    subplot(2,2,4)
    stats{4} = groupStats(AmpMatCentralID,stimIDCentral,'inAxis',true,'RepeatedMeasures',true,'color',colMat);
    hold on
%     plot(freqMatCentral','Color',[0.5 0.5 0.5]);
    ylabel('Peak theta power')
    title('Central arm');
    
    saveas(gcf,strcat(parentDir,'\Compiled\LFPFreq\',tag,'Freq.png'),'png');
    saveas(gcf,strcat(parentDir,'\Compiled\LFPFreq\',tag,'Freq.fig'),'fig');
    saveas(gcf,strcat(parentDir,'\Compiled\LFPFreq\',tag,'Freq.eps'),'epsc');
    save(strcat(parentDir,'\Compiled\LFPFreq\',tag,'Freq.mat'),'stats');
    
else
    colMat = [85/243 85/243 85/243;...
        224/243 163/243 46/243;... 
        8/243 133/243 161/243;...
        56/243 61/243 150/243];   
    figure
    set(gcf,'Position',[100 100 350 900])    
    set(gcf,'renderer','painters');  
    if strcmp(tag,'mEC')
        iiRange = 2;
    else
        iiRange = [1 2 3];
    end
    stimID =[];
    AmpMatID  = [];
    freqMatID = [];
    
    for ii = iiRange
        for jj = 1%:length(target)
%             freqMat = peakFreq{ii,1}(:,[1 4 2 5]); % Side base, Side stim, Center base, center stim
%             freqMatID = reshape(freqMat,[size(freqMat,1)*size(freqMat,2),1]);
%             AmpMat = peakAmp{ii,1}(:,[1 4 2 5]); % Side base, Side stim, Center base, center stim
%             AmpMatID = reshape(AmpMat,[size(AmpMat,1)*size(AmpMat,2),1]);            
%             zoneID = [ones(size(freqMat,1),1);ones(size(freqMat,1),1);ones(size(freqMat,1),1)*2;ones(size(freqMat,1),1)*2];
%             stimID = [ones(size(freqMat,1),1);ones(size(freqMat,1),1)*2;ones(size(freqMat,1),1);ones(size(freqMat,1),1)*2];            
%             subplot(2,3,ii)
%             stats{ii,1} = groupStats(freqMatID,[zoneID stimID],'inAxis',true,'color',[colMat(1,:);colMat(ii+1,:);colMat(1,:);colMat(ii+1,:)]);
%             hold on
%             plot([ones(1,size(freqMat,1)); ones(1,size(freqMat,1))*2], freqMat(:,1:2)','Color',[0.5 0.5 0.5]);
%             plot([ones(1,size(freqMat,1))*3.5; ones(1,size(freqMat,1))*4.5], freqMat(:,3:4)','Color',[0.5 0.5 0.5]);  
%             ylabel('Peak theta frequency')
%             ylim([6 10])
%             title(strcat(reg{ii},'  Side___Central'));
%             
%             subplot(2,3,ii+3)
%             stats{ii,2} = groupStats(AmpMatID,[zoneID stimID],'inAxis',true,'color',[colMat(1,:);colMat(ii+1,:);colMat(1,:);colMat(ii+1,:)]);
%             hold on
%             plot([ones(1,size(AmpMat,1)); ones(1,size(AmpMat,1))*2], AmpMat(:,1:2)','Color',[0.5 0.5 0.5]);
%             plot([ones(1,size(AmpMat,1))*3.5; ones(1,size(AmpMat,1))*4.5], AmpMat(:,3:4)','Color',[0.5 0.5 0.5]);  
%             ylabel('Peak theta amplitude')
%             ylim([4.4 5.5])
%             title(strcat(reg{ii},'  Side___Central')); 

            freqMat = peakFreq{ii,1}(:,[2 5]); % Side base, Side stim, Center base, center stim
            freqMatID = [freqMatID; reshape(freqMat,[size(freqMat,1)*size(freqMat,2),1])];
            AmpMat = peakAmp{ii,1}(:,[2 5]); % Side base, Side stim, Center base, center stim
            AmpMatID = [AmpMatID; reshape(AmpMat,[size(AmpMat,1)*size(AmpMat,2),1])];            
            %zoneID = [ones(size(freqMat,1),1);ones(size(freqMat,1),1);ones(size(freqMat,1),1)*2;ones(size(freqMat,1),1)*2];
            stimID = [stimID; ones(size(freqMat,1),1)*1;ones(size(freqMat,1),1)*(ii+1)];%;ones(size(freqMat,1),1);ones(size(freqMat,1),1)*2];   
        end
    end
    subplot(2,1,1)
    stats.peakFreq = groupStats(freqMatID, stimID,'inAxis',true,'color',colMat);
    hold on
    ylabel('Peak theta frequency')
    ylim([6 10])
    title(' Central arm');
    
    subplot(2,1,2)
    stats.peakAmp = groupStats(AmpMatID, stimID,'inAxis',true,'color',colMat);
    hold on
    ylabel('Peak theta amplitude')
    ylim([4.4 5.5])
    title(' Central arm');    
    
    saveas(gcf,strcat(parentDir,'\Compiled\LFPFreq\',tag,'Freq.png'),'png');
    saveas(gcf,strcat(parentDir,'\Compiled\LFPFreq\',tag,'Freq.fig'),'fig');
    saveas(gcf,strcat(parentDir,'\Compiled\LFPFreq\',tag,'Freq.eps'),'epsc');
    save(strcat(parentDir,'\Compiled\LFPFreq\',tag,'Freq.mat'),'stats');
end
end