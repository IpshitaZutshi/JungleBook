function compileMiceThetaCompressionVer2

parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
tag = 'mECBilateral';

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ17\Final','IZ18\Final','IZ15\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline',...
        'IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'}; 
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ33\Final','IZ32\Final','IZ34\Final','IZ27\Final','IZ28\Final','IZ29\Final'};
    reg = {'mEC','CA3','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Both'};
end

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for rr = 1:length(reg)
    for cc = 1:length(target)
        for zz = 1:2
            placefield_difference{rr,cc}{zz} = [];
            placefield_center{rr,cc}{zz} = [];
            ccgs_place_offset{rr,cc}{zz} = [];
            ccgs_time_offset{rr,cc}{zz} = [] ;
            ccgs_phase_offset{rr,cc}{zz} = [];      
            ccgs_time_peaks{rr,cc}{zz} = [] ;
            ccgs_phase_peaks{rr,cc}{zz} = [];
            compression{rr,cc}{zz} = [];
            compressionMouse{rr,cc}{zz} = [];
        end
        region{rr,cc} = [];
        putativeCellType{rr,cc} = [];                  
        sessCount{rr,cc} = [];
        mouseCount{rr,cc} = [];
    end
end

for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('ThetaCompression.mat','file')
        load('ThetaCompression.mat');
    else 
        disp(['ThetaCompression not computed for mouse' mice{m}])
        continue;
    end
    
    for rr = 1:length(reg)
        for cc = 1:length(target)
            for zz = 1:2      
                placefield_difference{rr,cc}{zz} = [placefield_difference{rr,cc}{zz}; thetaComp.placefield_difference{rr,cc}{zz}];
               % placefield_size{rr,cc}{zz} = [placefield_size{rr,cc}{zz}; thetaComp.placefield_size{rr,cc}{zz}];
                placefield_center{rr,cc}{zz} = [placefield_center{rr,cc}{zz}; thetaComp.placefield_center{rr,cc}{zz}];
                ccgs_place_offset{rr,cc}{zz} = [ccgs_place_offset{rr,cc}{zz}; thetaComp.ccgs_place_offset{rr,cc}{zz}];  
                ccgs_time_offset{rr,cc}{zz} = [ccgs_time_offset{rr,cc}{zz}; thetaComp.ccgs_time_offset{rr,cc}{zz}];
                ccgs_phase_offset{rr,cc}{zz} = [ccgs_phase_offset{rr,cc}{zz}; thetaComp.ccgs_phase_offset{rr,cc}{zz}];
                ccgs_time_peaks{rr,cc}{zz} = [ccgs_time_peaks{rr,cc}{zz}; thetaComp.ccgs_time_peaks{rr,cc}{zz}];
                ccgs_phase_peaks{rr,cc}{zz} = [ccgs_phase_peaks{rr,cc}{zz}; thetaComp.ccgs_phase_peaks{rr,cc}{zz}];                
                compression{rr,cc}{zz} = [compression{rr,cc}{zz}; thetaComp.compression{rr,cc}{zz}];
                % Calculate compression per mouse
%                 [fit_phase, fit_time] = getCompression(thetaComp.placefield_center{rr,cc}{zz},thetaComp.placefield_difference{rr,cc}{zz},...
%                     thetaComp.ccgs_phase_peaks{rr,cc}{zz},thetaComp.ccgs_time_peaks{rr,cc}{zz},thetaComp.region{rr,cc},thetaComp.putativeCellType{rr,cc});
% %                 [fit_phase, fit_time] = getCorrCoeff(thetaComp.placefield_center{rr,cc}{zz},thetaComp.ccgs_place_offset{rr,cc}{zz},...
% %                     thetaComp.ccgs_phase_peaks{rr,cc}{zz},thetaComp.ccgs_time_peaks{rr,cc}{zz},thetaComp.region{rr,cc},thetaComp.putativeCellType{rr,cc});
%                 compressionMouse{rr,cc}{zz} = [compressionMouse{rr,cc}{zz};[fit_time fit_phase]];
            end
            region{rr,cc} = [region{rr,cc};thetaComp.region{rr,cc}];
            putativeCellType{rr,cc} = [putativeCellType{rr,cc};thetaComp.putativeCellType{rr,cc}];                  
            sessCount{rr,cc} = [sessCount{rr,cc};thetaComp.sessCount{rr,cc}];
            mouseCount{rr,cc} = [mouseCount{rr,cc}; ones(length(thetaComp.region{rr,cc}),1)*m];
        end
    end
end

close all
pause(2)

if strcmp(tag,'CA3')==1 || strcmp(tag,'CA3Saline')==1

    colMat = [85/243 85/243 85/243;...
            8/243 133/243 161/243;...    
            224/243 163/243 46/243;...
            56/243 61/243 150/243];      
    figure(1)
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[20 20 1871 920])    
    timeSlope = [];  timePeriod = [];
    phaseSlope = []; phasePeriod = [];
    manipID = [];
    for ii = 2:3       
        for jj = 1
            for kk = 1:2

                idxtoKeep =  placefield_center{ii,jj}{kk}==1 & (ccgs_place_offset{ii,jj}{kk}>-600 & ccgs_place_offset{ii,jj}{kk}<600) &...                    
                    region{ii,jj} == 1 & putativeCellType{ii,jj} == 1; 
                [fit_phase, fit_time] = getCompression(placefield_center{ii,jj}{kk},placefield_difference{ii,jj}{kk},...
                    ccgs_phase_peaks{ii,jj}{kk},ccgs_time_peaks{ii,jj}{kk},region{ii,jj},putativeCellType{ii,jj});
                
                figure(1)
                subplot(2,7,2*(ii-2)+kk)   
                dataPlace = repmat(ccgs_place_offset{ii,jj}{kk}(idxtoKeep),[1 size(ccgs_time_peaks{ii,jj}{kk},2)]);
                dataTime = ccgs_time_peaks{ii,jj}{kk}(idxtoKeep,:);
                dataPlace = reshape(dataPlace,[size(dataPlace,1)*size(dataPlace,2),1]);
                dataTime = reshape(dataTime,[size(dataTime,1)*size(dataTime,2),1]);                   
                s = scatter(dataTime,dataPlace,2,'filled');
                s.MarkerFaceAlpha = 0.4;
                xlabel('Theta timescale (ms)')
                ylabel('Place field dist (ms)')
                x_bins1 = -400:8:400;
                hold on
                for i = -5:5
                    plot(x_bins1*fit_time(1)+i*fit_time(2),x_bins1,'--r')
                end
                ylim([-400,400]),xlim([-400,400])    
                title(strcat('Slope:',num2str(fit_time(1)),' Period:',num2str(fit_time(2))))

                figure(1)
                subplot(2,7,7+2*(ii-2)+kk)
                dataPlace = repmat(ccgs_place_offset{ii,jj}{kk}(idxtoKeep),[1 size(ccgs_phase_peaks{ii,jj}{kk},2)]);
                dataPhase = ccgs_phase_peaks{ii,jj}{kk}(idxtoKeep,:);
                dataPlace = reshape(dataPlace,[size(dataPlace,1)*size(dataPlace,2),1]);
                dataPhase = reshape(dataPhase,[size(dataPhase,1)*size(dataPhase,2),1]);
                dataPhase(dataPhase==0) = nan;
                s = scatter(dataPhase,dataPlace,2,'filled');
                s.MarkerFaceAlpha = 0.4;
                xlabel('Theta phase (rad)')
                ylabel('Place field dist (ms)')
                x_bins1 = -400:8:400;
                hold on
                for i = -5:5
                    plot(x_bins1*fit_phase(1)+i*fit_phase(2),x_bins1,'--r')
                end
                ylim([-400,400]),xlim([-16,16])    
                title(strcat('Slope:',num2str(fit_phase(1)),' Period:',num2str(fit_phase(2))))                

                timeSlope = [timeSlope;compressionMouse{ii,jj}{kk}(:,1)];  
                timePeriod = [timePeriod;compressionMouse{ii,jj}{kk}(:,2)];
                phaseSlope = [phaseSlope;compressionMouse{ii,jj}{kk}(:,8)];  
                phasePeriod = [phasePeriod;compressionMouse{ii,jj}{kk}(:,9)];
                if kk == 1 && ii == 2
                    manipID = [manipID;ones(size(compressionMouse{ii,jj}{kk},1),1)*1];
                elseif kk == 2 && ii == 2
                    manipID = [manipID;ones(size(compressionMouse{ii,jj}{kk},1),1)*2];
                elseif kk == 1 && ii == 3
                    manipID = [manipID;ones(size(compressionMouse{ii,jj}{kk},1),1)*3];
                elseif kk == 2 && ii == 3
                    manipID = [manipID;ones(size(compressionMouse{ii,jj}{kk},1),1)*4];
                end
%                 timeSlope = [timeSlope;compression{ii,jj}{kk}(:,1)];  
%                 timePeriod = [timePeriod;compression{ii,jj}{kk}(:,2)];
%                 phaseSlope = [phaseSlope;compression{ii,jj}{kk}(:,8)];  
%                 phasePeriod = [phasePeriod;compression{ii,jj}{kk}(:,9)];
%                 if kk == 1 && ii == 2
%                     manipID = [manipID;ones(size(compression{ii,jj}{kk},1),1)*1];
%                 elseif kk == 2 && ii == 2
%                     manipID = [manipID;ones(size(compression{ii,jj}{kk},1),1)*2];
%                 elseif kk == 1 && ii == 3
%                     manipID = [manipID;ones(size(compression{ii,jj}{kk},1),1)*3];
%                 elseif kk == 2 && ii == 3
%                     manipID = [manipID;ones(size(compression{ii,jj}{kk},1),1)*4];
%                 end

                if ii == 2 && kk == 1
                    figure(1)
                    subplot(2,7,7)
                    scatter(placefield_difference{ii,jj}{kk}(idxtoKeep)*1.75,ccgs_place_offset{ii,jj}{kk}(idxtoKeep),25,'.')
                    [R,pVal] = corr(placefield_difference{ii,jj}{kk}(idxtoKeep)*1.75,ccgs_place_offset{ii,jj}{kk}(idxtoKeep),'Rows','pairwise','type','Spearman');
                    %lsline
                    title(strcat('R:',num2str(R),' p:',num2str(pVal)))    
                    xlabel('Place field dist (cm)')
                    ylabel('Place field dist (ms)')                    
                end
            end
        end
    end
    figure(1)
    subplot(2,7,5)  
    stats.timePeriod = groupStats(timePeriod,manipID,'inAxis',true,'color',colMat);
    ylabel('time period')
    
    subplot(2,7,6)  
    stats.timeSlope = groupStats(timeSlope,manipID,'inAxis',true,'color',colMat);
    ylabel('time slope')    
    
    subplot(2,7,12)  
    stats.phasePeriod = groupStats(phasePeriod,manipID,'inAxis',true,'color',colMat);
    ylabel('phase period')    
    
    subplot(2,7,13)  
    stats.phaseSlope = groupStats(phaseSlope,manipID,'inAxis',true,'color',colMat);
    ylabel('phase slope')    

    saveas(gcf,strcat(parentDir,'Compiled\Theta compression\Compression',tag,'.png'));
    saveas(gcf,strcat(parentDir,'Compiled\Theta compression\Compression',tag,'.eps'),'epsc');
    saveas(gcf,strcat(parentDir,'Compiled\Theta compression\Compression',tag,'.fig'));     
    save(strcat(parentDir,'Compiled\Theta compression\Stats',tag,'.mat'),'stats');
   
else
   colMat = [85/243 85/243 85/243;...
            224/243 163/243 46/243;...       
            8/243 133/243 161/243;...                
            56/243 61/243 150/243];      
    figure(1)
    set(gcf,'Renderer','painters')
    set(gcf,'Position',[1 41 1220 563])    
    timeSlope = [];  timePeriod = [];
    phaseSlope = []; phasePeriod = [];
    manipID = [];
    
    for ii = 1:3       
        for jj = 1
            idxtoKeep =  (placefield_center{ii,jj}{1}==1 | placefield_center{ii,jj}{2}==1 ) & (abs(ccgs_place_offset{ii,jj}{1}) <300 & abs(ccgs_place_offset{ii,jj}{2}) <300) &...                    
                region{ii,jj} == 1 & putativeCellType{ii,jj} == 1 & (abs(ccgs_time_offset{ii,jj}{1}) <75 & abs(ccgs_time_offset{ii,jj}{2}) <75);
 
%                 [fit_phase, fit_time] = getCompression(placefield_center{ii,jj}{kk},ccgs_place_offset{ii,jj}{kk},...
%                     ccgs_phase_peaks{ii,jj}{kk},ccgs_time_peaks{ii,jj}{kk},region{ii,jj},putativeCellType{ii,jj});

            figure(1)
            subplot(3,4,4*(ii-1)+1)   
            scatter(ccgs_place_offset{ii,jj}{1}(idxtoKeep),ccgs_time_offset{ii,jj}{1}(idxtoKeep),'.');
            ylabel('Peak offset (ms)')
            xlabel('Theta offset (ms)')
            title('No stim')
            ylim([-50 50])
            
            subplot(3,4,4*(ii-1)+2)   
            scatter(ccgs_place_offset{ii,jj}{2}(idxtoKeep),ccgs_time_offset{ii,jj}{2}(idxtoKeep),'.');
            ylabel('Peak offset (ms)')
            xlabel('Theta offset (ms)')
            title(strcat(reg{ii},' stim'))
            ylim([-50 50])
            
            subplot(3,4,4*(ii-1)+3)   
            scatter(ccgs_place_offset{ii,jj}{1}(idxtoKeep),ccgs_time_offset{ii,jj}{2}(idxtoKeep),'.');
            ylabel('Peak offset (ms)')
            xlabel('Theta offset (ms)')
            title('stim theta vs no stim place')
            ylim([-50 50])            
        end
    end
    
    saveas(figure(1),strcat(parentDir,'Compiled\Theta compression\Compression',tag,'.png'));
    saveas(figure(1),strcat(parentDir,'Compiled\Theta compression\Compression',tag,'.eps'),'epsc');
    saveas(figure(1),strcat(parentDir,'Compiled\Theta compression\Compression',tag,'.fig'));     
    save(strcat(parentDir,'Compiled\Theta compression\Stats',tag,'.mat'),'stats');

end
end

function [fit_params_phase, fit_params_time] = getCompression(placefield_center,ccgs_place_offset,ccgs_phase_peaks,ccgs_time_peaks, region, putativeCellType)
    
% Fit for theta phase
idxKeep = placefield_center==1 & ...
        (ccgs_place_offset>-600 & ccgs_place_offset<600) &...
        region == 1 & putativeCellType == 1;
    
if sum(idxKeep)>10
    StartPoint = [0.001,5,400,20,0.00001,0,0];
    LowerLimits = [0.00005,0,100,5,0,-1,-10];
    UpperLimits = [0.003,10,700,200,0.01,1,100];
    data1 = repmat(ccgs_place_offset(idxKeep),[1 size(ccgs_phase_peaks,2)]);
    data2 = ccgs_phase_peaks(idxKeep,:);
    data1 = reshape(data1,[size(data1,1)*size(data1,2),1]);
    data2 = reshape(data2,[size(data2,1)*size(data2,2),1]);
    x_out = data1;
    y_out = data2;  
    x_bins = -400:8:400;
    y_bins = -20:0.2:20;
    [fit_params_phase,~] = customFit(x_out,y_out,x_bins,y_bins,StartPoint,LowerLimits,UpperLimits);

    % Fit for theta timescale
    StartPoint =  [0.2,   100,  400,   800,  0.00001,  0.6,     0];
    LowerLimits = [0.001,   50,  100,   400, 0.000001, 0.2,  -500];
    UpperLimits = [0.8,  150,  700,  1500,    0.001,  1.5,   500];
    data1 = repmat(ccgs_place_offset(idxKeep),[1 size(ccgs_time_peaks,2)]);
    data2 = ccgs_time_peaks(idxKeep,:);
    data1 = reshape(data1,[size(data1,1)*size(data1,2),1]);
    data2 = reshape(data2,[size(data2,1)*size(data2,2),1]);
    x_out = data1;
    y_out = data2;    
    x_bins = -400:8:400;
    y_bins = -490:10:490; 
    [fit_params_time,~] = customFit(x_out,y_out,x_bins,y_bins,StartPoint,LowerLimits,UpperLimits);    
else
    fit_params_phase(1,1:7) = nan;
    fit_params_time(1,1:7) = nan;
end
end
% 
% function [fit_phase,fit_time] = getCorrCoeff(thetaComp.placefield_center{rr,cc}{zz},thetaComp.ccgs_place_offset{rr,cc}{zz},...
%                     thetaComp.ccgs_phase_peaks{rr,cc}{zz},thetaComp.ccgs_time_peaks{rr,cc}{zz},thetaComp.region{rr,cc},thetaComp.putativeCellType{rr,cc});
% end