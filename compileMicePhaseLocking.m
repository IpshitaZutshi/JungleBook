function compileMicePhaseLocking

parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
tag = 'CA1';
taglfpName = {'pyrCh'};

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','CA1Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final'...
        'IZ21\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Final'};  % To add: IZ23, IZ24, IZ25, IZ26
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ32\Final','IZ33\Final','IZ34\Final','IZ27\Final','IZ28\Final','IZ29\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ32\Saline','IZ33\Saline','IZ34\Saline','IZ27\Saline','IZ28\Saline','IZ29\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Bilateral'};
end

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};


for tt = 1:length(taglfpName)
    taglfp = taglfpName{tt};
    for rr = 1:length(reg)
        for cc = 1:length(target)
            for zz = 1:length(zone)
               compiledPhaseData{rr,cc}{zz} = [];
               p{rr,cc}{zz} = [];
               r{rr,cc}{zz} = [];
               meanPhase{rr,cc}{zz} = [];
               region{rr,cc}{zz} = [];
               celltype{rr,cc}{zz} = [];
            end
        end
    end

    for m = 1:length(mice)

        cd(strcat(parentDir, mice{m},'\Summ'));
        if exist('PhaseLocking.mat','file')
            load('PhaseLocking.mat');
        else 
            disp(['Phase Locking not computed for mouse' mice{m}])
            continue;
        end

        for ii = 1:length(reg)
            for jj = 1:length(target)
                for kk = 1:length(zone)      
                    for pp = 2:size(phaseData.p.(taglfp){ii,jj}{kk},3) % First array is nans                
                        compiledPhaseData{ii,jj}{kk}  = [compiledPhaseData{ii,jj}{kk} phaseData.phasedistros.(taglfp){ii,jj}{kk}(:,:,pp)];
                        p{ii,jj}{kk} = [p{ii,jj}{kk} phaseData.p.(taglfp){ii,jj}{kk}(:,:,pp)];
                        r{ii,jj}{kk} = [r{ii,jj}{kk} phaseData.r.(taglfp){ii,jj}{kk}(:,:,pp)];
                        meanPhase{ii,jj}{kk} = [meanPhase{ii,jj}{kk} phaseData.m.(taglfp){ii,jj}{kk}(:,:,pp)];
                        region{ii,jj}{kk} = [region{ii,jj}{kk} phaseData.region{ii,jj}{kk}(:,:,pp)];
                        celltype{ii,jj}{kk} = [celltype{ii,jj}{kk} phaseData.putativeCellType{ii,jj}{kk}(:,:,pp)];
                    end
                end
            end
        end
    end

    % Divide into CA1 pyr cells and INs
    % Divide into DG/ CA3 pyr cells and INs

    % Angle bins
    binned = linspace(0,2*pi,180+1)';binned(end) = [];
    binSize = binned(2)-binned(1);
    binned = binned + binSize/2;
    leg = {'CA1 Pyr','CA1 IN','DG/CA3 Pyr','DG/CA3 IN'};
    colMat = [85/243 85/243 85/243;...
        224/243 163/243 46/243;... 
        8/243 133/243 161/243;...
        56/243 61/243 150/243];   

    if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
        figure(1)
        set(gcf,'Position',[100 100 1600 1200])    
        figure(2)
        set(gcf,'Position',[100 100 1200 1100])    
        set(gcf,'renderer','painters');    
        for ii = 2:length(reg)
            for jj = 1:length(target)
                % Extract cell id and region info
                b = cellfun(@char,region{ii,jj}{1},'UniformOutput',false);
                [a, ~, iR] = unique(b,'sorted');
                iR = reshape(iR, size(region{ii,jj}{1}));

                % Extract cell id and region info
                b = cellfun(@char,celltype{ii,jj}{1},'UniformOutput',false);
                [c, ~, iC] = unique(b,'sorted');
                iC = reshape(iC, size(celltype{ii,jj}{1}));        

                if length(a{1,1})==1
                    if length(c{1,1})==1
                        idx{1} = iR == 2 & iC == 3; %CA1 pyr
                        idx{2} = iR == 2 & iC == 2; % CA1 IN
                        idx{3} = (iR == 3 | iR == 4) & iC == 3; %DG Pyr
                        idx{4} = (iR == 3 | iR == 4) & iC == 2; %DG IN
                    else
                        idx{1} = iR == 2 & iC == 2; %CA1 pyr
                        idx{2} = iR == 2 & iC == 1; % CA1 IN
                        idx{3} = (iR == 3 | iR == 4) & iC == 2; %DG Pyr
                        idx{4} = (iR == 3 | iR == 4) & iC == 1; %DG IN    
                    end
                else
                    if length(c{1,1})==1
                        idx{1} = iR == 1 & iC == 3; %CA1 pyr
                        idx{2} = iR == 1 & iC == 2; % CA1 IN
                        idx{3} = (iR == 2 | iR == 3) & iC == 3; %DG Pyr
                        idx{4} = (iR == 2 | iR == 3) & iC == 2; %DG IN
                    else
                        idx{1} = iR == 1 & iC == 2; %CA1 pyr
                        idx{2} = iR == 1 & iC == 1; % CA1 IN
                        idx{3} = (iR == 2 | iR == 3) & iC == 2; %DG Pyr
                        idx{4} = (iR == 2 | iR == 3) & iC == 1; %DG IN    
                    end
                end
                for cc = 1:4             
                    for kk = 1:length(zone)
                        for bb = 1:20
                           phasedist(bb,:) = nanmean(compiledPhaseData{ii,jj}{kk}((9*(bb-1)+1):(9*bb),idx{cc}),1);
                        end
                        figure(1)
                        hax = subplot(8,12,24*(cc-1)+12*(jj-1)+6*(ii-2)+kk);
                        bins = [0:18:360];
                        nanRows = mean(~isnan(phasedist))>0;
                        imagesc(bins(2:end),1:size(phasedist(:,nanRows),2),zscore(phasedist(:,nanRows),0,1)')
                        %imagesc(bins(2:end),1:size(phasedist(:,nanRows),2),phasedist(:,nanRows)')
                        hold on
                        caxis([-3 3])            
                        %imagesc((binned*180/pi)+360,1:size(data(kk).phasedistros,2),zscore(data(kk).phasedistros,0,1)')
                        xlim([0 360])
                        set(hax,'XTick',[0 90 180 270 360]) 
                        plot([0:360],cos((pi/180)*[0:360])*0.05*size(phasedist(:,nanRows),2)+0.95*size(phasedist(:,nanRows),2),'color',[1 1 1],'LineWidth',1.5)
                        set(gca,'YDir','normal')
                        xlabel('Theta phase');
                        if kk == 1
                            ylabel(leg{cc},'fontweight','bold');            
                        end
                        title(strcat(target(jj),'.',zone(kk)));            

                        figure(2)
                        if kk ~= 3 && kk ~=6
                            if kk == 1 || kk == 2
                                zoneloc = kk;
                            elseif kk == 4 || kk == 5
                                zoneloc = kk-3;
                            end

                            subplot(4,4,4*(cc-1)+2*(jj-1)+zoneloc)                        

                            if ii==2 && (kk ==1 || kk==2)
                                col = [85/243 85/243 85/243];
                            elseif ii == 2
                                col = [8/243 133/243 161/243];
                            elseif ii==3 && (kk ==1 || kk==2)
                                col = [224/243 163/243 46/243]; 
                            elseif ii==3 
                                col = [56/243 61/243 150/243];                     
                            end

                            avgDist= nanmean(phasedist(:,nanRows),2);
                            stderr = nanstd(phasedist(:,nanRows),0,2)/sqrt(size(phasedist(:,nanRows),2));
                            fill([bins(2:end)'; fliplr(bins(2:end))'],[avgDist-stderr;flipud(avgDist+stderr)],col,'linestyle','none','FaceAlpha',0.2);
                            hold on            
                            plot(bins(2:end),avgDist,'Color',col,'LineWidth',1.5)
                            hold on
                            xlabel('Theta phase');
                            title(strcat(target(jj),' ',zone{zoneloc}(1:(end-1))));  
                            if jj == 1
                                ylabel(leg{cc},'fontweight','bold')
                            end
                        end
                        clear phasedist            
                    end
                end
            end
        end
        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLocking',taglfp,'.png'),'png');
        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLocking',taglfp,'.fig'),'fig');
        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLocking',taglfp,'.eps'),'epsc');

        saveas(figure(2),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\AvgPhaseLocking',taglfp,'.png'),'png');
        saveas(figure(2),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\AvgPhaseLocking',taglfp,'.fig'),'fig');
        saveas(figure(2),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\AvgPhaseLocking',taglfp,'.eps'),'epsc');

    else
        for ii = 1:length(reg)
            figure(ii)
            set(gcf,'Position',[100 100 800 1400])
            set(gcf,'renderer','painters');
            figure(ii+3)
            set(gcf,'Position',[100 100 800 1400])    
            set(gcf,'renderer','painters');
            for jj = 1%:length(target)

                % Extract cell id and region info
                b = cellfun(@char,region{ii,jj}{1},'UniformOutput',false);
                [a, ~, iR] = unique(b,'sorted');
                iR = reshape(iR, size(region{ii,jj}{1}));

                % Extract cell id and region info
                b = cellfun(@char,celltype{ii,jj}{1},'UniformOutput',false);
                [c, ~, iC] = unique(b,'sorted');
                iC = reshape(iC, size(celltype{ii,jj}{1}));        

                if length(a{1,1})==1
                    if length(c{1,1})==1
                        idx{1} = iR == 2 & iC == 3; %CA1 pyr
                        idx{2} = iR == 2 & iC == 2; % CA1 IN
                        idx{3} = (iR == 3 | iR == 4) & iC == 3; %DG Pyr
                        idx{4} = (iR == 3 | iR == 4) & iC == 2; %DG IN
                    else
                        idx{1} = iR == 2 & iC == 2; %CA1 pyr
                        idx{2} = iR == 2 & iC == 1; % CA1 IN
                        idx{3} = (iR == 3 | iR == 4) & iC == 2; %DG Pyr
                        idx{4} = (iR == 3 | iR == 4) & iC == 1; %DG IN    
                    end
                else
                    if length(c{1,1})==1
                        idx{1} = iR == 1 & iC == 3; %CA1 pyr
                        idx{2} = iR == 1 & iC == 2; % CA1 IN
                        idx{3} = (iR == 2 | iR == 3) & iC == 3; %DG Pyr
                        idx{4} = (iR == 2 | iR == 3) & iC == 2; %DG IN
                    else
                        idx{1} = iR == 1 & iC == 2; %CA1 pyr
                        idx{2} = iR == 1 & iC == 1; % CA1 IN
                        idx{3} = (iR == 2 | iR == 3) & iC == 2; %DG Pyr
                        idx{4} = (iR == 2 | iR == 3) & iC == 1; %DG IN    
                    end
                end
                for cc = 1:4             
                    for kk = 1:length(zone)
                        for bb = 1:20
                           phasedist(bb,:) = nanmean(compiledPhaseData{ii,jj}{kk}((9*(bb-1)+1):(9*bb),idx{cc}),1);
                        end
                        figure(ii)
                        set(gcf,'Renderer','painters')
                        hax = subplot(8,6,12*(cc-1)+6*(jj-1)+kk);
                        bins = [0:18:360];
                        nanRows = mean(~isnan(phasedist))>0;
                        imagesc(bins(2:end),1:size(phasedist(:,nanRows),2),zscore(phasedist(:,nanRows),0,1)')
                        %imagesc(bins(2:end),1:size(phasedist(:,nanRows),2),phasedist(:,nanRows)')
                        hold on
                        caxis([-3 3])            
                        %imagesc((binned*180/pi)+360,1:size(data(kk).phasedistros,2),zscore(data(kk).phasedistros,0,1)')
                        xlim([0 360])
                        set(hax,'XTick',[0 90 180 270 360]) 
                        plot([0:360],cos((pi/180)*[0:360])*0.05*size(phasedist(:,nanRows),2)+0.95*size(phasedist(:,nanRows),2),'color',[1 1 1],'LineWidth',1.5)
                        set(gca,'YDir','normal')
                        xlabel('Theta phase');
                        if kk == 1
                            ylabel(leg{cc},'fontweight','bold');            
                        end
                        title(strcat(target(jj),'.',zone(kk)));            

                        figure(ii+3)
                        subplot(4,2,2*(cc-1)+jj)
                        if kk ~= 3 && kk ~=6
                            switch kk
                                case 1
                                    col = colMat(1,:);
                                    lstyle = '-';
                                case 2
                                    col = colMat(1,:);
                                    lstyle = '--';
                                case 3
                                    col = colMat(1,:);
                                    lstyle = ':';
                                case 4
                                    col = colMat(ii+1,:);
                                    lstyle = '-';                    
                                case 5
                                    col = colMat(ii+1,:);
                                    lstyle = '--';                    
                                case 6
                                    col = colMat(ii+1,:);
                                    lstyle = ':';      
                            end
                            avgDist= nanmean(phasedist(:,nanRows),2);
                            stderr = nanstd(phasedist(:,nanRows),0,2)/sqrt(size(phasedist(:,nanRows),2));
                            fill([bins(2:end)'; fliplr(bins(2:end))'],[avgDist-stderr;flipud(avgDist+stderr)],col,'linestyle','none','FaceAlpha',0.2);
                            hold on            
                            plot(bins(2:end),avgDist,'Color',col,'LineWidth',1.5,'LineStyle',lstyle)
                            hold on
                            xlabel('Theta phase');
                            title(target(jj));  
                            if jj == 1
                                ylabel(leg{cc},'fontweight','bold')
                            end
                        end                            

                        clear phasedist            
                    end
                    if cc <3 % Only plot for CA1 cells
                        data{ii,jj,cc}.p = [];
                        data{ii,jj,cc}.m = [];
                        data{ii,jj,cc}.r = [];
                       
                        for kk  = 1:6
                            data{ii,jj,cc}.p = [data{ii,jj,cc}.p p{ii,jj}{kk}(idx{cc})'];
                            data{ii,jj,cc}.r = [data{ii,jj,cc}.r r{ii,jj}{kk}(idx{cc})'];
                            data{ii,jj,cc}.m = [data{ii,jj,cc}.m meanPhase{ii,jj}{kk}(idx{cc})'];
                        end
                    end
                end
            end
            saveas(figure(ii),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLocking',reg{ii},taglfp,'.png'),'png');
            saveas(figure(ii),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLocking',reg{ii},taglfp,'.fig'),'fig');
            saveas(figure(ii),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLocking',reg{ii},taglfp,'.eps'),'epsc');

            saveas(figure(ii+3),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\AvgPhaseLocking',reg{ii},taglfp,'.png'),'png');
            saveas(figure(ii+3),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\AvgPhaseLocking',reg{ii},taglfp,'.fig'),'fig');
            saveas(figure(ii+3),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\AvgPhaseLocking',reg{ii},taglfp,'.eps'),'epsc');    
        end
        close all
        figure
        set(gcf,'Renderer','painters')
        set(gcf,'Position',[60 520 1900 450])
        for ii = 1:3
            for cc = 1:2
                idxSig = (data{ii,1,cc}.p<0.05);
                idxSig = idxSig(:,[1 4 2 5]);

                idStim = zeros(size(idxSig));
                idReg= zeros(size(idxSig));
                idStim(:,[1 3]) = 1; idStim(:,[2 4]) = 2;
                idReg(:,[1 2]) = 1; idReg(:,[3 4]) = 3;            

                subplot(3,6,3*(cc-1)+ii)
                bar(sum(idxSig)./size(idxSig,1),'FaceColor',colMat(ii+1,:));
                ylabel('Proportion of significant cells')
                title(strcat(reg{ii},' ',leg{cc}))
                xticklabels({'RetB','RetS','StemB','StemS'})

                subplot(3,6,3*(cc-1)+ii+6)
                temp = data{ii,1,cc}.r(:,[1 4 2 5]);
                temp = temp(idxSig);
                tempidStim = idStim(idxSig);
                tempidReg = idReg(idxSig);
                groupStats(temp,[tempidStim tempidReg],'inAxis',true,'color',colMat([1 1 1 ii+1],:))
                ylabel('Mean resultant length')

                subplot(3,6,3*(cc-1)+ii+12)
                temp = data{ii,1,cc}.m(:,[1 4 2 5]);
                temp = rad2deg(temp(idxSig));
                tempidStim = idStim(idxSig);
                tempidReg = idReg(idxSig);
                groupStats(temp,[tempidReg tempidStim],'inAxis',true,'color',colMat([1 1 1 ii+1],:)) 
                ylabel('Mean angle')
            end
        end
        saveas(gcf,strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingMeasures',reg{ii},taglfp,'.png'),'png');
        saveas(gcf,strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingMeasures',reg{ii},taglfp,'.fig'),'fig');
        saveas(gcf,strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingMeasures',reg{ii},taglfp,'.eps'),'epsc');
    end
    
end
end