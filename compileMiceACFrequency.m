function compileMiceACFrequency

parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
tag = 'CA3Saline';

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline',...
        'IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Final'}; 
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ34\Final','IZ32\Final','IZ33\Final','IZ27\Final','IZ28\Final','IZ29\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Both'};
end

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for rr = 1:length(reg)
    for cc = 1:length(target)
        for zz = 1:length(zone)
           spectrum{rr,cc}{zz} = [];
           LFPfreq{rr,cc}{zz} = [];
           Unitfreq{rr,cc}{zz} = [];
           Unitpower{rr,cc}{zz} = [];
           freq{rr,cc}{zz} = [];
           region{rr,cc}{zz} = [];
           celltype{rr,cc}{zz} = [];     
        end
    end
end

for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('ACData.mat','file')
        load('ACData.mat');
    else 
        disp(['AC Data not computed for mouse' mice{m}])
        continue;
    end
    
    for ii = 1:length(reg)
        for jj = 1:length(target)
            for kk = 1:length(zone)      
                for pp = 2:size(ACData.Unitfreq{ii,jj}{kk},3) % First array is nans                
                    spectrum{ii,jj}{kk}  = [spectrum{ii,jj}{kk} ACData.spectrum{ii,jj}{kk}(:,:,pp)'];
                    freq{ii,jj}{kk} = ACData.freq{ii,jj}{kk}(:,:,2);
                    LFPfreq{ii,jj}{kk}  = [LFPfreq{ii,jj}{kk} ACData.LFPfreq{ii,jj}{kk}(:,:,pp)];
                    Unitfreq{ii,jj}{kk}  = [Unitfreq{ii,jj}{kk} ACData.Unitfreq{ii,jj}{kk}(:,:,pp)];
                    Unitpower{ii,jj}{kk}  = [Unitpower{ii,jj}{kk} ACData.Unitpower{ii,jj}{kk}(:,:,pp)];
                    region{ii,jj}{kk} = [region{ii,jj}{kk} ACData.region{ii,jj}{kk}(:,:,pp)];
                    celltype{ii,jj}{kk} = [celltype{ii,jj}{kk} ACData.putativeCellType{ii,jj}{kk}(:,:,pp)];
                end
            end
        end
    end
end

% Divide into CA1 pyr cells and INs
% Divide into DG/ CA3 pyr cells and INs

leg = {'CA1 Pyr','CA1 IN','DG/CA3 Pyr','DG/CA3 IN'};
freqRange = freq{2,1}{2};  

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    colMat = [85/243 85/243 85/243;...
        8/243 133/243 161/243;...
        224/243 163/243 46/243;...     
        56/243 61/243 150/243];     
    figure(1)
    set(gcf,'Position',[100 100 1600 1200])    
    figure(2)
    set(gcf,'Position',[100 100 1200 1100])    
    set(gcf,'renderer','painters');    
%     figure(3)
%     set(gcf,'Position',[100 100 1200 1100])
%     set(gcf,'renderer','painters');
    for jj = 1:length(target)
        for cc = 1:4
            datatemp{cc} = [];
            zoneid{cc} = [];
            stimid{cc} = [];
        end
        for ii = 2:length(reg)        
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
                    spectralPow = spectrum{ii,jj}{kk}(:,idx{cc});
                    lfp_freq = LFPfreq{ii,jj}{kk}(idx{cc});
                    unit_freq = Unitfreq{ii,jj}{kk}(idx{cc});
                    unit_power = Unitpower{ii,jj}{kk}(idx{cc});
                    
                    if kk <4
                        idx_pick{kk} = unit_power>0.2;
                    else
                        idx_pick{kk} = idx_pick{kk-3};
                    end
                    spectralPow = spectralPow(:,idx_pick{kk});
                    lfp_freq = lfp_freq(idx_pick{kk});
                    unit_freq = unit_freq(idx_pick{kk});
                    
                    [~,idxsort] = sort(unit_freq);
                    %end                    
                    figure(1)
                    hax = subplot(8,12,24*(cc-1)+12*(jj-1)+6*(ii-2)+kk);
                    %nanRows = mean(~isnan(spectralPow))>0;                    
                    imagesc(freqRange,1:size(spectralPow,2),zscore(spectralPow(:,idxsort),0,1)')
                    caxis([-2.5 2.5]) 
                    hold on
                    %Superimpose cell freq
                    plot(unit_freq(idxsort),1:1:numel(unit_freq),'LineStyle','none','Marker','o','MarkerFaceColor','w','MarkerEdgeColor','none','Markersize',1);                    
                    %Superimpose lfp freq
                    plot(lfp_freq(idxsort),1:1:numel(lfp_freq),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none','Markersize',1);
                    
                   % set(gca,'YDir','normal')
                    xlabel('AC frequency (Hz)');
                    if kk == 1
                        ylabel(leg{cc},'fontweight','bold');            
                    end
                    title(strcat(target(jj),'.',zone(kk)));            

                    
                    if kk ~= 3 && kk ~=6
                        if kk == 1 || kk == 2
                            zoneloc = kk;
                            if ii == 2
                                stimidx = 1;
                            else
                                stimidx = 3;
                            end
                        elseif kk == 4 || kk == 5
                            zoneloc = kk-3;
                            if ii == 2
                                stimidx = 2;
                            else
                                stimidx = 4;
                            end
                        end
%                         figure(2)
%                         subplot(4,4,4*(cc-1)+2*(jj-1)+zoneloc)                        
                        
%                         if ii==2 && (kk ==1 || kk==2)
%                             col = [85/243 85/243 85/243];
%                         elseif ii == 2
%                             col = [8/243 133/243 161/243];
%                         elseif ii==3 && (kk ==1 || kk==2)
%                             col = [224/243 163/243 46/243]; 
%                         elseif ii==3 
%                             col = [56/243 61/243 150/243];                     
%                         end
                        
                        data_temp = (unit_freq-lfp_freq);
                        %dataprop = sum(data_temp>=-2 & data_temp<=4)/sum(~isnan(data_temp));
                        data_temp = data_temp(data_temp>=-2 & data_temp<=4);
                        
                        datatemp{cc} = [datatemp{cc} data_temp];
                        zoneid{cc} = [zoneid{cc} ones(1,length(data_temp))*zoneloc];
                        stimid{cc} = [stimid{cc} ones(1,length(data_temp))*stimidx];
                        %data{jj,zoneloc}(cc,kk) = dataprop; 
%                         h = cdfplot(data_temp);
%                         set(h,'Color',col,'LineWidth',1.5);
%                         y = histcounts((unit_freq-lfp_freq),-4:0.2:4);
%                         plot(-4:0.2:3.8,y,'Color',col,'LineWidth',1.5)
%                         hold on 
%                         xlabel('Spike freq - LFP freq');
%                         title(strcat(target(jj),' ',zone{zoneloc}(1:(end-1))));  
%                         if jj == 1
%                             ylabel(leg{cc},'fontweight','bold')
%                         end
                        
%                         figure(3)
%                         subplot(4,4,4*(cc-1)+2*(jj-1)+zoneloc)  
%                         h = cdfplot(unit_power);
%                         set(h,'Color',col,'LineWidth',1.5);
%                         hold on 
%                         xlabel('AC power');  
%                         xlim([0 4])
%                         title(strcat(target(jj),' ',zone{zoneloc}(1:(end-1))));    
%                         if jj == 1
%                             ylabel(leg{cc},'fontweight','bold')
%                         end                        
                    end      
                end
            end
        end
        figure(2)
        for cc = 1:4
            subplot(4,2,2*(cc-1)+jj)    
            stats{jj,cc} = groupStats(datatemp{cc},[zoneid{cc}' stimid{cc}'],'inAxis',true,'color',[colMat;colMat]);
            title('Side___________________Central')
        end
    end
  
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC FrequencyHeatPlot.png'),'png');
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC FrequencyHeatPlot.fig'),'fig');
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC FrequencyHeatPlot.eps'),'epsc');

    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Frequency.png'),'png');
    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Frequency.fig'),'fig');
    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Frequency.eps'),'epsc');
    save(strcat(parentDir,'\Compiled\AC Frequency\',tag,'\Stats.mat'),'stats');

else
    colMat = [85/243 85/243 85/243;...       
        224/243 163/243 46/243;...    
        8/243 133/243 161/243;...           
        56/243 61/243 150/243];    
    for ii = 1:length(reg)
        figure(ii)
        set(gcf,'Position',[100 100 800 1400])
        set(gcf,'renderer','painters');
        figure(ii+3)
        set(gcf,'Position',[100 100 800 1400])    
        set(gcf,'renderer','painters');
%         figure(ii+6)
%         set(gcf,'Position',[100 100 800 1400])    
%         set(gcf,'renderer','painters');        
        for jj = 1:length(target)
            for cc = 1:4
                datatemp{cc} = [];
                zoneid{cc} = [];
                stimid{cc} = [];
            end
            % Extract cell id and region info
            if ~isempty(region{ii,jj}{1})
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
                        spectralPow = spectrum{ii,jj}{kk}(:,idx{cc});
                        lfp_freq = LFPfreq{ii,jj}{kk}(idx{cc});
                        unit_freq = Unitfreq{ii,jj}{kk}(idx{cc});
                        unit_power = Unitpower{ii,jj}{kk}(idx{cc});

                        if kk <4
                            idx_pick{kk} = unit_power>0.2;
                        else
                            idx_pick{kk} = idx_pick{kk-3};
                        end
                        spectralPow = spectralPow(:,idx_pick{kk});
                        lfp_freq = lfp_freq(idx_pick{kk});
                        unit_freq = unit_freq(idx_pick{kk});

                        [~,idxsort{kk}] = sort(unit_freq);

                        figure(ii)
                        hax = subplot(8,6,12*(cc-1)+6*(jj-1)+kk);
                        imagesc(freqRange,1:size(spectralPow,2),zscore(spectralPow(:,idxsort{kk}),0,1)')
                        caxis([-2.5 2.5]) 
                        hold on
                        %Superimpose cell freq
                        plot(unit_freq(idxsort{kk}),1:1:numel(unit_freq),'LineStyle','none','Marker','o','MarkerFaceColor','w','MarkerEdgeColor','none','Markersize',1);                    
                        %Superimpose lfp freq
                        plot(lfp_freq(idxsort{kk}),1:1:numel(lfp_freq),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none','Markersize',1);

                       % set(gca,'YDir','normal')
                        xlabel('AC frequency (Hz)');
                        if kk == 1
                            ylabel(leg{cc},'fontweight','bold');            
                        end
                        title(strcat(target(jj),'.',zone(kk)));            


                        if kk ~= 3 && kk ~=6
                            switch kk
                                case 1
                                    col = colMat(1,:);
                                    lstyle = '-';
                                    stimidx = 1;
                                    zoneloc = kk;
                                case 2
                                    col = colMat(1,:);
                                    lstyle = '--';
                                    stimidx = 1;
                                    zoneloc = kk;
                                case 4
                                    col = colMat(ii+1,:);
                                    lstyle = '-';    
                                    stimidx = 2;
                                    zoneloc = kk-3;
                                case 5
                                    col = colMat(ii+1,:);
                                    lstyle = '--';   
                                    stimidx = 2;
                                    zoneloc = kk-3;
                            end

%                             figure(ii+3)
%                             subplot(4,2,2*(cc-1)+jj);
                            data_temp = (unit_freq-lfp_freq);
                            data_temp = data_temp(data_temp>=-2 & data_temp<=4);
                            datatemp{cc} = [datatemp{cc} data_temp];
                            zoneid{cc} = [zoneid{cc} ones(1,length(data_temp))*zoneloc];
                            stimid{cc} = [stimid{cc} ones(1,length(data_temp))*stimidx];
%                             if ~isempty(data_temp)
%                                 h = cdfplot(data_temp);
%                                 set(h,'Color',col,'LineWidth',1.5,'LineStyle',lstyle);
%                                 hold on 
%                                 xlabel('Spike freq - LFP freq');
%                             end
%                             title(strcat(target(jj)));  
%                             if jj == 1
%                                 ylabel(leg{cc},'fontweight','bold')
%                             end

%                             figure(ii+6)
%                             subplot(4,2,2*(cc-1)+jj);
%                             if ~isempty(unit_power)
%                                 h = cdfplot(unit_power);
%                                 set(h,'Color',col,'LineWidth',1.5,'LineStyle',lstyle);
%                                 hold on 
%                                 xlabel('AC power');  
%                             end
%                             xlim([0 4])
%                             title(strcat(target(jj)));  
%                             if jj == 1
%                                 ylabel(leg{cc},'fontweight','bold')
%                             end

                        end
                        clear phasedist            
                    end
                    figure(ii+3)
                    subplot(4,2,2*(cc-1)+jj);
                    if cc<3
                        stats{ii,jj}{cc} = groupStats(datatemp{cc},[zoneid{cc}' stimid{cc}'],'inAxis',true,'color',[colMat(1,:);colMat(ii+1,:);colMat(1,:);colMat(ii+1,:)]);
                        title(strcat(leg{cc},'--',target{jj},'Side__Central'))
                        ylabel('Cell freq-lfp freq')
                    end
                end
            end
        end
    end
        
    saveas(figure(ii),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC FrequencyHeatPlot',reg{ii},'.png'),'png');
    saveas(figure(ii),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC FrequencyHeatPlot',reg{ii},'.fig'),'fig');
    saveas(figure(ii),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC FrequencyHeatPlot',reg{ii},'.eps'),'epsc');

    saveas(figure(ii+3),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Frequency',reg{ii},'.png'),'png');
    saveas(figure(ii+3),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Frequency',reg{ii},'.fig'),'fig');
    saveas(figure(ii+3),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Frequency',reg{ii},'.eps'),'epsc');  
    
    save(strcat(parentDir,'\Compiled\AC Frequency\',tag,'\Stats.mat'),'stats');    
%         
%         saveas(figure(ii+6),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Power',reg{ii},'.png'),'png');
%         saveas(figure(ii+6),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Power',reg{ii},'.fig'),'fig');
%         saveas(figure(ii+6),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\AC Power',reg{ii},'.eps'),'epsc');           
    end
end
