function compileMiceACFrequencyVer2

parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
tag = 'CA3';

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline',...
        'IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Final'}; 
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ34\Final','IZ33\Final','IZ32\Final','IZ27\Final','IZ28\Final','IZ29\Final'};
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
        for zz = 1:length(zone)
           spectrum{rr,cc}{zz}  = [];
           freq{rr,cc}{zz} = [];
           LFPfreq{rr,cc}{zz}  = [];
           Unitfreq{rr,cc}{zz}  = [];
           Unitpower{rr,cc}{zz}  = [];
           region{rr,cc}{zz} = [];
           celltype{rr,cc}{zz} = [];
           placefield{rr,cc}{zz} = [];
           location{rr,cc}{zz} = [];
           avgRate{rr,cc}{zz} = [];             
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
            for kk = 1:4
               for pp = 2:size(ACData.placefield{ii,jj}{kk},3) % First array is nans   
                   placefield{ii,jj}{kk} = [placefield{ii,jj}{kk}; ACData.placefield{ii,jj}{kk}(:,:,pp)];
                   location{ii,jj}{kk} = [location{ii,jj}{kk}; ACData.fieldMax{ii,jj}{kk}(:,:,pp)'];
                   avgRate{ii,jj}{kk} = [avgRate{ii,jj}{kk}; ACData.avgRate{ii,jj}{kk}(:,:,pp)'];                    
               end
            end
        end
    end
end

% Divide into CA1 pyr cells and INs

leg = {'CA1 Pyr','CA1 IN','DG/CA3 Pyr','DG/CA3 IN'};
freqRange = freq{2,1}{2};  
ret = round((100/175)*80);
stem = round((100/175)*110);

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    colMat = [85/243 85/243 85/243;...
        8/243 133/243 161/243;...
        224/243 163/243 46/243;...     
        56/243 61/243 150/243];     
    figure(1)
    set(gcf,'Position',[1081 30 1920 936])
    set(gcf,'renderer','painters');    
    figure(2)
    set(gcf,'Position',[100 100 800 800])    
    set(gcf,'renderer','painters');

    for jj = 1%:length(target)
        for cc = 1:2 % pyr/interneuron
            for nn = 1:2 % side/ central
                datatemp{nn,cc} = [];
                unitfreqtemp{nn,cc} = [];
                lfpfreqtemp{nn,cc} = [];
                stimid{nn,cc} = [];                
            end
        end
        % Extract cell id and region info
        if ~isempty(region{2,jj}{1})
            b = cellfun(@char,region{2,jj}{1},'UniformOutput',false);
            [a, ~, iR] = unique(b,'sorted');
            iR = reshape(iR, size(region{ii,jj}{1}));

            % Extract cell id and region info
            b = cellfun(@char,celltype{2,jj}{1},'UniformOutput',false);
            [c, ~, iC] = unique(b,'sorted');
            iC = reshape(iC, size(celltype{2,jj}{1}));        
            
            %Find stem and return place cells
            iCentral = (location{2,jj}{1}>=stem | location{2,jj}{2}>=stem | location{2,jj}{3}>=stem | location{2,jj}{4}>=stem |...
                location{3,jj}{1}>=stem | location{3,jj}{2}>=stem | location{3,jj}{3}>=stem | location{3,jj}{4}>=stem) & ...
                (location{2,jj}{1}<95 | location{2,jj}{2}<95 | location{2,jj}{3}<95| location{2,jj}{4}<95|...
                location{3,jj}{1}<95 | location{3,jj}{2}<95| location{3,jj}{3}<95 | location{3,jj}{4}<95);
            iSide = location{2,jj}{1}<=ret | location{2,jj}{2}<=ret| location{2,jj}{3}<=ret| location{2,jj}{4}<=ret|...
                location{3,jj}{1}<=ret| location{3,jj}{2}<=ret| location{3,jj}{3}<=ret| location{3,jj}{4}<=ret;
            
            if length(a{1,1})==1
                if length(c{1,1})==1              
                    idx{1} = iR == 2 & iC == 3; %CA1 pyr
                    idx{2} = iR == 2 & iC == 2; % CA1 IN                      
                else
                    idx{1} = iR == 2 & iC == 2; %CA1 pyr
                    idx{2} = iR == 2 & iC == 1; % CA1 IN
                end
            else
                if length(c{1,1})==1
                    idx{1} = iR == 1 & iC == 3; %CA1 pyr
                    idx{2} = iR == 1 & iC == 2; % CA1 IN
                else
                    idx{1} = iR == 1 & iC == 2; %CA1 pyr
                    idx{2} = iR == 1 & iC == 1; % CA1 IN
                end
            end
            for ii = 1:3
                for cc = 1:2   
                    for kk = 1:2
                        if ii == 1
                            lfp_freq_base = LFPfreq{2,jj}{kk}(idx{cc});
                            unit_freq_base = Unitfreq{2,jj}{kk}(idx{cc});
                            unit_power_base = Unitpower{2,jj}{kk}(idx{cc});
                            lfp_freq = LFPfreq{2,jj}{kk+3}(idx{cc});
                            unit_freq = Unitfreq{2,jj}{kk+3}(idx{cc});
                            unit_power = Unitpower{2,jj}{kk+3}(idx{cc});
                        elseif ii == 2
                            lfp_freq_base = LFPfreq{2,jj}{kk}(idx{cc});
                            unit_freq_base = Unitfreq{2,jj}{kk}(idx{cc});
                            unit_power_base = Unitpower{2,jj}{kk}(idx{cc});
                            lfp_freq = LFPfreq{3,jj}{kk}(idx{cc});
                            unit_freq = Unitfreq{3,jj}{kk}(idx{cc});
                            unit_power = Unitpower{3,jj}{kk}(idx{cc});
                        elseif ii == 3
                            lfp_freq_base = LFPfreq{2,jj}{kk}(idx{cc});
                            unit_freq_base = Unitfreq{2,jj}{kk}(idx{cc});
                            unit_power_base = Unitpower{2,jj}{kk}(idx{cc});
                            lfp_freq = LFPfreq{3,jj}{kk+3}(idx{cc});
                            unit_freq = Unitfreq{3,jj}{kk+3}(idx{cc});
                            unit_power = Unitpower{3,jj}{kk+3}(idx{cc});
                        end
                        idxSide = iSide(idx{cc});
                        idxCentral = iCentral(idx{cc});                        
%                       
                        if cc == 1 && kk == 1 % Only pick cells with fields in the return arm
                            idx_pick_base = idxSide' & ((unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11));
                            idx_pick = idxSide' & ((unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) |(unit_power>0.1 & unit_freq>=7 & unit_freq<=11));
                        elseif cc == 1 && kk == 2 %  Only pick cells with fields in the central arm
                            idx_pick_base = idxCentral' & ((unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11));
                            idx_pick = idxCentral' &((unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11));
                        else                    % If INs, pick all
                            idx_pick_base = (unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11);
                            idx_pick = (unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11);
                        end

                        lfp_freq_base = lfp_freq_base(idx_pick_base);
                        lfp_freq = lfp_freq(idx_pick);
                        unit_freq_base = unit_freq_base(idx_pick_base);
                        unit_freq = unit_freq(idx_pick);      
                        unitdiff_base = (unit_freq_base-lfp_freq_base);
                        unitdiff = (unit_freq-lfp_freq);

                        [~,idxsortbase] = sort(unit_freq_base,'descend');
                        [~,idxsort] = sort(unit_freq,'descend');
%                         [~,idxsortbase] = sort(unitdiff_base,'descend');
%                         [~,idxsort] = sort(unitdiff,'descend');

                        figure(1)
                        subplot(3,4,4*(ii-1)+2*(cc-1)+kk);
                        hold on
%                         h = cdfplot(unitdiff_base);
%                         set(h,'color',colMat(1,:),'LineWidth',1.5)
                        % Plot baseline
                        %Superimpose cell freq
                        plot(unit_freq_base(idxsortbase),1:1:numel(unit_freq_base),'LineStyle','none','Marker','o','MarkerFaceColor',colMat(1,:),'MarkerEdgeColor','none','Markersize',2);                    
                        %Superimpose lfp freq
                        plot(lfp_freq_base(idxsortbase),1:1:numel(lfp_freq_base),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none','Markersize',2);
                        line([nanmean(lfp_freq_base) nanmean(lfp_freq_base)],[0 numel(lfp_freq_base)],'Color','k')
                        line([nanmean(unit_freq_base) nanmean(unit_freq_base)],[0 numel(unit_freq_base)],'Color','k')                        
                        xlim([7 11])
                        ylim([0 numel(lfp_freq_base)])
                        xlabel('AC frequency (Hz)');
                        ylabel('Cell number')
                        if kk == 1
                            ylabel(leg{cc},'fontweight','bold');            
                        end                        
%                         h = cdfplot(unitdiff);
%                         set(h,'color',colMat(ii+1,:),'LineWidth',1.5)
                        % Plot stim
                        %Superimpose cell freq
                        plot(unit_freq(idxsort),1:1:numel(unit_freq),'LineStyle','none','Marker','o','MarkerFaceColor',colMat(ii+1,:),'MarkerEdgeColor','none','Markersize',2);                    
                        %Superimpose lfp freq
                        plot(lfp_freq(idxsort),1:1:numel(lfp_freq),'LineStyle','none','Marker','o','MarkerFaceColor','c','MarkerEdgeColor','none','Markersize',2);
                        line([nanmean(lfp_freq) nanmean(lfp_freq)],[0 numel(lfp_freq)],'Color','r')
                        line([nanmean(unit_freq) nanmean(unit_freq)],[0 numel(unit_freq)],'Color','r')                        
                        xlim([7 11])
                        %ylim([0 numel(lfp_freq)])
                        xlabel('AC frequency (Hz)');
                        ylabel('Cell number')                                                                       
                        if kk == 1
                            ylabel(leg{cc},'fontweight','bold');            
                        end
                        title(strcat(reg{ii},'.',zone(kk),num2str(length(unit_freq_base)),',',num2str(length(unit_freq))));            

                        if  ii == 1
                            data_temp = (unit_freq_base-lfp_freq_base);
                            datatemp{kk,cc} = [datatemp{kk,cc} data_temp];
                            unitfreqtemp{kk,cc} = [unitfreqtemp{kk,cc} unit_freq_base];
                            lfpfreqtemp{kk,cc} = [lfpfreqtemp{kk,cc} lfp_freq_base];                        
                            stimid{kk,cc} = [stimid{kk,cc} ones(1,length(data_temp))*1];
                        end
                        
                        data_temp = (unit_freq-lfp_freq);
                        datatemp{kk,cc} = [datatemp{kk,cc} data_temp];
                        unitfreqtemp{kk,cc} = [unitfreqtemp{kk,cc} unit_freq];
                        lfpfreqtemp{kk,cc} = [lfpfreqtemp{kk,cc} lfp_freq];         
                        if ii == 1
                            stimid{kk,cc} = [stimid{kk,cc} ones(1,length(data_temp))*2]; 
                        elseif  ii == 2
                            stimid{kk,cc} = [stimid{kk,cc} ones(1,length(data_temp))*3];
                        elseif ii == 3
                            stimid{kk,cc} = [stimid{kk,cc} ones(1,length(data_temp))*4];
                        end
                    end                        
                end                 
            end
        end                
    end                              
    figure(2)
    for kk = 1:2
        for cc = 1:2
            subplot(2,2,2*(kk-1)+cc);
            statsdiff{kk,cc} = groupStats(datatemp{kk,cc}(stimid{kk,cc}~=2),stimid{kk,cc}(stimid{kk,cc}~=2),'repeatedMeasures',true,'inAxis',true,'color',colMat);
            statsunitfreq{kk,cc} = groupStats(unitfreqtemp{kk,cc}(stimid{kk,cc}~=2),stimid{kk,cc}(stimid{kk,cc}~=2),'repeatedMeasures',true,'inAxis',true,'color',colMat,'doPlot',false);
            statslfpfreq{kk,cc} = groupStats(lfpfreqtemp{kk,cc}(stimid{kk,cc}~=2),stimid{kk,cc}(stimid{kk,cc}~=2),'repeatedMeasures',true,'inAxis',true,'color',colMat,'doPlot',false);
            title(strcat(leg{cc},'-',zone{kk}))
            ylabel('Cell freq-lfp freq')
        end
    end
    
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACSummary.png'),'png');
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACSummary.fig'),'fig');
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACSummary.eps'),'epsc');

    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACStats.png'),'png');
    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACStats.fig'),'fig');
    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACStats.eps'),'epsc');  

    save(strcat(parentDir,'\Compiled\AC Frequency\',tag,'\Stats.mat'),'statsdiff','statsunitfreq','statslfpfreq');    
else
    colMat = [85/243 85/243 85/243;...       
        224/243 163/243 46/243;...    
        8/243 133/243 161/243;...           
        56/243 61/243 150/243];  
    
    figure(1)
    set(gcf,'Position',[1081 30 1920 936])
    set(gcf,'renderer','painters');    
    figure(2)
    set(gcf,'Position',[100 100 800 800])    
    set(gcf,'renderer','painters');

    for ii = 1:length(reg)
        for jj = 1%:length(target)
            for cc = 1:2
                for nn = 1:2
                    datatemp{nn,cc} = [];
                    unitfreqtemp{nn,cc} = [];
                    lfpfreqtemp{nn,cc} = [];
                    zoneid{nn,cc} = [];
                    stimid{nn,cc} = [];                
                end
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

                %Find stem and return place cells
                iCentral = (location{ii,jj}{1}>=stem | location{ii,jj}{2}>=stem | location{ii,jj}{3}>=stem | location{ii,jj}{4}>=stem) & ...
                    (location{ii,jj}{1}<95| location{ii,jj}{2}<95 | location{ii,jj}{3}<95 | location{ii,jj}{4}<95);
                iSide = location{ii,jj}{1}<=ret | location{ii,jj}{2}<=ret| location{ii,jj}{3}<=ret| location{ii,jj}{4}<=ret;
                
                if length(a{1,1})==1
                    if length(c{1,1})==1
                        idx{1} = iR == 2 & iC == 3; %CA1 pyr
                        idx{2} = iR == 2 & iC == 2; % CA1 IN
                    else
                        idx{1} = iR == 2 & iC == 2; %CA1 pyr
                        idx{2} = iR == 2 & iC == 1; % CA1 IN
                    end
                else
                    if length(c{1,1})==1
                        idx{1} = iR == 1 & iC == 3; %CA1 pyr
                        idx{2} = iR == 1 & iC == 2; % CA1 IN
                    else
                        idx{1} = iR == 1 & iC == 2; %CA1 pyr
                        idx{2} = iR == 1 & iC == 1; % CA1 IN
                    end
                end
                for cc = 1:2             
                    for kk = [1 2]     
                        
                        lfp_freq_base = LFPfreq{ii,jj}{kk}(idx{cc});
                        unit_freq_base = Unitfreq{ii,jj}{kk}(idx{cc});
                        unit_power_base = Unitpower{ii,jj}{kk}(idx{cc});
                        lfp_freq = LFPfreq{ii,jj}{kk+3}(idx{cc});
                        unit_freq = Unitfreq{ii,jj}{kk+3}(idx{cc});
                        unit_power = Unitpower{ii,jj}{kk+3}(idx{cc});
                        idxSide = iSide(idx{cc});
                        idxCentral = iCentral(idx{cc});                               
                        
                        if cc == 1 && kk == 1 % Only pick cells with fields in the return arm
                            idx_pick_base = idxSide' & ((unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11));
                            idx_pick = idxSide' & ((unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) |(unit_power>0.1 & unit_freq>=7 & unit_freq<=11));
                        elseif cc == 1 && kk == 2 %  Only pick cells with fields in the central arm
                            idx_pick_base = idxCentral' & ((unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11));
                            idx_pick = idxCentral' &((unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11));
                        else                    % If INs, pick all
                            idx_pick_base = (unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11);
                            idx_pick = (unit_power_base>0.1 & unit_freq_base>=7 & unit_freq_base<=11) | (unit_power>0.1 & unit_freq>=7 & unit_freq<=11);
                        end

%                         idx_pick_base = (unit_freq_base>=7 & unit_freq_base<=11);%unit_power>0.1 & (unit_freq>=7 & unit_freq<=11) | 
%                         idx_pick = (unit_freq>=7 & unit_freq<=11);% &(unit_freq_base>=7 & unit_freq_base<=11);
                        
                        lfp_freq_base = lfp_freq_base(idx_pick_base);
                        unit_freq_base = unit_freq_base(idx_pick_base);
                        lfp_freq = lfp_freq(idx_pick);
                        unit_freq = unit_freq(idx_pick);

                        [~,idxsort_base] = sort(unit_freq_base,'descend');
                        [~,idxsort] = sort(unit_freq,'descend');
                        
                        figure(1)
                        subplot(3,4,4*(ii-1)+2*(cc-1)+kk);
                        hold on
                        %Superimpose cell freq
                        plot(unit_freq_base(idxsort_base),1:1:numel(unit_freq_base),'LineStyle','none','Marker','o','MarkerFaceColor',colMat(1,:),'MarkerEdgeColor','none','Markersize',2);                    
                        %Superimpose lfp freq
                        plot(lfp_freq_base(idxsort_base),1:1:numel(lfp_freq_base),'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','none','Markersize',2);
                        line([nanmean(lfp_freq_base) nanmean(lfp_freq_base)],[0 numel(lfp_freq_base)],'Color','m','LineWidth',1.5)
                        line([nanmean(unit_freq_base) nanmean(unit_freq_base)],[0 numel(unit_freq_base)],'Color','m','LineWidth',1.5) 
                        xlim([7 11])
                        ylim([0 numel(lfp_freq_base)])
                        
                         % Plot stim
                        %Superimpose cell freq
                        plot(unit_freq(idxsort),1:1:numel(unit_freq),'LineStyle','none','Marker','o','MarkerFaceColor',colMat(ii+1,:),'MarkerEdgeColor','none','Markersize',2);                    
                        %Superimpose lfp freq
                        plot(lfp_freq(idxsort),1:1:numel(lfp_freq),'LineStyle','none','Marker','o','MarkerFaceColor','c','MarkerEdgeColor','none','Markersize',2);
                        line([nanmean(lfp_freq) nanmean(lfp_freq)],[0 numel(lfp_freq)],'Color','r')
                        line([nanmean(unit_freq) nanmean(unit_freq)],[0 numel(unit_freq)],'Color','r')                        
                        xlim([7 11])
                        %ylim([0 numel(lfp_freq)])
                        xlabel('AC frequency (Hz)');
                        ylabel('Cell number')                                                                       
                        if kk == 1
                            ylabel(leg{cc},'fontweight','bold');            
                        end
                        title(strcat(reg{ii},'.',zone(kk),num2str(length(unit_freq_base)),',',num2str(length(unit_freq))));           

                        data_temp = (unit_freq_base-lfp_freq_base);
                        datatemp{kk,cc} = [datatemp{kk,cc} data_temp];
                        unitfreqtemp{kk,cc} = [unitfreqtemp{kk,cc} unit_freq_base];
                        lfpfreqtemp{kk,cc} = [lfpfreqtemp{kk,cc} lfp_freq_base];                        
                        stimid{kk,cc} = [stimid{kk,cc} ones(1,length(data_temp))*1];
                        
                        data_temp = (unit_freq-lfp_freq);
                        datatemp{kk,cc} = [datatemp{kk,cc} data_temp];
                        unitfreqtemp{kk,cc} = [unitfreqtemp{kk,cc} unit_freq];
                        lfpfreqtemp{kk,cc} = [lfpfreqtemp{kk,cc} lfp_freq];   
                        stimid{kk,cc} = [stimid{kk,cc} ones(1,length(data_temp))*2];
                    end
                    
                    figure(2)
                    subplot(3,4,4*(ii-1)+2*(cc-1)+1);
                    statsdiff{ii,jj}{1,cc} = groupStats(datatemp{1,cc},stimid{1,cc},'repeatedMeasures',true,'inAxis',true,'color',[colMat(1,:);colMat(ii+1,:)]);%,'repeatedMeasures',true);
                    statsunitfreq{ii,jj}{1,cc} = groupStats(unitfreqtemp{1,cc},stimid{1,cc},'repeatedMeasures',true,'inAxis',true,'doPlot',false);%,'repeatedMeasures',true);
                    statslfpfreq{ii,jj}{1,cc} = groupStats(lfpfreqtemp{1,cc},stimid{1,cc},'repeatedMeasures',true,'inAxis',true,'doPlot',false);%,'repeatedMeasures',true);
                    title(strcat(leg{cc},'--Side Arm'))
                    ylabel('Cell freq-lfp freq')
                    
                    subplot(3,4,4*(ii-1)+2*(cc-1)+2);
                    statsdiff{ii,jj}{2,cc} = groupStats(datatemp{2,cc},stimid{2,cc},'repeatedMeasures',true,'inAxis',true,'color',[colMat(1,:);colMat(ii+1,:)]);%,'repeatedMeasures',true);
                    statsunitfreq{ii,jj}{2,cc} = groupStats(unitfreqtemp{2,cc},stimid{2,cc},'repeatedMeasures',true,'inAxis',true,'doPlot',false);%,'repeatedMeasures',true);
                    statslfpfreq{ii,jj}{2,cc} = groupStats(lfpfreqtemp{2,cc},stimid{2,cc},'repeatedMeasures',true,'inAxis',true,'doPlot',false);%,'repeatedMeasures',true);
                    [statsdiff{ii,jj}{2,cc}.signrank.p,~,statsdiff{ii,jj}{2,cc}.signrank.stats] = signrank(datatemp{2,cc}(stimid{2,cc}==1),datatemp{2,cc}(stimid{2,cc}==2));
                    [statsunitfreq{ii,jj}{2,cc}.signrank.p,~,statsunitfreq{ii,jj}{2,cc}.signrank.stats] = signrank(unitfreqtemp{2,cc}(stimid{2,cc}==1),unitfreqtemp{2,cc}(stimid{2,cc}==2));
                    [statslfpfreq{ii,jj}{2,cc}.signrank.p,~,statslfpfreq{ii,jj}{2,cc}.signrank.stats] = signrank(lfpfreqtemp{2,cc}(stimid{2,cc}==1),lfpfreqtemp{2,cc}(stimid{2,cc}==2));
                    title(strcat(leg{cc},'--Central Arm'))
                    ylabel('Cell freq-lfp freq')
                    

                end                    
            end                
        end            
    end   
    
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACSummary.png'),'png');
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACSummary.fig'),'fig');
    saveas(figure(1),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACSummary.eps'),'epsc');

    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACStats.png'),'png');
    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACStats.fig'),'fig');
    saveas(figure(2),strcat(parentDir,'\Compiled\AC Frequency\',tag,'\ACStats.eps'),'epsc');  

    save(strcat(parentDir,'\Compiled\AC Frequency\',tag,'\Stats.mat'),'statsdiff','statsunitfreq','statslfpfreq');    
end    

end

