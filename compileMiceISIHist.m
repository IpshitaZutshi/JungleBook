function compileMiceISIHist

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
zonedef = {'Side','Center'};
target = {'STEM', 'RETURN'};
bins = -3:0.05:1.5;
freq = 1./(10.^bins(1:90));

for rr = 1:length(reg)
    for cc = 1:length(target)
        for zz = 1:length(zone)
           spectrum{rr,cc}{zz} = [];
           LFPfreq{rr,cc}{zz} = [];
           region{rr,cc}{zz} = [];
           celltype{rr,cc}{zz} = [];     
        end
    end
end

for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('ISIData.mat','file')
        load('ISIData.mat');
    else 
        disp(['ISI Data not computed for mouse' mice{m}])
        continue;
    end
    
    for ii = 1:length(reg)
        for jj = 1:length(target)
            for kk = [1 2 4 5]               
                spectrum{ii,jj}{kk}  = [spectrum{ii,jj}{kk}; ISIData.ISI{ii,jj}{kk}];
                LFPfreq{ii,jj}{kk}  = [LFPfreq{ii,jj}{kk}; ISIData.LFPfreq{ii,jj}{kk}];
                region{ii,jj}{kk} = [region{ii,jj}{kk}; ISIData.region{ii,jj}{kk}(2:end)];
                celltype{ii,jj}{kk} = [celltype{ii,jj}{kk}; ISIData.putativeCellType{ii,jj}{kk}(2:end)];
            end
        end
    end
end

% Divide into CA1 pyr cells and INs

leg = {'CA1 Pyr','CA1 IN','DG/CA3 Pyr','DG/CA3 IN'};

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    colMat = [85/243 85/243 85/243;...
        8/243 133/243 161/243;...
        224/243 163/243 46/243;...     
        56/243 61/243 150/243];     
    
    figure(1)
    set(gcf,'Position',[1081 30 1566 527])
    set(gcf,'renderer','painters');    
    figure(2)
    set(gcf,'Position',[1081 30 1566 527])
    set(gcf,'renderer','painters');        

    for jj = 1%:length(target)

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
                figure(cc)
                for kk = [1 2]     

                    ISIbase = spectrum{2,jj}{kk}(idx{cc},:)./sum(spectrum{2,jj}{kk}(idx{cc},:),2);
                    ISImEC = spectrum{2,jj}{kk+3}(idx{cc},:)./sum(spectrum{2,jj}{kk+3}(idx{cc},:),2);
                    ISICA3 = spectrum{3,jj}{kk}(idx{cc},:)./sum(spectrum{3,jj}{kk}(idx{cc},:),2);  
                    ISIBoth = spectrum{3,jj}{kk+3}(idx{cc},:)./sum(spectrum{3,jj}{kk+3}(idx{cc},:),2);     
                    
                    LFPbase = LFPfreq{2,jj}{kk}(idx{cc});
                    LFPmEC = LFPfreq{2,jj}{kk+3}(idx{cc});
                    LFPCA3 = LFPfreq{3,jj}{kk}(idx{cc});
                    LFPBoth = LFPfreq{3,jj}{kk+3}(idx{cc});

                    [~,idxsort_base] = sort(nanmean(ISIbase(:,10:30),2),'descend');
                    colMap = cbrewer('seq','Blues',20);
                    colormap(colMap)

                    subplot(2,5,5*(kk-1)+1);
                    imagesc(freq(10:80),1:1:size(ISIbase,1),zscore(ISIbase(idxsort_base,10:80),[],2),'AlphaData', .8)                    
                    caxis([-0.5 3])
                    set(gca, 'XScale', 'log')
                    set(gca, 'XDir', 'reverse')
                    hold on
                    line([nanmean(LFPbase) nanmean(LFPbase)],[1 size(ISIbase,1)],'Color','r','LineWidth',1)
                    xlabel('Frequency')
                    xticks([0.1 1 10 100]) 
                    ylabel('Cell number')
                    title(strcat('No stim ', zonedef{kk}))   

                    subplot(2,5,5*(kk-1)+2);
                    imagesc(freq(10:80),1:1:size(ISImEC,1),zscore(ISImEC(idxsort_base,10:80),[],2),'AlphaData', .8)
                    caxis([-0.5 3])
                    set(gca, 'XScale', 'log')
                    set(gca, 'XDir', 'reverse')
                    hold on
                    line([nanmean(LFPmEC) nanmean(LFPmEC)],[1 size(ISImEC,1)],'Color','r','LineWidth',1)
                    xlabel('Frequency')
                    xticks([0.1 1 10 100]) 
                    ylabel('Cell number')
                    title(strcat('ipsi mEC ', zonedef{kk}))   

                    subplot(2,5,5*(kk-1)+3);
                    imagesc(freq(10:80),1:1:size(ISICA3,1),zscore(ISICA3(idxsort_base,10:80),[],2),'AlphaData', .8)
                    caxis([-0.5 3])
                    set(gca, 'XScale', 'log')
                    set(gca, 'XDir', 'reverse')
                    hold on
                    line([nanmean(LFPCA3) nanmean(LFPCA3)],[1 size(ISICA3,1)],'Color','r','LineWidth',1)
                    xlabel('Frequency')
                    xticks([0.1 1 10 100]) 
                    ylabel('Cell number')
                    title(strcat('CA3 ', zonedef{kk}))                       
                    
                    subplot(2,5,5*(kk-1)+4);
                    imagesc(freq(10:80),1:1:size(ISIBoth,1),zscore(ISIBoth(idxsort_base,10:80),[],2),'AlphaData', .8)
                    caxis([-0.5 3])
                    set(gca, 'XScale', 'log')
                    set(gca, 'XDir', 'reverse')
                    hold on
                    line([nanmean(LFPBoth) nanmean(LFPBoth)],[1 size(ISIBoth,1)],'Color','r','LineWidth',1)
                    xlabel('Frequency')
                    xticks([0.1 1 10 100]) 
                    ylabel('Cell number')
                    title(strcat('Both ', zonedef{kk}))                                     

                    subplot(2,5,5*(kk-1)+5);
                    plot(freq(10:80),nanmean(ISIbase(:,10:80),1),'Color',colMat(1,:),'LineWidth',1.5)
                    hold on
                    plot(freq(10:80),nanmean(ISImEC(:,10:80),1),'Color',colMat(2,:),'LineWidth',1.5)
                    plot(freq(10:80),nanmean(ISICA3(:,10:80),1),'Color',colMat(3,:),'LineWidth',1.5)
                    plot(freq(10:80),nanmean(ISIBoth(:,10:80),1),'Color',colMat(4,:),'LineWidth',1.5)                    
                    set(gca, 'XScale', 'log')
                    set(gca, 'XDir', 'reverse')
                    ylabel('Avg distribution')
                    xlabel('Frequency')
                    xticks([0.1 1 10 100]) 
                    line([nanmean(LFPbase) nanmean(LFPbase)],[0 max(nanmean(ISIbase(:,10:80),1))],'Color',colMat(1,:),'LineWidth',1)
                    line([nanmean(LFPmEC) nanmean(LFPmEC)],[0 max(nanmean(ISIbase(:,10:80),1))],'Color',colMat(2,:),'LineWidth',1)
                    line([nanmean(LFPCA3) nanmean(LFPCA3)],[0 max(nanmean(ISIbase(:,10:80),1))],'Color',colMat(3,:),'LineWidth',1)
                    line([nanmean(LFPBoth) nanmean(LFPBoth)],[0 max(nanmean(ISIbase(:,10:80),1))],'Color',colMat(4,:),'LineWidth',1)                    
                end
            end                    
        end                
    end            

    saveas(figure(1),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryPyr.png'),'png');
    saveas(figure(1),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryPyr.fig'),'fig');
    saveas(figure(1),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryPyr.eps'),'epsc');  
    saveas(figure(2),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryIN.png'),'png');
    saveas(figure(2),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryIN.fig'),'fig');
    saveas(figure(2),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryIN.eps'),'epsc');      
    
else
    colMat = [85/243 85/243 85/243;...       
        224/243 163/243 46/243;...    
        8/243 133/243 161/243;...           
        56/243 61/243 150/243];  
    
    figure(1)
    set(gcf,'Position',[1968 213 1575 717])
    set(gcf,'renderer','painters');    
    figure(2)
    set(gcf,'Position',[1968 213 1575 717])
    set(gcf,'renderer','painters');        

    for ii = 1:length(reg)
        for jj = 1%:length(target)

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
                    figure(cc)
                    for kk = [1 2]     
                        
                        ISIbase = spectrum{ii,jj}{kk}(idx{cc},:)./sum(spectrum{ii,jj}{kk}(idx{cc},:),2);
                        ISIstim = spectrum{ii,jj}{kk+3}(idx{cc},:)./sum(spectrum{ii,jj}{kk+3}(idx{cc},:),2);
                        LFPbase = LFPfreq{ii,jj}{kk}(idx{cc});
                        LFPstim = LFPfreq{ii,jj}{kk+3}(idx{cc});
                        
                        [~,idxsort_base] = sort(nanmean(ISIbase(:,10:30),2),'descend');
                        colMap = cbrewer('seq','Blues',20);
                        colormap(colMap)
                        
                        subplot(3,6,6*(ii-1)+2*(kk-1)+1);
                        imagesc(freq(10:80),1:1:size(ISIbase,1),zscore(ISIbase(idxsort_base,10:80),[],2),'AlphaData', .8)
                        caxis([-0.5 3])
                        set(gca, 'XScale', 'log')
                        set(gca, 'XDir', 'reverse')
                        hold on
                        line([nanmean(LFPbase) nanmean(LFPbase)],[1 size(ISIbase,1)],'Color','r','LineWidth',1)
                        xlabel('Frequency')
                        xticks([0.1 1 10 100]) 
                        ylabel('Cell number')
                        title(strcat('No stim ', zonedef{kk}))   
                        
                        subplot(3,6,6*(ii-1)+2*(kk-1)+2);
                        imagesc(freq(10:80),1:1:size(ISIstim,1),zscore(ISIbase(idxsort_base,10:80),[],2),'AlphaData', .8)
                        caxis([-0.5 3])
                        set(gca, 'XScale', 'log')
                        set(gca, 'XDir', 'reverse')
                        hold on
                        line([nanmean(LFPstim) nanmean(LFPstim)],[1 size(ISIstim,1)],'Color','r','LineWidth',1)
                        xlabel('Frequency')
                        xticks([0.1 1 10 100]) 
                        ylabel('Cell number')
                        title(strcat(reg{ii},' ', zonedef{kk}))            
                        
                        subplot(3,6,6*(ii-1)+4+kk);
                        plot(freq(10:80),nanmean(ISIbase(:,10:80),1),'Color',colMat(1,:),'LineWidth',1.5)
                        hold on
                        plot(freq(10:80),nanmean(ISIstim(:,10:80),1),'Color',colMat(ii+1,:),'LineWidth',1.5)
                        set(gca, 'XScale', 'log')
                        set(gca, 'XDir', 'reverse')
                        ylabel('Avg distribution')
                        xlabel('Frequency')
                        xticks([0.1 1 10 100]) 
                        line([nanmean(LFPstim) nanmean(LFPstim)],[0 max(nanmean(ISIbase(:,10:80),1))],'Color',colMat(ii+1,:),'LineWidth',1)
                        line([nanmean(LFPbase) nanmean(LFPbase)],[0 max(nanmean(ISIbase(:,10:80),1))],'Color',colMat(1,:),'LineWidth',1)
                    end
                end                    
            end                
        end            
    end   
    
    saveas(figure(1),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryPyr.png'),'png');
    saveas(figure(1),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryPyr.fig'),'fig');
    saveas(figure(1),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryPyr.eps'),'epsc');  
    saveas(figure(2),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryIN.png'),'png');
    saveas(figure(2),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryIN.fig'),'fig');
    saveas(figure(2),strcat(parentDir,'Compiled\ISI\',tag,'ISISummaryIN.eps'),'epsc');      
end    

end

