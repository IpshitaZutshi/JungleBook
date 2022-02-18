% Compile data across all sessions

function compileMicePhasePrecession(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
addParameter(p,'numBins',100,@isnumeric); %Number of bins in the rateMap
addParameter(p,'savePlot',true,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;
numBins = p.Results.numBins;
savePlot = p.Results.savePlot;

tag = 'CA3'; % or mEC

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ21\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ11\Final','IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final'...
        'IZ21\Final','IZ27\Final','IZ28\Final','IZ29\Final','IZ30\Final','IZ31\Final'};  % To add: IZ23, IZ24, IZ25, IZ26
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ34\Final','IZ33\Final','IZ32\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Both'};
end

target = {'STEM', 'RETURN'};

if ~isempty(analogEv)
    for ii = 1:numAnalog
        analogCh(ii) = (analogEv-1)+ii;
    end
end

for rr = 1:3
    for cc = 1:2
        for zz = 1:4 %4 fields, left trajectories, right trajectories, stim left and right trajectories (correct)
            PPStats.slope{rr,cc}{zz} = [];
            PPStats.intercept{rr,cc}{zz} = [];
            PPStats.r2{rr,cc}{zz} = [];
            PPStats.p{rr,cc}{zz} = [];                
            PPStats.boundaries{rr,cc}{zz} = [];
            PPStats.region{rr,cc}{zz} = [];
            PPStats.putativeCellType{rr,cc}{zz} = [];      
        end
    end
end

%% Loop through the mice
for m = 1:length(mice)
   
    cd(strcat(parentDir, mice{m},'\Summ'));
    
    if exist('phasePrecession.mat','file')
        load('phasePrecession.mat');
    else 
        disp(['Phase precession not computed for mouse' mice{m}])
        continue;
    end
    
    for rr = 1:length(reg)
        for cc = 1:length(target)
            for zz = 1:4  
                PPStats.slope{rr,cc}{zz} = [PPStats.slope{rr,cc}{zz};PP.slope{rr,cc}{zz}];
                PPStats.intercept{rr,cc}{zz} = [PPStats.intercept{rr,cc}{zz};PP.intercept{rr,cc}{zz}];
                PPStats.p{rr,cc}{zz} = [PPStats.p{rr,cc}{zz};PP.p{rr,cc}{zz}];
                PPStats.r2{rr,cc}{zz} = [PPStats.r2{rr,cc}{zz};PP.r2{rr,cc}{zz}];                            
                PPStats.boundaries{rr,cc}{zz} = [PPStats.boundaries{rr,cc}{zz};PP.boundaries{rr,cc}{zz}];
                PPStats.region{rr,cc}{zz} = [PPStats.region{rr,cc}{zz};PP.region{rr,cc}{zz}];
                PPStats.putativeCellType{rr,cc}{zz} = [PPStats.putativeCellType{rr,cc}{zz};PP.putativeCellType{rr,cc}{zz}];     

            end
        end
    end
end

target = {'STEM', 'RETURN'};
zone = {'LeftOFF','LeftON','RightOFF','RightON'};%,'LeftOFFIN','LeftONIN','RightOFFIN','RightONIN'};

col = [85/243 85/243 85/243;...
    224/243 163/243 46/243;...
	8/243 133/243 161/243;...
    56/243 61/243 150/243];

ret = round((numBins/175)*80);
stem = round((numBins/175)*110);

%Plot differently depending on whether plotting mEC only, or manipulations
%with separate analog channels, or CA3
if strcmp(tag,'CA1') == 1 || strcmp(tag,'mECBilateral') == 1 || strcmp(tag,'mEC') == 1
    
    figure(1)
    set(gcf,'renderer','painters');
    set(gcf,'Position',[100 100 1800 1000])    
    
    for ii = 1:3
        for jj = 1:length(target)
            for zz = 1:4   
                
                if ~isempty(PPStats.boundaries{ii,jj}{zz})
                   % Separate CA1, combine DG & CA3, select principal cells
                    idxSelect = PPStats.putativeCellType{ii,jj}{1,zz} == 2 & PPStats.region{ii,jj}{1,zz} == 1;
                    fieldloc = ((PPStats.boundaries{ii,jj}{zz}(:,2)+PPStats.boundaries{ii,jj}{zz}(:,1))./2)*100;

                    if ii == 1 && zz == 1 
                        fieldloc_sel{jj}(:,zz) = fieldloc(idxSelect);
                        slope_sel{jj}(:,zz) =  PPStats.slope{ii,jj}{1,zz}(idxSelect);
                        intercept_sel{jj}(:,zz) =  PPStats.intercept{ii,jj}{1,zz}(idxSelect);
                        p_sel{jj}(:,zz) =  PPStats.p{ii,jj}{1,zz}(idxSelect);
                    else 
                        fieldloc_sel{jj} = catpad(2,fieldloc_sel{jj},fieldloc(idxSelect));
                        slope_sel{jj} =  catpad(2,slope_sel{jj},PPStats.slope{ii,jj}{1,zz}(idxSelect));
                        intercept_sel{jj} =  catpad(2,intercept_sel{jj},PPStats.intercept{ii,jj}{1,zz}(idxSelect));
                        p_sel{jj} =  catpad(2,p_sel{jj},PPStats.p{ii,jj}{1,zz}(idxSelect));
                    end
                end
            end
        end
    end
    
    for jj = 1:length(target)
        
        slope_combined = [slope_sel{jj}(:,[1 2 5 6 9 10]);slope_sel{jj}(:,[3 4 7 8 11 12])];  
        intercept_combined = [intercept_sel{jj}(:,[1 2 5 6 9 10]);intercept_sel{jj}(:,[3 4 7 8 11 12])];  
        field_combined = [fieldloc_sel{jj}(:,[1 2 5 6 9 10]);fieldloc_sel{jj}(:,[3 4 7 8 11 12])];  
        p_combined = [p_sel{jj}(:,[1 2 5 6 9 10]);p_sel{jj}(:,[3 4 7 8 11 12])];

        binsplit = 0:4:numBins;
        
        for kk = 1:6         
            for pp = 1:(length(binsplit)-1)
                idxF = field_combined(:,kk)>=binsplit(pp) & field_combined(:,kk)<binsplit(pp+1);
                medianSlope(pp) = nanmean(slope_combined(idxF,kk));
                medianInt(pp) = nanmean(intercept_combined(idxF,kk));                 
            end
            subplot(3,6,6*(jj-1)+kk)
            if kk == 1  || kk == 3 || kk == 5
                barcol = col(1,:);
            else
                barcol = col(((kk/2)+1),:);
            end
            bar(binsplit(2:end),medianSlope,'FaceColor',barcol,'EdgeColor','none','FaceAlpha',0.8);
            xlabel('Field location')
            ylabel('Slope')
            line([ret ret],[-1 1],'Color','red','LineWidth',1.5)
            line([stem stem],[-1 1],'Color','red','LineWidth',1.5)
            ylim([-1 1])
            
            numsig = sum(p_combined(:,kk)<0.05)./sum(~isnan(p_combined(:,kk)));
            title(strcat('Prop sig = ', num2str(numsig)));                
            
            fieldret = field_combined(:,kk)>=0 & field_combined(:,kk)<=ret; 
            fieldstem = field_combined(:,kk)>=stem+10 & field_combined(:,kk)<=(numBins-15); 
            if kk == 1 
                MeanSlope(:,1) = slope_combined(fieldret,kk);
                MeanInt(:,1) = intercept_combined(fieldret,kk);
            else 
                MeanSlope = catpad(2,MeanSlope,slope_combined(fieldret,kk));
                MeanInt = catpad(2,MeanInt, intercept_combined(fieldret,kk));                    
            end   
            MeanSlope = catpad(2,MeanSlope,slope_combined(fieldstem,kk));
            MeanInt = catpad(2,MeanInt, intercept_combined(fieldstem,kk));                 
        end
        
        subplot(3,6,12+jj)
        meanR = nanmean(MeanSlope(:,[1 3 5 7 9 11]),1);
        stdR = nanstd(MeanSlope(:,[1 3 5 7 9 11]),1)./sqrt(sum(~isnan(MeanSlope(:,[1 3 5 7 9 11]))));
        errorbar(meanR,stdR,'Color',col(1,:),'LineWidth',1.5,'LineStyle','-')
        hold on
        meanS = nanmean(MeanSlope(:,[2 4 6 8 10 12]),1);
        stdS = nanstd(MeanSlope(:,[2 4 6 8 10 12]),1)./sqrt(sum(~isnan(MeanSlope(:,[2 4 6 8 10 12]))));
        errorbar(meanS,stdS,'Color',col(1,:),'LineWidth',1.5,'LineStyle','--')           
        xlim([0 7])
        ylim([-0.75 0.25])
        line([0 7],[0 0],'Color',[200/243 200/243 200/243],'LineStyle','--')
        ylabel('Slope')
        title(strcat('Slope ',target{jj}))     
        xticks(0:1:7)
        xticklabels({'','B','CA1','B','mEC','B','Both',''})

        subplot(3,6,14+jj) 
        a = MeanInt(:,[1 3 5 7 9 11]);
        a(a<0) = a(a<0)+2*pi;
        a = rad2deg(a);
        meanR = nanmean(a);
        stdR = nanstd(a)./sqrt(sum(~isnan(a)));
        errorbar(meanR,stdR,'Color',col(1,:),'LineWidth',1.5,'LineStyle','-')
        hold on
        
        a = MeanInt(:,[2 4 6 8 10 12]);
        a(a<0) = a(a<0)+2*pi; 
        a = rad2deg(a);
        meanS = nanmean(a);
        stdS = nanstd(a)./sqrt(sum(~isnan(a)));
        errorbar(meanS,stdS,'Color',col(1,:),'LineWidth',1.5,'LineStyle','--')           
        xlim([0 7])
        ylim([90 240])
        line([0 7],[0 0],'Color',[200/243 200/243 200/243],'LineStyle','--') %0 is the theta trough
        ylabel('Intercept (180=peak)')
        title(strcat('Intercept ',target{jj}))
        xticks(0:1:7)
        xticklabels({'','B','CA1','B','mEC','B','Both',''})
        
        clear MeanSlope MeanInt
    end

    if savePlot
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.png'));
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.eps'),'epsc');
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.fig'));
    end
    
elseif strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    
    col = [85/243 85/243 85/243;...
        8/243 133/243 161/243;...
        224/243 163/243 46/243;...
        56/243 61/243 150/243];

    figure(1)
    set(gcf,'renderer','painters');
    set(gcf,'Position',[100 100 1800 1000])
    for ii = 2:3
        for jj = 1:length(target)
            for zz = 1:4   
                
               % Separate CA1, combine DG & CA3, select principal cells
                idxSelect = PPStats.putativeCellType{ii,jj}{1,zz} == 2 & PPStats.region{ii,jj}{1,zz} == 1;
                fieldloc = ((PPStats.boundaries{ii,jj}{zz}(:,2)+PPStats.boundaries{ii,jj}{zz}(:,1))./2)*100;

                if ii == 2 && zz == 1 
                    fieldloc_sel{jj}(:,zz) = fieldloc(idxSelect);
                    slope_sel{jj}(:,zz) =  PPStats.slope{ii,jj}{1,zz}(idxSelect);
                    intercept_sel{jj}(:,zz) =  PPStats.intercept{ii,jj}{1,zz}(idxSelect);
                    p_sel{jj}(:,zz) =  PPStats.p{ii,jj}{1,zz}(idxSelect);
                else 
                    fieldloc_sel{jj} = catpad(2,fieldloc_sel{jj},fieldloc(idxSelect));
                    slope_sel{jj} =  catpad(2,slope_sel{jj},PPStats.slope{ii,jj}{1,zz}(idxSelect));
                    intercept_sel{jj} =  catpad(2,intercept_sel{jj},PPStats.intercept{ii,jj}{1,zz}(idxSelect));
                    p_sel{jj} =  catpad(2,p_sel{jj},PPStats.p{ii,jj}{1,zz}(idxSelect));
                end   
            end
        end
    end
        
    for jj = 1:length(target)
        
        slope_combined = [slope_sel{jj}(:,[1 2 5 6]);slope_sel{jj}(:,[3 4 7 8])];  
        intercept_combined = [intercept_sel{jj}(:,[1 2 5 6]);intercept_sel{jj}(:,[3 4 7 8])];  
        field_combined = [fieldloc_sel{jj}(:,[1 2 5 6]);fieldloc_sel{jj}(:,[3 4 7 8])];  
        p_combined = [p_sel{jj}(:,[1 2 5 6]);p_sel{jj}(:,[3 4 7 8])];

        binsplit = 0:5:numBins;
        
        for kk = 1:4         
            for pp = 1:(length(binsplit)-1)
                idxF = field_combined(:,kk)>=binsplit(pp) & field_combined(:,kk)<binsplit(pp+1);
                medianSlope(pp) = nanmean(slope_combined(idxF,kk));
                medianInt(pp) = nanmean(intercept_combined(idxF,kk));                 
            end
            subplot(3,4,4*(jj-1)+kk)
            bar(binsplit(2:end),medianSlope,'FaceColor',col(kk,:),'EdgeColor','none','FaceAlpha',0.8);
            xlabel('Field location')
            ylabel('Slope')
            line([ret ret],[-1 1],'Color','red','LineWidth',1.5)
            line([stem stem],[-1 1],'Color','red','LineWidth',1.5)
            numsig = sum(p_combined(:,kk)<0.05)./sum(~isnan(p_combined(:,kk)));
            title(strcat('Prop sig = ', num2str(numsig)));                  
            
            fieldret = field_combined(:,kk)>=0 & field_combined(:,kk)<=ret & p_combined(:,kk)<=0.05; 
            fieldstem = field_combined(:,kk)>=stem & field_combined(:,kk)<=(numBins-15) & p_combined(:,kk)<=0.05; 
            if kk == 1 
                MeanSlope(:,1) = slope_combined(fieldret,kk);
                MeanInt(:,1) = intercept_combined(fieldret,kk);
            else 
                MeanSlope = catpad(2,MeanSlope,slope_combined(fieldret,kk));
                MeanInt = catpad(2,MeanInt, intercept_combined(fieldret,kk));                    
            end   
            MeanSlope = catpad(2,MeanSlope,slope_combined(fieldstem,kk));
            MeanInt = catpad(2,MeanInt, intercept_combined(fieldstem,kk));                 
        end
        
        subplot(3,4,8+jj)
        meanR = nanmean(MeanSlope(:,[1 3 5 7]),1);
        stdR = nanstd(MeanSlope(:,[1 3 5 7]),1)./sqrt(sum(~isnan(MeanSlope(:,[1 3 5 7]))));
        errorbar(meanR,stdR,'Color',col(1,:),'LineWidth',1.5,'LineStyle','-')
        hold on
        meanS = nanmean(MeanSlope(:,[2 4 6 8]),1);
        stdS = nanstd(MeanSlope(:,[2 4 6 8]),1)./sqrt(sum(~isnan(MeanSlope(:,[2 4 6 8]))));
        errorbar(meanS,stdS,'Color',col(1,:),'LineWidth',1.5,'LineStyle','--')           
        xlim([0 5])
        ylim([-0.75 0.25])
        line([0 5],[0 0],'Color',[200/243 200/243 200/243],'LineStyle','--')
        ylabel('Slope')
        title(strcat('Slope ',target{jj}))            
        xticklabels({'','B','mEC','CA3','Both',''})

        subplot(3,4,10+jj) 
        
        a = MeanInt(:,[1 3 5 7]);
        a(a<0) = a(a<0)+2*pi;
        a = rad2deg(a);
        meanR = nanmean(a);
        stdR = nanstd(a)./sqrt(sum(~isnan(a)));
        errorbar(meanR,stdR,'Color',col(1,:),'LineWidth',1.5,'LineStyle','-')
        hold on
        
        a = MeanInt(:,[2 4 6 8]);
        a(a<0) = a(a<0)+2*pi;
        a = rad2deg(a);
        meanS = nanmean(a);
        stdS = nanstd(a)./sqrt(sum(~isnan(a)));
        errorbar(meanS,stdS,'Color',col(1,:),'LineWidth',1.5,'LineStyle','--')           
        xlim([0 5])
        ylim([90 240])
        line([0 5],[0 0],'Color',[200/243 200/243 200/243],'LineStyle','--')
        ylabel('Intercept (180=peak)')
        title(strcat('Intercept ',target{jj}))
        xticklabels({'','B','mEC','CA3','Both',''})
        
        clear MeanSlope MeanInt
    end
    
    if savePlot
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.png'));
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.eps'),'epsc');
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.fig'));
    end
end
end
