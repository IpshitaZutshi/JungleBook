% Compile data across all sessions

function compileMicePhasePrecessionVer3(varargin)

%% Defaults and Parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'analogEv',64,@isnumeric);
addParameter(p,'numAnalog',2,@isnumeric);
addParameter(p,'sigVal',0.05,@isnumeric);
addParameter(p,'numBins',100,@isnumeric); %Number of bins in the rateMap
addParameter(p,'savePlot',false,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
analogEv = p.Results.analogEv;
numAnalog = p.Results.numAnalog;
sigVal = p.Results.sigVal;
numBins = p.Results.numBins;
savePlot = p.Results.savePlot;

tag = 'CA3'; % or mEC

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final'...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline',...
        'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Final'};  % To add: IZ34,'IZ29\Saline'
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ32\Saline','IZ33\Saline','IZ34\Saline','IZ27\Saline','IZ28\Saline','IZ29\Saline'};
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
            PPStats.mice{rr,cc}{zz} = [];   
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
                PPStats.mice{rr,cc}{zz} = [PPStats.mice{rr,cc}{zz};ones(length(PP.p{rr,cc}{zz}),1)*m]; 
            end
        end
    end
end

target = {'STEM', 'RETURN'};
zone = {'LeftOFF','LeftON','RightOFF','RightON'};%,'LeftOFFIN','LeftONIN','RightOFFIN','RightONIN'};

colMat = [85/243 85/243 85/243;...
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
        for jj = 1%:length(target)
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
                        miceId{jj}(:,zz) = PPStats.mice{ii,jj}{1,zz}(idxSelect);
                    else 
                        fieldloc_sel{jj} = catpad(2,fieldloc_sel{jj},fieldloc(idxSelect));
                        slope_sel{jj} =  catpad(2,slope_sel{jj},PPStats.slope{ii,jj}{1,zz}(idxSelect));
                        intercept_sel{jj} =  catpad(2,intercept_sel{jj},PPStats.intercept{ii,jj}{1,zz}(idxSelect));
                        p_sel{jj} =  catpad(2,p_sel{jj},PPStats.p{ii,jj}{1,zz}(idxSelect));
                        miceId{jj} = catpad(2,miceId{jj},PPStats.mice{ii,jj}{1,zz}(idxSelect));
                    end
                end
            end
        end
    end
    
    for jj = 1%:length(target)
        
        slope_combined = [slope_sel{jj}(:,[1 2 5 6 9 10]);slope_sel{jj}(:,[3 4 7 8 11 12])];  
        intercept_combined = [intercept_sel{jj}(:,[1 2 5 6 9 10]);intercept_sel{jj}(:,[3 4 7 8 11 12])];  
        field_combined = [fieldloc_sel{jj}(:,[1 2 5 6 9 10]);fieldloc_sel{jj}(:,[3 4 7 8 11 12])];  
        p_combined = [p_sel{jj}(:,[1 2 5 6 9 10]);p_sel{jj}(:,[3 4 7 8 11 12])];
        mice_combined = [miceId{jj}(:,[1 2 5 6 9 10]);miceId{jj}(:,[3 4 7 8 11 12])];
        
        for kk = 1:3       
           kkIdx = 2*(kk-1)+1;
           %idxtoKeepAll = (field_combined(:,kkIdx) >63 | field_combined(:,kkIdx+1) >63) & (field_combined(:,kkIdx) <95 | field_combined(:,kkIdx+1) < 95);
           idxtoKeepAll = (field_combined(:,kkIdx:kkIdx+1) >63) & (field_combined(:,kkIdx:kkIdx+1) <95);
           idxtoKeepSig = field_combined(:,kkIdx:kkIdx+1) >63 & field_combined(:,kkIdx:kkIdx+1) < 95 & p_combined(:,kkIdx:kkIdx+1) < sigVal;
           
           slopesBaseAll = slope_combined(idxtoKeepAll(:,1),kkIdx);
           slopesStimAll = slope_combined(idxtoKeepAll(:,2),kkIdx+1);          
           interceptBaseAll = intercept_combined(idxtoKeepAll(:,1),kkIdx)+pi;
           interceptStimAll = intercept_combined(idxtoKeepAll(:,2),kkIdx+1)+pi;
           
           miceBase = mice_combined(idxtoKeepSig(:,1),kkIdx);
           miceStim = mice_combined(idxtoKeepSig(:,2),kkIdx+1); 
           
           slopesBase = slope_combined(idxtoKeepSig(:,1),kkIdx);
           slopesStim = slope_combined(idxtoKeepSig(:,2),kkIdx+1);          
           interceptBase = intercept_combined(idxtoKeepSig(:,1),kkIdx)+pi;
           interceptStim = intercept_combined(idxtoKeepSig(:,2),kkIdx+1)+pi;
           
           for mm = unique(miceBase)'
               idxBase = miceBase == mm;
               avgSlopeBase(mm,1) = nanmean(slopesBase(idxBase));

               idxStim = miceStim == mm;
               avgSlopeStim(mm,1) = nanmean(slopesStim(idxStim));               
           end
           slopRat = abs(avgSlopeStim-1)./abs(avgSlopeBase-1);
           
           subplot(3,8,8*(kk-1)+1)
           [N,edges] = histcounts(slopesBaseAll,[-2:0.2:2]);
           h1 = bar(edges(1:end-1),N); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.4;
           hold on
           [N1,edges] = histcounts(slopesBase,[-2:0.2:2]);
           h1 = bar(edges(1:end-1),N1); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.9;  
           line([nanmedian(slopesBaseAll) nanmedian(slopesBaseAll)],[0 max(N)],'Color','k','LineWidth',1.5)
           line([nanmedian(slopesBase) nanmedian(slopesBase)],[0 max(N)],'Color','r','LineWidth',1.5)
           sigPercent = (sum(idxtoKeepSig(:,1))./sum(idxtoKeepAll(:,1))*100);
           title(strcat('Base ',reg{kk},' %Sig: ',num2str(sigPercent)))
           xlabel('Precession slope')
           ylabel('Cell count')

           subplot(3,8,8*(kk-1)+2)
           [N,edges] = histcounts(slopesStimAll,[-2:0.2:2]);
           h1 = bar(edges(1:end-1),N); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.4;
           hold on
           [N1,edges] = histcounts(slopesStim,[-2:0.2:2]);
           h1 = bar(edges(1:end-1),N1); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.9;  
           line([nanmedian(slopesStimAll) nanmedian(slopesStimAll)],[0 max(N)],'Color','k','LineWidth',1.5)
           line([nanmedian(slopesStim) nanmedian(slopesStim)],[0 max(N)],'Color','r','LineWidth',1.5)      
           sigPercent = (sum(idxtoKeepSig(:,2))./sum(idxtoKeepAll(:,2))*100);
           title(strcat('Stim ',reg{kk},' %Sig: ',num2str(sigPercent)))
           xlabel('Precession slope')
           ylabel('Cell count')
           
           subplot(3,8,8*(kk-1)+3)
           stats{kk}.Slopesall = groupStats([slopesBaseAll;slopesStimAll],[ones(length(slopesBaseAll),1);ones(length(slopesStimAll),1)*2],'inAxis',true,'Color',[colMat(1,:);colMat(kk+1,:)]);
           title('All slopes')
           ylim([-1.5 1.5])
           
           subplot(3,8,8*(kk-1)+4)
           stats{kk}.SlopesSig = groupStats([slopesBase;slopesStim],[ones(length(slopesBase),1);ones(length(slopesStim),1)*2],'inAxis',true,'Color',[colMat(1,:);colMat(kk+1,:)]);
           [stats{kk}.SlopesSig.ranksum.p, ~,stats{kk}.SlopesSig.ranksum.stats] = ranksum(slopesBase,slopesStim);
           [stats{kk}.SlopesSig.signrank.p(1), ~,stats{kk}.SlopesSig.signrank.stats{1}] = signrank(slopesBase);
           [stats{kk}.SlopesSig.signrank.p(2), ~,stats{kk}.SlopesSig.signrank.stats{2}] = signrank(slopesStim);
           title('Significant slopes')
           ylim([-1.5 0.5])
           
           subplot(3,8,8*(kk-1)+5)
           c1 = polarhistogram(interceptBaseAll,25,'Normalization','probability');
           c1.DisplayStyle = 'stairs';   hold on
           meanAngle = circ_mean(interceptBaseAll);
           polarplot([meanAngle meanAngle],[0 .15])  
           title('No stim all')
           [stats{kk}.InterceptAll.p stats{kk}.InterceptAll.table] = circ_wwtest([interceptBaseAll;interceptStimAll],[ones(length(interceptBaseAll),1);ones(length(interceptStimAll),1)*2]);
            
           subplot(3,8,8*(kk-1)+6)
           c1 = polarhistogram(interceptStimAll,25,'Normalization','probability');
           c1.DisplayStyle = 'stairs';   hold on
           meanAngle = circ_mean(interceptStimAll);
           polarplot([meanAngle meanAngle],[0 .15])  
           title(strcat('Stim all ,p:',num2str(stats{kk}.InterceptAll.p)))

           subplot(3,8,8*(kk-1)+7)
           c1 = polarhistogram(interceptBase,25,'Normalization','probability');
           c1.DisplayStyle = 'stairs';   hold on
           meanAngle = circ_mean(interceptBase);
           polarplot([meanAngle meanAngle],[0 .15])  
           title('No stim sig')
           [stats{kk}.InterceptSig.p stats{kk}.InterceptSig.table] = circ_wwtest([interceptBase;interceptStim],[ones(length(interceptBase),1);ones(length(interceptStim),1)*2]);
            
           subplot(3,8,8*(kk-1)+8)
           c1 = polarhistogram(interceptStim,25,'Normalization','probability');
           c1.DisplayStyle = 'stairs';   hold on
           meanAngle = circ_mean(interceptStim);
           polarplot([meanAngle meanAngle],[0 .15])  
           title(strcat('Stim sig ,p:',num2str(stats{kk}.InterceptSig.p)))           
        end       
    end

    if savePlot
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession2.png'));
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession2.eps'),'epsc');
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession2.fig'));
        save(strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecessionStats2.mat'),'stats')
    end
    
elseif strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    
    colMat = [85/243 85/243 85/243;...
        8/243 133/243 161/243;...
        224/243 163/243 46/243;...
        56/243 61/243 150/243];

    figure(1)
    set(gcf,'renderer','painters');
    set(gcf,'Position',[100 100 1800 1000])
    for ii = 2:3
        for jj = 1%:length(target)
            for zz = 1:4   
                
               % Separate CA1, combine DG & CA3, select principal cells
                idxSelect = PPStats.putativeCellType{ii,jj}{1,zz} == 2 & PPStats.region{ii,jj}{1,zz} == 1;
                fieldloc = ((PPStats.boundaries{ii,jj}{zz}(:,2)+PPStats.boundaries{ii,jj}{zz}(:,1))./2)*100;

                if ii == 2 && zz == 1 
                    fieldloc_sel{jj}(:,zz) = fieldloc(idxSelect);
                    slope_sel{jj}(:,zz) =  PPStats.slope{ii,jj}{1,zz}(idxSelect);
                    intercept_sel{jj}(:,zz) =  PPStats.intercept{ii,jj}{1,zz}(idxSelect);
                    p_sel{jj}(:,zz) =  PPStats.p{ii,jj}{1,zz}(idxSelect);
                    r_sel{jj}(:,zz) =  PPStats.r2{ii,jj}{1,zz}(idxSelect);
                    miceId{jj}(:,zz) = PPStats.mice{ii,jj}{1,zz}(idxSelect);
                else 
                    fieldloc_sel{jj} = catpad(2,fieldloc_sel{jj},fieldloc(idxSelect));
                    slope_sel{jj} =  catpad(2,slope_sel{jj},PPStats.slope{ii,jj}{1,zz}(idxSelect));
                    intercept_sel{jj} =  catpad(2,intercept_sel{jj},PPStats.intercept{ii,jj}{1,zz}(idxSelect));
                    p_sel{jj} =  catpad(2,p_sel{jj},PPStats.p{ii,jj}{1,zz}(idxSelect));
                    r_sel{jj} = catpad(2,r_sel{jj},PPStats.r2{ii,jj}{1,zz}(idxSelect));
                    miceId{jj} = catpad(2,miceId{jj},PPStats.mice{ii,jj}{1,zz}(idxSelect));
                end   
            end
        end
    end
        
    for jj = 1%:length(target)
        
       slope_combined = [slope_sel{jj}(:,[1 2 5 6]);slope_sel{jj}(:,[3 4 7 8])];  
       intercept_combined = [intercept_sel{jj}(:,[1 2 5 6]);intercept_sel{jj}(:,[3 4 7 8])];  
       field_combined = [fieldloc_sel{jj}(:,[1 2 5 6]);fieldloc_sel{jj}(:,[3 4 7 8])];  
       p_combined = [p_sel{jj}(:,[1 2 5 6]);p_sel{jj}(:,[3 4 7 8])];
       r_combined = [r_sel{jj}(:,[1 2 5 6]);r_sel{jj}(:,[3 4 7 8])];
       mice_combined = [miceId{jj}(:,[1 2 5 6]);miceId{jj}(:,[3 4 7 8])];
        
       idxtoKeepAll = field_combined(:,:) >63 & field_combined(:,:) < 95;% & slope_combined(:,:)<=1;
       idxtoKeepSig = field_combined(:,:) >63 & field_combined(:,:) < 95 & p_combined(:,:) <sigVal;% & slope_combined(:,:)<=1;
       
       slopesBaseAll = slope_combined(idxtoKeepAll(:,1),1);
       slopesmECAll = slope_combined(idxtoKeepAll(:,2),2);     
       slopesCA3All = slope_combined(idxtoKeepAll(:,3),3);
       slopesBothAll = slope_combined(idxtoKeepAll(:,4),4);     
       
       interceptBaseAll = intercept_combined(idxtoKeepAll(:,1),1)+pi;
       interceptmECAll = intercept_combined(idxtoKeepAll(:,2),2)+pi;
       interceptCA3All = intercept_combined(idxtoKeepAll(:,3),3)+pi;
       interceptBothAll = intercept_combined(idxtoKeepAll(:,4),4)+pi;
 
       slopesBase = slope_combined(idxtoKeepSig(:,1),1);
       slopesmEC = slope_combined(idxtoKeepSig(:,2),2);     
       slopesCA3 = slope_combined(idxtoKeepSig(:,3),3);
       slopesBoth = slope_combined(idxtoKeepSig(:,4),4);     
       
       interceptBase = intercept_combined(idxtoKeepSig(:,1),1)+pi;
       interceptmEC = intercept_combined(idxtoKeepSig(:,2),2)+pi;
       interceptCA3 = intercept_combined(idxtoKeepSig(:,3),3)+pi;
       interceptBoth = intercept_combined(idxtoKeepSig(:,4),4)+pi;
           
       miceBase = mice_combined(idxtoKeepSig(:,1),1);
       micemEC = mice_combined(idxtoKeepSig(:,2),2); 
       miceCA3 = mice_combined(idxtoKeepSig(:,3),3);
       miceBoth = mice_combined(idxtoKeepSig(:,4),4);  
       
       for mm = unique(miceBase)'
           idxBase = miceBase == mm;
           avgSlopeBase(mm,1) = nanmedian(slopesBaseAll(idxBase));
           idxmEC = micemEC == mm;
           avgSlopemEC(mm,1) = nanmedian(slopesmECAll(idxmEC));                       
           idxCA3 = miceCA3 == mm;
           avgSlopeCA3(mm,1) = nanmedian(slopesCA3All(idxCA3));                          
           idxStim = miceBoth == mm;
           avgSlopeBoth(mm,1) = nanmedian(slopesBothAll(idxStim));           
       end
       a = [avgSlopeBase avgSlopemEC avgSlopeCA3 avgSlopeBoth];
%        slopRatCA3 = abs(avgSlopeCA3-1)./abs(avgSlopeBase-1);
%        slopRatBoth = abs(avgSlopeBoth-1)./abs(avgSlopeBase-1);       
       
       subplot(3,6,1)
       [N,edges] = histcounts(slopesBaseAll,[-2:0.2:2]);
       h1 = bar(edges(1:end-1),N); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.4;
       hold on
       [N1,edges] = histcounts(slopesBase,[-2:0.2:2]);
       h1 = bar(edges(1:end-1),N1); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.9;  
       line([nanmedian(slopesBaseAll) nanmedian(slopesBaseAll)],[0 max(N)],'Color','k','LineWidth',1.5)
       line([nanmedian(slopesBase) nanmedian(slopesBase)],[0 max(N)],'Color','r','LineWidth',1.5)
       sigPercent = (sum(idxtoKeepSig(:,1))./sum(idxtoKeepAll(:,1)))*100;
       title(strcat('Base %Sig: ',num2str(sigPercent)))
       xlabel('Precession slope')
       ylabel('Cell count')

       subplot(3,6,2)
       [N,edges] = histcounts(slopesmECAll,[-2:0.2:2]);
       h1 = bar(edges(1:end-1),N); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.4;
       hold on
       [N1,edges] = histcounts(slopesmEC,[-2:0.2:2]);
       h1 = bar(edges(1:end-1),N1); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.9;  
       line([nanmedian(slopesmECAll) nanmedian(slopesmECAll)],[0 max(N)],'Color','k','LineWidth',1.5)
       line([nanmedian(slopesmEC) nanmedian(slopesmEC)],[0 max(N)],'Color','r','LineWidth',1.5)      
       sigPercent = (sum(idxtoKeepSig(:,2))./sum(idxtoKeepAll(:,2)))*100;
       title(strcat('mEC %Sig: ',num2str(sigPercent)))
       xlabel('Precession slope')
       ylabel('Cell count')
       
       subplot(3,6,3)
       [N,edges] = histcounts(slopesCA3All,[-2:0.2:2]);
       h1 = bar(edges(1:end-1),N); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.4;
       hold on
       [N1,edges] = histcounts(slopesCA3,[-2:0.2:2]);
       h1 = bar(edges(1:end-1),N1); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.9;  
       line([nanmedian(slopesCA3All) nanmedian(slopesCA3All)],[0 max(N)],'Color','k','LineWidth',1.5)
       line([nanmedian(slopesCA3) nanmedian(slopesCA3)],[0 max(N)],'Color','r','LineWidth',1.5)      
       sigPercent = (sum(idxtoKeepSig(:,3))./sum(idxtoKeepAll(:,3)))*100;
       title(strcat('CA3 %Sig: ',num2str(sigPercent)))
       xlabel('Precession slope')
       ylabel('Cell count')

       subplot(3,6,4)
       [N,edges] = histcounts(slopesBothAll,[-2:0.2:2]);
       h1 = bar(edges(1:end-1),N); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.4;
       hold on
       [N1,edges] = histcounts(slopesBoth,[-2:0.2:2]);
       h1 = bar(edges(1:end-1),N1); h1.FaceColor = [0.5 0.5 0.5]; h1.FaceAlpha = 0.9;  
       line([nanmedian(slopesBothAll) nanmedian(slopesBothAll)],[0 max(N)],'Color','k','LineWidth',1.5)
       line([nanmedian(slopesBoth) nanmedian(slopesBoth)],[0 max(N)],'Color','r','LineWidth',1.5)      
       sigPercent = (sum(idxtoKeepSig(:,4))./sum(idxtoKeepAll(:,4)))*100;
       title(strcat('Both %Sig: ',num2str(sigPercent)))
       xlabel('Precession slope')
       ylabel('Cell count')
       

       subplot(3,6,5)
       stats.Slopesall = groupStats([slopesBaseAll;slopesmECAll;slopesCA3All;slopesBothAll],...
           [ones(length(slopesBaseAll),1);ones(length(slopesmECAll),1)*2;ones(length(slopesCA3All),1)*3;ones(length(slopesBothAll),1)*4],'inAxis',true,'Color',colMat);
       [stats.Slopesall.signrank.p(1),~,stats.Slopesall.signrank.stats{1}] = signrank(slopesBaseAll);
       [stats.Slopesall.signrank.p(2),~,stats.Slopesall.signrank.stats{2}] = signrank(slopesmECAll);
       [stats.Slopesall.signrank.p(3),~,stats.Slopesall.signrank.stats{3}] = signrank(slopesCA3All);
       [stats.Slopesall.signrank.p(4),~,stats.Slopesall.signrank.stats{4}] = signrank(slopesBothAll);       
       title(strcat(num2str(stats.Slopesall.signrank.p(1:4))))

       subplot(3,6,6)
       stats.SlopesSig = groupStats([slopesBase;slopesmEC;slopesCA3;slopesBoth],...
          [ones(length(slopesBase),1);ones(length(slopesmEC),1)*2;ones(length(slopesCA3),1)*3;ones(length(slopesBoth),1)*4],'inAxis',true,'Color',colMat);
%        stats.SlopesSig = groupStats([slopesBase;slopesCA3;slopesBoth],...
%            [ones(length(slopesBase),1);ones(length(slopesCA3),1)*3;ones(length(slopesBoth),1)*4],'inAxis',true,'Color',colMat);       
%        stats.SlopesSig = groupStats([slopesBase;slopesBoth],...
%            [ones(length(slopesBase),1);ones(length(slopesBoth),1)*4],'inAxis',true,'Color',colMat); 
       [stats.SlopesSig.signrank.p(1),~,stats.SlopesSig.signrank.stats{1}] = signrank(slopesBase);
       [stats.SlopesSig.signrank.p(2),~,stats.SlopesSig.signrank.stats{2}] = signrank(slopesmEC);
       [stats.SlopesSig.signrank.p(3),~,stats.SlopesSig.signrank.stats{3}] = signrank(slopesCA3);
       [stats.SlopesSig.signrank.p(4),~,stats.SlopesSig.signrank.stats{4}] = signrank(slopesBoth);
              
       title(strcat(num2str(stats.SlopesSig.signrank.p(1:4))))
       %title('Significant slopes')    
       
       %ylim([-2 2])

       subplot(3,6,7)
       c1 = polarhistogram(interceptBaseAll,25,'Normalization','probability');
       c1.DisplayStyle = 'stairs';   hold on
       meanAngle = circ_mean(interceptBaseAll);
       polarplot([meanAngle meanAngle],[0 .15])  
       title('No stim all')
       [stats.InterceptAll.p stats.InterceptAll.table] = circ_wwtest([interceptBaseAll;interceptmECAll;interceptCA3All;interceptBothAll],...
           [ones(length(interceptBaseAll),1);ones(length(interceptmECAll),1)*2;ones(length(interceptCA3All),1)*3;ones(length(interceptBothAll),1)*4]);

       subplot(3,6,8)
       c1 = polarhistogram(interceptmECAll,25,'Normalization','probability');
       c1.DisplayStyle = 'stairs';   hold on
       meanAngle = circ_mean(interceptmECAll);
       polarplot([meanAngle meanAngle],[0 .15])  
       title(strcat('mEC all ,p:',num2str(stats.InterceptAll.p)))

       subplot(3,6,9)
       c1 = polarhistogram(interceptCA3All,25,'Normalization','probability');
       c1.DisplayStyle = 'stairs';   hold on
       meanAngle = circ_mean(interceptCA3All);
       polarplot([meanAngle meanAngle],[0 .15])  
       title(strcat('CA3 all'))
       
       subplot(3,6,10)
       c1 = polarhistogram(interceptBothAll,25,'Normalization','probability');
       c1.DisplayStyle = 'stairs';   hold on
       meanAngle = circ_mean(interceptBothAll);
       polarplot([meanAngle meanAngle],[0 .15])  
       title(strcat('Stim both'))
       
       subplot(3,6,13)
       c1 = polarhistogram(interceptBase,25,'Normalization','probability');
       c1.DisplayStyle = 'stairs';   hold on
       meanAngle = circ_mean(interceptBase);
       polarplot([meanAngle meanAngle],[0 .15])  
       title('No stim sig')
       [stats.InterceptSig.p, stats.InterceptSig.table] = circ_wwtest([interceptBase;interceptmEC;interceptCA3;interceptBoth],...
           [ones(length(interceptBase),1);ones(length(interceptmEC),1)*2;ones(length(interceptCA3),1)*3;ones(length(interceptBoth),1)*4]);

       subplot(3,6,14)
       c1 = polarhistogram(interceptmEC,25,'Normalization','probability');
       c1.DisplayStyle = 'stairs';   hold on
       meanAngle = circ_mean(interceptmEC);
       polarplot([meanAngle meanAngle],[0 .15])  
       title(strcat('mEC sig ,p:',num2str(stats.InterceptSig.p)))

       subplot(3,6,15)
       c1 = polarhistogram(interceptCA3,25,'Normalization','probability');
       c1.DisplayStyle = 'stairs';   hold on
       meanAngle = circ_mean(interceptCA3);
       polarplot([meanAngle meanAngle],[0 .15])  
       title(strcat('CA3 sig'))
       
       subplot(3,6,16)
       c1 = polarhistogram(interceptBoth,25,'Normalization','probability');
       c1.DisplayStyle = 'stairs';   hold on
       meanAngle = circ_mean(interceptBoth);
       polarplot([meanAngle meanAngle],[0 .15])  
       title(strcat('Stim sig'))
     
    end

    if savePlot
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.png'));
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.eps'),'epsc');
        saveas(figure(1),strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecession.fig'));
        save(strcat(parentDir,'Compiled\Phase Precession\',tag,'\PhasePrecessionStats.mat'),'stats')
    end
end
end
