%% Plot phase precession distribution histograms
function data = plotPhasePrecession_hist(data)

%subsets = {'Repeated'}; %'Unique'
saveType = 'PDF';
folderpath = 'W:\Ipshita\Recordings\Wfs1 Project\E32';
%folderpath = 'Z:\Ipshita\Recordings\E41';
%cellTypes = {'CA1' 'DG' 'IN'}; %'HeadDirection'
cellTypes = {'MEC'};
methods = {'Pooled'}; % 'Single' 'Averaged'
[colors] = defineCellTypeColors_2015Poster;
load('C:\Users\Ipshita\Documents\PhD documents\Wfs1 project compiled data\E32\E32Data.mat');

%load('C:\Users\Ipshita\Documents\PhD documents\Wfs1 project compiled data\E41\E41Data.mat');
%load('C:\Users\Ipshita\Documents\PhD documents\Wfs1 project compiled data\E10\E10spikeData.mat');

axisMin = -1.5;
axisMax = 1.5;
binWidth = .1;
axis1 = axisMin:binWidth:axisMax;
axis2 = axis1-binWidth/2;
axis2(1) = -inf; axis2(end+1) = inf;
%data=[];

     
        for i = 1:length(session) 
            cd(fullfile(folderpath,path{i})) %Select the Proper directory
            saveDir = fullfile(folderpath,path{i},'Analysis');
            %load([saveDir,'\trainData'])
            %load([saveDir,'\rateMapData'])
           
            
          % if strcmp(session{i},'8 Hz')== 1 || strcmp(session{i},'12 Hz')== 1 ||strcmp(session{i},'50 Hz')== 1 
           if strcmp(session{i},'Baseline')== 1
                load([saveDir,'\phasePrecessionStats_Train_LFP6to10'])
            if isempty(data)
                
                data = phasePrecessionStats_Train.dSp.pooledTrial;
            else
                data = [data phasePrecessionStats_Train.dSp.pooledTrial];
            end
            else
                continue
            end
        end   
            
            
        [nSess nCells] = size(data)
        slope = nan(nSess,nCells);
        pVal = nan(nSess,nCells);
        slopeMean = nan(1,nCells);
        pValMean = nan(1,nCells);
        
        for c = 1:nCells
            for ss = 1:nSess
                if ~isempty(data{ss,c})
                    slope(ss,c) = data{ss,c}.s;
                    pVal(ss,c) = data{ss,c}.p;
                end
            end
            
            slopeMean(1,c) = nanmean(slope(:,c));
            nSigN = length(find(pVal(:,c)<0.05 & slope(:,c) < 0));
            nSigP = length(find(pVal(:,c)<0.05 & slope(:,c) > 0));
            nSig = max([nSigP nSigN]);
            nValues = length(find(~isnan(pVal(:,c))));
            if nValues == 0;
                pValMean(1,c) = nan;
            else
                if nSig/nValues >= 0.5
                    pValMean(1,c) = 0.025;
                else
                    pValMean(1,c) = 0.1;
                end
            end
        end
%slopeData = slope;    
%pValData = pVal;
        
        for c = 1:length(cellTypes)
           % ID = ismember(generalizedDefinition_Population,cellTypes{c});
            for m = 1:length(methods)
                switch methods{m}
                    case 'Pooled'
                        %slopeData = slope(:,ID);
                        slopeData = slope;
                        %pValData = pVal(:,ID);
                        pValData = pVal;
                    case 'Single'
                        slopeData = slope(2,ID);
                        pValData = pVal(2,ID);
                    case 'Averaged'
                        slopeData = slopeMean(ID);
                        pValData = pValMean(ID);
                end
                
                Sig = pValData<0.05;
                if length(slopeData(:)) < 5
                    continue
                end
                
                countSig = histc(slopeData(Sig),axis2);
                countNS = histc(slopeData(~Sig),axis2);
                medianSig = median(slopeData(Sig));
                medianAll = nanmedian(slopeData(:));
                
                figure
                total = sum(countNS(1:end-1)+countSig(1:end-1));
                bar(axis1,(countNS(1:end-1)+countSig(1:end-1))/total,'FaceColor','w','EdgeColor',[0 0 0],'LineWidth',1.5)
                hold on
                bar(axis1,countSig(1:end-1)/total,'FaceColor',colors.(cellTypes{c}),'EdgeColor',[0 0 0],'LineWidth',1.5)
                legend({'All Slopes' 'Sig Slopes'},'Location','Best')
                
               title([cellTypes{c},' ',methods{m},' Sessions'])
                xlabel('Relative Slope')
                ylabel('Proportion of Cells')
                plotVertLine(medianSig,{'--','LineWidth',1.5,'Color',colors.(cellTypes{c})});
                plotVertLine(medianAll,{'--k','LineWidth',1.5});
                xlim([axisMin-binWidth axisMax+binWidth])
                ylim([0 .35])
%                 if strcmp(cellTypes{c},'NonSpatial')
%                     ylim([0 .5])
%                 end
                set(gca, 'Box', 'off')
                set(gcf, 'Position', get(0,'Screensize'));
                
                if ~exist(folderpath,'dir')
                    mkdir(folderpath)
                end
                switch saveType
                    case 'AI'
                        saveas(figure(1),[folderpath,'\Hist_',cellTypes{c},'.ai'],'epsc');
                    case 'BMP'
                        saveas(figure(1),[folderpath,'\',methods{m},'Hist_',cellTypes{c},'.bmp'],'bmp16');
                    case 'PDF'
                        saveas(figure(1),[folderpath,'\Hist_',cellTypes{c},'_Baseline6to10.pdf']);
                end
            end
              %  close all
        end
end
           
      
           
