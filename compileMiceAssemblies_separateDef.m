function compileMiceAssemblies_separateDef(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;

mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
    'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ30\Final','IZ31\Final','IZ27\Saline','IZ28\Saline',...
    'IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
mice{2} = {'IZ24\Final','IZ25\Final','IZ26\Final'};  
mice{3} = {'IZ18\Final','IZ31\Final','IZ20\Final','IZ30\Final','IZ15\Final'};
mice{4} = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
mice{5} = {'IZ32\Saline','IZ33\Saline','IZ34\Saline'};

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

% First determine baseline variability
animalID = mice{1};
col1 = cbrewer('seq','BuPu',11);
col = col1(11:-1:1,:);

for ii = 1:2
    baseManip{ii} = [];
    mecManip{ii} = [];
    contraBaseManip{ii} = [];
    contraManip{ii} = [];
    biBaseManip{ii} = [];
    biManip{ii} = [];
    ca1baseManip{ii} = [];
    ca1Manip{ii} = [];
    mecca1baseManip{ii} = [];
    mecca1Manip{ii} = [];
end

nummecbase = [];
nummec = [];
numcontrabase= [];
numcontra= [];
numbibase = [];
numbi = [];
numca1base = [];
numca1 = [];
nummecca1base = [];
nummecca1 = [];
numca3base = [];
numca3 = [];
nummecca3 = [];

nummecbaseA = [];
nummecA = [];
numcontrabaseA= [];
numcontraA= [];
numbibaseA = [];
numbiA = [];
numca1baseA = [];
numca1A = [];
nummecca1baseA = [];
nummecca1A = [];
numca3baseA = [];
numca3A = [];
nummecca3A = [];


for ii = 1:3
    ca3baseManip{ii} = [];
    ca3Manip{ii} = [];
    mecca3Manip{ii} =  [];
    ca3salbaseManip{ii} = [];
    ca3salManip{ii} = [];
    mecca3salManip{ii} =  [];    
end


%% Now cycle through to get all the other stimulations
for tag = 1:5
    animalID = mice{tag};
    
    for aa = 1:length(animalID)    
        cd(strcat(parentDir, animalID{aa},'\Summ'));
        if exist(strcat('Summ\AssembliesSeparate.mat'),'file')
            load(strcat('Summ\AssembliesSeparate.mat'));
        else 
            disp(['Assembly not computed for mouse' animalID{aa}])
            continue;
        end
        
        if tag == 1
            if aa < 13 % non CA3 mice
                if ~isempty(Assemblies.expressionBase{2,1})
                    baseManip{1} = [baseManip{1}; -(Assemblies.expressionBase{2,1}(:,3)- Assemblies.expressionBase{2,1}(:,4))./(Assemblies.expressionBase{2,1}(:,3)+ Assemblies.expressionBase{2,1}(:,4))];
                    mecManip{1} = [mecManip{1}; -(Assemblies.expressionBase{2,1}(:,1)-Assemblies.expressionBase{2,1}(:,2))./(Assemblies.expressionBase{2,1}(:,1)+Assemblies.expressionBase{2,1}(:,2))];                    
                end
                if ~isempty(Assemblies.expressionStim{2,1})
                    baseManip{2} = [baseManip{2}; -(Assemblies.expressionStim{2,1}(:,2)- Assemblies.expressionStim{2,1}(:,1))./(Assemblies.expressionStim{2,1}(:,2)+ Assemblies.expressionStim{2,1}(:,1))];
                    mecManip{2} = [mecManip{2}; -(Assemblies.expressionStim{2,1}(:,3)-Assemblies.expressionStim{2,1}(:,4))./(Assemblies.expressionStim{2,1}(:,3)+Assemblies.expressionStim{2,1}(:,4))];
                end  
                nummecbase = [nummecbase Assemblies.numberBase{2,1}];
                nummec = [nummec Assemblies.numberStim{2,1}];

                nummecbaseA = [nummecbaseA Assemblies.numberBaseA{2,1}];
                nummecA = [nummecA Assemblies.numberStimA{2,1}];                
            else
               baseManip{1} = [baseManip{1}; -(Assemblies.expressionBase(:,5)- Assemblies.expressionBase(:,6))./(Assemblies.expressionBase(:,5)+Assemblies.expressionBase(:,6))]; 
               mecManip{1} = [mecManip{1}; -(Assemblies.expressionBase(:,1)- Assemblies.expressionBase(:,2))./(Assemblies.expressionBase(:,1)+Assemblies.expressionBase(:,2))]; 

               baseManip{2} = [baseManip{2}; (Assemblies.expressionmEC(:,2)- Assemblies.expressionmEC(:,1))./(Assemblies.expressionmEC(:,2)+Assemblies.expressionmEC(:,1))]; 
               mecManip{2} = [mecManip{2}; (Assemblies.expressionmEC(:,5)- Assemblies.expressionmEC(:,6))./(Assemblies.expressionmEC(:,5)+Assemblies.expressionmEC(:,6))];                
               
               nummecbase = [nummecbase Assemblies.numberBase];
               nummec = [nummec Assemblies.numbermEC];         
               
               nummecbaseA = [nummecbaseA Assemblies.numberBaseA];
               nummecA = [nummecA Assemblies.numbermECA];                   
            end            
        end
            
        if tag == 2
            if ~isempty(Assemblies.expressionBase{1,1})
                contraBaseManip{1} = [contraBaseManip{1}; -(Assemblies.expressionBase{1,1}(:,3)-Assemblies.expressionBase{1,1}(:,4))./(Assemblies.expressionBase{1,1}(:,3)+Assemblies.expressionBase{1,1}(:,4))];
                contraBaseManip{2} = [contraBaseManip{2}; -(Assemblies.expressionStim{1,1}(:,2)-Assemblies.expressionStim{1,1}(:,1))./(Assemblies.expressionStim{1,1}(:,2)+Assemblies.expressionStim{1,1}(:,1))];
                contraManip{1} = [contraManip{1}; -(Assemblies.expressionBase{1,1}(:,1)-Assemblies.expressionBase{1,1}(:,2))./(Assemblies.expressionBase{1,1}(:,1)+Assemblies.expressionBase{1,1}(:,2))];
                contraManip{2} = [contraManip{2}; -(Assemblies.expressionStim{1,1}(:,3)-Assemblies.expressionStim{1,1}(:,4))./(Assemblies.expressionStim{1,1}(:,3)+Assemblies.expressionStim{1,1}(:,4))];                
                
                numcontrabase= [numcontrabase Assemblies.numberBase{1,1}];
                numcontra= [numcontra Assemblies.numberStim{1,1}];
                
                numcontrabaseA= [numcontrabaseA Assemblies.numberBaseA{1,1}];
                numcontraA= [numcontraA Assemblies.numberStimA{1,1}];                
            end
            if ~isempty(Assemblies.expressionBase{3,1})
                biBaseManip{1} = [biBaseManip{1}; -(Assemblies.expressionBase{3,1}(:,3)-Assemblies.expressionBase{3,1}(:,4))./(Assemblies.expressionBase{3,1}(:,3)+Assemblies.expressionBase{3,1}(:,4))];
                biBaseManip{2} = [biBaseManip{2}; -(Assemblies.expressionStim{3,1}(:,2)-Assemblies.expressionStim{3,1}(:,1))./(Assemblies.expressionStim{3,1}(:,2)+Assemblies.expressionStim{3,1}(:,1))];
                biManip{1} = [biManip{1}; -(Assemblies.expressionBase{3,1}(:,1)-Assemblies.expressionBase{3,1}(:,2))./(Assemblies.expressionBase{3,1}(:,1)+Assemblies.expressionBase{3,1}(:,2))];
                biManip{2} = [biManip{2}; -(Assemblies.expressionStim{3,1}(:,3)-Assemblies.expressionStim{3,1}(:,4))./(Assemblies.expressionStim{3,1}(:,3)+Assemblies.expressionStim{3,1}(:,4))];                
                
                numbibase = [numbibase Assemblies.numberBase{3,1}];
                numbi = [numbi Assemblies.numberStim{3,1}];

                numbibaseA = [numbibaseA Assemblies.numberBaseA{3,1}];
                numbiA = [numbiA Assemblies.numberStimA{3,1}];                
            end         
        end
        
        if tag == 3
            if ~isempty(Assemblies.expressionBase{1,1})
                ca1baseManip{1} = [ca1baseManip{1}; -(Assemblies.expressionBase{1,1}(:,3)-Assemblies.expressionBase{1,1}(:,4))./(Assemblies.expressionBase{1,1}(:,3)+Assemblies.expressionBase{1,1}(:,4))];
                ca1baseManip{2} = [ca1baseManip{2}; -(Assemblies.expressionStim{1,1}(:,2)-Assemblies.expressionStim{1,1}(:,1))./(Assemblies.expressionStim{1,1}(:,2)+Assemblies.expressionStim{1,1}(:,1))];
                ca1Manip{1} = [ca1Manip{1}; -(Assemblies.expressionBase{1,1}(:,1)-Assemblies.expressionBase{1,1}(:,2))./(Assemblies.expressionBase{1,1}(:,1)+Assemblies.expressionBase{1,1}(:,2))];
                ca1Manip{2} = [ca1Manip{2}; -(Assemblies.expressionStim{1,1}(:,3)-Assemblies.expressionStim{1,1}(:,4))./(Assemblies.expressionStim{1,1}(:,3)+Assemblies.expressionStim{1,1}(:,4))];                
                
                numca1base = [numca1base Assemblies.numberBase{1,1}];
                numca1 = [numca1 Assemblies.numberStim{1,1}];   

                numca1baseA = [numca1baseA Assemblies.numberBaseA{1,1}];
                numca1A = [numca1A Assemblies.numberStimA{1,1}];                   
            end
            if ~isempty(Assemblies.expressionBase{3,1})
                mecca1baseManip{1} = [mecca1baseManip{1}; -(Assemblies.expressionBase{3,1}(:,3)-Assemblies.expressionBase{3,1}(:,4))./(Assemblies.expressionBase{3,1}(:,3)+Assemblies.expressionBase{3,1}(:,4))];
                mecca1baseManip{2} = [mecca1baseManip{2}; -(Assemblies.expressionStim{3,1}(:,2)-Assemblies.expressionStim{3,1}(:,1))./(Assemblies.expressionStim{3,1}(:,2)+Assemblies.expressionStim{3,1}(:,1))];
                mecca1Manip{1} = [mecca1Manip{1}; -(Assemblies.expressionBase{3,1}(:,1)-Assemblies.expressionBase{3,1}(:,2))./(Assemblies.expressionBase{3,1}(:,1)+Assemblies.expressionBase{3,1}(:,2))];
                mecca1Manip{2} = [mecca1Manip{2}; -(Assemblies.expressionStim{3,1}(:,3)-Assemblies.expressionStim{3,1}(:,4))./(Assemblies.expressionStim{3,1}(:,3)+Assemblies.expressionStim{3,1}(:,4))];                                       
                
                nummecca1base = [nummecca1base Assemblies.numberBase{3,1}];
                nummecca1 = [nummecca1 Assemblies.numberStim{3,1}];   
                
                nummecca1baseA = [nummecca1baseA Assemblies.numberBaseA{3,1}];
                nummecca1A = [nummecca1A Assemblies.numberStimA{3,1}];                  
            end               
        end        
        
        if tag == 4
            ca3baseManip{1} = [ca3baseManip{1}; -(Assemblies.expressionBase(:,5)- Assemblies.expressionBase(:,6))./(Assemblies.expressionBase(:,5)+Assemblies.expressionBase(:,6))];
            ca3Manip{1} = [ca3Manip{1}; -(Assemblies.expressionBase(:,1)- Assemblies.expressionBase(:,3))./(Assemblies.expressionBase(:,1)+Assemblies.expressionBase(:,3))];
            mecca3Manip{1} = [mecca3Manip{1}; -(Assemblies.expressionBase(:,1)- Assemblies.expressionBase(:,4))./(Assemblies.expressionBase(:,1)+Assemblies.expressionBase(:,4))];
            
            ca3baseManip{2} = [ca3baseManip{2}; -(Assemblies.expressionCA3(:,3)- Assemblies.expressionCA3(:,1))./(Assemblies.expressionCA3(:,3)+Assemblies.expressionCA3(:,1))];
            ca3Manip{2} = [ca3Manip{2}; -(Assemblies.expressionCA3(:,5)- Assemblies.expressionCA3(:,6))./(Assemblies.expressionCA3(:,5)+Assemblies.expressionCA3(:,6))];
            mecca3Manip{2}= [mecca3Manip{2}; -(Assemblies.expressionCA3(:,3)- Assemblies.expressionCA3(:,4))./(Assemblies.expressionCA3(:,3)+Assemblies.expressionCA3(:,4))];
            
            ca3baseManip{3} = [ca3baseManip{3}; -(Assemblies.expressionBoth(:,4)- Assemblies.expressionBoth(:,1))./(Assemblies.expressionBoth(:,4)+Assemblies.expressionBoth(:,1))];
            ca3Manip{3}  = [ca3Manip{3}; -(Assemblies.expressionBoth(:,4)- Assemblies.expressionBoth(:,3))./(Assemblies.expressionBoth(:,4)+Assemblies.expressionBoth(:,3))];
            mecca3Manip{3} = [mecca3Manip{3}; -(Assemblies.expressionBoth(:,5)- Assemblies.expressionBoth(:,6))./(Assemblies.expressionBoth(:,5)+Assemblies.expressionBoth(:,6))];
            
            numca3base = [numca3base Assemblies.numberBase];
            numca3 = [numca3 Assemblies.numberCA3];
            nummecca3 = [nummecca3 Assemblies.numberBoth];
            
            numca3baseA = [numca3baseA Assemblies.numberBaseA];
            numca3A = [numca3A Assemblies.numberCA3A];
            nummecca3A = [nummecca3A Assemblies.numberBothA];            
        end            
        
        if tag == 5
            ca3salbaseManip{1} = [ca3salbaseManip{1}; -(Assemblies.expressionBase(:,5)- Assemblies.expressionBase(:,6))./(Assemblies.expressionBase(:,5)+Assemblies.expressionBase(:,6))];
            ca3salManip{1} = [ca3salManip{1}; -(Assemblies.expressionBase(:,1)- Assemblies.expressionBase(:,3))./(Assemblies.expressionBase(:,1)+Assemblies.expressionBase(:,3))];
            mecca3salManip{1} = [mecca3salManip{1}; -(Assemblies.expressionBase(:,1)- Assemblies.expressionBase(:,4))./(Assemblies.expressionBase(:,1)+Assemblies.expressionBase(:,4))];
            
            ca3salbaseManip{2} = [ca3salbaseManip{2}; -(Assemblies.expressionCA3(:,3)- Assemblies.expressionCA3(:,1))./(Assemblies.expressionCA3(:,3)+Assemblies.expressionCA3(:,1))];
            ca3salManip{2} = [ca3salManip{2}; -(Assemblies.expressionCA3(:,5)- Assemblies.expressionCA3(:,6))./(Assemblies.expressionCA3(:,5)+Assemblies.expressionCA3(:,6))];
            mecca3salManip{2}= [mecca3salManip{2}; -(Assemblies.expressionCA3(:,3)- Assemblies.expressionCA3(:,4))./(Assemblies.expressionCA3(:,3)+Assemblies.expressionCA3(:,4))];
            
            ca3salbaseManip{3} = [ca3salbaseManip{3}; -(Assemblies.expressionBoth(:,4)- Assemblies.expressionBoth(:,1))./(Assemblies.expressionBoth(:,4)+Assemblies.expressionBoth(:,1))];
            ca3salManip{3}  = [ca3salManip{3}; -(Assemblies.expressionBoth(:,4)- Assemblies.expressionBoth(:,3))./(Assemblies.expressionBoth(:,4)+Assemblies.expressionBoth(:,3))];
            mecca3salManip{3} = [mecca3salManip{3}; -(Assemblies.expressionBoth(:,5)- Assemblies.expressionBoth(:,6))./(Assemblies.expressionBoth(:,5)+Assemblies.expressionBoth(:,6))];
        end          
    
    end
end

matrix = [];
for ii = 1:2
    matrix(ii,:) = [nanmedian(baseManip{ii}) nanmedian(mecManip{ii}) nanmedian(contraBaseManip{ii}) nanmedian(contraManip{ii}) ...
        nanmedian(biBaseManip{ii}) nanmedian(biManip{ii}) nanmedian(ca1baseManip{ii}) nanmedian(ca1Manip{ii}) nanmedian(mecca1baseManip{ii}) nanmedian(mecca1Manip{ii})]; 
end

[stats.mEC.Stability.p, ~,stats.mEC.Stability.stats] = ranksum(baseManip{1},mecManip{2});
[stats.contra.Stability.p, ~,stats.contra.Stability.stats] = ranksum(contraBaseManip{1},contraManip{2});
[stats.bi.Stability.p, ~,stats.bi.Stability.stats] = ranksum(biBaseManip{1},biManip{2});
[stats.CA1.Stability.p, ~,stats.CA1.Stability.stats] = ranksum(ca1baseManip{1},ca1Manip{2});
[stats.CA1mec.Stability.p, ~,stats.CA1mec.Stability.stats] = ranksum(mecca1baseManip{1},mecca1Manip{2});

figure
set(gcf,'renderer','painters')
subplot(2,3,1)
groupStats([baseManip(1) mecManip(2)],[],'inAxis',true,'plotType','boxplot');
title('ipsi')
ylim([-1 1])

subplot(2,3,2)
groupStats([contraBaseManip(1) contraManip(2)],[],'inAxis',true,'plotType','boxplot');
title('contra')
ylim([-1 1])

subplot(2,3,3)
groupStats([biBaseManip(1) biManip(2)],[],'inAxis',true,'plotType','boxplot');
title('bi')
ylim([-1 1])

subplot(2,3,4)
groupStats([ca1baseManip(1) ca1Manip(2)],[],'inAxis',true,'plotType','boxplot');
title('ca1')
ylim([-1 1])

subplot(2,3,5)
groupStats([mecca1baseManip(1) mecca1Manip(2)],[],'inAxis',true,'plotType','boxplot');
title('ca1mec')
ylim([-1 1])

subplot(2,3,6)
groupStats([ca3baseManip(1) ca3Manip(2) mecca3Manip(3)],[],'inAxis',true,'plotType','boxplot');
title('CA3')
ylim([-1 1])


saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyStability.png'))
saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyStability.fig'),'fig')    
saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyStability.eps'),'epsc') 


figure
set(gcf,'renderer','painters')
subplot(3,5,1)
imagesc(matrix(:,1:2))
colormap(col)
caxis([-0.8 0.2])
title(num2str(stats.mEC.Stability.p));

subplot(3,5,2)
imagesc(matrix(:,3:4))
colormap(col)
caxis([-0.8 0.2])
title(num2str(stats.contra.Stability.p));

subplot(3,5,3)
imagesc(matrix(:,5:6))
colormap(col)
caxis([-0.8 0.2])
title(num2str(stats.bi.Stability.p));

subplot(3,5,4)
imagesc(matrix(:,7:8))
colormap(col)
caxis([-0.8 0.2])
title(num2str(stats.CA1.Stability.p));

subplot(3,5,5)
imagesc(matrix(:,9:10))
colormap(col)
caxis([-0.8 0.2])
title(num2str(stats.CA1mec.Stability.p));


matrix = [];
matrix = [nanmedian(ca3baseManip{1}) nanmedian(ca3Manip{1}) nanmedian(mecca3Manip{1});...
    nanmedian(ca3baseManip{2}) nanmedian(ca3Manip{2}) nanmedian(mecca3Manip{2});...
    nanmedian(ca3baseManip{3}) nanmedian(ca3Manip{3}) nanmedian(mecca3Manip{3})];
subplot(3,5,6:10)
imagesc(matrix)
colormap(col)
colorbar
caxis([-0.8 0.2])

[stats.CA3.Stability] = groupStats([ca3baseManip(1) ca3Manip(2) mecca3Manip(3)],[],'doPlot',false);
title(num2str(stats.CA3.Stability.kruskalWallis.p))

subplot(3,5,11:15)
[stats.Expression] = groupStats([mecManip(1) contraManip(1) biManip(1) ca1Manip(1) mecca1Manip(1) ca3Manip(1) mecca3Manip(1) ca3salManip(1)],[],'doPlot',true,'inAxis',true,'plotType','boxplot');

saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyExpression.png'))
saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyExpression.fig'),'fig')    
saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyExpression.eps'),'epsc') 

figure
set(gcf,'renderer','painters')
subplot(2,6,1)
groupStats([{nummecbase} {nummec}],[],'inAxis',true,'plotType','boxplot');
[stats.numcells.mEC.p,~,stats.numcells.mEC.stats]  = ranksum(nummecbase,nummec);
[stats.numcells.mEC.num]  = [length(nummecbase) length(nummec)];
title(num2str(stats.numcells.mEC.p));
ylim([0 15])

subplot(2,6,2)
groupStats([{numcontrabase} {numcontra}],[],'inAxis',true,'plotType','boxplot');
[stats.numcells.contra.p,~,stats.numcells.contra.stats]  = ranksum(numcontrabase,numcontra);
[stats.numcells.contra.num]  = [length(numcontrabase) length(numcontra)];
title(num2str(stats.numcells.contra.p));
ylim([0 15])

subplot(2,6,3)
groupStats([{numbibase} {numbi}],[],'inAxis',true,'plotType','boxplot');
[stats.numcells.bi.p,~,stats.numcells.bi.stats]  = ranksum(numbibase,numbi);
[stats.numcells.bi.num]  = [length(numbibase) length(numbi)];
title(num2str(stats.numcells.bi.p));
ylim([0 15])

subplot(2,6,4)
groupStats([{numca1base} {numca1}],[],'inAxis',true,'plotType','boxplot');
[stats.numcells.ca1.p,~,stats.numcells.ca1.stats]  = ranksum(numca1base,numca1);
[stats.numcells.ca1.num]  = [length(numca1base) length(numca1)];
title(num2str(stats.numcells.ca1.p));
ylim([0 15])

subplot(2,6,5)
groupStats([{nummecca1base} {nummecca1}],[],'inAxis',true,'plotType','boxplot');
[stats.numcells.mecca1.p,~,stats.numcells.mecca1.stats]  = ranksum(nummecca1base,nummecca1);
[stats.numcells.mecca1.num]  = [length(nummecca1base) length(nummecca1)];
title(num2str(stats.numcells.mecca1.p));
ylim([0 15])

subplot(2,6,6)
stats.numcells.ca3 = groupStats([{numca3base} {numca3} {nummecca3}],[],'inAxis',true,'plotType','boxplot');
title(num2str(stats.numcells.ca3.kruskalWallis.p));
ylim([0 15])

%%%%%%%%%%%%%%%%
subplot(2,6,7)
groupStats([{nummecbaseA} {nummecA}],[],'inAxis',true,'plotType','boxplot');
[stats.numAss.mEC.p,~,stats.numAss.mEC.stats]  = ranksum(nummecbaseA,nummecA);
[stats.numAss.mEC.num]  = [length(nummecbaseA) length(nummecA)];
title(num2str(stats.numAss.mEC.p));
ylim([0 15])

subplot(2,6,8)
groupStats([{numcontrabaseA} {numcontraA}],[],'inAxis',true,'plotType','boxplot');
[stats.numAss.contra.p,~,stats.numAss.contra.stats]  = ranksum(numcontrabaseA,numcontraA);
[stats.numAss.contra.num]  = [length(numcontrabaseA) length(numcontraA)];
title(num2str(stats.numAss.contra.p));
ylim([0 15])

subplot(2,6,9)
groupStats([{numbibaseA} {numbiA}],[],'inAxis',true,'plotType','boxplot');
[stats.numAss.bi.p,~,stats.numAss.bi.stats]  = ranksum(numbibaseA,numbiA);
[stats.numAss.bi.num]  = [length(numbibaseA) length(numbiA)];
title(num2str(stats.numAss.bi.p));
ylim([0 15])

subplot(2,6,10)
groupStats([{numca1base} {numca1}],[],'inAxis',true,'plotType','boxplot');
[stats.numAss.ca1.p,~,stats.numAss.ca1.stats]  = ranksum(numca1baseA,numca1A);
[stats.numAss.ca1.num]  = [length(numca1baseA) length(numca1A)];
title(num2str(stats.numAss.ca1.p));
ylim([0 15])

subplot(2,6,11)
groupStats([{nummecca1baseA} {nummecca1A}],[],'inAxis',true,'plotType','boxplot');
[stats.numAss.mecca1.p,~,stats.numAss.mecca1.stats]  = ranksum(nummecca1baseA,nummecca1A);
[stats.numAss.mecca1.num]  = [length(nummecca1baseA) length(nummecca1A)];
title(num2str(stats.numAss.mecca1.p));
ylim([0 15])

subplot(2,6,12)
stats.numAss.ca3 = groupStats([{numca3baseA} {numca3A} {nummecca3A}],[],'inAxis',true,'plotType','boxplot');
title(num2str(stats.numAss.ca3.kruskalWallis.p));
ylim([0 15])

saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyNum.png'))
saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyNum.fig'),'fig')    
saveas(gcf,strcat(parentDir,'Compiled\Assemblies\AssemblyNum.eps'),'epsc') 

save(strcat(parentDir,'Compiled\Assemblies\separateStats.mat'),'stats') 
end