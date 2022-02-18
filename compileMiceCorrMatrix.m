function compileMiceCorrMatrix(varargin)

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

mecManip = [];
contraManip = [];
biManip = [];
ca1Manip = [];
ca1mecManip = [];
ca3Manip = [];
ca3salManip = [];

%% Now cycle through to get all the other stimulations
for tag = 1:5
    animalID = mice{tag};
    
    for aa = 1:length(animalID)    
        cd(strcat(parentDir, animalID{aa},'\Summ'));
        if exist(strcat('Summ\CorrMatrix.mat'),'file')
            load(strcat('Summ\CorrMatrix.mat'));
        else 
            disp(['Assembly not computed for mouse' animalID{aa}])
            continue;
        end
        
        if tag == 1
            if aa < 13 % non CA3 mice
                if ~isempty(CorrMatrix{2,1})
                    mecManip = [mecManip; CorrMatrix{2,1}];                 
                end              
            else
                if ~isempty(CorrMatrix)
                    mecManip = [mecManip; CorrMatrix(:,1:3)];                 
                end                 
            end            
        end
            
        if tag == 2
            if ~isempty(CorrMatrix{1,1})
                contraManip = [contraManip; CorrMatrix{1,1}];                 
            end 
            if ~isempty(CorrMatrix{3,1})
                biManip = [biManip; CorrMatrix{3,1}];                 
            end       
        end
        
        if tag == 3
            if ~isempty(CorrMatrix{1,1})
                ca1Manip = [ca1Manip; CorrMatrix{1,1}];                 
            end 
            if ~isempty(CorrMatrix{3,1})
                ca1mecManip = [ca1mecManip; CorrMatrix{3,1}];                 
            end                 
        end        
        
        if tag == 4
            if ~isempty(CorrMatrix)
                ca3Manip = [ca3Manip; CorrMatrix(:,[1,4:8])];                 
            end        
        end            
        
        if tag == 5
            if ~isempty(CorrMatrix)
                ca3salManip = [ca3salManip; CorrMatrix(:,[1,4:5])];                 
            end
        end
    
    end
end

figure
set(gcf,'Position',[680 48 580 962])
set(gcf,'Renderer','painters')

subplot(3,3,1)
stats.mec = groupStats([{mecManip(:,1)} {mecManip(:,2)} {mecManip(:,3)}],[],'RepeatedMeasures',true,'inAxis',true,'plotType','boxplot');
title('ipsi mEC')
ylim([0 1])

subplot(3,3,2)
stats.contra = groupStats([{contraManip(:,1)} {contraManip(:,2)} {contraManip(:,3)}],[],'RepeatedMeasures',true,'inAxis',true,'plotType','boxplot');
title('contra mEC')
ylim([0 1])

subplot(3,3,3)
stats.bi = groupStats([{biManip(:,1)} {biManip(:,2)} {biManip(:,3)}],[],'RepeatedMeasures',true,'inAxis',true,'plotType','boxplot');
title('bilateral mEC')
ylim([0 1])

subplot(3,3,4)
stats.ca1 = groupStats([{ca1Manip(:,1)} {ca1Manip(:,2)} {ca1Manip(:,3)}],[],'RepeatedMeasures',true,'inAxis',true,'plotType','boxplot');
title('ca1')
ylim([0 1])

subplot(3,3,5)
stats.ca1mec = groupStats([{ca1mecManip(:,1)} {ca1mecManip(:,2)} {ca1mecManip(:,3)}],[],'RepeatedMeasures',true,'inAxis',true,'plotType','boxplot');
title('ca1 mEC')
ylim([0 1])

subplot(3,3,7)
stats.ca3 = groupStats([{ca3Manip(:,1)} {ca3Manip(:,2)} {ca3Manip(:,3)} {ca3Manip(:,4)} {ca3Manip(:,5)}],[],'RepeatedMeasures',true,'inAxis',true,'plotType','boxplot');
title('ca3')
ylim([0 1])

subplot(3,3,8)
stats.ca3sal = groupStats([{ca3salManip(:,1)} {ca3salManip(:,2)} {ca3salManip(:,3)}],[],'inAxis',true,'plotType','boxplot');
title('ca3 saline')
ylim([0 1])

saveas(gcf,strcat(parentDir,'Compiled\Assemblies\CorrMatrix.png'))
saveas(gcf,strcat(parentDir,'Compiled\Assemblies\CorrMatrix.fig'),'fig')    
saveas(gcf,strcat(parentDir,'Compiled\Assemblies\CorrMatrix.eps'),'epsc') 

save(strcat(parentDir,'Compiled\Assemblies\CorrMatrixStats.mat'),'stats') 
end