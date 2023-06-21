function compileMiceAssembliesMazeRipples(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;

mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
    'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ30\Final','IZ31\Final','IZ27\Saline','IZ28\Saline',...
    'IZ29\Saline'};
mice{2} = {'IZ18\Final','IZ31\Final','IZ21\Final','IZ30\Final','IZ15\Final'};


for ii = 1:3
    mecManip{ii} = [];
    ca1Manip{ii} = [];
end

target = 2;

%% Now cycle through to get all the other stimulations
for tag = 1:2
    animalID = mice{tag};
    
    for aa = 1:length(animalID)    
        cd(strcat(parentDir, animalID{aa},'\Summ'));
        if exist(strcat('Summ\AssembliesEntireTrack.mat'),'file') %AssembliesEntireTrackLowThresh
            load(strcat('Summ\AssembliesEntireTrack.mat'));
        else 
            disp(['Assembly not computed for mouse' animalID{aa}])
            continue;
        end
        
        if tag == 1
            if ~isempty(Assemblies{2,2})                    
                mecManip{1} = [mecManip{1}; -(Assemblies{2,target}(:,1)-Assemblies{2,target}(:,2))./(Assemblies{2,target}(:,1)+Assemblies{2,target}(:,2))];                    
                mecManip{2} = [mecManip{2}; -(Assemblies{2,target}(:,3)-Assemblies{2,target}(:,4))./(Assemblies{2,target}(:,3)+Assemblies{2,target}(:,4))];
                mecManip{3} = [mecManip{3}; -(Assemblies{2,target}(:,5)-Assemblies{2,target}(:,6))./(Assemblies{2,target}(:,5)+Assemblies{2,target}(:,6))];
            end                       
        end
            
        if tag == 2
             if ~isempty(Assemblies{1,2})                    
                ca1Manip{1} = [ca1Manip{1}; -(Assemblies{1,target}(:,1)-Assemblies{1,target}(:,2))./(Assemblies{1,target}(:,1)+Assemblies{1,target}(:,2))];                    
                ca1Manip{2} = [ca1Manip{2}; -(Assemblies{1,target}(:,3)-Assemblies{1,target}(:,4))./(Assemblies{1,target}(:,3)+Assemblies{1,target}(:,4))];
                ca1Manip{3} = [ca1Manip{3}; -(Assemblies{1,target}(:,5)-Assemblies{1,target}(:,6))./(Assemblies{1,target}(:,5)+Assemblies{1,target}(:,6))];
             end  
        end
    end
end


figure
set(gcf,'Renderer','painters')
subplot(1,2,1)
data{1} = abs(mecManip{1});
data{2} = abs(mecManip{2});
data{3} = abs(mecManip{3});
dataFin{1} = data{1}(data{1}<2 & data{2}<2 & data{3}<2);
dataFin{2} = data{2}(data{1}<2 & data{2}<2 & data{3}<2);
dataFin{3} = data{3}(data{1}<2 & data{2}<2 & data{3}<2);

stats.mEC = groupStats(dataFin,{},'inAxis',true,'repeatedMeasures',true);

subplot(1,2,2)
data{1} = abs(ca1Manip{1});
data{2} = abs(ca1Manip{2});
data{3} = abs(ca1Manip{3});
dataFin{1} = data{1}(data{1}<2 & data{2}<2 & data{3}<2);
dataFin{2} = data{2}(data{1}<2 & data{2}<2 & data{3}<2);
dataFin{3} = data{3}(data{1}<2 & data{2}<2 & data{3}<2);


stats.CA1 = groupStats(dataFin,{},'inAxis',true,'repeatedMeasures',true);

parentdir = 'C:\Users\ipshi\Dropbox (NYU Langone Health)\RippleManuscript\Revision 1\Revisions\';
saveas(gcf,strcat(parentDir,'AssembliesMazeRipples.png'))
saveas(gcf,strcat(parentDir,'AssembliesMazeRipples.fig'),'fig')    
saveas(gcf,strcat(parentDir,'AssembliesMazeRipples.eps'),'epsc') 

save(strcat(parentDir,'AssembliesMazeRipplesStats.mat'),'stats') 
end