function compileMiceAssemblies_proportionchange(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;

mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
    'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ30\Final','IZ31\Final','IZ27\Saline','IZ28\Saline',...
    'IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
mice{2} = {'IZ25\Final','IZ26\Final','IZ24\Final'};  
mice{3} = {'IZ18\Final','IZ31\Final','IZ20\Final','IZ30\Final','IZ15\Final'};
mice{4} = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
mice{5} = {'IZ32\Saline','IZ33\Saline','IZ34\Saline'};

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

% First determine baseline variability
animalID = mice{1};

baseManip = [];
mecManip = [];

for aa = 1:length(animalID)
    
    cd(strcat(parentDir, animalID{aa},'\Summ'));
    
    if exist(strcat('Summ\Assemblies.mat'),'file')
        load(strcat('Summ\Assemblies.mat'));
    else 
        disp(['Assembly not computed for mouse' animalID{aa}])
        continue;
    end
    
    if aa < 13 % non CA3 mice
        for ii = 1:3
            if ~isempty(Assemblies{ii,1})
                baseManip = [baseManip; Assemblies{ii,1}(:,3:4)];
                if ii == 2
                    mecManip = [mecManip; (Assemblies{ii,1}(:,3)-Assemblies{ii,1}(:,6))./(Assemblies{ii,1}(:,3)+Assemblies{ii,1}(:,6))];
                end
            end
        end
    else
       baseManip = [baseManip; Assemblies{1,1}(:,5:6)]; 
       mecManip = [mecManip;(Assemblies{1,1}(:,5)-Assemblies{1,1}(:,8))./(Assemblies{1,1}(:,5)+Assemblies{1,1}(:,8))];
    end
end

manip = (baseManip(:,1)-baseManip(:,2))./(baseManip(:,1)+baseManip(:,2));
threshPos = mean(manip)+2*std(manip);
threshNeg = mean(manip)-2*std(manip);

contraManip = [];
biManip = [];
ca1Manip = [];
mecca1Manip = [];
ca3Manip = [];
mecca3Manip =  [];
ca3salManip =[];

%% Now cycle through to get all the other stimulations
for tag = 2:5
    animalID = mice{tag};
    
    for aa = 1:length(animalID)    
        cd(strcat(parentDir, animalID{aa},'\Summ'));
        if exist(strcat('Summ\Assemblies.mat'),'file')
            load(strcat('Summ\Assemblies.mat'));
        else 
            disp(['Assembly not computed for mouse' animalID{aa}])
            continue;
        end
            
        if tag == 2
            if ~isempty(Assemblies{1,1})
                contraManip = [contraManip; (Assemblies{1,1}(:,3)-Assemblies{1,1}(:,6))./(Assemblies{1,1}(:,3)+Assemblies{1,1}(:,6))];
            end
            if ~isempty(Assemblies{3,1})
                biManip = [biManip; (Assemblies{3,1}(:,3)-Assemblies{3,1}(:,6))./(Assemblies{3,1}(:,3)+Assemblies{3,1}(:,6))];
            end            
        end
        
        if tag == 3
            if ~isempty(Assemblies{1,1})
                ca1Manip = [ca1Manip; ((Assemblies{1,1}(:,3)-Assemblies{1,1}(:,6)))./((Assemblies{1,1}(:,3)+Assemblies{1,1}(:,6)))];
            end
            if ~isempty(Assemblies{3,1})
                mecca1Manip = [mecca1Manip; (Assemblies{3,1}(:,3)-Assemblies{3,1}(:,6))./(Assemblies{3,1}(:,3)+Assemblies{3,1}(:,6))];
            end            
        end        
        
        if tag == 4
            if ~isempty(Assemblies{1,1})
                ca3Manip = [ca3Manip; (Assemblies{1,1}(:,5)-Assemblies{1,1}(:,9))./(Assemblies{1,1}(:,5)+Assemblies{1,1}(:,9))];
            end
            if ~isempty(Assemblies{1,1})
                mecca3Manip = [mecca3Manip; (Assemblies{1,1}(:,5)-Assemblies{1,1}(:,11))./(Assemblies{1,1}(:,5)+Assemblies{1,1}(:,11))];
            end            
        end         
        
        if tag == 5
            if ~isempty(Assemblies{1,1})
                ca3salManip = [ca3salManip; (Assemblies{1,1}(:,5)-Assemblies{1,1}(:,9))./(Assemblies{1,1}(:,5)+Assemblies{1,1}(:,9))];
            end          
        end           
    
    end
end

% find percent above, no change and below - column order - baseline, mec,
% contra, bi, ca1, both, ca3, both, ca3 saline
ChangeIdx = zeros(3,9);
ChangeIdx(:,1) = [sum(manip<=threshPos & manip>=threshNeg);sum(manip>threshPos);sum(manip<threshNeg)]./length(manip);
ChangeIdx(:,2) = [sum(mecManip<=threshPos & mecManip>=threshNeg);sum(mecManip>threshPos);sum(mecManip<threshNeg)]./length(mecManip);
ChangeIdx(:,3) = [sum(contraManip<=threshPos & contraManip>=threshNeg);sum(contraManip>threshPos);sum(contraManip<threshNeg)]./length(contraManip);
ChangeIdx(:,4) = [sum(biManip<=threshPos & biManip>=threshNeg);sum(biManip>threshPos);sum(biManip<threshNeg)]./length(biManip);
ChangeIdx(:,5) = [sum(ca1Manip<=threshPos & ca1Manip>=threshNeg);sum(ca1Manip>threshPos);sum(ca1Manip<threshNeg)]./length(ca1Manip);
ChangeIdx(:,6) = [sum(mecca1Manip<=threshPos & mecca1Manip>=threshNeg);sum(mecca1Manip>threshPos);sum(mecca1Manip<threshNeg)]./length(mecca1Manip);
ChangeIdx(:,7) = [sum(ca3Manip<=threshPos & ca3Manip>=threshNeg);sum(ca3Manip>threshPos);sum(ca3Manip<threshNeg)]./length(ca3Manip);
ChangeIdx(:,8) = [sum(mecca3Manip<=threshPos & mecca3Manip>=threshNeg);sum(mecca3Manip>threshPos);sum(mecca3Manip<threshNeg)]./length(mecca3Manip);
ChangeIdx(:,9) = [sum(ca3salManip<=threshPos & ca3salManip>=threshNeg);sum(ca3salManip>threshPos);sum(ca3salManip<threshNeg)]./length(ca3salManip);

figure
set(gcf,'Renderer','painters')
bar(ChangeIdx','stacked')

ChangeNum = zeros(3,9);
ChangeNum(:,1) = [sum(manip<=threshPos & manip>=threshNeg);sum(manip>threshPos);sum(manip<threshNeg)];
ChangeNum(:,2) = [sum(mecManip<=threshPos & mecManip>=threshNeg);sum(mecManip>threshPos);sum(mecManip<threshNeg)];
ChangeNum(:,3) = [sum(contraManip<=threshPos & contraManip>=threshNeg);sum(contraManip>threshPos);sum(contraManip<threshNeg)];
ChangeNum(:,4) = [sum(biManip<=threshPos & biManip>=threshNeg);sum(biManip>threshPos);sum(biManip<threshNeg)];
ChangeNum(:,5) = [sum(ca1Manip<=threshPos & ca1Manip>=threshNeg);sum(ca1Manip>threshPos);sum(ca1Manip<threshNeg)];
ChangeNum(:,6) = [sum(mecca1Manip<=threshPos & mecca1Manip>=threshNeg);sum(mecca1Manip>threshPos);sum(mecca1Manip<threshNeg)];
ChangeNum(:,7) = [sum(ca3Manip<=threshPos & ca3Manip>=threshNeg);sum(ca3Manip>threshPos);sum(ca3Manip<threshNeg)];
ChangeNum(:,8) = [sum(mecca3Manip<=threshPos & mecca3Manip>=threshNeg);sum(mecca3Manip>threshPos);sum(mecca3Manip<threshNeg)];
ChangeNum(:,9) = [sum(ca3salManip<=threshPos & ca3salManip>=threshNeg);sum(ca3salManip>threshPos);sum(ca3salManip<threshNeg)];

saveas(gcf,strcat(parentDir,'Compiled\Assembliesproportion.png'));
saveas(gcf,strcat(parentDir,'Compiled\Assembliesproportion.eps'),'epsc');
saveas(gcf,strcat(parentDir,'Compiled\Assembliesproportion.fig'));

end