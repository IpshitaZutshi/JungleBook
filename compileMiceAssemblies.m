function compileMiceAssemblies(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
   
% if strcmp(tag,'CA1') == 1
%     mice{3} = {'IZ18\Final','IZ31\Final','IZ20\Final','IZ30\Final','IZ15\Final'};
%     reg = {'CA1','mEC','Both'};
% elseif strcmp(tag,'mEC') == 1
%     mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
%         'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Saline','IZ28\Saline',...
%         'IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
%     reg = {'CA1','mEC','Both'};
% elseif strcmp(tag,'CA3') == 1
%     mice{4} = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
%     reg = {'CA3','mEC','Both'};
% elseif strcmp(tag,'CA3Saline') == 1
%     mice{5} = {'IZ32\Saline','IZ33\Saline','IZ34\Saline'};
%     reg = {'CA3','mEC','Both'};
% elseif strcmp(tag,'mECBilateral') == 1 
%     mice{2} = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
%     reg = {'contramEC','ipsimEC','Both'};
% end
% 
mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
    'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Saline','IZ28\Saline',...
    'IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
mice{2} = {'IZ25\Final','IZ26\Final','IZ24\Final'};  
mice{3} = {'IZ18\Final','IZ31\Final','IZ20\Final','IZ30\Final','IZ15\Final'};
mice{4} = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'};
mice{5} = {'IZ32\Saline','IZ33\Saline','IZ34\Saline'};

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for ii = 1:18
    compiledAssemblies{ii} = [];
end

for tag = 2:5
    if tag == 1
        range = 2;
    elseif tag == 2 || tag == 3
        %range = [1 3];        
        range = [1 2 3];        
    elseif tag == 4 || tag == 5
        %range = [3 4];
        range = [1 2 3 4];
    end
        
    for m = 1:length(mice{tag})
        
        cd(strcat(parentDir, mice{tag}{m},'\Summ'));
        if exist(strcat('Summ\Assemblies.mat'),'file')
            load(strcat('Summ\Assemblies.mat'));
        else 
            disp(['Assembly not computed for mouse' mice{m}])
            continue;
        end

        for ii = range
            for jj = 1
               if tag <4
%                    if tag == 1
%                        currIdx = 2;
%                    elseif tag == 2 && ii == 1
%                        currIdx = 3;
%                    elseif tag == 2 && ii == 3
%                        currIdx = 4;
%                    elseif tag == 3 && ii == 1
%                        currIdx = 5;
%                    elseif tag == 3 && ii == 3
%                        currIdx = 6;   
%                    end
                   if tag == 1
                       currIdx = 2;
                   elseif tag == 2 && ii == 1
                       currIdx = 3;
                   elseif tag == 2 && ii == 2
                       currIdx = 2;
                   elseif tag == 2 && ii == 3
                       currIdx = 4;
                   elseif tag == 3 && ii == 1
                       currIdx = 7;
                   elseif tag == 3 && ii == 2
                       currIdx = 6;                         
                   elseif tag == 3 && ii == 3
                       currIdx = 8;   
                   end
                       
                   if strcmp(mice{tag}{m},'IZ27\Saline')==1 || strcmp(mice{tag}{m},'IZ28\Saline')==1 || strcmp(mice{tag}{m},'IZ29\Saline')==1 ||...
                        strcmp(mice{tag}{m},'IZ32\Saline')==1 || strcmp(mice{tag}{m},'IZ33\Saline')==1 || strcmp(mice{tag}{m},'IZ34\Saline')==1  
                       
                       currManip = Assemblies{jj}(:,[1 2]);
                       if tag == 1
                           baseManip = Assemblies{jj}(:,[5 6]);
                       end
                   else 
                       if ~isempty(Assemblies{ii,jj})
                           currManip = Assemblies{ii,jj}(:,[1 2]);
                           baseManip = Assemblies{ii,jj}(:,[3 4]);                           
                       end
                   end
                   if ~isempty(currManip)
                       assemblyrat = abs(currManip(:,2)-currManip(:,1))./abs(currManip(:,2)+currManip(:,1));
                       compiledAssemblies{currIdx} = [compiledAssemblies{currIdx};assemblyrat];
                       if tag == 2
                           assemblyrat = abs(baseManip(:,2)-baseManip(:,1))./abs(baseManip(:,2)+baseManip(:,1));
                           compiledAssemblies{1} = [compiledAssemblies{1};assemblyrat];
                       elseif tag == 3
                           assemblyrat = abs(baseManip(:,2)-baseManip(:,1))./abs(baseManip(:,2)+baseManip(:,1));
                           compiledAssemblies{5} = [compiledAssemblies{5};assemblyrat];
                       end
                   end
               else
                   if tag == 4 && ii == 1
                       currIdx = 9;
                   elseif tag == 4 && ii == 2
                       currIdx = 10;
                   elseif tag == 4 && ii == 3
                       currIdx = 11;
                   elseif tag == 4 && ii == 4
                       currIdx = 12;                       
                   elseif tag == 5 && ii == 1
                       currIdx = 14;
                   elseif tag == 5 && ii == 2
                       currIdx = 15;
                   elseif tag == 5 && ii == 3
                       currIdx = 16;
                   elseif tag == 5 && ii == 4
                       currIdx = 17;    
                   end

                   if ii ==1
                       currManip = Assemblies{jj}(:,[5 6]);
                       assemblyrat = abs(currManip(:,2)-currManip(:,1))./abs(currManip(:,2)+currManip(:,1));
                       compiledAssemblies{currIdx} = [compiledAssemblies{currIdx};assemblyrat];                            
                   else
                       currManip = Assemblies{jj}(:,[1 ii]);
                       if ~isempty(currManip)                   
                           assemblyrat = abs(currManip(:,2)-currManip(:,1))./abs(currManip(:,2)+currManip(:,1));
                           compiledAssemblies{currIdx} = [compiledAssemblies{currIdx};assemblyrat];       
                       end                
                   end
                   
                   %Add CA3 vs CA3&mEC plot
                   if tag == 4 && ii == 4
                       currManip = Assemblies{jj}(:,[3 4]);
                       if ~isempty(currManip)                       
                           assemblyrat = abs(currManip(:,2)-currManip(:,1))./abs(currManip(:,2)+currManip(:,1));
                           compiledAssemblies{13} = [compiledAssemblies{13};assemblyrat];      
                       end
                   elseif tag == 5 && ii == 4
                       currManip = Assemblies{jj}(:,[3 4]);
                       if ~isempty(currManip)                        
                           assemblyrat = abs(currManip(:,2)-currManip(:,1))./abs(currManip(:,2)+currManip(:,1));
                           compiledAssemblies{18} = [compiledAssemblies{18};assemblyrat];
                       end
                   end                      
               end
            end
        end
    end
end

assemblyData = [];
for ii = 1:18
    assemblyData = catpad(2,assemblyData,compiledAssemblies{ii});    
end

colMat = [0.5 0.5 0.5;8/243 133/243 161/243;103/243 189/243 170/243;56/243 61/243 150/243];
figure
set(gcf,'Position',[50 50 1000 850])
set(gcf,'Renderer','painters')
subplot(2,2,1)
stimdata = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,1:4),[size(assemblyData,1)*4,1]);
stats.mECBilateral = groupStats(assembly,stimdata,'inAxis',true,'color',colMat);
title('baseline ipsimEC   contramEC   bilateralmEC')
ylabel('assembly (stim-base)/(stim+base)')
ylim([0 1])

colMat = [0.5 0.5 0.5;8/243 133/243 161/243;224/243 163/243 46/243;115/243 82/243 68/243];
subplot(2,2,2)
stimdata = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3;ones(size(assemblyData,1),1)*4];
assembly = reshape(assemblyData(:,[5 6 7 8]),[size(assemblyData,1)*4,1]);
stats.CA1 = groupStats(assembly,stimdata,'inAxis',true,'color',colMat);
title('baseline ipsimEC   CA1   ipsimECCA1')
ylabel('assembly (stim-base)/(stim+base)')
%ylim([0 1])

colMat = [0.5 0.5 0.5;8/243 133/243 161/243;187/243 86/243 149/243;94/243 60/243 108/243;52/243 52/243 52/243];%;193/243 90/243 99/243; 133/243 128/243 177/243];
subplot(2,2,3)
stimdata = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3;ones(size(assemblyData,1),1)*4;...
    ones(size(assemblyData,1),1)*5];%;ones(size(assemblyData,1),1)*6];
assembly = reshape(assemblyData(:,[9 10 11 12 13]),[size(assemblyData,1)*5,1]);
stats.CA3 = groupStats(assembly,stimdata,'inAxis',true,'color',colMat);
title('baseline ipsimEC   CA3 ipsimECCA3 CA3vsipsimECCA3')
ylabel('assembly (stim-base)/(stim+base)')
%ylim([0 1])

colMat = [0.5 0.5 0.5;8/243 133/243 161/243;193/243 90/243 99/243; 133/243 128/243 177/243;52/243 52/243 52/243];
subplot(2,2,4)
stimdata = [ones(size(assemblyData,1),1)*1;ones(size(assemblyData,1),1)*2;ones(size(assemblyData,1),1)*3;ones(size(assemblyData,1),1)*4;ones(size(assemblyData,1),1)*5];
assembly = reshape(assemblyData(:,[14 15 16 17 18]),[size(assemblyData,1)*5,1]);
stats.CA3Saline = groupStats(assembly,stimdata,'inAxis',true,'color',colMat);
title('Saline   ipsimEC   CA3 ipsimECCA3 CA3vsipsimECCA3')
ylabel('assembly (stim-base)/(stim+base)')
%ylim([0 1])

saveas(gcf,strcat(parentDir,'Compiled\Assemblies.png'));
saveas(gcf,strcat(parentDir,'Compiled\Assemblies.eps'),'epsc');
saveas(gcf,strcat(parentDir,'Compiled\Assemblies.fig'));
save(strcat(parentDir,'Compiled\Assemblies.mat'),'stats');
end
