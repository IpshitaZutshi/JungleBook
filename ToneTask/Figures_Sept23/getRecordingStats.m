function getRecordingStats

Ctrlsess= {'IZ41\Final\IZ41_220624_sess5','IZ41\Final\IZ41_220701_sess9',...
    'IZ41\Final\IZ41_220704_sess11','IZ41\Final\IZ41_220708_sess13',...
    'IZ41\Final\IZ41_220629_sess7','IZ41\Final\IZ41_220714_sess14',...
    'IZ45\Final\IZ45_230410_sess11','IZ45\Final\IZ45_230420_sess19',...
    'IZ45\Final\IZ45_230414_sess15','IZ45\Final\IZ45_230417_sess16','IZ45\Final\IZ45_230424_sess21',... 
    'IZ46\Final\IZ46_230406_sess9','IZ46\Final\IZ46_230407_sess10',...
    'IZ46\Final\IZ46_230410_sess11','IZ46\Final\IZ46_230413_sess14','IZ46\Final\IZ46_230420_sess19'};

Tasksess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...27
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'};  

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

for ss = 1:length(Ctrlsess)
    cd(strcat(expPath,Ctrlsess{ss}))  
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);
    possible = {'Pyramidal Cell', 'Narrow Interneuron', 'Wide Interneuron'};
    vals = [0 1 1];
    list = cell_metrics.putativeCellType;
    [tf, idx] = ismember(list, possible);
    list2 = [];
    list2(tf) = vals(idx(tf));
    Ctrlnum(ss,:) = [length(list2) sum(list2) sum(~list2)];
end

for ss = 1:length(Tasksess)
    cd(strcat(expPath,Tasksess{ss}))  
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);
    possible = {'Pyramidal Cell', 'Narrow Interneuron', 'Wide Interneuron'};
    vals = [0 1 1];
    list = cell_metrics.putativeCellType;
    [tf, idx] = ismember(list, possible);
    list2 = [];
    list2(tf) = vals(idx(tf));
    Tasknum(ss,:) = [length(list2) sum(list2) sum(~list2)];    
end


end
