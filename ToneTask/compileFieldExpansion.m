function fieldData = compileFieldExpansion

plotfig = 0;

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
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

fieldData.placefield = [];
fieldData.avgfield = [];
fieldData.avgfieldEdges = [];
fieldData.fieldStart = [];
fieldData.fieldEnd = [];
fieldData.cellnum = [];
fieldData.sessnum = [];
fieldData.avgCOM = [];
fieldData.trialCOM = [];
fieldData.fieldStartRate = [];
fieldData.fieldEndRate = [];

for ii = 1:length(sess)
    data = sessionFieldExpansion('expPath',strcat(expPath,sess{ii}),'plotfig',false);
    idx = find(~isnan(data.avgfield));
    fieldData.placefield = [fieldData.placefield; data.placefield(idx,:)];
    fieldData.avgfield = [fieldData.avgfield; data.avgfield(idx)'];
    fieldData.avgfieldEdges = [fieldData.avgfieldEdges; data.avgfieldEdges(idx,:)];
    fieldData.fieldStart = [fieldData.fieldStart; data.fieldStart(idx,:)];
    fieldData.fieldEnd = [fieldData.fieldEnd; data.fieldEnd(idx,:)];
    fieldData.avgCOM = [fieldData.avgCOM; data.avgCOM(idx)'];
    fieldData.trialCOM = [fieldData.trialCOM; data.trialCOM(idx,:)];
    fieldData.fieldStartRate = [fieldData.fieldStartRate; data.fieldStartRate(idx,:)];
    fieldData.fieldEndRate = [fieldData.fieldEndRate; data.fieldEndRate(idx,:)];

    fieldData.cellnum = [fieldData.cellnum;idx'];
    fieldData.sessnum = [fieldData.sessnum;ones(length(idx),1)*ii];
end

if plotfig

    col = [238/255 67/255 69/255;...
        241/255 114/255 42/255;...
        247/255 149/255 33/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];

    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')
    subplot(2,2,1)
    
    groupStats({fieldData.placefield(:,1),fieldData.placefield(:,2),fieldData.placefield(:,3),...
       fieldData.placefield(:,4),fieldData.placefield(:,5),fieldData.placefield(:,6)},1:6,'inAxis',true,'color',col)
    title('Place field peak') 
    
    subplot(2,2,2)
    subMat  = (fieldData.placefield-fieldData.avgfield);
    groupStats({subMat(:,1),subMat(:,2),subMat(:,3),subMat(:,4),subMat(:,5),subMat(:,6)},1:6,'inAxis',true,'color',col)
    title('Peak loc - avg loc') 
    
    subplot(2,2,3)
    groupStats({fieldData.fieldStart(:,1),fieldData.fieldStart(:,2),fieldData.fieldStart(:,3),...
       fieldData.fieldStart(:,4),fieldData.fieldStart(:,5),fieldData.fieldStart(:,6)},1:6,'inAxis',true,'color',col)
    title('Place field start') 
    
    subplot(2,2,4)
    groupStats({fieldData.fieldEnd(:,1),fieldData.fieldEnd(:,2),fieldData.fieldEnd(:,3),...
       fieldData.fieldEnd(:,4),fieldData.fieldEnd(:,5),fieldData.fieldEnd(:,6)},1:6,'inAxis',true,'color',col)
    title('Place field end') 
end

end