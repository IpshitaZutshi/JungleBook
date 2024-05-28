function precessData = compilePhasePrecession(varargin)

%% defaults
p = inputParser;
addParameter(p,'plotfig',true,@islogical)
addParameter(p,'tonecell',false,@islogical)
addParameter(p,'saveMat',true,@islogical)

parse(p,varargin{:})
plotfig = p.Results.plotfig;
tonecell = p.Results.tonecell;
saveMat = p.Results.saveMat;

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   7
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...11
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    20
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 29
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...33
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'}; %37

% Unless broken or below CA1, keep the lfp channel as 63, the top-most
% channel on shank 2
lfpChan = [63 63 63 63 63 63 63,...
    64 64 67 79,...
    63 63 63 63 63 63 63 63 63,...
    63 63 63 63 63 63 63 63 63,...
    63 63 63 63,...
    63 63 63 63];

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

precessData.slope = [];
precessData.offset = [];
precessData.corrT = [];
precessData.sigP = [];
precessData.cellnum = [];
precessData.sessnum = [];

for ii = 1:length(sess)
    try
        data = sessionPhasePrecession('expPath',strcat(expPath,sess{ii}),'tonecell',tonecell,'lfpChan',lfpChan(ii),'plotfig',false);
        precessData.slope = [precessData.slope; data.slope];
        precessData.offset = [precessData.offset; data.offset];
        precessData.corrT = [precessData.corrT; data.corrT];
        precessData.sigP = [precessData.sigP; data.sigP];
        precessData.cellnum = [precessData.cellnum; data.cellnum'];
        precessData.sessnum = [precessData.sessnum;ones(length(data.cellnum),1)*ii];
    catch
        disp(strcat('Error on sess ',sess{ii}))
        continue
    end
end

% if saveMat
%    save(strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\PhasePrecess',num2str(tonecell),'.mat'), 'precessData'); 
% end

if plotfig
   figure
   set(gcf,'Renderer','painters')
   set(gcf,'Color','w')
   col = [241/255 114/255 42/255;...
        247/255 149/255 33/255;...
        249/255 197/255 81/255;...
        143/255 189/255 107/255;...
        87/255 116/255 144/255];

   for ii = 1:5
        subplot(2,5,ii)
        histogram(precessData.slope(:,ii),-2:0.2:2,'FaceColor',[0.7 0.7 0.7])
        hold on
        line([nanmedian(precessData.slope(:,ii)) nanmedian(precessData.slope(:,ii))],[0 5],'Color',[0.7 0.7 0.7],'LineWidth',1.5)
        histogram(precessData.slope(precessData.sigP(:,ii)<0.05,ii),-2:0.2:2,'FaceColor',col(ii,:))
        line([nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii)) nanmedian(precessData.slope(precessData.sigP(:,ii)<0.05,ii))],[0 5],'Color',col(ii,:),'LineWidth',1.5)
        idxSig{ii} = precessData.sigP(:,ii)<0.05;
   end

   subplot(2,5,6)
   groupStats({precessData.slope(:,1),precessData.slope(:,2),precessData.slope(:,3),...
       precessData.slope(:,4),precessData.slope(:,5)},1:5,'inAxis',true,'color',col)
   title('Slope all') 

   subplot(2,5,7)
   groupStats({precessData.offset(:,1),precessData.offset(:,2),precessData.offset(:,3),...
       precessData.offset(:,4),precessData.offset(:,5)},1:5,'inAxis',true,'color',col)
   title('Intercept all') 
   
   subplot(2,5,8)
   groupStats({precessData.slope(idxSig{1},1),precessData.slope(idxSig{2},2),precessData.slope(idxSig{3},3),...
       precessData.slope(idxSig{4},4),precessData.slope(idxSig{5},5)},1:5,'inAxis',true,'color',col)
   title('Slope sig') 
   
   subplot(2,5,9)
   groupStats({precessData.offset(idxSig{1},1),precessData.offset(idxSig{2},2),precessData.offset(idxSig{3},3),...
       precessData.offset(idxSig{4},4),precessData.offset(idxSig{5},5)},1:5,'inAxis',true,'color',col)
   title('Intercept sig') 
end

end