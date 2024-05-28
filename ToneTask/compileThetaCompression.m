function compileThetaCompression

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


for jj = 1:2
    compiledThetaComp.placefield{jj} = [];
    compiledThetaComp.tonefield{jj} = [];
    compiledThetaComp.thetaDT{jj} = [];
    compiledThetaComp.placeDT{jj} = [];
    compiledThetaComp.placeDist{jj} = [];
end

for ii = 1:length(sess)
    thetaComp = sessionThetaCompression('expPath',strcat(expPath,sess{ii}),'plotfig',false);
    for jj = 1:2
        compiledThetaComp.placefield{jj} = [compiledThetaComp.placefield{jj}; thetaComp.placefield{jj}];
        compiledThetaComp.tonefield{jj} = [compiledThetaComp.tonefield{jj}; thetaComp.tonefield{jj}];
        compiledThetaComp.thetaDT{jj} = [compiledThetaComp.thetaDT{jj}; thetaComp.thetaDT{jj}];
        compiledThetaComp.placeDT{jj} = [compiledThetaComp.placeDT{jj}; thetaComp.placeDT{jj}];
        compiledThetaComp.placeDist{jj} = [compiledThetaComp.placeDist{jj}; thetaComp.placeDist{jj}];
    end
end

saveas('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\thetaCompression.mat',compiledThetaComp);

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')

selectidx{1} = (abs(compiledThetaComp.placeDist{1})<10) & abs(compiledThetaComp.thetaDT{1})<250;
selectidx{2} = (abs(compiledThetaComp.placeDist{2})<10) & abs(compiledThetaComp.thetaDT{2}) < 250;

% 
% selectidx_place = compiledThetaComp.placefield{1}==1 & compiledThetaComp.placefield{2}==1 & compiledThetaComp.thetaDT{1} <250 & compiledThetaComp.thetaDT{1} >-250 ...
%         & compiledThetaComp.thetaDT{2} <250 & compiledThetaComp.thetaDT{2} >-250 & ...
%         compiledThetaComp.placeDT{1} <800 & compiledThetaComp.placeDT{1} >-800 ...
%         & compiledThetaComp.placeDT{2} <800 & compiledThetaComp.placeDT{2} >-800;
% 
% selectidx_tone = compiledThetaComp.tonefield{1}==1 & compiledThetaComp.tonefield{2}==1 & compiledThetaComp.thetaDT{1} <250 & compiledThetaComp.thetaDT{1} >-250 ...
%         & compiledThetaComp.thetaDT{2} <250 & compiledThetaComp.thetaDT{2} >-250 & ...
%         compiledThetaComp.placeDT{1} <800 & compiledThetaComp.placeDT{1} >-800 ...
%         & compiledThetaComp.placeDT{2} <800 & compiledThetaComp.placeDT{2} >-800;

for ii = 1:2
    subplot(2,4,ii)
    %scatter(compiledThetaComp.thetaDT{ii}(selectidx_place),compiledThetaComp.placeDT{ii}(selectidx_place),'.')
    scatter(compiledThetaComp.placeDist{ii}(selectidx{ii}),compiledThetaComp.thetaDT{ii}(selectidx{ii}),'.')
    hold on
    lsline
    [r,p]= corrcoef(compiledThetaComp.placeDist{ii}(selectidx{ii}),compiledThetaComp.thetaDT{ii}(selectidx{ii}),'Rows','complete');
    xlabel('delta Position bins')
    ylabel('Theta CCG(ms)')
    if ii == 1
        title(strcat('Place: Short (1-3), r=', num2str(r(1,2)),', p= ',num2str(p(1,2))));
    else
        title(strcat('Place: Long (5-6), r=', num2str(r(1,2)),', p= ',num2str(p(1,2))));  
    end
end

subplot(2,4,5)
groupStats([{compiledThetaComp.thetaDT{1}(selectidx{1})} {compiledThetaComp.thetaDT{2}(selectidx{2})}],[],'inAxis',true);
xlabel({'Short','Long'})
title('Theta timescale')

subplot(2,4,6)
groupStats([{compiledThetaComp.placeDist{1}(selectidx{1})} {compiledThetaComp.placeDist{2}(selectidx{2})}],[],'inAxis',true);
xlabel({'Short','Long'})
title('Place timescale')
% 
% for ii = 1:2
%     subplot(2,4,2+ii)
%     scatter(compiledThetaComp.thetaDT{ii}(selectidx_tone),compiledThetaComp.placeDT{ii}(selectidx_tone),'.')
%     hold on
%     lsline
%     [r,p]= corrcoef(compiledThetaComp.thetaDT{ii}(selectidx_tone),compiledThetaComp.placeDT{ii}(selectidx_tone));
%     xlim([-250 250])
%     ylim([-800 800])
%     xlabel('Theta timescale(ms)')
%     ylabel('Place field timecsale(ms)')
%     if ii == 1
%         title(strcat('Tone: Short (1-3), r=', num2str(r(1,2)),', p= ',num2str(p(1,2))));
%     else
%         title(strcat('Tone: Long (5-6), r=', num2str(r(1,2)),', p= ',num2str(p(1,2))));  
%     end
% end

end
