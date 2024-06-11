function compileThetaCompression

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
    'IZ40\Final\IZ40_220707_sess16','IZ40\Final\IZ40_220714_sess18',...27
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

plotfig = 0;
compiledthetaComp.place_offset = [];
compiledthetaComp.ccgs_time_offset = [];
compiledthetaComp.pairIDs = [];
compiledthetaComp.tone_offset = [];
compiledthetaComp.tone_ccgs_time_offset = [];
compiledthetaComp.pairIDs_tone = [];


for ii = 1:length(sess)
    thetaComp = sessionThetaCompression('expPath',strcat(expPath,sess{ii}),'plotfig',false);
    compiledthetaComp.place_offset = [compiledthetaComp.place_offset; thetaComp.place_offset'];
    compiledthetaComp.ccgs_time_offset = [compiledthetaComp.ccgs_time_offset; thetaComp.ccgs_time_offset'];
    compiledthetaComp.pairIDs = [compiledthetaComp.pairIDs; thetaComp.pairIDs];
    compiledthetaComp.tone_offset = [compiledthetaComp.tone_offset; thetaComp.tone_offset'];
    compiledthetaComp.tone_ccgs_time_offset = [compiledthetaComp.tone_ccgs_time_offset; thetaComp.tone_ccgs_time_offset'];
    compiledthetaComp.pairIDs_tone = [compiledthetaComp.pairIDs_tone; thetaComp.pairIDs_tone];
end

save('Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\thetaCompression.mat','compiledthetaComp');

if plotfig
    figure
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')
    
    numBins = 30;
    subplot(1,2,1)
    % % Create a 2D histogram
    % [counts, centers] = hist3([compiledthetaComp.place_offset,compiledthetaComp.ccgs_time_offset], 'Nbins', [numBins numBins]);
    % imagesc(centers{1}, centers{2}, counts.');
    % set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
    % colorbar;
    [R,P] = corrcoef(compiledthetaComp.place_offset,compiledthetaComp.ccgs_time_offset,'rows','complete');
    scatter(compiledthetaComp.place_offset,compiledthetaComp.ccgs_time_offset,'.')
    hold on
    lsline
    title(num2str(R));
    
    subplot(1,2,2)
    scatter(compiledthetaComp.tone_offset,compiledthetaComp.tone_ccgs_time_offset,'.')
    hold on
    lsline
    [R,P] = corrcoef(compiledthetaComp.tone_offset,compiledthetaComp.tone_ccgs_time_offset,'rows','complete');
    title(num2str(R));
end
end
