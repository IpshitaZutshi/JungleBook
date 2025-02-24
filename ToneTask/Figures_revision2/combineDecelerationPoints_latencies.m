function combDec = combineDecelerationPoints_latencies(varargin)

p = inputParser;
addParameter(p,'plotfig',true,@islogical);

parse(p,varargin{:});
plotfig = p.Results.plotfig;

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

combDec.velPSTH = [];
combDec.accPSTH = [];
combDec.decType = [];
combDec.timetoLick1 = [];

combDec2.velPSTH = [];
combDec2.accPSTH = [];
combDec2.decType = [];
combDec2.timetoLick1 = [];

for ii = 1:length(sess)

    cd(strcat(expPath,sess{ii}))   
    Dec = findDeceleration_lickLatencies('plotfig',false,'error',1);
    combDec.velPSTH = [combDec.velPSTH;Dec.velPSTH];
    combDec.accPSTH = [combDec.accPSTH;Dec.accPSTH];
    combDec.timetoLick1 = [combDec.timetoLick1 Dec.timetoLick1];
    combDec.decType = [combDec.decType Dec.decType];

    Dec = findDeceleration_lickLatencies('plotfig',false,'error',0);
    combDec2.velPSTH = [combDec2.velPSTH;Dec.velPSTH];
    combDec2.accPSTH = [combDec2.accPSTH;Dec.accPSTH];
    combDec2.timetoLick1 = [combDec2.timetoLick1 Dec.timetoLick1];
    combDec2.decType = [combDec2.decType Dec.decType];

end

if plotfig
    fig2= figure;
    set(gcf,'Renderer','painters')
    tAxis = linspace(-2,5,211);
    col = {'b','k','m'};
    for ii = 1
        plotAvgStd(combDec.velPSTH(combDec.decType==ii,:),2,2,1,fig2,tAxis',col{ii})
        hold on
        plotAvgStd(combDec2.velPSTH(combDec2.decType==ii,:),2,2,1,fig2,tAxis',col{ii+1})
        title('Speed')    
        xlim([-2 5])
    
        plotAvgStd(combDec.accPSTH(combDec.decType==ii,:),2,2,3,fig2,tAxis',col{ii})  
        hold on
        plotAvgStd(combDec2.accPSTH(combDec2.decType==ii,:),2,2,3,fig2,tAxis',col{ii+1})
        xlim([-2 5])
        title('Acceleration')
    end
    
    subplot(2,2,2)
    histogram(combDec.timetoLick1(combDec.decType==1))
    % scatter(combDec.timetoLick1(combDec.decType==1), combDec.posY(combDec.decType==1), 5,'b.')
    % hold on
    % scatter(combDec.timetoLick1(combDec.decType==2), combDec.posY(combDec.decType==2), 5,'k.')
    % scatter(combDec.timetoLick1(combDec.decType==3), combDec.posY(combDec.decType==3), 5,'m.')
    xlabel('Latency between lick and max deceleration')
    xlim([-2 4])
    % ylabel('Position')
    
    % subplot(2,2,4)
    % scatter(combDec.posX(combDec.decType==1), combDec.posY(combDec.decType==1), 5,'b.')
    % hold on
    % scatter(combDec.posX(combDec.decType==2), combDec.posY(combDec.decType==2), 5,'k.')
    % scatter(combDec.posX(combDec.decType==3), combDec.posY(combDec.decType==3), 5,'m.')
    % xlabel('x position')
    % ylabel('y position')
end

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\Revision2\';
saveas(gcf,strcat(expPath,'decelerationLatency.png'));
saveas(gcf,strcat(expPath,'decelerationLatency.eps'),'epsc');
saveas(gcf,strcat(expPath,'decelerationLatency.fig'));

end

function plotAvgStd(array,numrows,numcol,subplotlocation,figureHandle,xAxis,col)

    subplot(numrows, numcol, subplotlocation, 'Parent', figureHandle);

    meanpsth = nanmean(array,1);
    stdpsth = nanstd(array,1)./sqrt(size(array,1));
    lArr  = meanpsth-stdpsth;
    uArr = meanpsth+stdpsth;

    fill([xAxis; flipud(xAxis)],[lArr'; flipud(uArr')],col,'linestyle','none','FaceAlpha',0.5);                    
    hold on
    hi = line(xAxis,meanpsth,'LineWidth',1,'Color',col);

end