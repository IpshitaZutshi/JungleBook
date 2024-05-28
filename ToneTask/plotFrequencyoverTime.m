function plotFrequencyoverTime(varargin)

p = inputParser;
addParameter(p,'fighandle',[]);
addParameter(p,'numrows',1,@isnumeric);
addParameter(p,'numcol',1,@isnumeric);
addParameter(p,'rowloc',1,@isnumeric);
addParameter(p,'colloc',1,@isnumeric);
addParameter(p,'correct',true,@islogical);

parse(p,varargin{:});
fighandle = p.Results.fighandle;
numrows = p.Results.numrows;
numcol = p.Results.numcol;
rowloc = p.Results.rowloc;
colloc = p.Results.colloc;
correct = p.Results.correct;

sess = {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 
    'IZ47\Final\IZ47_230626_sess15','IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28',... 
    };

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

gain = [120/11.6, 120/32.27 120/55.53 120/79.62 120/102.79 120/120];
freqExp = log10(22000/1000);

if isempty(fighandle)
    figure
else
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc], 'Parent', fighandle);
end
avgTone = [];
for ss = 1:length(sess)
    cd(strcat(expPath,sess{ss}))

    file = dir(['*.TrialBehavior.Behavior.mat']);
    load(file.name);
    file = dir(['*.Tracking.Behavior.mat']);
    load(file.name);    
    positions = [];
    for pf = 1:(size(behavTrials.timestamps,1)-1)    
        [idx] = InIntervals(tracking.timestamps,(behavTrials.timestamps(pf,:)-0.03));
        if sum(idx)>0
            positions.forward{pf} = [tracking.timestamps(idx) tracking.position.x(idx) tracking.position.y(idx)];
            if behavTrials.linTrial(pf)==1
                positions.tone{pf} = [tracking.timestamps(idx) tracking.position.x(idx)*nan tracking.position.y(idx)*nan];
            else
                y = tracking.position.y(idx);
                tonepos = [];
                for ii = 1:length(y)
                    freq = (y(ii)*gain(behavTrials.toneGain(pf)+1))/122;
                    tonepos(ii) = 1000*(10.^(freqExp*freq));
                end
                %tonepos(tonepos>25000) = nan;
                positions.tone{pf} = [tracking.timestamps(idx) tracking.position.x(idx) tonepos'];
            end
            [idx] = InIntervals(tracking.timestamps,[behavTrials.timestamps(pf,2) behavTrials.timestamps(pf+1,1)]);   
            positions.reverse{pf} = [tracking.timestamps(idx) tracking.position.x(idx) tracking.position.y(idx)];
        else
            position.forward{pf} = [];
            positions.reverse{pf} = [];
        end
    end    
    
    for ii = 1:length(positions.tone)
        if behavTrials.lickLoc(ii)>=3 && behavTrials.linTrial(ii)==0 && behavTrials.correct(ii)==correct
            plot(positions.tone{ii}(end:-1:(end-160)),'LineWidth',0.5,'Color',[0.5 0.5 0.5])
            hold on
            avgTone = [avgTone;positions.tone{ii}(end:-1:(end-160))];
        end
    end
end

plot(nanmedian(avgTone,1),'LineWidth',1.5,'Color','k')
% line([0 180],[16000 16000])
% line([0 180],[8000 8000])
% line([0 180],[4000 4000])
line([30 30],[0 25000])
line([60 60],[0 25000])
line([90 90],[0 25000])
end