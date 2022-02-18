function compileMiceBehavVars(varargin)

%% Defaults and parms
p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
addParameter(p,'savePlot',true,@islogical);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
savePlot = p.Results.savePlot;

tag = 'CA3Saline'; 

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ21\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ11\Final','IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final',...
        'IZ21\Final','IZ24\Final', 'IZ25\Final', 'IZ26\Final','IZ27\Final','IZ28\Final',...
        'IZ29\Final','IZ30\Final','IZ31\Final','IZ32\Final','IZ33\Final','IZ34\Final'};  % To add: IZ23
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final'}; % To add: IZ35
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'}; % To add: IZ35
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    %
    reg = {'contramEC','ipsimEC','Both'};
end

tar = {'STEM', 'RETURN'};

for ii = 1:length(reg)
    for jj = 1:length(tar)
        behavVars{ii,jj}.time = [];
        behavVars{ii,jj}.mouse = [];
    end
end


for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m}));
    allSess = dir('*_sess*');
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        file = dir(['*.SessionArmChoice.Events.mat']);
        load(file.name);    
        file = dir(['*.SessionPulses.Events.mat']);
        load(file.name);
        efields = fieldnames(sessionPulses);
        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region;
            target = sessionPulses.(efields{jj}).target;
            timeStamps = diff(sessionArmChoice.(efields{jj}).timestamps);
            stimidx = logical(sessionPulses.(efields{jj}).stim(1:(end-1)));            
            % Average time of non stim trials
            TSnonStim = timeStamps(~stimidx);           
            % Average time of stim strials
            TSStim = timeStamps(stimidx);                
            behavVars{region,target}.time = [behavVars{region,target}.time;[median(TSnonStim) median(TSStim)]];            
            behavVars{region,target}.mouse = [behavVars{region,target}.mouse;m];
        end
    end
end


reg = {'CA1','mEC','Both'};
zone = {'Stem','Return'};


figure
set(gcf,'Renderer','painters')

if strcmp(tag,'CA3')==1 || strcmp(tag,'CA3Saline')==1
    set(gcf,'Position',[100 100 400 200])
    for jj = 1:size(behavVars,2)
        if ~isempty(behavVars{ii,jj}.time)
            subplot(1,2,jj)
            data = [];
            cond = [];        
            data = [behavVars{2,jj}.time(1:end,1);behavVars{2,jj}.time(1:end,2);...
                behavVars{3,jj}.time(1:end,1);behavVars{3,jj}.time(1:end,2)];
            cond(1:size(behavVars{2,jj}.time,1)) = 1;
            cond(size(behavVars{2,jj}.time,1)+1:2*size(behavVars{2,jj}.time,1)) = 2;
            cond(2*size(behavVars{2,jj}.time,1)+1:3*size(behavVars{2,jj}.time,1)) = 3;
            cond(3*size(behavVars{2,jj}.time,1)+1:4*size(behavVars{2,jj}.time,1)) = 4;
            
            if size(data,1)>3
                boxchart([behavVars{2,jj}.time behavVars{3,jj}.time])
                hold on
                plot([behavVars{2,jj}.time behavVars{3,jj}.time]','Color',[0.5 0.5 0.5]);  
                stats{jj} = groupStats(data,cond,'doPlot',false,'repeatedMeasures',true);
                title(strcat(zone{jj},' ',num2str(stats{jj}.friedman.p)));

            end
        end

    end    
else
    set(gcf,'Position',[10 10 400 800])
    for ii = 1:size(behavVars,1)
        for jj = 1:size(behavVars,2)
            if ~isempty(behavVars{ii,jj}.time)
                subplot(3,2,2*(ii-1)+jj)
                data = [];
                cond = [];        
                data = [behavVars{ii,jj}.time(1:end,1);behavVars{ii,jj}.time(1:end,2)];
                cond(1:size(behavVars{ii,jj}.time,1)) = 1;
                cond(size(behavVars{ii,jj}.time,1)+1:2*size(behavVars{ii,jj}.time,1)) = 2;
                if size(data,1)>3
                    boxchart(behavVars{ii,jj}.time)
                    hold on
                    plot(behavVars{ii,jj}.time','Color',[0.5 0.5 0.5]);  
                    stats{ii,jj} = groupStats(data,cond,'doPlot',false);
                    [stats{ii,jj}.signrank.p,~,stats{ii,jj}.signrank.stats]  = signrank(behavVars{ii,jj}.time(1:end,1),behavVars{ii,jj}.time(1:end,2));
                    title(strcat(reg{ii},' ',zone{jj},' ',num2str(stats{ii,jj}.signrank.p)));

                end
            end

        end
    end
end
if savePlot
    saveas(gcf,strcat(parentDir,'\Compiled\Behavior\BehaviorVars\',tag,'.png'));
    saveas(gcf,strcat(parentDir,'\Compiled\Behavior\BehaviorVars\',tag,'.fig'));
    saveas(gcf,strcat(parentDir,'\Compiled\Behavior\BehaviorVars\',tag,'.eps'),'epsc');
    save(strcat(parentDir,'\Compiled\Behavior\BehaviorVars\',tag,'Stats.mat'),'stats');
end

end