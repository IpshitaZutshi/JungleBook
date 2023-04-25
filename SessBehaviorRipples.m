function ripples = SessBehaviorRipples(varargin)

p = inputParser;
addParameter(p,'expPath',[],@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'plotfig',false,@islogical);
parse(p,varargin{:});

expPath = p.Results.expPath;
saveMat = p.Results.saveMat;
force = p.Results.force;
plotfig = p.Results.plotfig;

if ~exist('expPath') || isempty(expPath)
    expPath = uigetdir; % select folder
end

allpath = strsplit(genpath(expPath),';'); % all folders
cd(allpath{1});
allSess = dir('*_sess*');


if exist('Summ\mazeRipples.mat','file') && ~force 
    disp('Ripples already computed! Loading file.');
    load('Summ\mazeRipples.mat');
else
    for rr = 1:3
        for cc = 1:2
            for zz = 1:6
                mazeRipples.number{rr,cc}{zz} = [];
                mazeRipples.duration{rr,cc}{zz} = [];
                mazeRipples.frequency{rr,cc}{zz} = [];
                mazeRipples.peakPower{rr,cc}{zz}  = [];
                mazeRipples.timestamps{rr,cc}{zz} = [];                
            end
        end
    end
    for ii = 1:size(allSess,1)
        fprintf(' ** Examining session %3.i of %3.i... \n',ii, size(allSess,1));
        cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
        [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true);
        file = dir(('*.SessionPulses.Events.mat'));
        load(file.name);
        file = dir(('*.SessionArmChoice.Events.mat'));
        load(file.name);       
        file = dir(('*.ripples.events.mat'));
        load(file.name);    

        efields = fieldnames(sessionPulses);    

        for jj = 1:length(efields)
            region = sessionPulses.(efields{jj}).region; %1 is CA1/CA3, 2 is mec, 3 is both
            target = sessionPulses.(efields{jj}).target; %1 is stem, 2 is return

            rewardTS = sessionArmChoice.(efields{jj}).timestamps;
            startDelay = sessionArmChoice.(efields{jj}).delay.timestamps(1,:)';     
            endDelay = sessionArmChoice.(efields{jj}).delay.timestamps(2,:)';  
            
            for zz = 1:6
                %Extract relevant intervals for ripples
                switch zz
                    case 1  %First, no stim trials, return        
                        startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);
                        events = [startTS'; endTS'];
                    case 2  %No stim, stem
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
                        endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==0)+1);
                        events = [startTS';endTS'];
                    case 3 %No stim, delay
                        startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0);        
                        endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==0); 
                        events = [startTS';endTS'];  
                    case 4  % Stim, return
                        startTS = rewardTS(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        endTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);
                        events = [startTS';endTS'];                    
                    case 5   % Stim, stem
                        startTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
                        endTS = rewardTS(find(sessionPulses.(efields{jj}).stim(1:(end-1))==1)+1);
                        events = [startTS';endTS'];                      
                    case 6    %stim, delay
                        startTS = startDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1);        
                        endTS = endDelay(sessionPulses.(efields{jj}).stim(1:(end-1))==1); 
                        events = [startTS';endTS'];
                end

                if (zz == 3 || zz == 6) && sessionArmChoice.(efields{jj}).delay.dur < 1 
                    ripples.number = nan;
                    ripples.duration = nan;
                    ripples.peakPower = nan;
                    ripples.timestamps = nan;
                else
                    keepIdx = InIntervals(ripples.peaks,events');    
                end
                mazeRipples.number{region,target}{zz} = catpad(3,mazeRipples.number{region,target}{zz},sum(keepIdx));    
                mazeRipples.duration{region,target}{zz} = catpad(3,mazeRipples.duration{region,target}{zz},ripples.data.duration(keepIdx,:)); 
                mazeRipples.peakPower{region,target}{zz} = catpad(3,mazeRipples.peakPower{region,target}{zz},ripples.data.peakAmplitude(keepIdx,:)); 
                mazeRipples.frequency{region,target}{zz} = catpad(3,mazeRipples.frequency{region,target}{zz},ripples.data.peakFrequency(keepIdx,:)); 
                mazeRipples.timestamps{region,target}{zz} = catpad(3,mazeRipples.timestamps{region,target}{zz},ripples.timestamps(keepIdx,:)); 
            end        
            clear rewardTS startDelay events
        end

    end

    if saveMat
        save([expPath '\Summ\' 'mazeRipples.mat'], 'mazeRipples');
    end
end

if plotfig
    reg = {'CA3','mEC','Both'};
    zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
    target = {'STEM', 'RETURN'};

    for ii = 1:length(reg)
        figure(ii)
         for jj = 1%:length(target)
             for kk = 1:length(zone)      
                dat.number = [];
                dat.duration = [];
                dat.peakPower = [];
                for pp = 2:size(mazeRipples.number{ii,jj}{kk},3) % First array is nans
                    dat.number = [dat.number mazeRipples.number{ii,jj}{kk}(pp)]; 
                    if mazeRipples.number{ii,jj}{kk}(pp)> 0
                        dat.duration = [dat.duration mazeRipples.duration{ii,jj}{kk}(:,:,pp)']; 
                        dat.peakPower = [dat.peakPower mazeRipples.peakPower{ii,jj}{kk}(:,:,pp)']; 
                    end
                end  
                dat.duration = dat.duration(~isnan(dat.duration));
                dat.peakPower = dat.peakPower(~isnan(dat.peakPower));
                if kk < 4
                    loc = kk;
                    col = 'black';
                else
                    loc = kk-3;
                    col = 'blue';
                end
                subplot(2,9,9*(jj-1)+loc)  
                hold on
                if sum(dat.number)>0
                    histogram(dat.duration,0:0.005:0.1,'FaceColor',col,'FaceAlpha',0.4);
                    xlabel('Ripple duration')
                end

                subplot(2,9,9*(jj-1)+loc+3)   
                hold on
                if sum(dat.number)>0
                    histogram(dat.peakPower,0:0.5:15,'FaceColor',col,'FaceAlpha',0.4);
                    xlabel('Ripple power')       
                end

                subplot(2,9,9*(jj-1)+loc+6)   
                hold on       
                h2 = plot(dat.number,'o','MarkerFaceColor',col);
                h2.Color = col;
                title ('Number of ripples')

             end

        end
    %     saveas(figure(ii),strcat(expPath,'\Summ\Ripples',reg{ii},'.png'));
    %     saveas(figure(ii),strcat(expPath,'\Summ\Ripples',reg{ii},'.eps'));
    %     saveas(figure(ii),strcat(expPath,'\Summ\Ripples',reg{ii},'.fig'));

    end
end
end