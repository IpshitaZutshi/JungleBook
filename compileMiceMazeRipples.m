function compiledMiceRipples = compileMiceMazeRipples

tag = 'mEC';% mEC, CA3, Bilateral mEC

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','CA1Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final',...
         'IZ18\Final','IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
         'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'}; % To add, IZ16, IZ23
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ29\Final','IZ32\Final','IZ33\Final','IZ34\Final','IZ27\Final','IZ28\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ33\Saline','IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ34\Saline'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Bilateral'};
end


parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

numAnalog = 2;

cmap = cbrewer('qual','Pastel2',length(mice));

compiledMiceRipples.mice = mice;

for ii = 1:(numAnalog+1)
    for jj = 1:2   
        for kk = 1:6
            compiledMiceRipples.number{ii,jj}{kk} = [];
            compiledMiceRipples.duration{ii,jj}{kk} = [];
            compiledMiceRipples.frequency{ii,jj}{kk} = [];
            compiledMiceRipples.peakPower{ii,jj}{kk} = [];
        end
    end
end

for m = 1:length(mice)
    cd(strcat(parentDir, mice{m},'\Summ'));  
    if exist('mazeRipples.mat','file') 
        disp(['Loading ripples for mouse' mice{m}])
        load('mazeRipples.mat');
    else 
        continue
    end
    
    for ii = 1:(numAnalog+1)
        for jj = 1:2
            for kk = 1:6 
                for p = 1:size(mazeRipples.frequency{ii,jj}{kk},3)
                    compiledMiceRipples.number{ii,jj}{kk} = [compiledMiceRipples.number{ii,jj}{kk};mazeRipples.number{ii,jj}{kk}(:,:,p)];
                    compiledMiceRipples.duration{ii,jj}{kk} = [compiledMiceRipples.duration{ii,jj}{kk};mazeRipples.duration{ii,jj}{kk}(:,:,p)];
                    compiledMiceRipples.frequency{ii,jj}{kk} = [compiledMiceRipples.frequency{ii,jj}{kk};mazeRipples.frequency{ii,jj}{kk}(:,:,p)];
                    compiledMiceRipples.peakPower{ii,jj}{kk} = [compiledMiceRipples.peakPower{ii,jj}{kk};mazeRipples.peakPower{ii,jj}{kk}(:,:,p)*0.195];
                end                    
            end          
        end
    end
end
 

colMat = [0.5 0.5 0.5;
    8/243 133/243 161/243];
target = {'Stem','Return'};
zone = {'Side','Center','Delay'};

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
    figure    
    set(gcf,'Renderer','painters')
    set(gcf,'Color','w')
    colMat = [0.5 0.5 0.5;8/243 133/243 161/243;193/243 90/243 99/243; 133/243 128/243 177/243];
        
    subplot(2,3,1)
    dataRipB = compiledMiceRipples.number{2,1}{1};
    dataRipmEC = compiledMiceRipples.number{2,1}{4};
    dataRipCA3 = compiledMiceRipples.number{3,1}{1};
    dataRipBoth = compiledMiceRipples.number{3,1}{4};
    dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipmEC,1),1)*2;ones(size(dataRipCA3,1),1)*3;ones(size(dataRipBoth,1),1)*4];
    assembly = [dataRipB;dataRipmEC;dataRipCA3;dataRipBoth];
    stats{1}.side = groupStats(assembly, dataID,'color',colMat,'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesSEM');
    title('Center arm mEC, side arm ripple')
    
    subplot(2,3,2)
    dataRipB = compiledMiceRipples.number{2,1}{2};
    dataRipmEC = compiledMiceRipples.number{2,1}{5};
    dataRipCA3 = compiledMiceRipples.number{3,1}{2};
    dataRipBoth = compiledMiceRipples.number{3,1}{5};
    dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipmEC,1),1)*2;ones(size(dataRipCA3,1),1)*3;ones(size(dataRipBoth,1),1)*4];
    assembly = [dataRipB;dataRipmEC;dataRipCA3;dataRipBoth];
    stats{1}.center = groupStats(assembly, dataID,'color',colMat,'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesSEM');
    title('Center arm mEC, center arm ripple')
    
    subplot(2,3,3)
    dataRipB = compiledMiceRipples.number{2,1}{3};
    dataRipmEC = compiledMiceRipples.number{2,1}{6};
    dataRipCA3 = compiledMiceRipples.number{3,1}{3};
    dataRipBoth = compiledMiceRipples.number{3,1}{6};
    dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipmEC,1),1)*2;ones(size(dataRipCA3,1),1)*3;ones(size(dataRipBoth,1),1)*4];
    assembly = [dataRipB;dataRipmEC;dataRipCA3;dataRipBoth];
    stats{1}.delay = groupStats(assembly, dataID,'color',colMat,'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesSEM');    
    title('Center arm mEC, delay ripples')
 
     subplot(2,3,4)
    dataRipB = compiledMiceRipples.number{2,2}{1};
    dataRipmEC = compiledMiceRipples.number{2,2}{4};
    dataRipCA3 = compiledMiceRipples.number{3,2}{1};
    dataRipBoth = compiledMiceRipples.number{3,2}{4};
    dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipmEC,1),1)*2;ones(size(dataRipCA3,1),1)*3;ones(size(dataRipBoth,1),1)*4];
    assembly = [dataRipB;dataRipmEC;dataRipCA3;dataRipBoth];
    stats{2}.side = groupStats(assembly, dataID,'color',colMat,'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesSEM');
    title('Side arm mEC, side arm ripple')
    
    subplot(2,3,5)
    dataRipB = compiledMiceRipples.number{2,2}{2};
    dataRipmEC = compiledMiceRipples.number{2,2}{5};
    dataRipCA3 = compiledMiceRipples.number{3,2}{2};
    dataRipBoth = compiledMiceRipples.number{3,2}{5};
    dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipmEC,1),1)*2;ones(size(dataRipCA3,1),1)*3;ones(size(dataRipBoth,1),1)*4];
    assembly = [dataRipB;dataRipmEC;dataRipCA3;dataRipBoth];
    stats{2}.center = groupStats(assembly, dataID,'color',colMat,'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesSEM');
    title('Side arm mEC, center arm ripple')
    
    subplot(2,3,6)
    dataRipB = compiledMiceRipples.number{2,2}{3};
    dataRipmEC = compiledMiceRipples.number{2,2}{6};
    dataRipCA3 = compiledMiceRipples.number{3,2}{3};
    dataRipBoth = compiledMiceRipples.number{3,2}{6};
    dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipmEC,1),1)*2;ones(size(dataRipCA3,1),1)*3;ones(size(dataRipBoth,1),1)*4];
    assembly = [dataRipB;dataRipmEC;dataRipCA3;dataRipBoth];
    stats{2}.delay = groupStats(assembly, dataID,'color',colMat,'inAxis',true,'repeatedMeasures',true,'plotType','BoxLinesSEM');    
    title('Side arm mEC, delay ripples')
    
else
    for ii = 2 
        figure    
        set(gcf,'Renderer','painters')
        set(gcf,'Color','w')
        set(gcf,'Position',[1900 40 1900 970])
        
        for jj = 1:2 % Stim at stem, Stim at return
            for kk = 1:3
                subplot(4,6,3*(jj-1)+kk)
                % kk= 1, 2, 3 - Return, Stem, Delay zone. 
                dataRipB = compiledMiceRipples.number{ii,jj}{kk}; 
                dataRipS = compiledMiceRipples.number{ii,jj}{kk+3}; 
                dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipS,1),1)*2];
                assembly = [dataRipB;dataRipS];
                stats.number{ii,jj}{kk} = groupStats(assembly, dataID,'inAxis',true,'color',colMat,'repeatedMeasures',true,'plotType','BoxLinesSEM');
                [stats.number{ii,jj}{kk}.signrank.p,~,stats.number{ii,jj}{kk}.signrank.stats] = signrank(dataRipB,dataRipS); 
                stats.number{ii,jj}{kk}.n = [sum(~isnan(dataRipB)) sum(~isnan(dataRipS))];
                ylabel('Number of ripples')
                ylim([0 15])
                title(strcat(target{jj},' zone:',zone{kk}))           
                
                subplot(4,6,3*(jj-1)+kk+6)
                dataRipB = compiledMiceRipples.duration{ii,jj}{kk}; 
                dataRipS = compiledMiceRipples.duration{ii,jj}{kk+3}; 
                dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipS,1),1)*2];
                assembly = [dataRipB;dataRipS];
                if ~isempty(dataRipB) && ~isempty(dataRipB)
                    stats.duration{ii,jj}{kk} = groupStats(assembly, dataID,'inAxis',true,'color',colMat);
                    [stats.duration{ii,jj}{kk}.signrank.p,~,stats.duration{ii,jj}{kk}.signrank.stats] = ranksum(dataRipB,dataRipS); 
                    stats.duration{ii,jj}{kk}.n = [sum(~isnan(dataRipB)) sum(~isnan(dataRipS))];
                    ylabel('Ripple duration')
                    title(strcat(target{jj},' zone:',zone{kk}))           

                    subplot(4,6,3*(jj-1)+kk+12)
                    dataRipB = compiledMiceRipples.frequency{ii,jj}{kk}; 
                    dataRipS = compiledMiceRipples.frequency{ii,jj}{kk+3}; 
                    dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipS,1),1)*2];
                    assembly = [dataRipB;dataRipS];
                    stats.frequency{ii,jj}{kk} = groupStats(assembly, dataID,'inAxis',true,'color',colMat);
                    [stats.frequency{ii,jj}{kk}.signrank.p,~,stats.frequency{ii,jj}{kk}.signrank.stats] = ranksum(dataRipB,dataRipS); 
                    stats.frequency{ii,jj}{kk}.n = [sum(~isnan(dataRipB)) sum(~isnan(dataRipS))];
                    ylabel('Ripple frequency')
                    title(strcat(target{jj},' zone:',zone{kk}))       

                    subplot(4,6,3*(jj-1)+kk+18)
                    dataRipB = compiledMiceRipples.peakPower{ii,jj}{kk}; 
                    dataRipS = compiledMiceRipples.peakPower{ii,jj}{kk+3}; 
                    dataID = [ones(size(dataRipB,1),1)*1;ones(size(dataRipS,1),1)*2];
                    assembly = [dataRipB;dataRipS];
                    stats.amplitude{ii,jj}{kk} = groupStats(assembly, dataID,'inAxis',true,'color',colMat);
                    [stats.amplitude{ii,jj}{kk}.signrank.p,~,stats.amplitude{ii,jj}{kk}.signrank.stats] = ranksum(dataRipB,dataRipS); 
                    stats.amplitude{ii,jj}{kk}.n = [sum(~isnan(dataRipB)) sum(~isnan(dataRipS))];
                    ylabel('Ripple amplitude')
                    title(strcat(target{jj},' zone:',zone{kk}))                       
                end
                
            end
        end
    end
end

saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesMaze',tag,'.png'));
saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesMaze',tag,'.fig'));
saveas(gca,strcat(parentDir,'Compiled\Ripples\compiledRipplesMaze',tag,'.eps'),'epsc');
save(strcat(parentDir,'Compiled\Ripples\compiledRipplesMaze',tag,'.mat'),'stats')   

end