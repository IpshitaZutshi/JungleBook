function compileMiceDSRipples(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
 
%mice{1} = {'IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
mice{1} = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final'...
         'IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
         'IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};
mice{2} = {'IZ25\Final','IZ26\Final','IZ24\Final'};  
%mice{2} = {'IZ15\Final','IZ30\Final','IZ31\Final'};  
for ii = 1:4
    compiledDS{ii,1} = [];
    compiledDS{ii,2} = [];
end
for tag = 1:2
    if tag==1
        range = 2;
        rangetag = 1;
    elseif tag==2
        range = [1 2 3];
        rangetag = [2 3 4];
    end

    for m = 1:length(mice{tag})       
        cd(strcat(parentDir, mice{tag}{m},'\Summ'));
        if exist(strcat('Summ\DSRipples.mat'),'file')
            load(strcat('Summ\DSRipples.mat'));
        else 
            disp(['Assembly not computed for mouse' mice{tag}{m}])
            continue;
        end

        for ii = 1:length(range)
            compiledDS{rangetag(ii),1} = [compiledDS{rangetag(ii),1}; DSRipples{range(ii),1}];
            compiledDS{rangetag(ii),2} = [compiledDS{rangetag(ii),2}; DSRipples{range(ii),2}];
        end
    end
end

colMat = [0.5 0.5 0.5;
    8/243 133/243 161/243];

figure
set(gcf,'Renderer','painters')
set(gcf,'Color','w')
set(gcf,'Position',[2000 500 1500 345])

for ii = 1:4
    subplot(1,4,ii)
    avgerr = nanmean(compiledDS{ii,1});
    stderr = nanstd(compiledDS{ii,1})./sqrt(size(compiledDS{ii,1},2));
    fill([(1:1:length(avgerr))'; (length(avgerr):-1:1)'],[avgerr'-stderr';flipud(avgerr'+stderr')],colMat(1,:),'linestyle','none','FaceAlpha',0.2);
    hold on
    plot(avgerr,'Color',colMat(1,:),'LineWidth',1.5)    

    avgerr = nanmean(compiledDS{ii,2});
    stderr = nanstd(compiledDS{ii,2})./sqrt(size(compiledDS{ii,2},2));
    fill([(1:1:length(avgerr))'; (length(avgerr):-1:1)'],[avgerr'-stderr';flipud(avgerr'+stderr')],colMat(2,:),'linestyle','none','FaceAlpha',0.2);
    hold on
    plot(avgerr,'Color',colMat(2,:),'LineWidth',1.5)     
    ylim([0.01 0.07])
end
    
saveas(gcf,strcat(parentDir,'Compiled\Ripples\DownState\RippleHpcDS.png'));
saveas(gcf,strcat(parentDir,'Compiled\Ripples\DownState\RippleHpcDS.eps'),'epsc');
saveas(gcf,strcat(parentDir,'Compiled\Ripples\DownState\RippleHpcDS.fig'));
%save(strcat(parentDir,'Compiled\Ripples\Assemblies\RippleHpcDS.mat'),'stats');
end
