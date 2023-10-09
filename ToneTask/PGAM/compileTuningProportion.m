function [Summary] = compileTuningProportion(varargin)

p = inputParser;
addParameter(p,'plotfig',true,@islogical);

parse(p,varargin{:});
plotfig = p.Results.plotfig;

sess= {'IZ39_220622_sess8','IZ39_220624_sess10','IZ39_220629_sess12',...
    'IZ39_220702_sess14','IZ39_220714_sess18',...
    'IZ39_220705_sess16','IZ39_220707_sess17',...   
    'IZ40_220705_sess15','IZ40_220707_sess16',...
    'IZ40_220708_sess17','IZ40_220714_sess18',...
    'IZ43_220826_sess2','IZ43_220828_sess4',...
    'IZ43_220830_sess6','IZ43_220901_sess8',...
    'IZ43_220911_sess9','IZ43_220913_sess11','IZ43_220919_sess14',...
    'IZ43_220915_sess13','IZ43_220920_sess15',...    
    'IZ44_220827_sess4', 'IZ44_220828_sess5',...
    'IZ44_220829_sess6','IZ44_220830_sess7',...
    'IZ44_220912_sess10','IZ44_220913_sess11','IZ44_220919_sess14',...
    'IZ44_220915_sess13','IZ44_220920_sess15',...
    'IZ47_230707_sess24','IZ47_230710_sess25','IZ47_230712_sess27',...
    'IZ48_230628_sess17','IZ48_230703_sess21',...
    'IZ48_230705_sess22','IZ48_230714_sess28',...
    }; 

basepath = 'C:\Data\PGAMAnalysis\processedData\';
Summary.sigAll = [];
Summary.mutInfoAll = [];
Summary.kernelStrengthAll = [];
Summary.sessIDAll = [];
Summary.mouseIDAll = [];

tuningDist = [];
tuningSpace = [];
tuningInfoSpace = [];
tuningInfoDist = [];

for ss = 1:length(sess)
    cd(strcat(basepath,sess{ss}))
    
    load('results_struct.mat');

    load(strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\',sess{ss}(1:4),'\Final\',sess{ss},'\',sess{ss},'.cell_metrics.cellinfo.mat'))
    
    % Get neuronID 
    for aa = 1:size(results,2)
        neuronID(aa) = results(aa).neuron+1;
    end

    numCells = max(neuronID);

    sigMat = zeros(numCells,10);
    mutInfo = nan(numCells,10);
    kernelStrength = nan(numCells,10);
    sessID(1:numCells) = ss;
    pyrID = zeros(numCells,1);
    
    if strcmp(sess{ss}(1:4),'IZ39')==1
        mouseID(1:numCells) = 1;
        Summary.mouseIDProp(ss) = 1;
    elseif strcmp(sess{ss}(1:4),'IZ40')==1
        mouseID(1:numCells) = 2;
        Summary.mouseIDProp(ss) = 2;
    elseif strcmp(sess{ss}(1:4),'IZ43')==1
        mouseID(1:numCells) = 3;
        Summary.mouseIDProp(ss) = 3;
    elseif strcmp(sess{ss}(1:4),'IZ44')==1
        mouseID(1:numCells) = 4;
        Summary.mouseIDProp(ss) = 4;
    elseif strcmp(sess{ss}(1:4),'IZ47')==1
        mouseID(1:numCells) = 5;
        Summary.mouseIDProp(ss) = 5;
    elseif strcmp(sess{ss}(1:4),'IZ48')==1
        mouseID(1:numCells) = 6;
        Summary.mouseIDProp(ss) = 6;
    end
        
    for aa = 1:numCells
        if ~strcmp(cell_metrics.putativeCellType(aa),'Pyramidal Cell')
            continue
        end
        pyrID(aa) = 1;
        idx = find(neuronID == aa);

        for id = 1:length(idx)     
            if (results(idx(id)).pval < 0.001 && ~isnan(results(idx(id)).mutual_info))
                if id == 10 % If licks, take only thosse with a positive kernel
                    if results(idx(id)).signed_kernel_strength>=0
                        sigMat(aa,id) = 1;
                    else
                        sigMat(aa,id) = 0;
                    end
                else
                    sigMat(aa,id) = 1;
                end
                mutInfo(aa,id) = results(idx(id)).mutual_info;
                kernelStrength(aa,id) = results(idx(id)).signed_kernel_strength;

                if id == 1 % if Yfwd
                    tuningSpace = [tuningSpace; results(idx(id)).raw_rate_Hz];
                    tuningInfoSpace = [tuningInfoSpace;log2((results(idx(1)).mutual_info./results(idx(4)).mutual_info))];
                elseif id == 4 % if Freq
                    tuningDist = [tuningDist; results(idx(id)).raw_rate_Hz];
                    tuningInfoDist = [tuningInfoDist;log2((results(idx(1)).mutual_info./results(idx(4)).mutual_info))];
                end
            else
                sigMat(aa,id) = 0;
                mutInfo(aa,id) = results(idx(id)).mutual_info;
                kernelStrength(aa,id) = results(idx(id)).signed_kernel_strength;
            end
        end
    end    
    
    Summary.propSigAll(ss,1:5) = nansum(sigMat(:,[ 1 2 3 4 10]))./sum(pyrID);
    
    %Also add combinations
    % y and dist
    comb1  = sigMat(:,1) & sigMat(:,4);
    Summary.propSigAll(ss,6) = nansum(comb1)./sum(pyrID);
    
    % y and licks
    comb1  = sigMat(:,1) & sigMat(:,10);
    Summary.propSigAll(ss,7) = nansum(comb1)./sum(pyrID);
    
    % dist and licks
    comb1  = sigMat(:,4) & sigMat(:,10);
    Summary.propSigAll(ss,8) = nansum(comb1)./sum(pyrID);
    
    % y, licks and dist
    comb1  = sigMat(:,1) & sigMat(:,4) & sigMat(:,10);
    Summary.propSigAll(ss,9) = nansum(comb1)./sum(pyrID);    
    
    %% Calculate the proportion of cells tuning to distance and y-position
    
    idxBoth  =  sigMat(:,1) & sigMat(:,4);
    
    idxSpaceOnly  =  sigMat(:,1) & ~sigMat(:,4);     
    idxDistOnly  =  ~sigMat(:,1) & sigMat(:,4);   
    idxSpace  = mutInfo(idxBoth,1)>mutInfo(idxBoth,4);
    idxDist  = mutInfo(idxBoth,4)>mutInfo(idxBoth,1);
    Summary.propSigBoth(ss,:) = [sum(idxSpaceOnly)./sum(pyrID) sum(idxDistOnly)./sum(pyrID) ...
        sum(idxSpace)./sum(idxBoth) sum(idxDist)./sum(idxBoth)];
    
    Summary.sigAll = [Summary.sigAll; sigMat];
    Summary.mutInfoAll = [Summary.mutInfoAll; mutInfo];
    Summary.kernelStrengthAll = [Summary.kernelStrengthAll;kernelStrength];
    Summary.sessIDAll = [Summary.sessIDAll;sessID'];
    Summary.mouseIDAll = [Summary.mouseIDAll;mouseID'];

    clear neuronID sessID mouseID    
end

if plotfig
    %% Plot A: tuning to each variables
    figure('Position', [100 100 1300 900]);
    set(gcf,'Color','w')
    set(gcf,'Renderer','painters')
    subplot(3,5,[1:4])
    hold on

    % Scatter plots
    colors = [...
        1.0000    0.7333    0.6196; ... 
        0.8510    0.3255    0.0980; ... 
        1.0000    0.8902    0.6314; ... 
        0.9294    0.6941    0.1255; ... 
        0.6667    0.5098    0.7020; ... 
        0.4941    0.1843    0.5569; ... 
        0.6980    0.8196    0.5412; ... 
        0.4667    0.6745    0.1882; ... 
        0.5216    0.6824    0.7882; ... 
        0    0.4471    0.7412; ...
        0.5 0.5 0.5;...
        0 0 0];

    c = 1;

    % Extract tuning proportions per mouse and sessions
    names = {'IZ39','IZ40','IZ43','IZ44','IZ47','IZ48'};

    for mm = 1:6
        x = 0.25;
        propSigMouse = Summary.propSigAll(Summary.mouseIDProp==mm,:);
        propSigAvg = nanmean(propSigMouse);

        color = colors(c, :);
        colorL = colors(c+1, :);
        for v = 1:size(propSigMouse,2)
            s(mm) = scatter(ones(size(propSigMouse(:,v)))*x, propSigMouse(:,v), 14, color, 'filled');
            xvals = [x-0.03, x+0.03];
            yvals = repmat(propSigAvg(v),1,2);
            plot(xvals, yvals, "-", "Color", colorL, "LineWidth", 2)
            x = x + 0.25;   
        end
        legend(s, names, 'AutoUpdate', 'off', 'FontSize', 12)
        c = c + 2;
    end

    % Boxplots
    positions = [0.35, 0.6, 0.85, 1.1, 1.35, 1.6, 1.85, 2.1 2.35];
    for v = 1:size(Summary.propSigAll, 2)
       b = boxplot(Summary.propSigAll(:,v), 'Positions', positions(v), 'PlotStyle', 'traditional', ...
            'Colors', [0.5216 0.6824 0.7882], 'Symbol', ".", 'Widths', 0.09);
        set(b, {'LineWidth'}, {1.3})
    end

    h = findobj('LineStyle','--'); 
    set(h, 'LineStyle','-');
    g = findobj(gca,'Tag','Box');
    for j = 1:length(g)
        patch(get(g(j),'XData'), get(g(j),'YData'), [0.7451 0.8392 0.9020],'FaceAlpha',.5);
    end

    xticks([0.3, 0.55, 0.8, 1.05, 1.3, 1.55, 1.8, 2.05 2.3])
    xticklabels(["y", "yRev","yNoTone","relDistStop", "licks", "y relDistStop", ...
            "y licks","licks relDistStop","y licks relDistStop"])
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel', a,'FontSize', 11);
    xlim([0.2, 2.5])
    ylim([0,1])
    box off

    %% Plot B: tuning to both relDistStop and y 
    c = 1;
    subplot(3,5,5)
    hold on
    for mm = 1:6
        x = 0.25;
        propSigMouse = Summary.propSigBoth(Summary.mouseIDProp==mm,:);
        propSigAvg = nanmean(propSigMouse);

        color = colors(c, :);
        colorL = colors(c+1, :);
        for v = 1:size(propSigMouse,2)
            s(mm) = scatter(ones(size(propSigMouse(:,v)))*x, propSigMouse(:,v), 14, color, 'filled');
            xvals = [x-0.03, x+0.03];
            yvals = repmat(propSigAvg(v),1,2);
            plot(xvals, yvals, "-", "Color", colorL, "LineWidth", 2)
            x = x + 0.25;   
        end
        c = c + 2;
    end

    % Boxplots
    positions = [0.35, 0.6, 0.85, 1.1];
    for v = 1:size(Summary.propSigBoth, 2)
       b = boxplot(Summary.propSigBoth(:,v), 'Positions', positions(v), 'PlotStyle', 'traditional', ...
            'Colors', [0.5216 0.6824 0.7882], 'Symbol', ".", 'Widths', 0.09);
        set(b, {'LineWidth'}, {1.3})
    end

    h = findobj('LineStyle','--'); 
    set(h, 'LineStyle','-');
    g = findobj(gca,'Tag','Box');
    for j = 1:length(g)
        patch(get(g(j),'XData'), get(g(j),'YData'), [0.7451 0.8392 0.9020],'FaceAlpha',.5);
    end

    xticks([0.3, 0.55, 0.8, 1.05])
    xticklabels(["y Only", "Dist only", "y > relDistStop", "relDistStop > y"])
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel', a,'FontSize', 11);
    xlim([0.2, 1.2])
    ylim([0,1])
    box off

    %% Plot C: Distribution of mutual information

    idxBoth  =  Summary.sigAll(:,1) & Summary.sigAll(:,4);
    idxSpaceOnly  = Summary.sigAll(:,1) & ~Summary.sigAll(:,4); 
    idxDistOnly = ~Summary.sigAll(:,1) & Summary.sigAll(:,4); 

    infoExtract = Summary.mutInfoAll(idxBoth,:);
    InfoBothSpace = infoExtract(:,1)>infoExtract(:,4);
    infoBothDist  = infoExtract(:,4)>infoExtract(:,1);

    subplot(3,5,6)
    dist{1} = log2(Summary.mutInfoAll(idxSpaceOnly,1)./Summary.mutInfoAll(idxSpaceOnly,4));
    dist{2} = log2(Summary.mutInfoAll(idxDistOnly,1)./Summary.mutInfoAll(idxDistOnly,4));
    nhist(dist,'probability','samebins')

    subplot(3,5,7)
    scatter(Summary.mutInfoAll(idxBoth,1),Summary.mutInfoAll(idxBoth,4),10,[0.2 0.2 0.2],'filled','MarkerFaceAlpha',0.7)
    ylim([0 3])
    xlim([0 3])
    refline(1)
    xlabel('Mut Info, y')
    ylabel('Mut Info, Dist to stop')

    subplot(3,5,8)
    data = log2((infoExtract(:,1)./infoExtract(:,4)));%./(infoExtract(:,1)+infoExtract(:,4));
    histogram(data,-3:0.1:5,'Normalization','probability')
    xlabel('Log(2) of ratio of mutual info')
    ylabel('Proportion')

    subplot(3,5,9)
    tuningSelected = tuningSpace(tuningInfoSpace>0.5,:);
    [~,idx] = max(tuningSelected,[],2);
    [~,sortidx] = sort(idx);
    hm = zscore(tuningSelected,[],2);
    imagesc(results(1).x,1:1:length(sortidx),hm(sortidx,:))
    title('Tuning to Space')

    subplot(3,5,10)
    tuningSelected = tuningDist(tuningInfoDist<-0.5,:);
    [~,idx] = max(tuningSelected,[],2);
    [~,sortidx] = sort(idx);
    hm = zscore(tuningSelected,[],2);
    imagesc(results(4).x,1:1:length(sortidx),hm(sortidx,:))
    title('Tuning to Dist')

    %% Look at lick kernels

    idxLick  =  Summary.sigAll(:,10);
    kerStrength = Summary.kernelStrengthAll(idxLick==1,5:10);

    strengthChoice = kerStrength(:,1);
    strengthSpont = mean([kerStrength(:,2) kerStrength(:,4)],2);
    strengthHome = kerStrength(:,5);


    subplot(3,5,11)
    scatter(strengthChoice, strengthHome,10,[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.7)
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])
    lsline
    xlabel('Ker Strength, Choice licks')
    ylabel('Ker Strength, Home licks')
    [R,p] = corrcoef(strengthChoice, strengthHome);
    title(strcat('R=',num2str(R(1,2)),',p=',num2str(p(1,2))))

    subplot(3,5,12)
    scatter(strengthChoice, strengthSpont,10,[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.7)
    xlim([-1 1])
    ylim([-1 1])
    lsline
    xlabel('Ker Strength, Choice licks')
    ylabel('Ker Strength, Spont licks')
    [R,p] = corrcoef(strengthChoice, strengthSpont);
    title(strcat('R=',num2str(R(1,2)),',p=',num2str(p(1,2))))


    subplot(3,5,13)
    scatter(strengthSpont, strengthHome,10,[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.7)
    xlim([-1 1])
    ylim([-1 1])
    lsline
    xlabel('Ker Strength, Spont licks')
    ylabel('Ker Strength, Home licks')
    [R,p] = corrcoef(strengthSpont, strengthHome);
    title(strcat('R=',num2str(R(1,2)),',p=',num2str(p(1,2))))
end

end