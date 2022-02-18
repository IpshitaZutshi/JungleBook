function plotCoherence(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;

mice = {'IZ12\Final','IZ13\Final','IZ18\Final','IZ20\Final','IZ24\Final','IZ25\Final','IZ26\Final',...
    'IZ21\Final','IZ31\Final','IZ27\Final','IZ28\Final','IZ33\Final'};  

figure(1)
figure(2)
idxloc = 1;
for m = 1:length(mice)

    cd(strcat(parentDir, mice{m},'\Summ'));
    load('csdcohData.mat');

    if strcmp(mice{m},'IZ27\Final')==1 || strcmp(mice{m},'IZ28\Final')==1 || strcmp(mice{m},'IZ33\Final')==1
        figure(1)
        subplot(3,4,m)
        plot(nanmean(csdcohData.cohTheta{2,1}{2},3),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'k','LineWidth',1.5)
        hold on
        plot(nanmean(csdcohData.cohTheta{2,1}{5},3),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'b','LineWidth',1.5)
        plot(nanmean(csdcohData.cohTheta{3,1}{2},3),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'m','LineWidth',1.5)
        plot(nanmean(csdcohData.cohTheta{3,1}{5},3),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'c','LineWidth',1.5)        
        ylim([1 64])
        xlim([0.5 1])
        title(mice{m})               
      
        for pp = 1:size(csdcohData.cohTheta{2,1}{2},3)
            [~,idxBase(idxloc,1)] = min(csdcohData.cohTheta{2,1}{2}(1,1:25,pp));
            [~,idxBase(idxloc,2)] = min(csdcohData.cohTheta{2,1}{5}(1,1:25,pp));
            [~,idxBase(idxloc,3)] = min(csdcohData.cohTheta{3,1}{2}(1,1:25,pp));
            idxloc = idxloc+1;
        end
        
        figure(2)
        subplot(3,4,m)
        plot((rad2deg(csdcohData.cohPhaseTheta{2,1}{2}(:,:,1))),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'k.')
        hold on
        plot((rad2deg(csdcohData.cohPhaseTheta{2,1}{2}(:,:,1)))+180,size(csdcohData.cohTheta{2,1}{2},2):-1:1,'k.')
        plot((rad2deg(csdcohData.cohPhaseTheta{2,1}{5}(:,:,1))),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'b.')
        plot((rad2deg(csdcohData.cohPhaseTheta{2,1}{5}(:,:,1)))+180,size(csdcohData.cohTheta{2,1}{2},2):-1:1,'b.')
        plot((rad2deg(csdcohData.cohPhaseTheta{3,1}{2}(:,:,1))),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'m.')
        plot((rad2deg(csdcohData.cohPhaseTheta{3,1}{2}(:,:,1)))+180,size(csdcohData.cohTheta{2,1}{2},2):-1:1,'m.')
        plot((rad2deg(csdcohData.cohPhaseTheta{3,1}{5}(:,:,1))),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'c.')      
        plot((rad2deg(csdcohData.cohPhaseTheta{3,1}{5}(:,:,1)))+180,size(csdcohData.cohTheta{2,1}{2},2):-1:1,'c.')    
        ylim([1 size(csdcohData.cohTheta{2,1}{2},2)])
        xlim([-200 400])
        title(mice{m})        
    else
        figure(1)
        subplot(3,4,m)        
        plot(nanmean(csdcohData.cohTheta{2,1}{2},3),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'k','LineWidth',1.5)
        hold on
        plot(nanmean(csdcohData.cohTheta{2,1}{5},3),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'b','LineWidth',1.5)
        ylim([1 size(csdcohData.cohTheta{2,1}{2},2)])
        xlim([0.5 1])
        title(mice{m})
        
        figure(2)
        subplot(3,4,m)
        plot((rad2deg(csdcohData.cohPhaseTheta{2,1}{2}(:,:,1))),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'k.')
        hold on
        plot((rad2deg(csdcohData.cohPhaseTheta{2,1}{2}(:,:,1)))+180,size(csdcohData.cohTheta{2,1}{2},2):-1:1,'k.')
        plot((rad2deg(csdcohData.cohPhaseTheta{2,1}{5}(:,:,1))),size(csdcohData.cohTheta{2,1}{2},2):-1:1,'b.')   
        plot((rad2deg(csdcohData.cohPhaseTheta{2,1}{5}(:,:,1)))+180,size(csdcohData.cohTheta{2,1}{2},2):-1:1,'b.')  
        ylim([1 size(csdcohData.cohTheta{2,1}{2},2)])
        xlim([-200 400])
        title(mice{m})   
        
        for pp = 1:size(csdcohData.cohTheta{2,1}{2},3)
            [~,idxBase(idxloc,1)] = min(csdcohData.cohTheta{2,1}{2}(1,1:25,pp));
            [~,idxBase(idxloc,2)] = min(csdcohData.cohTheta{2,1}{5}(1,1:25,pp));
            idxBase(idxloc,3) = nan;
            idxloc = idxloc+1;
        end
    end
end

idxStim(:,1) = idxBase(:,1)-idxBase(:,2);
idxStim(:,2) = idxBase(:,1)-idxBase(:,3);
idxStim = idxStim*20;
figure
boxplot(idxStim)
hold on
scatter(ones*0.8, idxStim(:,1),'k','.','jitter','on')
[stats.mEC.p,~,stats.mEC.stats] = signrank(idxStim(:,1),0);
stats.mEC.n = size(idxStim,1);

scatter(ones*1.8, idxStim(:,2),'k','.','jitter','on')
[stats.CA3.p,~,stats.CA3.stats] = signrank(idxStim(:,2),0);
stats.CA3.n = sum(~isnan(idxStim(:,2)));

title(strcat('mEC  ',num2str(stats.mEC.p),'  CA3   ',num2str(stats.CA3.p)))

end