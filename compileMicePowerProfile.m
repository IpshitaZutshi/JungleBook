function compileMicePowerProfile

%mice = {'IZ8\Final','IZ9\Final','IZ11\Final','IZ15','IZ12\Final','IZ13',};%If in a 'Final' Folder, include the Final in the string
mice = {'IZ27\Final','IZ28\Final','IZ29',};%If in a 'Final' Folder, include the Final in the string
%shankNum = [2,2,3,3,1,1];
shankNum = [1,1,3];
parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

nf = 1;

cmap = cbrewer('qual','Paired',6);

for m = 1:length(mice)
    
    expPath = strcat(parentDir, mice{m});
    allpath = strsplit(genpath(expPath),';'); % all folders
    cd(allpath{1});
    allSess = dir('*_sess*');
    load(strcat(allSess(1).folder,'\',allSess(1).name,'\',allSess(1).name,'.sessionInfo.mat'));
    
    channels = [1:sessionInfo.nChannels]-1;
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('PowerProfile.mat','file')
        load(strcat('PowerProfile.mat'));
    else 
        disp(['Power Profile not computed for mouse' mice{m}])
        continue;
    end
    
    for ii = 2:length(reg)
        for jj = 1:length(target)
            figure((2*(ii-1))+jj)
            set(gcf, 'Position', get(0, 'Screensize'));
            for kk = 1:length(zone)      
                if kk~=3 && kk~=6 
                    meanTheta = mean(PowerProfile.theta{ii,jj}{kk}(:,:,2:end),3);
                    meansg  = mean(PowerProfile.sg{ii,jj}{kk}(:,:,2:end),3);
                    meanmg  = mean(PowerProfile.mg{ii,jj}{kk}(:,:,2:end),3);
                    meanhfo  = mean(PowerProfile.hfo{ii,jj}{kk}(:,:,2:end),3);
                    stdTheta = (std(PowerProfile.theta{ii,jj}{kk}(:,:,2:end),[],3))/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                    stdsg  = std(PowerProfile.sg{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                    stdmg  = std(PowerProfile.mg{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                    stdhfo  = std(PowerProfile.hfo{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);

                    if kk < 4 
                        col = [0 0 0];
                    elseif (jj==1 && kk == 4) || (jj==2 && kk == 5)
                        col = [85/243 85/243 85/243];
                    elseif (jj==1 && kk == 5) || (jj==2 && kk == 4)
                        col = [8/243 133/243 161/243];
                    end
                    
                    [Lia] = ismember(sessionInfo.AnatGrps(shankNum(m)).Channels, channels);
                    nC = 1:length(sessionInfo.AnatGrps(shankNum(m)).Channels);
                    nC = nC(Lia);
                    chanImp = sessionInfo.AnatGrps(shankNum(m)).Channels(Lia)+1;

                    subplot(6,4,4*(floor((m-1)/2))+1)
                    hold on
                    dev1 = meanTheta(chanImp)-stdTheta(chanImp);
                    dev2 = meanTheta(chanImp)+stdTheta(chanImp);
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(meanTheta(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('theta power')
                    set(gca,'YDir','reverse');
%                     if m==3 || m==4
%                         ylim([0 11])
%                     end

                    subplot(6,4,4*(floor((m-1)/2))+2)
                    hold on
                    dev1 = meansg(chanImp)-stdsg(chanImp);
                    dev2 = meansg(chanImp)+stdsg(chanImp);
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(meansg(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('sg power')
                    set(gca,'YDir','reverse');
%                     if m==3 || m==4
%                         ylim([0 11])
%                     end

                    subplot(6,4,4*(floor((m-1)/2))+3)
                    hold on
                    dev1 = meanmg(chanImp)-stdmg(chanImp);
                    dev2 = meanmg(chanImp)+stdmg(chanImp);
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(meanmg(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('mg power')
                    set(gca,'YDir','reverse');
%                     if m==3 || m==4
%                         ylim([0 11])
%                     end

                    subplot(6,4,4*(floor((m-1)/2))+4)
                    hold on
                    dev1 = meanhfo(chanImp)-stdhfo(chanImp);
                    dev2 = meanhfo(chanImp)+stdhfo(chanImp);
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(meanhfo(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('hfo power')
                    set(gca,'YDir','reverse');
%                     if m==3 || m==4
%                         ylim([0 11])
%                     end
                    
                    if kk == 1 || kk ==2
                          
                        diffTheta = PowerProfile.theta{ii,jj}{kk} - PowerProfile.theta{ii,jj}{kk+3};
                        meanTheta = mean(diffTheta(:,:,2:end),3);
                        stdTheta = (std(diffTheta(:,:,2:end),[],3))/sqrt(size(diffTheta,3)-1);
                        
                        diffsg = PowerProfile.sg{ii,jj}{kk} - PowerProfile.sg{ii,jj}{kk+3};
                        meansg= mean(diffsg(:,:,2:end),3);
                        stdsg = (std(diffsg(:,:,2:end),[],3))/sqrt(size(diffsg,3)-1);
                        
                        diffmg = PowerProfile.mg{ii,jj}{kk} - PowerProfile.mg{ii,jj}{kk+3};
                        meanmg = mean(diffmg(:,:,2:end),3);
                        stdmg = (std(diffmg(:,:,2:end),[],3))/sqrt(size(diffmg,3)-1);                        
                        
                        diffhfo = PowerProfile.hfo{ii,jj}{kk} - PowerProfile.hfo{ii,jj}{kk+3};
                        meanhfo = mean(diffhfo(:,:,2:end),3);
                        stdhfo = (std(diffhfo(:,:,2:end),[],3))/sqrt(size(diffhfo,3)-1);     
                        
                        if (jj==1 && kk == 1) || (jj==2 && kk == 2)
                            col = [85/243 85/243 85/243];
                        elseif (jj==1 && kk == 2) || (jj==2 && kk == 1)
                            col = [8/243 133/243 161/243];
                        end
                        
                        subplot(6,4,4*(floor((m-1)/2))+13)
                        hold on
                        dev1 = meanTheta(chanImp)-stdTheta(chanImp);
                        dev2 = meanTheta(chanImp)+stdTheta(chanImp);
                        fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                        plot(meanTheta(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                        ylabel('Channels'); xlabel('Power'); title('theta power')
                        set(gca,'YDir','reverse');
                        if m==3 || m==4
                            ylim([0 11])
                        end                        

                        subplot(6,4,4*(floor((m-1)/2))+14)
                        hold on
                        dev1 = meansg(chanImp)-stdsg(chanImp);
                        dev2 = meansg(chanImp)+stdsg(chanImp);
                        fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                        plot(meansg(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                        ylabel('Channels'); xlabel('Power'); title('sg power')
                        set(gca,'YDir','reverse');
                        if m==3 || m==4
                            ylim([0 11])
                        end                        

                        subplot(6,4,4*(floor((m-1)/2))+15)
                        hold on
                        dev1 = meanmg(chanImp)-stdmg(chanImp);
                        dev2 = meanmg(chanImp)+stdmg(chanImp);
                        fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                        plot(meanmg(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                        ylabel('Channels'); xlabel('Power'); title('mg power')
                        set(gca,'YDir','reverse');
                        if m==3 || m==4
                            ylim([0 11])
                        end                        

                        subplot(6,4,4*(floor((m-1)/2))+16)
                        hold on
                        dev1 = meanhfo(chanImp)-stdhfo(chanImp);
                        dev2 = meanhfo(chanImp)+stdhfo(chanImp);
                        fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                        plot(meanhfo(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                        ylabel('Channels'); xlabel('Power'); title('hfo power')
                        set(gca,'YDir','reverse');
                        if m==3 || m==4
                            ylim([0 11])
                        end                        
                           
                    end
                end
            end
        end
    end
end

for ii = 2:length(reg)
    for jj = 1:length(target)
        saveas(figure((2*(ii-1))+jj),strcat(parentDir,'\Compiled\compiledPowerProfile_',reg{ii},'_',target{jj},'.png'),'png');
        saveas(figure((2*(ii-1))+jj),strcat(parentDir,'\Compiled\compiledPowerProfile_',reg{ii},'_',target{jj},'.fig'),'fig');
        saveas(figure((2*(ii-1))+jj),strcat(parentDir,'\Compiled\compiledPowerProfile_',reg{ii},'_',target{jj},'.eps'),'epsc');
    end
end
close all
end