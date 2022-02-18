function compileMicePowerProfileVer2

mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};%If in a 'Final' Folder, include the Final in the string
shankNum = [1,1,2];
tag = 'Bilateral';
parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM'};%, 'RETURN'};

nf = 1;

for m = 1:length(mice)
    
    expPath = strcat(parentDir, mice{m});
    allpath = strsplit(genpath(expPath),';'); % all folders
    cd(allpath{1});
    allSess = dir('*_sess*');
    load(strcat(allSess(1).folder,'\',allSess(1).name,'\',allSess(1).name,'.sessionInfo.mat'));
    load(strcat(allSess(1).folder,'\',allSess(1).name,'\',allSess(1).name,'.region.mat'));   
    Chstart  = [];
    %Find index of pyramidal channel
    pyrCh = region.CA1sp;
    for ch = 1:size(sessionInfo.AnatGrps,2)
        if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
            Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
        end
    end     
        
    channels = [1:sessionInfo.nChannels]-1;
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('PowerProfile.mat','file')
        load(strcat('PowerProfile.mat'));
    else 
        disp(['Power Profile not computed for mouse' mice{m}])
        continue;
    end
    
    for jj = 1%:length(target)
        for ii = 1:length(reg)
            figure(jj) % One plot for each manipulation
            set(gcf, 'Position', get(0, 'Screensize'));
            set(gcf,'renderer','painters')
            for kk = 1:length(zone)      

                meanTheta = nanmean(PowerProfile.theta{ii,jj}{kk}(:,:,2:end),3);
                meansg  = nanmean(PowerProfile.sg{ii,jj}{kk}(:,:,2:end),3);
                meanmg  = nanmean(PowerProfile.mg{ii,jj}{kk}(:,:,2:end),3);
                meanhfo  = nanmean(PowerProfile.hfo{ii,jj}{kk}(:,:,2:end),3);
                stdTheta = (nanstd(PowerProfile.theta{ii,jj}{kk}(:,:,2:end),[],3))/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                stdsg  = nanstd(PowerProfile.sg{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                stdmg  = nanstd(PowerProfile.mg{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                stdhfo  = nanstd(PowerProfile.hfo{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);

                if kk < 4 && ii == 1 % If CA1 baseline
                    col = [0 0 0];
                elseif kk < 4 && ii == 2 % If mEC baseline
                    col = [85/243 85/243 85/243];                   
                elseif kk< 4 && ii == 3 % If Both baseline
                    col = [160/243 160/243 160/243];       
                elseif ii==1 && kk >=4 %If CA1 stim
                    col = [224/243 163/243 46/243];                    
                elseif ii==2 && kk >=4 %If mEC stim 
                    col = [8/243 133/243 161/243];
                elseif ii==3 && kk >=4 %If both stim
                    col = [187/243 86/243 149/243];
                end

                    
                [Lia] = ismember(sessionInfo.AnatGrps(shankNum(m)).Channels, channels);
                nC = 1:length(sessionInfo.AnatGrps(shankNum(m)).Channels);
                nC = nC(Lia);
                chanImp = sessionInfo.AnatGrps(shankNum(m)).Channels(Lia)+1;
                
                chanImp = chanImp(Chstart-6:Chstart+16);
                nC = 1:1:length(chanImp);     
%                 x = 1:68/length(chanImp):68; %Determine current sampling                     
                
                subplot(6,4,4*(rem(kk,3))+1)
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meanTheta(chanImp);
                stdProfile = stdTheta(chanImp);           
%                 Meanprofile = interp1(x,Meanprofile,nC);
%                 stdProfile = interp1(x,stdProfile,nC);  
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel('Channels'); xlabel('Power'); title('theta power')
                set(gca,'YDir','reverse');

                subplot(6,4,4*(rem(kk,3))+2)
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meansg(chanImp);
                stdProfile = stdsg(chanImp);       
%                 Meanprofile = interp1(x,Meanprofile,nC);
%                 stdProfile = interp1(x,stdProfile,nC);
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel('Channels'); xlabel('Power'); title('sg power')
                set(gca,'YDir','reverse');

                subplot(6,4,4*(rem(kk,3))+3)
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meanmg(chanImp);
                stdProfile = stdmg(chanImp);           
%                 Meanprofile = interp1(x,Meanprofile,nC);
%                 stdProfile = interp1(x,stdProfile,nC);                
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel('Channels'); xlabel('Power'); title('mg power')
                set(gca,'YDir','reverse');

                subplot(6,4,4*(rem(kk,3))+4)
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meanhfo(chanImp);
                stdProfile = stdhfo(chanImp);              
%                 Meanprofile = interp1(x,Meanprofile,nC);
%                 stdProfile = interp1(x,stdProfile,nC);                
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel('Channels'); xlabel('Power'); title('hfo power')
                set(gca,'YDir','reverse');
                
                % Now plot the difference
               if kk == 1 || kk ==2

                    diffTheta = -(PowerProfile.theta{ii,jj}{kk} - PowerProfile.theta{ii,jj}{kk+3});
                    meanTheta = mean(diffTheta(:,:,2:end),3);
                    stdTheta = (std(diffTheta(:,:,2:end),[],3))/sqrt(size(diffTheta,3)-1);

                    diffsg = -(PowerProfile.sg{ii,jj}{kk} - PowerProfile.sg{ii,jj}{kk+3});
                    meansg= mean(diffsg(:,:,2:end),3);
                    stdsg = (std(diffsg(:,:,2:end),[],3))/sqrt(size(diffsg,3)-1);

                    diffmg = -(PowerProfile.mg{ii,jj}{kk} - PowerProfile.mg{ii,jj}{kk+3});
                    meanmg = mean(diffmg(:,:,2:end),3);
                    stdmg = (std(diffmg(:,:,2:end),[],3))/sqrt(size(diffmg,3)-1);                        

                    diffhfo = -(PowerProfile.hfo{ii,jj}{kk} - PowerProfile.hfo{ii,jj}{kk+3});
                    meanhfo = mean(diffhfo(:,:,2:end),3);
                    stdhfo = (std(diffhfo(:,:,2:end),[],3))/sqrt(size(diffhfo,3)-1);     

                    if ii == 1 % If CA1 baseline
                        col = [224/243 163/243 46/243];  
                    elseif ii == 2 % If mEC baseline
                        col = [8/243 133/243 161/243];                  
                    elseif ii == 3 % If Both baseline
                        col = [187/243 86/243 149/243];      
                    end
                                                    
                    subplot(6,4,4*(rem(kk,3))+13)
                    hold on
                    dev1 = meanTheta(chanImp)-stdTheta(chanImp);
                    dev2 = meanTheta(chanImp)+stdTheta(chanImp);
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(meanTheta(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('theta power')
                    set(gca,'YDir','reverse');
                    xlim([-5*10^15 5*10^15])
                    
                    subplot(6,4,4*(rem(kk,3))+14)
                    hold on
                    dev1 = meansg(chanImp)-stdsg(chanImp);
                    dev2 = meansg(chanImp)+stdsg(chanImp);
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(meansg(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('sg power')
                    set(gca,'YDir','reverse');                     
                    xlim([-7*10^14 7*10^14])

                    subplot(6,4,4*(rem(kk,3))+15)
                    hold on
                    dev1 = meanmg(chanImp)-stdmg(chanImp);
                    dev2 = meanmg(chanImp)+stdmg(chanImp);
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(meanmg(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('mg power')
                    set(gca,'YDir','reverse');  
                    xlim([-5*10^13 5*10^13])

                    subplot(6,4,4*(rem(kk,3))+16)
                    hold on
                    dev1 = meanhfo(chanImp)-stdhfo(chanImp);
                    dev2 = meanhfo(chanImp)+stdhfo(chanImp);
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(meanhfo(chanImp),nC(Lia),'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('hfo power')
                    set(gca,'YDir','reverse');  
                    xlim([-3*10^12 3*10^12])

                end


            end
        end
    end
end


for jj = 1:length(target)
    saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile',tag,'_',target{jj},'.png'),'png');
    saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile',tag,'_',target{jj},'.fig'),'fig');
    saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile',tag,'_',target{jj},'.eps'),'epsc');
end

%close all
end