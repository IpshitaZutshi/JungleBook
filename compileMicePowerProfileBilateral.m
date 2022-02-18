function PD = compileMicePowerProfileBilateral

mice = {'IZ24\Final'}%,'IZ25\Final','IZ26\Final','IZ26\Final'};%If in a 'Final' Folder, include the Final in the string
shankNum = [2,1,1,3];
tag = 'Bilateral';
parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM'};%, 'RETURN'};

nf = 1;
for ii = 1:length(reg)
    for kk = 1:2
        PD.Theta{ii,kk} = [];
        PD.sg{ii,kk} = [];
        PD.mg{ii,kk} = [];
        PD.hfo{ii,kk} = [];
    end
end

for m = 1:length(mice)
    
    expPath = strcat(parentDir, mice{m});
    allpath = strsplit(genpath(expPath),';'); % all folders
    cd(allpath{1});
    allSess = dir('*_sess*');
    load(strcat(allSess(1).folder,'\',allSess(1).name,'\',allSess(1).name,'.sessionInfo.mat'));
    load(strcat(allSess(1).folder,'\',allSess(1).name,'\',allSess(1).name,'.region.mat'));   
    Chstart  = [];
    %Find index of pyramidal channel
    if m ==3
        pyrCh = region.CA1so;
    else
        pyrCh = region.CA1sp;
    end
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
        meanThetaBase = [];
        meansgBase  = [];
        meanmgBase  = [];
        meanhfoBase = [];
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

                   if strcmp(mice{m},'IZ24\Final') == 1 && ii ==3
                        diffTheta = -(PowerProfile.theta{ii,jj}{kk} - PowerProfile.theta{ii,jj}{kk+3})./PowerProfile.theta{ii,jj}{kk};
                        meanTheta = nanmean(diffTheta(:,:,2),3);
                        stdTheta = (nanstd(diffTheta(:,:,2),[],3))/sqrt(size(diffTheta,3)-2);

                        diffsg = -(PowerProfile.sg{ii,jj}{kk} - PowerProfile.sg{ii,jj}{kk+3})./PowerProfile.sg{ii,jj}{kk};
                        meansg= nanmean(diffsg(:,:,2),3);
                        stdsg = (nanstd(diffsg(:,:,2),[],3))/sqrt(size(diffsg,3)-2);

                        diffmg = -(PowerProfile.mg{ii,jj}{kk} - PowerProfile.mg{ii,jj}{kk+3})./PowerProfile.mg{ii,jj}{kk};
                        meanmg = nanmean(diffmg(:,:,2),3);
                        stdmg = (nanstd(diffmg(:,:,2),[],3))/sqrt(size(diffmg,3)-2);                        

                        diffhfo = -(PowerProfile.hfo{ii,jj}{kk} - PowerProfile.hfo{ii,jj}{kk+3})./PowerProfile.hfo{ii,jj}{kk};
                        meanhfo = nanmean(diffhfo(:,:,2),3);
                        stdhfo = (nanstd(diffhfo(:,:,2),[],3))/sqrt(size(diffhfo,3)-2); 
                   else
                        diffTheta = -(PowerProfile.theta{ii,jj}{kk} - PowerProfile.theta{ii,jj}{kk+3})./PowerProfile.theta{ii,jj}{kk};
                        meanTheta = nanmean(diffTheta(:,:,2:end),3);
                        stdTheta = (nanstd(diffTheta(:,:,2:end),[],3))/sqrt(size(diffTheta,3)-1);

                        diffsg = -(PowerProfile.sg{ii,jj}{kk} - PowerProfile.sg{ii,jj}{kk+3})./PowerProfile.sg{ii,jj}{kk};
                        meansg= nanmean(diffsg(:,:,2:end),3);
                        stdsg = (nanstd(diffsg(:,:,2:end),[],3))/sqrt(size(diffsg,3)-1);

                        diffmg = -(PowerProfile.mg{ii,jj}{kk} - PowerProfile.mg{ii,jj}{kk+3})./PowerProfile.mg{ii,jj}{kk};
                        meanmg = nanmean(diffmg(:,:,2:end),3);
                        stdmg = (nanstd(diffmg(:,:,2:end),[],3))/sqrt(size(diffmg,3)-1);                        

                        diffhfo = -(PowerProfile.hfo{ii,jj}{kk} - PowerProfile.hfo{ii,jj}{kk+3})./PowerProfile.hfo{ii,jj}{kk};
                        meanhfo = nanmean(diffhfo(:,:,2:end),3);
                        stdhfo = (nanstd(diffhfo(:,:,2:end),[],3))/sqrt(size(diffhfo,3)-1);     
                   end

                   if m < 4
                        if ii == 1 % If CA1 baseline
                            col = [224/243 163/243 46/243];  
                        elseif ii == 2 % If mEC baseline
                            col = [8/243 133/243 161/243];                  
                        elseif ii == 3 % If Both baseline
                            col = [187/243 86/243 149/243];      
                        end
                        stimloc = ii;
                    else
                        if ii == 2 % If CA1 baseline
                            col = [224/243 163/243 46/243]; 
                            stimloc = 1;
                        elseif ii == 1 % If mEC baseline
                            col = [8/243 133/243 161/243]; 
                            stimloc = 2;
                        elseif ii == 3 % If Both baseline
                            col = [187/243 86/243 149/243];    
                            stimloc = 3;
                        end
                    end
                    chanImp = chanImp(Chstart-6:Chstart+16);
                    nC = 1:1:length(chanImp);  
                
                    subplot(6,4,4*(rem(kk,3))+13)
                    hold on
                    Meanprofile = meanTheta(chanImp);
                    Meanprofile = smooth(Meanprofile);
                    stdProfile = stdTheta(chanImp); 
                    stdProfile = smooth(stdProfile);
                    PD.Theta{stimloc,kk} = catpad(2,PD.Theta{stimloc,kk},Meanprofile);                    
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('theta power')
                    set(gca,'YDir','reverse');
                    xlim([-1 1])
                    
                    subplot(6,4,4*(rem(kk,3))+14)
                    hold on
                    Meanprofile = meansg(chanImp);
                    Meanprofile = smooth(Meanprofile);
                    stdProfile = stdsg(chanImp); 
                    stdProfile = smooth(stdProfile);
                    PD.sg{stimloc,kk} = catpad(2,PD.sg{stimloc,kk},Meanprofile);
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('sg power')
                    set(gca,'YDir','reverse');                     
                    xlim([-1 1])

                    subplot(6,4,4*(rem(kk,3))+15)
                    hold on
                    Meanprofile = meanmg(chanImp);
                    Meanprofile = smooth(Meanprofile);
                    stdProfile = stdmg(chanImp); 
                    stdProfile = smooth(stdProfile);
                    PD.mg{stimloc,kk} = catpad(2,PD.mg{stimloc,kk},Meanprofile); 
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('mg power')
                    set(gca,'YDir','reverse');  
                    xlim([-1 1])

                    subplot(6,4,4*(rem(kk,3))+16)
                    hold on
                    Meanprofile = meanhfo(chanImp);
                    Meanprofile = smooth(Meanprofile);
                    stdProfile = stdhfo(chanImp); 
                    stdProfile = smooth(stdProfile);
                    PD.hfo{stimloc,kk} = catpad(2,PD.hfo{stimloc,kk},Meanprofile); 
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    ylabel('Channels'); xlabel('Power'); title('hfo power')
                    set(gca,'YDir','reverse');  
                    xlim([-1 1])

                end


            end
        end
    end
end

figure(jj)
col = [115/243 82/243 66/243; 56/243 61/243 150/243; 0/243 0/243 0/243];
nC = 1:1:23;

for ii = 1:length(reg)
    for kk = 1:2
        subplot(6,4,4*(rem(kk,3))+13)
        stdProfile = nanstd(PD.Theta{ii,kk},[],2)./sqrt(length(mice));
        Meanprofile = nanmean(PD.Theta{ii,kk},2);
        dev1 = Meanprofile-stdProfile;
        dev2 = Meanprofile+stdProfile;
        fill([dev1' flip(dev2')],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
        plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5);  
        
        subplot(6,4,4*(rem(kk,3))+14)
        stdProfile = nanstd(PD.sg{ii,kk},[],2)./sqrt(length(mice));
        Meanprofile = nanmean(PD.sg{ii,kk},2);
        dev1 = Meanprofile-stdProfile;
        dev2 = Meanprofile+stdProfile;
        fill([dev1' flip(dev2')],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
        plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5);  

        subplot(6,4,4*(rem(kk,3))+15)
        stdProfile = nanstd(PD.mg{ii,kk},[],2)./sqrt(length(mice));
        Meanprofile = nanmean(PD.mg{ii,kk},2);
        dev1 = Meanprofile-stdProfile;
        dev2 = Meanprofile+stdProfile;
        fill([dev1' flip(dev2')],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
        plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5); 

        subplot(6,4,4*(rem(kk,3))+16)
        stdProfile = nanstd(PD.hfo{ii,kk},[],2)./sqrt(length(mice));
        Meanprofile = nanmean(PD.hfo{ii,kk},2);
        dev1 = Meanprofile-stdProfile;
        dev2 = Meanprofile+stdProfile;
        fill([dev1' flip(dev2')],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
        plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5);  

    end
end


for jj = 1%:length(target)
    saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile',tag,'_',target{jj},'.png'),'png');
    saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile',tag,'_',target{jj},'.fig'),'fig');
    saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile',tag,'_',target{jj},'.eps'),'epsc');
end

%close all
end