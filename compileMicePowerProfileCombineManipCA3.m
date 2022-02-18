function PD = compileMicePowerProfileCombineManipCA3

mice = {'IZ27\Final'};%If in a 'Final' Folder, include the Final in the string
shankNum = [1];

parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

reg = {'mEC','Both'};
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

    channels = [1:sessionInfo.nChannels]-1;
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('PowerProfile.mat','file')
        load(strcat('PowerProfile.mat'));
    else 
        disp(['Power Profile not computed for mouse' mice{m}])
        continue;
    end
    
    for jj = 1%:length(target)
        for kk = 1:3
            baseTheta{kk} = [];
            basesg{kk} = [];
            basemg{kk} = [];        
            basehfo{kk} = [];
        end
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
                
                subplot(3,4,4*(rem(kk,3))+1)
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

                subplot(3,4,4*(rem(kk,3))+2)
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

                subplot(3,4,4*(rem(kk,3))+3)
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

                subplot(3,4,4*(rem(kk,3))+4)
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
            end
        end
        
        for kk = 1:3

            meanTheta = nanmean(baseTheta{kk}(:,:,2:end),3);
            meansg = nanmean(basesg{kk}(:,:,2:end),3);
            meanmg = nanmean(basemg{kk}(:,:,2:end),3);
            meanhfo = nanmean(basehfo{kk}(:,:,2:end),3);

            subplot(3,4,4*(rem(kk,3))+1)
            stdProfile = nanstd(baseTheta{kk},[],3)./sqrt(size(baseTheta{kk},3)-1);
            stdProfile = stdProfile(chanImp);
            Meanprofile = meanTheta(chanImp);
            dev1 = Meanprofile-stdProfile;
            dev2 = Meanprofile+stdProfile;
            fill([dev1' flip(dev2')],[nC flip(nC)],[85/243 85/243 85/243],'FaceAlpha',.2,'EdgeColor','none')
            plot(Meanprofile,nC,'color',[85/243 85/243 85/243],'LineWidth',1.5);  

            subplot(3,4,4*(rem(kk,3))+2)
            stdProfile = nanstd(basesg{kk},[],3)./sqrt(size(basesg{kk},3)-1);
            stdProfile = stdProfile(chanImp);
            Meanprofile = meansg(chanImp);
            dev1 = Meanprofile-stdProfile;
            dev2 = Meanprofile+stdProfile;
            fill([dev1' flip(dev2')],[nC flip(nC)],[85/243 85/243 85/243],'FaceAlpha',.2,'EdgeColor','none')
            plot(Meanprofile,nC,'color',[85/243 85/243 85/243],'LineWidth',1.5);  

            subplot(3,4,4*(rem(kk,3))+3)
            stdProfile = nanstd(basemg{kk},[],3)./sqrt(size(basemg{kk},3)-1);
            stdProfile = stdProfile(chanImp);    
            Meanprofile = meanmg(chanImp);
            dev1 = Meanprofile-stdProfile;
            dev2 = Meanprofile+stdProfile;
            fill([dev1' flip(dev2')],[nC flip(nC)],[85/243 85/243 85/243],'FaceAlpha',.2,'EdgeColor','none')
            plot(Meanprofile,nC,'color',[85/243 85/243 85/243],'LineWidth',1.5); 

            subplot(3,4,4*(rem(kk,3))+4)
            stdProfile = nanstd(basehfo{kk},[],3)./sqrt(size(basehfo{kk},3)-1);
            stdProfile = stdProfile(chanImp);      
            Meanprofile = meanhfo(chanImp);
            dev1 = Meanprofile-stdProfile;
            dev2 = Meanprofile+stdProfile;
            fill([dev1' flip(dev2')],[nC flip(nC)],[85/243 85/243 85/243],'FaceAlpha',.2,'EdgeColor','none')
            plot(Meanprofile,nC,'color',[85/243 85/243 85/243],'LineWidth',1.5);  

        end

    end
end


for jj = 1%:length(target)
    saveas(figure(jj),strcat(parentDir,mice{1},'\Summ\','compiledPowerProfile.png'),'png');
    saveas(figure(jj),strcat(parentDir,mice{1},'\Summ\','compiledPowerProfile.fig'),'fig');
    saveas(figure(jj),strcat(parentDir,mice{1},'\Summ\','compiledPowerProfile.eps'),'epsc');
end

%close all
end