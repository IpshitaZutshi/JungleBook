function compileMicePowerProfileCA3

mice = {'IZ28\Final'};%,'IZ28\Final','IZ29\Final'};%If in a 'Final' Folder, include the Final in the string
shankNum = [1,1,3];
parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';

reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

nf = 1;

for m = 1:length(mice)
    
    expPath = strcat(parentDir, mice{m});
    allpath = strsplit(genpath(expPath),';'); % all folders
    cd(allpath{1});
    allSess = dir('*_sess*');
    load(strcat(allSess(1).folder,'\',allSess(1).name,'\',allSess(1).name,'.sessionInfo.mat'));
    
    channels = [1:sessionInfo.nChannels]-1;
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('PowerProfileCSD.mat','file')
        load(strcat('PowerProfileCSD.mat'));
    else 
        disp(['Power Profile not computed for mouse' mice{m}])
        continue;
    end
    
    for jj = 1:length(target)
        for ii = 2:length(reg)
            figure(jj) % One plot for each manipulation
            set(gcf, 'Position', get(0, 'Screensize'));
            for kk = 1:length(zone)      

                meanTheta = nanmean(PowerProfile.theta{ii,jj}{kk}(:,:,2:end),3);
                meansg  = nanmean(PowerProfile.sg{ii,jj}{kk}(:,:,2:end),3);
                meanmg  = nanmean(PowerProfile.mg{ii,jj}{kk}(:,:,2:end),3);
                meanhfo  = nanmean(PowerProfile.hfo{ii,jj}{kk}(:,:,2:end),3);
                stdTheta = (nanstd(PowerProfile.theta{ii,jj}{kk}(:,:,2:end),[],3))/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                stdsg  = nanstd(PowerProfile.sg{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                stdmg  = nanstd(PowerProfile.mg{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);
                stdhfo  = nanstd(PowerProfile.hfo{ii,jj}{kk}(:,:,2:end),[],3)/sqrt(size(PowerProfile.theta{ii,jj}{kk},3)-1);

                if kk < 4 && ii == 2 % If baseline
                    col = [0 0 0];
                elseif kk< 4 && ii == 3 % If CA3 baseline
                    col = [56/243 61/243 150/243];                    
                elseif ii==2 && kk >=4 %If mEC 
                    col = [85/243 85/243 85/243];
                elseif ii==3 && kk >=4 %If mEC & CA3
                    col = [8/243 133/243 161/243];
                end

                    
                [Lia] = ismember(sessionInfo.AnatGrps(shankNum(m)).Channels, channels);
                nC = 1:length(sessionInfo.AnatGrps(shankNum(m)).Channels);
                nC = nC(Lia);
                chanImp = sessionInfo.AnatGrps(shankNum(m)).Channels(Lia)+1;
                
                if strcmp(mice{m},'IZ28\Final') == 1 %Remove noisy channels on top
                    chanImp = chanImp(6:end);
                elseif strcmp(mice{m},'IZ29\Final') == 1 %Remove noisy channels on top
                    chanImp = chanImp(3:end);
                end
                chanImp = chanImp(chanImp<=62);
                nC = 1:1:64;     
                x = 1:(64/length(chanImp)):64; %Determine current sampling                     
                
                subplot(3,4,4*(rem(kk,3))+1)
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meanTheta(chanImp);
                stdProfile = stdTheta(chanImp);           
                Meanprofile = interp1(x,Meanprofile,nC);
                stdProfile = interp1(x,stdProfile,nC);  
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1 flip(dev2)],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel('Channels'); xlabel('Power'); title('theta power')
                set(gca,'YDir','reverse');

                subplot(3,4,4*(rem(kk,3))+2)
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meansg(chanImp);
                stdProfile = stdsg(chanImp);       
                Meanprofile = interp1(x,Meanprofile,nC);
                stdProfile = interp1(x,stdProfile,nC);
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1 flip(dev2)],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel('Channels'); xlabel('Power'); title('sg power')
                set(gca,'YDir','reverse');

                subplot(3,4,4*(rem(kk,3))+3)
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meanmg(chanImp);
                stdProfile = stdmg(chanImp);           
                Meanprofile = interp1(x,Meanprofile,nC);
                stdProfile = interp1(x,stdProfile,nC);                
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1 flip(dev2)],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel('Channels'); xlabel('Power'); title('mg power')
                set(gca,'YDir','reverse');

                subplot(3,4,4*(rem(kk,3))+4)
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meanhfo(chanImp);
                stdProfile = stdhfo(chanImp);              
                Meanprofile = interp1(x,Meanprofile,nC);
                stdProfile = interp1(x,stdProfile,nC);                
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1 flip(dev2)],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel('Channels'); xlabel('Power'); title('hfo power')
                set(gca,'YDir','reverse');

            end
        end
    end
end


for jj = 1:length(target)
%     saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile_',target{jj},'.png'),'png');
%     saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile_',target{jj},'.fig'),'fig');
%     saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfile_',target{jj},'.eps'),'epsc');
end

close all
end