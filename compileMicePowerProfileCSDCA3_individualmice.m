function PD = compileMicePowerProfileCSDCA3_individualmice

shankNum = [1,1,3];
parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
tag = 'CA3Saline';

if strcmp(tag(4:end),'Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ29\Saline'};
else
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final'};
end
reg = {'CA3','mEC','Both'};
zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for ii = 1:3
    for kk = 1:3
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
    pyrCh = region.CA1sp;
    Chstart  = [];
    %Find index of pyramidal channel
    for ch = 1:size(sessionInfo.AnatGrps,2)
        if ismember(pyrCh, sessionInfo.AnatGrps(ch).Channels)
            Chstart = find(sessionInfo.AnatGrps(ch).Channels==pyrCh);         
        end
    end     
            
    channels = [1:sessionInfo.nChannels]-1;
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist('PowerProfileCSD.mat','file')
        load(strcat('PowerProfileCSD.mat'));
    else 
        disp(['Power Profile not computed for mouse' mice{m}])
        continue;
    end
    
    for jj = 1%:length(target)
        for ii = 2:length(reg)
            figure(jj) % One plot for each manipulation 
            set(gcf, 'Position', [100 100 1600 800]);
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

                if kk < 4 && ii == 2 % If baseline
                    col = [0 0 0];
                elseif kk< 4 && ii == 3 % If CA3 baseline
                    col = [224/243 163/243 46/243];                  
                elseif ii==2 && kk >=4 %If mEC 
                    col = [8/243 133/243 161/243];
                elseif ii==3 && kk >=4 %If mEC & CA3
                    col = [56/243 61/243 150/243];  
                end

                chanImp = 1:(length(sessionInfo.AnatGrps(shankNum(m)).Channels))-2;
                
                if strcmp(mice{m},'IZ29\Final') == 1 || strcmp(mice{m},'IZ29\Saline') == 1
                    chanImp = chanImp(3:(end-2));
                elseif strcmp(mice{m},'IZ28\Final') == 1 || strcmp(mice{m},'IZ28\Saline') == 1
                    chanImp = chanImp((Chstart):(end-8));
                else
                    chanImp = chanImp(Chstart:(end-6));
                end
   
                nC = 1:1:length(chanImp);                 
                
                subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+2*(m-1)+1))
                hold on
                Meanprofile = meanTheta(chanImp)';
                stdProfile = stdTheta(chanImp)';            
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1 flip(dev2)],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                ylabel(zone{kk}(1:end-1),'fontweight','bold'); 
                xlabel('Power'); 
                if m == 3                         
                    yticks([1:1:length(chanImp)])                     
                else
                    yticks([2:3:length(chanImp)])                    
                end
                yticklabels({})
                grid on
                if kk == 3
                    title('theta power')
                end
                set(gca,'YDir','reverse');

                subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+2*(m-1)+2))
                hold on
                % Interpolate so that each profile has 64 datapoints
                Meanprofile = meansg(chanImp)';
                stdProfile = stdsg(chanImp)';       
%                 Meanprofile = interp1(x,Meanprofile,nC);
%                 stdProfile = interp1(x,stdProfile,nC);
                dev1 = Meanprofile-stdProfile;
                dev2 = Meanprofile+stdProfile;
                fill([dev1 flip(dev2)],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                xlabel('Power'); 
                if m == 3                         
                    yticks([1:1:length(chanImp)])                     
                else
                    yticks([2:3:length(chanImp)])                    
                end
                yticklabels({})
                grid on                
                if kk == 3                
                    title('sg power')
                end
                set(gca,'YDir','reverse');

%                 subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+2*(m-1)+3))
%                 hold on
%                 % Interpolate so that each profile has 64 datapoints
%                 Meanprofile = meanmg(chanImp)';
%                 stdProfile = stdmg(chanImp)';           
% %                 Meanprofile = interp1(x,Meanprofile,nC);
% %                 stdProfile = interp1(x,stdProfile,nC);                
%                 dev1 = Meanprofile-stdProfile;
%                 dev2 = Meanprofile+stdProfile;
%                 fill([dev1 flip(dev2)],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
%                 plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
%                 xlabel('Power'); 
%                 if m == 3                         
%                     yticks([1:1:length(chanImp)])                     
%                 else
%                     yticks([2:3:length(chanImp)])                    
%                 end
%                 yticklabels({})
%                 grid on                
%                 if kk == 3
%                     title('mg power')
%                 end
%                 set(gca,'YDir','reverse');
% 
%                 subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+2*(m-1)+4))
%                 hold on
%                 % Interpolate so that each profile has 64 datapoints
%                 Meanprofile = meanhfo(chanImp)';
%                 stdProfile = stdhfo(chanImp)';              
% %                 Meanprofile = interp1(x,Meanprofile,nC);
% %                 stdProfile = interp1(x,stdProfile,nC);                
%                 dev1 = Meanprofile-stdProfile;
%                 dev2 = Meanprofile+stdProfile;
%                 fill([dev1 flip(dev2)],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
%                 plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
%                 xlabel('Power'); 
%                 if m == 3                         
%                     yticks([1:1:length(chanImp)])                     
%                 else
%                     yticks([2:3:length(chanImp)])                    
%                 end
%                 yticklabels({})
%                 grid on                
%                 if kk == 3                
%                     title('hfo power')
%                 end
%                 set(gca,'YDir','reverse');
                
                % Now plot the difference
                if ii == 2 && kk<4
                    
                    diffTheta = -(PowerProfile.theta{ii,jj}{kk} - PowerProfile.theta{ii,jj}{kk+3});%./PowerProfile.theta{2,jj}{kk};
                    meanTheta = nanmean(diffTheta(:,:,2:end),3);
                    stdTheta = (nanstd(diffTheta(:,:,2:end),[],3))/sqrt(size(diffTheta,3)-1);

                    diffsg = -(PowerProfile.sg{ii,jj}{kk} - PowerProfile.sg{ii,jj}{kk+3});%./PowerProfile.sg{2,jj}{kk};
                    meansg= nanmean(diffsg(:,:,2:end),3);
                    stdsg = (nanstd(diffsg(:,:,2:end),[],3))/sqrt(size(diffsg,3)-1);

                    diffmg = -(PowerProfile.mg{ii,jj}{kk} - PowerProfile.mg{ii,jj}{kk+3});%./PowerProfile.mg{2,jj}{kk};
                    meanmg = nanmean(diffmg(:,:,2:end),3);
                    stdmg = (nanstd(diffmg(:,:,2:end),[],3))/sqrt(size(diffmg,3)-1);                        

                    diffhfo = -(PowerProfile.hfo{ii,jj}{kk} - PowerProfile.hfo{ii,jj}{kk+3});%./PowerProfile.hfo{2,jj}{kk};
                    meanhfo = nanmean(diffhfo(:,:,2:end),3);
                    stdhfo = (nanstd(diffhfo(:,:,2:end),[],3))/sqrt(size(diffhfo,3)-1); 


                    col = [8/243 133/243 161/243];

                    subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+1))
                    hold on
                    Meanprofile = meanTheta(chanImp)';
                    stdProfile = stdTheta(chanImp)';                      
                    Meanprofile = smooth(Meanprofile);
                    stdProfile = smooth(stdProfile);
                    PD.Theta{1,kk} = catpad(2,PD.Theta{1,kk},Meanprofile);
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    ylabel(zone{kk}(1:end-1),'fontweight','bold'); 
                    xlabel('Power');  
                    if m == 3                         
                        yticks([1:1:length(chanImp)])                     
                    else
                        yticks([2:3:length(chanImp)])                    
                    end
                    yticklabels({})
                    grid on                    
                    if kk == 3
                        title('Diff theta')
                    end
                    set(gca,'YDir','reverse');
                    
                    subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+2))
                    hold on
                    Meanprofile = meansg(chanImp)';
                    stdProfile = stdsg(chanImp)';                     
                    Meanprofile = smooth(Meanprofile); 
                    stdProfile = smooth(stdProfile);
                    PD.sg{1,kk} = catpad(2,PD.sg{1,kk},Meanprofile);
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    xlabel('Power'); 
                    if m == 3                         
                        yticks([1:1:length(chanImp)])                     
                    else
                        yticks([2:3:length(chanImp)])                    
                    end
                    yticklabels({})
                    grid on                                        
                    if kk == 3
                        title('Diff sg')
                    end                    
                    set(gca,'YDir','reverse');                     

%                     subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+3))
%                     hold on
%                     Meanprofile = meanmg(chanImp)';
%                     stdProfile = stdmg(chanImp)';                                      
%                     Meanprofile = smooth(Meanprofile); 
%                     stdProfile = smooth(stdProfile);
%                     PD.mg{1,kk} = catpad(2,PD.mg{1,kk},Meanprofile);
%                     dev1 = Meanprofile-stdProfile;
%                     dev2 = Meanprofile+stdProfile;
%                     fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
%                     plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
%                     xlabel('Power'); 
%                     if m == 3                         
%                         yticks([1:1:length(chanImp)])                     
%                     else
%                         yticks([2:3:length(chanImp)])                    
%                     end
%                     yticklabels({})
%                     grid on                                        
%                     if kk == 3
%                         title('Diff mg')
%                     end                      
%                     set(gca,'YDir','reverse');  
% 
%                     subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+4))
%                     hold on
%                     Meanprofile = meanhfo(chanImp)';
%                     stdProfile = stdhfo(chanImp)';                      
%                     Meanprofile = smooth(Meanprofile);                    
%                     stdProfile = smooth(stdProfile);
%                     PD.hfo{1,kk} = catpad(2,PD.hfo{1,kk},Meanprofile);
%                     dev1 = Meanprofile-stdProfile;
%                     dev2 = Meanprofile+stdProfile;
%                     fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
%                     plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
%                     xlabel('Power'); 
%                     if m == 3                         
%                         yticks([1:1:length(chanImp)])                     
%                     else
%                         yticks([2:3:length(chanImp)])                    
%                     end
%                     yticklabels({})
%                     grid on                      
%                     if kk == 3
%                         title('Diff hfo')
%                     end                      
%                     set(gca,'YDir','reverse');  
    
                    diffTheta = -(PowerProfile.theta{2,jj}{kk} - PowerProfile.theta{3,jj}{kk});%./PowerProfile.theta{2,jj}{kk};
                    meanTheta = nanmean(diffTheta(:,:,2:end),3);
                    stdTheta = (nanstd(diffTheta(:,:,2:end),[],3))/sqrt(size(diffTheta,3)-1);

                    diffsg = -(PowerProfile.sg{2,jj}{kk} - PowerProfile.sg{3,jj}{kk});%./PowerProfile.sg{2,jj}{kk};
                    meansg= nanmean(diffsg(:,:,2:end),3);
                    stdsg = (nanstd(diffsg(:,:,2:end),[],3))/sqrt(size(diffsg,3)-1);

                    diffmg = -(PowerProfile.mg{2,jj}{kk} - PowerProfile.mg{3,jj}{kk});%./PowerProfile.mg{2,jj}{kk};
                    meanmg = nanmean(diffmg(:,:,2:end),3);
                    stdmg = (nanstd(diffmg(:,:,2:end),[],3))/sqrt(size(diffmg,3)-1);                        

                    diffhfo = -(PowerProfile.hfo{2,jj}{kk} - PowerProfile.hfo{3,jj}{kk});%./PowerProfile.hfo{2,jj}{kk};
                    meanhfo = nanmean(diffhfo(:,:,2:end),3);
                    stdhfo = (nanstd(diffhfo(:,:,2:end),[],3))/sqrt(size(diffhfo,3)-2); 
                    
                    col = [224/243 163/243 46/243];

                    subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+1))
                    hold on
                    Meanprofile = meanTheta(chanImp)';
                    stdProfile = stdTheta(chanImp)';                     
                    Meanprofile = smooth(Meanprofile);                    
                    stdProfile = smooth(stdProfile);
                    PD.Theta{2,kk} = catpad(2,PD.Theta{2,kk},Meanprofile);
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;                    
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    xlabel('Power');
                    if m == 3                         
                        yticks([1:1:length(chanImp)])                     
                    else
                        yticks([2:3:length(chanImp)])                    
                    end
                    yticklabels({})
                    grid on                      
                    set(gca,'YDir','reverse');
                    
                    subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+2))
                    hold on
                    Meanprofile = meansg(chanImp)';
                    stdProfile = stdsg(chanImp)';                     
                    Meanprofile = smooth(Meanprofile);                    
                    stdProfile = smooth(stdProfile);
                    PD.sg{2,kk} = catpad(2,PD.sg{2,kk},Meanprofile);
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;                      
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    xlabel('Power');
                    if m == 3                         
                        yticks([1:1:length(chanImp)])                     
                    else
                        yticks([2:3:length(chanImp)])                    
                    end
                    yticklabels({})
                    grid on                      
                    set(gca,'YDir','reverse');                     

%                     subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+3))
%                     hold on
%                     Meanprofile = meanmg(chanImp)';
%                     stdProfile = stdmg(chanImp)';                      
%                     Meanprofile = smooth(Meanprofile);                    
%                     stdProfile = smooth(stdProfile);
%                     PD.mg{2,kk} = catpad(2,PD.mg{2,kk},Meanprofile);
%                     dev1 = Meanprofile-stdProfile;
%                     dev2 = Meanprofile+stdProfile;
%                     fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
%                     plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
%                     xlabel('Power'); 
%                     if m == 3                         
%                         yticks([1:1:length(chanImp)])                     
%                     else
%                         yticks([2:3:length(chanImp)])                    
%                     end
%                     yticklabels({})
%                     grid on                      
%                     set(gca,'YDir','reverse');  
% 
%                     subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+4))
%                     hold on
%                     Meanprofile = meanhfo(chanImp)';
%                     stdProfile = stdhfo(chanImp)';                     
%                     Meanprofile = smooth(Meanprofile);                    
%                     stdProfile = smooth(stdProfile);
%                     PD.hfo{2,kk} = catpad(2,PD.hfo{2,kk},Meanprofile);
%                     dev1 = Meanprofile-stdProfile;
%                     dev2 = Meanprofile+stdProfile;
%                     fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
%                     plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
%                     xlabel('Power'); 
%                     if m == 3                         
%                         yticks([1:1:length(chanImp)])                     
%                     else
%                         yticks([2:3:length(chanImp)])                    
%                     end
%                     yticklabels({})
%                     grid on                      
%                     set(gca,'YDir','reverse');  
%                     

                    diffTheta = -(PowerProfile.theta{2,jj}{kk} - PowerProfile.theta{3,jj}{kk+3});%./PowerProfile.theta{2,jj}{kk};
                    meanTheta = nanmean(diffTheta(:,:,2:end),3);
                    stdTheta = (nanstd(diffTheta(:,:,2:end),[],3))/sqrt(size(diffTheta,3)-1);

                    diffsg = -(PowerProfile.sg{2,jj}{kk} - PowerProfile.sg{3,jj}{kk+3});%./PowerProfile.sg{2,jj}{kk};
                    meansg= nanmean(diffsg(:,:,2:end),3);
                    stdsg = (nanstd(diffsg(:,:,2:end),[],3))/sqrt(size(diffsg,3)-1);

                    diffmg = -(PowerProfile.mg{2,jj}{kk} - PowerProfile.mg{3,jj}{kk+3});%./PowerProfile.mg{2,jj}{kk};
                    meanmg = nanmean(diffmg(:,:,2:end),3);
                    stdmg = (nanstd(diffmg(:,:,2:end),[],3))/sqrt(size(diffmg,3)-1);                        

                    diffhfo = -(PowerProfile.hfo{2,jj}{kk} - PowerProfile.hfo{3,jj}{kk+3});%./PowerProfile.hfo{2,jj}{kk};
                    meanhfo = nanmean(diffhfo(:,:,2:end),3);
                    stdhfo = (nanstd(diffhfo(:,:,2:end),[],3))/sqrt(size(diffhfo,3)-1); 

                    col = [56/243 61/243 150/243];  
 
                    subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+1))
                    hold on
                    Meanprofile = meanTheta(chanImp)';
                    stdProfile = stdTheta(chanImp)';                     
                    Meanprofile = smooth(Meanprofile);                    
                    stdProfile = smooth(stdProfile);
                    PD.Theta{3,kk} = catpad(2,PD.Theta{3,kk},Meanprofile);
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;                    
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    xlabel('Power');
                    if m == 3                         
                        yticks([1:1:length(chanImp)])                     
                    else
                        yticks([2:3:length(chanImp)])                    
                    end
                    yticklabels({})
%                     xlim([-2 2])    
%                     xticks([-1 -0.5 0 0.5 1])
                    grid on                      
                    set(gca,'YDir','reverse');
                    
                    subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+2))
                    hold on
                    Meanprofile = meansg(chanImp)';
                    stdProfile = stdsg(chanImp)';                     
                    Meanprofile = smooth(Meanprofile);                    
                    stdProfile = smooth(stdProfile);
                    PD.sg{3,kk} = catpad(2,PD.sg{3,kk},Meanprofile);
                    dev1 = Meanprofile-stdProfile;
                    dev2 = Meanprofile+stdProfile;                      
                    fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
                    plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
                    xlabel('Power');
                    if m == 3                         
                        yticks([1:1:length(chanImp)])                     
                    else
                        yticks([2:3:length(chanImp)])                    
                    end
                    yticklabels({})
%                     xlim([-2 2])    
%                     xticks([-1 -0.5 0 0.5 2])                    
                    grid on                      
                    set(gca,'YDir','reverse');                     


%                     subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+3))
%                     hold on
%                     Meanprofile = meanmg(chanImp)';
%                     stdProfile = stdmg(chanImp)';                        
%                     Meanprofile = smooth(Meanprofile);                    
%                     stdProfile = smooth(stdProfile);
%                     PD.mg{3,kk} = catpad(2,PD.mg{3,kk},Meanprofile);
%                     dev1 = Meanprofile-stdProfile;
%                     dev2 = Meanprofile+stdProfile;
%                     fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
%                     plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
%                     xlabel('Power'); 
%                     if m == 3                         
%                         yticks([1:1:length(chanImp)])                     
%                     else
%                         yticks([2:3:length(chanImp)])                    
%                     end
%                     yticklabels({})
% %                     xlim([-1 1])    
% %                     xticks([-1 -0.5 0 0.5 1])
%                     grid on                 
%                     set(gca,'YDir','reverse');  
% 
%                     subplot(6,2*length(mice),(2*length(mice)*(rem(kk,3))+6*length(mice)+2*(m-1)+4))
%                     hold on
%                     Meanprofile = meanhfo(chanImp)';
%                     stdProfile = stdhfo(chanImp)';                     
%                     Meanprofile = smooth(Meanprofile);                    
%                     stdProfile = smooth(stdProfile);
%                     PD.hfo{3,kk} = catpad(2,PD.hfo{3,kk},Meanprofile);
%                     dev1 = Meanprofile-stdProfile;
%                     dev2 = Meanprofile+stdProfile;
%                     fill([dev1' flip(dev2')],[nC flip(nC)],col,'FaceAlpha',.2,'EdgeColor','none')
%                     plot(Meanprofile,nC,'color',col,'LineWidth',1.5); 
%                     xlabel('Power');
%                     if m == 3
%                         yticks([1:1:length(chanImp)])
%                     else
%                         yticks([2:3:length(chanImp)])
%                     end
%                     yticklabels({})
% %                     xlim([-1 1])    
% %                     xticks([-1 -0.5 0 0.5 1])
%                     grid on                      
%                     set(gca,'YDir','reverse');  
                         
                end             
            end             
        end
    end
end

% figure(jj)
% col = [187/243 86/243 149/243; 115/243 82/243 66/243; 0/243 0/243 0/243];
% nC = 1:1:64;
% 
% for ii = 1:3
%     for kk = 1:2
%         subplot(6,4,4*(rem(kk,3))+13)
%         stdProfile = nanstd(PD.Theta{ii,kk},[],2)./sqrt(length(mice));
%         Meanprofile = nanmean(PD.Theta{ii,kk},2);
%         dev1 = Meanprofile-stdProfile;
%         dev2 = Meanprofile+stdProfile;
%         fill([dev1' flip(dev2')],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
%         plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5);  
%         
%         subplot(6,4,4*(rem(kk,3))+14)
%         stdProfile = nanstd(PD.sg{ii,kk},[],2)./sqrt(length(mice));
%         Meanprofile = nanmean(PD.sg{ii,kk},2);
%         dev1 = Meanprofile-stdProfile;
%         dev2 = Meanprofile+stdProfile;
%         fill([dev1' flip(dev2')],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
%         plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5);  
% 
%         subplot(6,4,4*(rem(kk,3))+15)
%         stdProfile = nanstd(PD.mg{ii,kk},[],2)./sqrt(length(mice));
%         Meanprofile = nanmean(PD.mg{ii,kk},2);
%         dev1 = Meanprofile-stdProfile;
%         dev2 = Meanprofile+stdProfile;
%         fill([dev1' flip(dev2')],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
%         plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5); 
% 
%         subplot(6,4,4*(rem(kk,3))+16)
%         stdProfile = nanstd(PD.hfo{ii,kk},[],2)./sqrt(length(mice));
%         Meanprofile = nanmean(PD.hfo{ii,kk},2);
%         dev1 = Meanprofile-stdProfile;
%         dev2 = Meanprofile+stdProfile;
%         fill([dev1' flip(dev2')],[nC flip(nC)],col(ii,:),'FaceAlpha',.2,'EdgeColor','none')
%         plot(Meanprofile,nC,'color',col(ii,:),'LineWidth',1.5);  
% 
%     end
% end
% 
% 
% for jj = 1:length(target)
%     saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfileIndividualMice',tag,'_',target{jj},'.png'),'png');
%     saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfileIndividualMice',tag,'_',target{jj},'.fig'),'fig');
%     saveas(figure(jj),strcat(parentDir,'\Compiled\compiledPowerProfileIndividualMice',tag,'_',target{jj},'.eps'),'epsc');
% end

%close all
end