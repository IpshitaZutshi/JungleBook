function compileMiceModIdxICA(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
tag = 'CA1';
chName =  {'pyramidal','radiatum','slm'};
   
if strcmp(tag,'CA1') == 1
    mice = {'IZ18\Final','IZ31\Final','IZ21\Final',};%,'IZ20\Final','IZ21\Final','IZ31\Final'};%,'IZ15\Final',,'IZ30\Final'
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final','IZ18\Final','IZ20\Final'...
        'IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline',...
        'IZ29\Saline','IZ30\Final','IZ31\Final','IZ32\Saline','IZ33\Saline','IZ34\Saline'};  
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ33\Final','IZ34\Final'};%'IZ32\Final','IZ29\Final',
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ33\Saline','IZ34\Saline'};%%'IZ29\Saline','IZ32\Saline',
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Both'};
end

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for ch = 1:3 
    for rr = 1:length(reg)
        for cc = 1:length(target)
            for zz = 1:length(zone)
               compiledModData.(chName{ch}){rr,cc}{zz} = [];
            end
        end
    end
end

for m = 1:length(mice)
    
    cd(strcat(parentDir, mice{m},'\Summ'));
    if exist(strcat('Summ\ModIdxICA.mat'),'file')
        load(strcat('Summ\ModIdxICA.mat'));
    else 
        disp(['Mod Idx not computed for mouse' mice{m}])
        continue;
    end
    
    for ch = 1:length(chName)
        for ii = 1:length(reg)
            for jj = 1%:length(target)
                for kk = 1:length(zone)      
%                     meancomod = nanmean(ModData.(chName{ch}){ii,jj}{kk},3);
%                     compiledModData.(chName{ch}){ii,jj}{kk}(:,:,m) = meancomod;
                   compiledModData.(chName{ch}){ii,jj}{kk} = catpad(3,compiledModData.(chName{ch}){ii,jj}{kk},ModData.(chName{ch}){ii,jj}{kk});
                end
            end
        end
    end
end

phaserange = 2:0.2:20;
amprange = 25:1:300;

figure(1)
set(gcf,'Position',[100 100 900 600])

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
   for ch = 1:3
        for jj = 1%:length(target)
             subplot(3,4,4*(ch-1)+1)  
             meancomod = nanmean(compiledModData.(chName{ch}){2,jj}{2},3);
             imagesc(phaserange(2:end),amprange(2:end),meancomod)  
             set(gca,'YDir','normal')
             colormap jet
             if ch == 1
                caxis([0 0.003])
             else
                caxis([0 0.004])
             end
             colorbar
             xlabel('Frequency phase');
             ylabel('Frequency amplitude');
             ylim([25 250])
             title('Baseline');
             
             subplot(3,4,4*(ch-1)+2)  
             meancomod = nanmean(compiledModData.(chName{ch}){2,jj}{5},3);
             imagesc(phaserange(2:end),amprange(2:end),meancomod)  
             set(gca,'YDir','normal')
             colormap jet
             if ch == 1
                caxis([0 0.003])
             else
                caxis([0 0.004])
             end
             colorbar
             xlabel('Frequency phase');
             ylabel('Frequency amplitude');
             ylim([25 250])
             title('mEC');
             
             subplot(3,4,4*(ch-1)+3)  
             meancomod = nanmean(compiledModData.(chName{ch}){3,jj}{2},3);
             imagesc(phaserange(2:end),amprange(2:end),meancomod)  
             set(gca,'YDir','normal')
             colormap jet
             if ch == 1
                caxis([0 0.003])
             else
                caxis([0 0.004])
             end
             colorbar
             xlabel('Frequency phase');
             ylabel('Frequency amplitude');
             ylim([25 250])
             title('CA3');
             
             subplot(3,4,4*(ch-1)+4)  
             meancomod = nanmean(compiledModData.(chName{ch}){3,jj}{5},3);
             imagesc(phaserange(2:end),amprange(2:end),meancomod)  
             set(gca,'YDir','normal')
             colormap jet
             if ch == 1
                caxis([0 0.003])
             else
                caxis([0 0.004])
             end
             colorbar
             xlabel('Frequency phase');
             ylabel('Frequency amplitude');
             ylim([25 250])
             title('mEC & CA3');
        end
   end
   saveas(gcf,strcat(parentDir,'\Compiled\',tag,'ICAThetaGammaCoupling.png'),'png');
   saveas(gcf,strcat(parentDir,'\Compiled\',tag,'ICAThetaGammaCoupling.fig'),'fig');
   saveas(gcf,strcat(parentDir,'\Compiled\',tag,'ICAThetaGammaCoupling.eps'),'epsc');
else
    for ch = 1:3
        modIdxBase = [];
        for ii = 1:length(reg)        
            for jj = 1%:length(target)
                 for kk = [2 5]
                     subplot(3,4,4*(ch-1)+ii+1)  
                     meancomod = nanmean(compiledModData.(chName{ch}){ii,jj}{kk},3);
                     if kk == 2
                         modIdxBase = catpad(3,modIdxBase,meancomod);
                     else
                         imagesc(phaserange(2:end),amprange(2:end),meancomod)  
                         %imagesc(phaserange(2:end),log2(amprange(2:end)),meancomod)
                         hold on
                         set(gca,'YDir','normal')
                         colormap jet
                         if ch == 1
                            caxis([0 0.004])
                         elseif ch==2
                            caxis([0 0.0015])
                         else
                            caxis([0 0.010])
                         end
                         colorbar
                         xlabel('Frequency phase');
                         ylabel('Frequency amplitude');
                         %ylim([25 200])
                         title(reg{ii});
                     end
                 end
            end
        end
        subplot(3,4,4*(ch-1)+1)  
        %imagesc(phaserange(2:end),log2(amprange(2:end)),nanmean(modIdxBase,3))
        imagesc(phaserange(2:end),amprange(2:end),nanmean(modIdxBase,3))
        set(gca,'YDir','normal')
        colormap jet
        colorbar
         if ch == 1
            caxis([0 0.004])
         elseif ch==2
            caxis([0 0.0015])
         else
            caxis([0 0.010])
         end
        %ylim([25 200])
        xlabel('Frequency phase');
        ylabel('Frequency amplitude');
        title('Baseline');
    end
    saveas(gcf,strcat(parentDir,'\Compiled\',tag,'ICAThetaGammaCoupling.png'),'png');
    saveas(gcf,strcat(parentDir,'\Compiled\',tag,'ICAThetaGammaCoupling.fig'),'fig');
    saveas(gcf,strcat(parentDir,'\Compiled\',tag,'ICAThetaGammaCoupling.eps'),'epsc');
end
end