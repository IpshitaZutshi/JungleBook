function compileMicePhaseAmpCSD(varargin)

p = inputParser;
addParameter(p,'parentDir','Z:\Homes\zutshi01\Recordings\CA1_silencing\',@isfolder);
parse(p,varargin{:});

parentDir = p.Results.parentDir;
tag = 'CA1';
chName =  {'pyramidal','radiatum','slm','ml'};
   
if strcmp(tag,'CA1') == 1
    mice = {'IZ18\Final','IZ31\Final','IZ20\Final','IZ21\Final'}; 
    %mice = {'IZ31\Final'};%,'IZ20\Final'};%,'IZ30\Final','IZ15\Final'};%,'IZ30\Final'};%
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ18\Final','IZ20\Final','IZ21\Final','IZ31\Final','IZ27\Final','IZ28\Final','IZ33\Final'};  
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ33\Final'};%,'IZ32\Final','IZ29\Final'};
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ27\Saline','IZ28\Saline','IZ33\Saline'};%%'IZ29\Saline','IZ32\Saline',
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Both'};
end

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for ch = 1:4 
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
    if exist(strcat('Summ\PhaseAmpCSD.mat'),'file')
        load(strcat('Summ\PhaseAmpCSD.mat'));
    else 
        disp(['Phase Amp not computed for mouse' mice{m}])
        continue;
    end
    
    for ch = 1:length(chName)
        for ii = 1:length(reg)
            for jj = 1%:length(target)
                for kk = [1 2 4 5]      
%                     meancomod = nanmean(ModData.(chName{ch}){ii,jj}{kk},3);
%                     compiledModData.(chName{ch}){ii,jj}{kk}(:,:,m) = meancomod;
                   if isfield(ModData,chName{ch})
                       if (strcmp(mice{m},'IZ21\Final')==1 || strcmp(mice{m},'IZ18\Final')==1) && (ii == 1 || ii==3) && ch == 1
                           continue
                       elseif ((strcmp(mice{m},'IZ15\Final')==1 )|| (strcmp(mice{m},'IZ30\Final')==1 )) && (ii ~= 1 || ii~=3 || ch ~= 1)
                           continue                           
                       else
                           compiledModData.(chName{ch}){ii,jj}{kk} = catpad(3,compiledModData.(chName{ch}){ii,jj}{kk},ModData.(chName{ch}){ii,jj}{kk});
                      end
                   end
                end
            end
        end
    end
end

phaserange = phasecenters;
amprange = amp_freq;

figure
set(gcf,'Position',[100 100 900 600])
colMat = [85/243 85/243 85/243;...
    224/243 163/243 46/243;... 
    8/243 133/243 161/243;...
    56/243 61/243 150/243];  
% 
% freqRange = {[150 250],[25 50],[50 120],[100 220]};
% thetafreqRange = [6 11];

if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
   for ch = 1:4
        data = [];
        manip = [];
        for jj = 1%:length(target)
             subplot(4,5,5*(ch-1)+1)  
             meancomod = nanmean(compiledModData.(chName{ch}){2,jj}{2},3);
             imagesc(phasecenters,log2(amp_freq),meancomod')
             hold on
             imagesc(phasecenters+2*pi,log2(amp_freq),meancomod')
             plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(amp_freq(end/2)),'k')
             ylabel('Frequency (Hz)')
             LogScale('y',2)       
             xlim([phasecenters(1) phasecenters(end)+2*pi]);
             ylim([log2(25) log2(150)])
             caxis([0.99 1.15])                          
             colorbar
             axis xy             
             colormap jet
             title('Baseline');
             
             subplot(4,5,5*(ch-1)+2)  
             meancomod = nanmean(compiledModData.(chName{ch}){2,jj}{5},3);
             imagesc(phasecenters,log2(amp_freq),meancomod')
             hold on
             imagesc(phasecenters+2*pi,log2(amp_freq),meancomod')
             plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(amp_freq(end/2)),'k')
             ylabel('Frequency (Hz)')
             LogScale('y',2)       
             xlim([phasecenters(1) phasecenters(end)+2*pi]);
             ylim([log2(25) log2(150)])
             caxis([0.99 1.15])                          
             colorbar
             axis xy             
             colormap jet
             title('mEC');
             
             subplot(4,5,5*(ch-1)+3)  
             meancomod = nanmean(compiledModData.(chName{ch}){3,jj}{2},3);
             imagesc(phasecenters,log2(amp_freq),meancomod')
             hold on
             imagesc(phasecenters+2*pi,log2(amp_freq),meancomod')
             plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(amp_freq(end/2)),'k')
             ylabel('Frequency (Hz)')
             LogScale('y',2)       
             xlim([phasecenters(1) phasecenters(end)+2*pi]);
             ylim([log2(25) log2(150)])
             caxis([0.99 1.15])                          
             colorbar
             axis xy             
             colormap jet
             title('CA3');
             
             subplot(4,5,5*(ch-1)+4)  
             meancomod = nanmean(compiledModData.(chName{ch}){3,jj}{5},3);
             imagesc(phasecenters,log2(amp_freq),meancomod')
             hold on
             imagesc(phasecenters+2*pi,log2(amp_freq),meancomod')
             plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(amp_freq(end/2)),'k')
             ylabel('Frequency (Hz)')
             LogScale('y',2)       
             xlim([phasecenters(1) phasecenters(end)+2*pi]);
             ylim([log2(25) log2(150)])
             caxis([0.99 1.1])             
             colorbar
             axis xy             
             colormap jet
             title('mEC & CA3');
             
%                                      
%              selidxY = phaserange(2:end)>=thetafreqRange(1) & phaserange(2:end)<=thetafreqRange(2);
%              selidxX = amprange(2:end) >= freqRange{ch}(1) & amprange(2:end) <= freqRange{ch}(2);
%              data1 = compiledModData.(chName{ch}){2,jj}{2}(selidxX,selidxY,:);
%              data2 = compiledModData.(chName{ch}){2,jj}{5}(selidxX,selidxY,:);                      
%              data3 = compiledModData.(chName{ch}){3,jj}{2}(selidxX,selidxY,:); 
%              data4 = compiledModData.(chName{ch}){3,jj}{5}(selidxX,selidxY,:); 
%              d1 = reshape(nanmean(data1,[1 2]),[1,size(data1,3)]);
%              d2 = reshape(nanmean(data2,[1 2]),[1,size(data2,3)]);   
%              d3 = reshape(nanmean(data3,[1 2]),[1,size(data3,3)]);
%              d4 = reshape(nanmean(data4,[1 2]),[1,size(data4,3)]);            
%              
%              data = [(d3-d1)./d1 (d2-d1)./d1 (d4-d1)./d1];
%              manip = [ones(1,length(d3))*1 ones(1,length(d2))*2 ones(1,length(d4))*3];  
%              
%              subplot(4,5,5*(ch-1)+5)  
%              stats{ch} = groupStats(data,manip,'inAxis',true,'plotData',false,'labelSummary',false,'color',colMat(2:end,:));
%              %ylim([-1 1])
%              for ss = 1:3
%                  [stats{ch}.signrank.p(ss),~,stats{ch}.signrank.stats{ss}] = signrank(data(manip==ss));
%              end
        end
   end
%    saveas(gcf,strcat(parentDir,'\Compiled\',tag,'CSDThetaGammaCoupling.png'),'png');
%    saveas(gcf,strcat(parentDir,'\Compiled\',tag,'CSDThetaGammaCoupling.fig'),'fig');
%    saveas(gcf,strcat(parentDir,'\Compiled\',tag,'CSDThetaGammaCoupling.eps'),'epsc');
%    save(strcat(parentDir,'\Compiled\',tag,'CSDThetaGammaCoupling.mat'),'stats');
else
    for ch = 1:4
        data = [];
        manip = [];
        modIdxBase = [];
        for ii = 1:length(reg)        
            for jj = 1%:length(target)
                 for kk = [2 5]
                     subplot(4,5,5*(ch-1)+ii+1)  
                     meancomod = nanmean(compiledModData.(chName{ch}){ii,jj}{kk},3);                
                     if kk == 2
                         modIdxBase = catpad(3,modIdxBase,compiledModData.(chName{ch}){ii,jj}{kk});
%                          selidxY = phaserange(2:end)>=thetafreqRange(1) & phaserange(2:end)<=thetafreqRange(2);
%                          selidxX = amprange(2:end) >= freqRange{ch}(1) & amprange(2:end) <= freqRange{ch}(2);
%                          data1 = compiledModData.(chName{ch}){ii,jj}{2}(selidxX,selidxY,:);
%                          data2 = compiledModData.(chName{ch}){ii,jj}{5}(selidxX,selidxY,:);
% %                          d1 = reshape(max(data1,[],[1 2]),[1,size(data1,3)]);
% %                          d2 = reshape(max(data2,[],[1 2]),[1,size(data2,3)]);                         
%                          d1 = reshape(nanmean(data1,[1 2]),[1,size(data1,3)]);
%                          d2 = reshape(nanmean(data2,[1 2]),[1,size(data2,3)]);    
%                          data = [data (d2-d1)./d1];
%                          manip = [manip ones(1,length(d1))*ii];                              
                     else
                         imagesc(phasecenters,log2(amp_freq),meancomod')
                         hold on
                         imagesc(phasecenters+2*pi,log2(amp_freq),meancomod')
                         plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(amp_freq(end/2)),'k')
                         ylabel('Amplitude Frequency (Hz)')
                         LogScale('y',2)       
                         xlim([phasecenters(1) phasecenters(end)+2*pi]);
                         ylim([log2(25) log2(150)])
                         caxis([0.99 1.1])
                         colorbar
                         axis xy             
                         colormap jet                         
                         title(reg{ii});
                     end
                 end
            end
        end
        subplot(4,5,5*(ch-1)+1)  
        imagesc(phasecenters,log2(amp_freq),nanmean(modIdxBase,3)')
        hold on
        imagesc(phasecenters+2*pi,log2(amp_freq),nanmean(modIdxBase,3)')
        plot(linspace(-pi,3*pi),cos(linspace(-pi,3*pi))+log2(amp_freq(end/2)),'k')
        ylabel('Amplitude Frequency (Hz)')
        LogScale('y',2)       
        xlim([phasecenters(1) phasecenters(end)+2*pi]);
        colorbar
        ylim([log2(25) log2(150)])
        caxis([0.99 1.1])        
        axis xy             
        colormap jet        
        title('Baseline');
%         
%         subplot(4,5,5*(ch-1)+5)  
%         stats{ch} = groupStats(data,manip,'inAxis',true,'plotData',false,'labelSummary',false,'color',colMat(2:end,:));
%         %ylim([-1 1])
%         for ss = 1:3
%             [stats{ch}.signrank.p(ss),~,stats{ch}.signrank.stats{ss}] = signrank(data(manip==ss));
%         end
    end
%     saveas(gcf,strcat(parentDir,'\Compiled\',tag,'CSDThetaGammaCoupling.png'),'png');
%     saveas(gcf,strcat(parentDir,'\Compiled\',tag,'CSDThetaGammaCoupling.fig'),'fig');
%     saveas(gcf,strcat(parentDir,'\Compiled\',tag,'CSDThetaGammaCoupling.eps'),'epsc');
%     save(strcat(parentDir,'\Compiled\',tag,'CSDThetaGammaCoupling.mat'),'stats');
end
end