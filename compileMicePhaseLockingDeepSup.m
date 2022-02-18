function compileMicePhaseLockingDeepSup

parentDir = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
tag = 'CA3';
taglfpName = {'pyrCh'};

if strcmp(tag,'CA1') == 1
    mice = {'IZ15\Final','IZ18\Final','IZ20\Final','IZ30\Final','IZ31\Final'};
    reg = {'CA1','mEC','CA1Both'};
elseif strcmp(tag,'mEC') == 1
    mice = {'IZ12\Final','IZ13\Final','IZ15\Final','IZ17\Final',...
        'IZ18\Final','IZ20\Final','IZ21\Final','IZ24\Final','IZ25\Final','IZ26\Final','IZ27\Saline','IZ28\Saline','IZ29\Saline',...
        'IZ30\Final','IZ31\Final','IZ32\Final','IZ33\Final'}; %
    reg = {'CA1','mEC','Both'};
elseif strcmp(tag,'CA3') == 1
    mice = {'IZ27\Final','IZ28\Final','IZ29\Final','IZ32\Final','IZ33\Final'};%%,'IZ34\Final'
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'CA3Saline') == 1
    mice = {'IZ33\Saline','IZ27\Saline','IZ28\Saline','IZ29\Saline','IZ32\Saline'};%,'IZ33\Saline','IZ27\Saline','IZ28\Saline','IZ29\Saline'};%'IZ34\Saline',
    reg = {'CA3','mEC','Both'};
elseif strcmp(tag,'mECBilateral') == 1 
    mice = {'IZ24\Final','IZ25\Final','IZ26\Final'};    
    reg = {'contramEC','ipsimEC','Bilateral'};
end

zone = {'returnB','stemB','delayB','returnS','stemS','delayS'};
target = {'STEM', 'RETURN'};

for tt = 1:length(taglfpName)
    taglfp = taglfpName{tt};
    for rr = 1:length(reg)
        for cc = 1:length(target)
            for zz = 1:length(zone)
               compiledPhaseData{rr,cc}{zz} = [];
               p{rr,cc}{zz} = [];
               r{rr,cc}{zz} = [];
               meanPhase{rr,cc}{zz} = [];
               region{rr,cc}{zz} = [];
               celltype{rr,cc}{zz} = [];
               deepSup{rr,cc}{zz} = [];
               miceID{rr,cc}{zz} = [];
            end
        end
    end

    for m = 1:length(mice)

        cd(strcat(parentDir, mice{m},'\Summ'));
        if exist('PhaseLocking.mat','file')
            load('PhaseLocking.mat');
        else 
            disp(['Phase Locking not computed for mouse' mice{m}])
            continue;
        end

        for ii = 1:length(reg)
            for jj = 1:length(target)
                for kk = 1:length(zone)      
                    for pp = 2:size(phaseData.p.(taglfp){ii,jj}{kk},3) % First array is nans                
                        compiledPhaseData{ii,jj}{kk}  = [compiledPhaseData{ii,jj}{kk} phaseData.phasedistros.(taglfp){ii,jj}{kk}(:,:,pp)];
                        p{ii,jj}{kk} = [p{ii,jj}{kk} phaseData.p.(taglfp){ii,jj}{kk}(:,:,pp)];
                        r{ii,jj}{kk} = [r{ii,jj}{kk} phaseData.r.(taglfp){ii,jj}{kk}(:,:,pp)];
                        meanPhase{ii,jj}{kk} = [meanPhase{ii,jj}{kk} phaseData.m.(taglfp){ii,jj}{kk}(:,:,pp)];
                        region{ii,jj}{kk} = [region{ii,jj}{kk} phaseData.region{ii,jj}{kk}(:,:,pp)];
                        celltype{ii,jj}{kk} = [celltype{ii,jj}{kk} phaseData.putativeCellType{ii,jj}{kk}(:,:,pp)];
                        deepSup{ii,jj}{kk} = [deepSup{ii,jj}{kk} phaseData.deepSuperficial{ii,jj}{kk}(:,:,pp)];
                        miceID{ii,jj}{kk} = [miceID{ii,jj}{kk} ones(1,length(phaseData.p.(taglfp){ii,jj}{kk}(:,:,pp)))*m];
                    end
                end
            end
        end
    end

    % Divide into CA1 pyr cells and INs
    % Divide into DG/ CA3 pyr cells and INs

    % Angle bins
    leg = {'CA1 Pyr','CA1 IN','DG/CA3 Pyr','DG/CA3 IN'};
    colMat = [85/243 85/243 85/243;...
        224/243 163/243 46/243;... 
        8/243 133/243 161/243;...
        56/243 61/243 150/243];  
    edges = deg2rad([0:10:360]);
    centers = edges(2:end) - mean(diff(edges));
    centers2 = rad2deg([centers centers+2*pi]);

    if strcmp(tag,'CA3') == 1 || strcmp(tag,'CA3Saline') == 1
        for ii = 2:length(reg)
            for jj = 1%:length(target)
                % Extract cell id and region info
                b = cellfun(@char,region{ii,jj}{1},'UniformOutput',false);
                [a, ~, iR] = unique(b,'sorted');
                iR = reshape(iR, size(region{ii,jj}{1}));

                % Extract cell id and region info
                b = cellfun(@char,celltype{ii,jj}{1},'UniformOutput',false);
                [c, ~, iC] = unique(b,'sorted');
                iC = reshape(iC, size(celltype{ii,jj}{1}));        

                if length(a{1,1})==1
                    if length(c{1,1})==1
                        idx{1} = iR == 2 & iC == 3 & (deepSup{ii,jj}{1}==1); %CA1 pyr superficial
                        idx{2} = iR == 2 & iC == 3 & (deepSup{ii,jj}{1}==2); % CA1 pyr deep
                        idx{3} = (iR == 3 | iR == 4) & iC == 3; %DG Pyr
                        idx{4} = (iR == 3 | iR == 4) & iC == 2; %DG IN
                    else
                        idx{1} = iR == 2 & iC == 2 & (deepSup{ii,jj}{1}==1); %CA1 pyr superficial
                        idx{2} = iR == 2 & iC == 2 & (deepSup{ii,jj}{1}==2); % CA1 pyr deep
                        idx{3} = (iR == 3 | iR == 4) & iC == 2; %DG Pyr
                        idx{4} = (iR == 3 | iR == 4) & iC == 1; %DG IN    
                    end
                else
                    if length(c{1,1})==1
                        idx{1} = iR == 1 & iC == 3 & (deepSup{ii,jj}{1}==1); %CA1 pyr superficial
                        idx{2} = iR == 1 & iC == 3 & (deepSup{ii,jj}{1}==2); %CA1 pyr deep
                        idx{3} = (iR == 2 | iR == 3) & iC == 3; %DG Pyr
                        idx{4} = (iR == 2 | iR == 3) & iC == 2; %DG IN
                    else
                        idx{1} = iR == 1 & iC == 2 & (deepSup{ii,jj}{1}==1); %CA1 pyr superficial
                        idx{2} = iR == 1 & iC == 2 & (deepSup{ii,jj}{1}==2); %CA1 pyr deep
                        idx{3} = (iR == 2 | iR == 3) & iC == 2; %DG Pyr
                        idx{4} = (iR == 2 | iR == 3) & iC == 1; %DG IN    
                    end
                end
                for cc = 1:4        
                    data{ii,jj,cc}.p = [];
                    data{ii,jj,cc}.m = [];
                    data{ii,jj,cc}.r = [];
                    data{ii,jj,cc}.mice = [];

                    for kk  = [1 4 2 5]
                        data{ii,jj,cc}.p = [data{ii,jj,cc}.p p{ii,jj}{kk}(idx{cc})'];
                        data{ii,jj,cc}.r = [data{ii,jj,cc}.r r{ii,jj}{kk}(idx{cc})'];
                        data{ii,jj,cc}.m = [data{ii,jj,cc}.m meanPhase{ii,jj}{kk}(idx{cc})'];
                        data{ii,jj,cc}.mice = [data{ii,jj,cc}.mice miceID{ii,jj}{kk}(idx{cc})'];
                    end
                end
            end
        end
        
        figure(1)
        set(gcf,'Position',[100 100 1400 600])
        set(gcf,'renderer','painters');
            
        for cc = 1:2
            idxSig = [data{2,1,cc}.p(:,[3 4])<0.05 data{3,1,cc}.p(:,[3 4])<0.05]; 
            idStim = zeros(size(idxSig));
            %idStim(:,[1 3]) = 1; idStim(:,[2 4]) = 2;
            idStim(:,1) = 1; idStim(:,2) = 2; idStim(:,3) = 3; idStim(:,4) = 4;      

            subplot(2,6,6*(cc-1)+1)
            temp = [data{2,1,cc}.r(:,[3 4]) data{3,1,cc}.r(:,[3 4])];
            temp = temp(idxSig);
            tempidStim = idStim(idxSig);
            stats.MRL{cc} = groupStats(temp(tempidStim~=2),tempidStim(tempidStim~=2),'inAxis',true,'color',colMat([1 3 2 4],:));
            title(num2str(sum(idxSig)./size(idxSig,1)))

            temp = [data{2,1,cc}.m(:,[3 4]) data{3,1,cc}.m(:,[3 4])];
            figure(2)
            c1 = polarhistogram(temp(idxSig(:,1),1),edges,'Normalization','probability');
            data1 = c1.Values;
            c2 = polarhistogram(temp(idxSig(:,2),2),edges,'Normalization','probability'); 
            data2 = c2.Values;
            c3 = polarhistogram(temp(idxSig(:,3),3),edges,'Normalization','probability');
            data3 = c3.Values;
            c4 = polarhistogram(temp(idxSig(:,4),4),edges,'Normalization','probability'); 
            data4 = c4.Values;          
            tempdata = temp(idxSig);
            tempId = idStim(idxSig);
            [stats.wwtest{cc}.angle.p stats.wwtest{cc}.angle.table] = circ_wwtest(tempdata(tempId~=2),tempId(tempId~=2));
            
            figure(1)
            subplot(2,6,6*(cc-1)+2)
            hold on        
            plot(centers2, smooth([data1 data1],5), 'color', [85/243 85/243 85/243],'LineWidth',2);
            plot(centers2, smooth([data2 data2],5), 'color', [8/243 133/243 161/243],'LineWidth',2);
            plot(centers2, smooth([data3 data3],5), 'color', [224/243 163/243 46/243],'LineWidth',2);
            plot(centers2, smooth([data4 data4],5), 'color', [56/243 61/243 150/243],'LineWidth',2);     
            set(gca,'YScale','linear','TickDir','out','XTick',[0:180:720]); xlim([0 720]); ylim([0 0.15]);
            ax = axis; xlim([0 720]);
            x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave)/2 +.5) * (ax(4)-ax(3))/5) + ax(3); x_wave = rad2deg(x_wave);
            plot(x_wave, y_wave,'color',[.3 .3 .3]);
            xlabel('\theta phase oriens [deg]'); ylabel('Proportion of cells');
            title(num2str(stats.wwtest{cc}.angle.p))
            %groupStats(temp,[tempidReg tempidStim],'inAxis',true,'color',colMat(ii+1,:)) 
           
            tempMice = [data{2,1,cc}.mice(:,[3 4]) data{3,1,cc}.mice(:,[3 4])];                
            tempBase = tempMice(idxSig(:,1),1);
            tempCA3 = tempMice(idxSig(:,3),3);
            tempBoth = tempMice(idxSig(:,4),4);
            tempAngle = temp(idxSig(:,1),1);
            tempAngleCA3 = temp(idxSig(:,3),3);
            tempAngleBoth = temp(idxSig(:,4),4);
            angleBase = [];
            angleCA3 = [];
            angleBoth = [];
            for mm = 1:length(mice)%unique(tempBase)'
                angleBase(mm,1) = circ_mean(tempAngle(tempBase==mm));
                angleCA3(mm,1) = circ_mean(tempAngleCA3(tempCA3==mm));
                angleBoth(mm,1) = circ_mean(tempAngleBoth(tempBoth==mm));
            end
            angleRat = abs(angleCA3./angleBase);
            angleRat2 = abs(angleBoth./angleBase);
                
            subplot(2,6,6*(cc-1)+3)
            temp = [data{2,1,cc}.m(:,[3 4]) data{3,1,cc}.m(:,[3 4])];
            temp2 = temp(idxSig);
            tempidStim = idStim(idxSig);
            %temp = rad2deg(temp);
            stats.Angle.Center{cc} = groupStats(temp2(tempidStim~=2),tempidStim(tempidStim~=2),'doPlot',false,'inAxis',true,'color',colMat([1 3 2 4],:));
            title(num2str(sum(idxSig)./size(idxSig,1)))            
%             tempdiff(:,1) = min(360 - (temp(:,2)-temp(:,1)),(temp(:,2)-temp(:,1)));
%             tempdiff(:,2) = min(360 - (temp(:,3)-temp(:,1)),(temp(:,3)-temp(:,1)));
%             tempdiff(:,3) = min(360 - (temp(:,4)-temp(:,1)),(temp(:,4)-temp(:,1)));
            c1 = polarhistogram(temp(idxSig(:,1),1),25,'Normalization','probability');
            c1.DisplayStyle = 'stairs';   hold on
            meanAngle = circ_mean(temp(idxSig(:,1),1));
            polarplot([meanAngle meanAngle],[0 .15])  
            title('No stim')
            
            subplot(2,6,6*(cc-1)+4)
            c2 = polarhistogram(temp(idxSig(:,2),2),25,'Normalization','probability');
            c2.DisplayStyle = 'stairs';   hold on
            meanAngle = circ_mean(temp(idxSig(:,2),2));
            polarplot([meanAngle meanAngle],[0 .15])  
            title('mEC')

            subplot(2,6,6*(cc-1)+5)
            c3 = polarhistogram(temp(idxSig(:,3),3),25,'Normalization','probability');
            c3.DisplayStyle = 'stairs';   hold on
            meanAngle = circ_mean(temp(idxSig(:,3),3));
            polarplot([meanAngle meanAngle],[0 .15])  
            title('CA3')
            
            subplot(2,6,6*(cc-1)+6)
            c4 = polarhistogram(temp(idxSig(:,4),4),25,'Normalization','probability');
            c4.DisplayStyle = 'stairs';   hold on
            meanAngle = circ_mean(temp(idxSig(:,4),4));
            polarplot([meanAngle meanAngle],[0 .15])  
            title('Both')            
            
%             subplot(2,7,7*(cc-1)+7)
%             temp = [data{2,1,cc}.m(:,[1 2]) data{3,1,cc}.m(:,[1 2])];
%             temp = temp(idxSig);
%             tempidStim = idStim(idxSig);
%             temp = rad2deg(temp);
%             stats.Angle.Side{cc} = groupStats(temp,tempidStim,'inAxis',true,'color',colMat([1 3 2 4],:));
%             title(num2str(sum(idxSig)./size(idxSig,1)))   
%             title('Side arm')
          
        end        
        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingStats',taglfp,'_DeepSup.png'),'png');
        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingStats',taglfp,'_DeepSup.fig'),'fig');
        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingStats',taglfp,'_DeepSup.eps'),'epsc');
        save(strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingStats',taglfp,'_DeepSup.mat'),'stats');
    else
        for ii = 1:length(reg)
            for jj = 1%:length(target)
                % Extract cell id and region info
                b = cellfun(@char,region{ii,jj}{1},'UniformOutput',false);
                [a, ~, iR] = unique(b,'sorted');
                iR = reshape(iR, size(region{ii,jj}{1}));

                % Extract cell id and region info
                b = cellfun(@char,celltype{ii,jj}{1},'UniformOutput',false);
                [c, ~, iC] = unique(b,'sorted');
                iC = reshape(iC, size(celltype{ii,jj}{1}));        

                if length(a{1,1})==1
                    if length(c{1,1})==1
                        idx{1} = iR == 2 & iC == 3 & (deepSup{ii,jj}{1}==1); %CA1 pyr superficial
                        idx{2} = iR == 2 & iC == 3 & (deepSup{ii,jj}{1}==2); %CA1 pyr deep
                        idx{3} = (iR == 3 | iR == 4) & iC == 3; %DG Pyr
                        idx{4} = (iR == 3 | iR == 4) & iC == 2; %DG IN
                    else
                        idx{1} = iR == 2 & iC == 2 & (deepSup{ii,jj}{1}==1); %CA1 pyr superficial
                        idx{2} = iR == 2 & iC == 2 & (deepSup{ii,jj}{1}==2); %CA1 pyr deep
                        idx{3} = (iR == 3 | iR == 4) & iC == 2; %DG Pyr
                        idx{4} = (iR == 3 | iR == 4) & iC == 1; %DG IN    
                    end
                else
                    if length(c{1,1})==1
                        idx{1} = iR == 1 & iC == 3 & (deepSup{ii,jj}{1}==1); %CA1 pyr superficial
                        idx{2} = iR == 1 & iC == 3 & (deepSup{ii,jj}{1}==2); %CA1 pyr deep
                        idx{3} = (iR == 2 | iR == 3) & iC == 3; %DG Pyr
                        idx{4} = (iR == 2 | iR == 3) & iC == 2; %DG IN
                    else
                        idx{1} = iR == 1 & iC == 2 & (deepSup{ii,jj}{1}==1); %CA1 pyr superficial
                        idx{2} = iR == 1 & iC == 2 & (deepSup{ii,jj}{1}==2); %CA1 pyr deep
                        idx{3} = (iR == 2 | iR == 3) & iC == 2; %DG Pyr
                        idx{4} = (iR == 2 | iR == 3) & iC == 1; %DG IN    
                    end
                end
                for cc = 1:4             
                    data{ii,jj,cc}.p = [];
                    data{ii,jj,cc}.m = [];
                    data{ii,jj,cc}.r = [];
                    data{ii,jj,cc}.mice = [];

                    for kk  = [1 4 2 5]
                        data{ii,jj,cc}.p = [data{ii,jj,cc}.p p{ii,jj}{kk}(idx{cc})'];
                        data{ii,jj,cc}.r = [data{ii,jj,cc}.r r{ii,jj}{kk}(idx{cc})'];
                        data{ii,jj,cc}.m = [data{ii,jj,cc}.m meanPhase{ii,jj}{kk}(idx{cc})'];
                        data{ii,jj,cc}.mice = [data{ii,jj,cc}.mice miceID{ii,jj}{kk}(idx{cc})'];
                    end
                end
            end
        end
        
        figure(1)
        set(gcf,'Position',[100 50 1300 1000])
        set(gcf,'renderer','painters');
            
        for ii = 1:3
            for cc = 1:2
                idxSig = (data{ii,1,cc}.p(:,[3 4])<0.05);
                idStim = zeros(size(idxSig));
                idReg= zeros(size(idxSig));
                %idStim(:,[1 3]) = 1; idStim(:,[2 4]) = 2;
                idStim(:,1) = 1; idStim(:,2) = 2;
                %idReg(:,[1 2]) = 1; idReg(:,[3 4]) = 3;            
                
                subplot(4,6,3*(cc-1)+ii)
                temp = data{ii,1,cc}.r(:,[3 4]);
                temp = temp(idxSig);
                tempidStim = idStim(idxSig);
                stats.MRL{ii,cc} = groupStats(temp,tempidStim,'inAxis',true,'color',colMat([1 ii+1],:));
                [stats.MRL{ii,cc}.ranksum.p,~,stats.MRL{ii,cc}.ranksum.stats] = ranksum(temp(tempidStim==1),temp(tempidStim==2)); 
                title(num2str(sum(idxSig)./size(idxSig,1)))
                ylabel('Mean vector length')
                xticklabels({'No stim', 'Stim'})
                        
                temp = data{ii,1,cc}.m(:,[3 4]);
                tempdata = temp(idxSig);
                tempId = idStim(idxSig);
                [stats.WWTEST{ii,cc}.angle.p stats.WWTEST{ii,cc}.angle.table] = circ_wwtest(tempdata,tempId);        
                
                figure(2)
                c1 = polarhistogram(temp(idxSig(:,1),1),edges,'Normalization','probability');
                data1 = c1.Values;
                c2 = polarhistogram(temp(idxSig(:,2),2),edges,'Normalization','probability'); 
                data2 = c2.Values;
                figure(1)
                subplot(4,6,3*(cc-1)+ii+6)
                hold on        
                plot(centers2, smooth([data1 data1],5), 'color', [0.5 0.5 0.5],'LineWidth',1.5);
                plot(centers2, smooth([data2 data2],5), 'color', colMat(ii+1,:),'LineWidth',1.5);
                set(gca,'YScale','linear','TickDir','out','XTick',[0:180:720]); xlim([0 720]); ylim([0 0.15]);
                ax = axis; xlim([0 720]);
                x_wave = deg2rad([0:0.5:720]); y_wave = ((cos(x_wave)/2 +.5) * (ax(4)-ax(3))/5) + ax(3); x_wave = rad2deg(x_wave);
                plot(x_wave, y_wave,'color',[.3 .3 .3]);
                xlabel('\theta phase slm [deg]'); ylabel('Proportion of cells');
                title(num2str(stats.WWTEST{ii,cc}.angle.p))
                %groupStats(temp,[tempidReg tempidStim],'inAxis',true,'color',colMat(ii+1,:))   
                
                subplot(4,6,3*(cc-1)+12+ii)
                temp = data{ii,1,cc}.m(:,[3 4]);
                
                tempMice = data{ii,1,cc}.mice(:,[3 4]);                
                tempBase = tempMice(idxSig(:,1),1);
                tempStim = tempMice(idxSig(:,2),2);
                tempAngle = temp(idxSig(:,1),1);
                tempAngleStim = temp(idxSig(:,2),2);
                angleBase = [];
                angleStim = [];
                for mm = 1:length(mice)%unique(tempBase)'
                    angleBase(mm,1) = circ_mean(tempAngle(tempBase==mm));
                    angleStim(mm,1) = circ_mean(tempAngleStim(tempStim==mm));
                end
                angleRat = abs(angleStim./angleBase);
                
                %temp = rad2deg(temp);                
                temp2 = temp(idxSig);
                
                tempidStim = idStim(idxSig);
                stats.Angle.Center{ii,cc} = groupStats(temp2,tempidStim,'inAxis',true,'doPlot',false,'color',colMat([1 ii+1],:));  
                [stats.Angle.Center{ii,cc}.ranksum.p, ~,stats.Angle.Center{ii,cc}.ranksum.stats] = ranksum(temp2(tempidStim==1),temp2(tempidStim==2));
               
                c1 = polarhistogram(temp(idxSig(:,1),1),25,'Normalization','probability');
                c1.DisplayStyle = 'stairs';   hold on
                meanAngle = circ_mean(temp(idxSig(:,1),1));
                polarplot([meanAngle meanAngle],[0 .15])  
                title('No stim')
                %hold on
                subplot(4,6,3*(cc-1)+18+ii)
                c2 = polarhistogram(temp(idxSig(:,2),2),25,'Normalization','probability');
                c2.DisplayStyle = 'stairs';   hold on
                meanAngle = circ_mean(temp(idxSig(:,2),2));
                polarplot([meanAngle meanAngle],[0 .15])                
                title('Stim')
                %ylabel('Preferred angle')
                %[stats.Angle.Center{ii,cc}.ranksum.p,~,stats.Angle.Center{ii,cc}.ranksum.stats] = ranksum(temp(tempidStim==1),temp(tempidStim==2)); 
%               
%                 subplot(4,6,3*(cc-1)+18+ii)
%                 temp = data{ii,1,cc}.m(:,[1 2]);
%                 temp = rad2deg(temp);                
%                 temp = temp(idxSig);
%                 tempidStim = idStim(idxSig);
%                 stats.Angle.Side{ii,cc} = groupStats(temp,tempidStim,'inAxis',true,'color',colMat([1 ii+1],:));  
%                 title('Side arm')
%                 ylabel('Preferred angle')
%                 [stats.Angle.Side{ii,cc}.ranksum.p,~,stats.Angle.Side{ii,cc}.ranksum.stats] = ranksum(temp(tempidStim==1),temp(tempidStim==2)); 

            end
        end

        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingStats',taglfp,'_DeepSup.png'),'png');
        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingStats',taglfp,'_DeepSup.fig'),'fig');
        saveas(figure(1),strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingStats',taglfp,'_DeepSup.eps'),'epsc');
        save(strcat(parentDir,'\Compiled\Phase Locking\',tag,'\PhaseLockingStats',taglfp,'_DeepSup.mat'),'stats');
    end
    
end
end