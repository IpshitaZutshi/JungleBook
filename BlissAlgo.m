
comparisons = {[1 2 3],[1 4 5],[1 6 7]};
X(X==0) = nan;
colnum = [4:16];    
colTitle = {'CSD pyr','CSD rad','CSD slm','Gamma pyr','Gamma rad','Gamma slm','rateCA1Pyr','rateCA1IN','rateCA3Pyr','rateCA3IN','PVcorr','phaseprecess','phasepref'};
diffTable = [];
diffTable(length(comparisons),length(colnum)) = nan;
pTable(length(comparisons),length(colnum)) = nan;
figure
set(gcf,'Position',[10 10 1500 800])
plotnum = 1;
for comp = 1:length(comparisons)
    compTable = comparisons{comp};
    for col = colnum
       %load data
        A=X((X(:,1)==compTable(1)& X(:,3)==comp),col);
        B=X((X(:,1)==compTable(2)& X(:,3)==comp),col);
        D=X((X(:,1)==compTable(3)& X(:,3)==comp),col);     
%         lA = log(A);
%         lB = log(B);
%         lD = log(D);
        n1=length(A);n2=length(B);n3=length(D);
        m=catpad(2,A,B,D);
        indep =[];
        k=0;
        indep = A.*B;
%         for i1 = 1:n1
%             for i2 = 1:n2
%                k=k+1;
%                indep(k)=m(i1,1)+m(i2,2);
%             end
%         end
        %indep = exp(indep);
        datam = catpad(2,A,B,D,indep);
        subplot(length(comparisons),length(colnum),plotnum)
        boxplot(datam)
        xticklabels({'Stim1','Stim2','Observed','Expected'});
        title(colTitle{col-3})
        plotnum = plotnum+1;
        if nanmedian(D) < 1 && nanmedian(indep) < 1
            diffTable(comp,col-3) = nanmedian(indep-D);
        else
            diffTable(comp,col-3) = nan;
        end
        %diffTable(comp,col-3) = nanmean(indep./D);
        if sum(~isnan(indep))>0
            [pTable(comp,col-3)] = signrank(indep,D);
        else
            pTable(comp,col-3) = nan;
        end
    end
end

figure
subplot(1,3,1)
imagesc(diffTable(:,1:6),'AlphaData',~isnan(diffTable(:,1:6)))
YlGnBu=cbrewer('div', 'PiYG', 40);
YlGnBu(YlGnBu<0) = 0;
YlGnBu = YlGnBu(40:-1:1,:);
colormap(YlGnBu)
caxis([-0.3 0.33])
colorbar

subplot(1,3,2)
imagesc(diffTable(:,7:10),'AlphaData',~isnan(diffTable(:,7:10)))
YlGnBu=cbrewer('div', 'PiYG', 11);
YlGnBu(YlGnBu<0) = 0;
YlGnBu = YlGnBu(11:-1:1,:);
colormap(YlGnBu)
caxis([-0.3 0.3])
colorbar

subplot(1,3,3)
imagesc(diffTable(:,11:13),'AlphaData',~isnan(diffTable(:,11:13)))
YlGnBu=cbrewer('div', 'PiYG', 30);
YlGnBu(YlGnBu<0) = 0;
YlGnBu = YlGnBu(30:-1:1,:);
colormap(YlGnBu)
caxis([-0.3 0.3])
colorbar
% 
% subplot(2,1,2)
% imagesc(pTable)
% caxis([0 0.2])

