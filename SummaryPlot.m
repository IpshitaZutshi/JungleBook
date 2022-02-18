comparisons = [1 2 3 4 5 6 7];

X(X==0) = nan;
colnum = [4:16];    
colTitle = {'CSD pyr','CSD rad','CSD slm','Gamma pyr','Gamma rad','Gamma slm','rateCA1Pyr','rateCA1IN','rateCA3Pyr','rateCA3IN','PVcorr','phaseangle','phaseprecess'};
diffTable(length(comparisons),length(colnum)) = nan;

for comp = 1:length(comparisons)
    for col = colnum
       %load data
        A=X((X(:,1)==comp),col);
        diffTable(comp,col-3) = nanmean(A);
    end
end

figure
subplot(1,3,1)
imagesc(diffTable(:,1:6),'AlphaData',~isnan(diffTable(:,1:6)))
YlGnBu=cbrewer('seq', 'YlGnBu', 25);
YlGnBu = YlGnBu(25:-1:1,:);
colormap(YlGnBu)
caxis([0.4 1.1])

subplot(1,3,2)
imagesc(diffTable(:,7:10),'AlphaData',~isnan(diffTable(:,7:10)))
YlGnBu=cbrewer('seq', 'BuPu', 25);
YlGnBu = YlGnBu(25:-1:1,:);
colormap(YlGnBu)
caxis([0.2 1.1])

subplot(1,3,3)
imagesc(diffTable(:,11:13),'AlphaData',~isnan(diffTable(:,11:13)))
YlGnBu=cbrewer('seq', 'YlGnBu', 25);
YlGnBu = YlGnBu(25:-1:1,:);
colormap(YlGnBu)
caxis([0.4 1.1])