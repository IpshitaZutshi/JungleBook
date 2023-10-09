function plotPGAMExampleFits(rowloc,colloc,results,fighandle,numrows,numcol, cellNum)

allCN = [results.neuron];
neuIdx = find(allCN==(cellNum-1));

idx = [3 1 4 10];

%% Loop to plot ylin, yfwd, disttoStop, licks, distribution of licksKernels
for ii = 1:length(idx)
    %Plot model rates
    subplot(numrows, numcol, [(rowloc-1)*numcol+colloc+(ii-1)*numcol (rowloc-1)*numcol+colloc+1+(ii-1)*numcol], 'Parent', fighandle);
    plot(results(neuIdx(idx(ii))).x,results(neuIdx(idx(ii))).raw_rate_Hz,'Color',[0.5 0.5 0.5])
    hold on
    plot(results(neuIdx(idx(ii))).x,results(neuIdx(idx(ii))).model_rate_Hz,'Color','r','LineWidth',1)
    title(strcat(num2str(results(neuIdx(idx(ii))).pval)))
    box off
    xlim([results(neuIdx(idx(ii))).x(1) results(neuIdx(idx(ii))).x(end)])
    xticks([])
    %Plot kernels
    if ii <4
        subplot(numrows, numcol, [(rowloc-1)*numcol+2+colloc+(ii-1)*numcol (rowloc-1)*numcol+colloc+3+(ii-1)*numcol], 'Parent', fighandle);
        a = results(neuIdx(idx(ii))).kernel_x';        
        posfill = results(neuIdx(idx(ii))).kernel_pCI';
        negfill = results(neuIdx(idx(ii))).kernel_mCI';
        fill([a;flipud(a)],[negfill;flipud(posfill)],[0.5 0.5 0.5],'linestyle','none','FaceAlpha',0.5);                    
        hold on
        plot(results(neuIdx(idx(ii))).kernel_x,results(neuIdx(idx(ii))).kernel,'Color',[0.5 0.5 0.5],'LineWidth',1)
        title(strcat(num2str(results(neuIdx(idx(ii))).mutual_info)))
        xlim([results(neuIdx(idx(ii))).kernel_x(1) results(neuIdx(idx(ii))).kernel_x(end)])
        xticks([])
        box off
    else
        subplot(numrows, numcol, [(rowloc-1)*numcol+2+colloc+(ii-1)*numcol (rowloc-1)*numcol+colloc+3+(ii-1)*numcol], 'Parent', fighandle);
        for kk = 5:9
            data(kk-4) = results(neuIdx(kk)).signed_kernel_strength;
        end
        bar(data)
        box off
        title(strcat(num2str(results(neuIdx(10)).signed_kernel_strength)))
    end
end

end