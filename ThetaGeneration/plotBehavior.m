%script to combine animal behavior summaries
% for ii = 1:3
%     for jj = 1:2
%         arr{ii,jj} = BehavDataSumm15{ii,jj};
%         arr{ii,jj} = [arr{ii,jj}; BehavDataSumm12{ii,jj}];
%     end
% end




% figure
% for ii = 1:size(arr,1)
%     for jj = 1:size(arr,2)
%         subplot(3,2,2*(ii-1)+jj)
%         err = std(arr{ii,jj}(1:3,:))/sqrt(3);
%         errorbar(nanmean(arr{ii,jj}(1:3,:)),err,'LineWidth',2,'Color',colormap(3,:));
%         xlabel('Trial Blocks')
%         ylabel('Performance')
%         hold on
%         f1 = scatter([1:1:8],arr{ii,jj}(1,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));
%         f1.MarkerFaceAlpha = 0.5;
%         f1 = scatter([1:1:8],arr{ii,jj}(2,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));
%         f1.MarkerFaceAlpha = 0.5;
%         f1 = scatter([1:1:8],arr{ii,jj}(3,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));        
%         f1.MarkerFaceAlpha = 0.5;
%         
%         err = std(arr{ii,jj}(4:6,:))/sqrt(3);
%         errorbar(nanmean(arr{ii,jj}(4:6,:)),err,'LineWidth',2,'Color',colormap(4,:));
%         xlabel('Trial Blocks')
%         ylabel('Performance')
%         hold on
%         f1 = scatter([1:1:8],arr{ii,jj}(4,:),'MarkerFaceColor',colormap(4,:),'MarkerEdgeColor',colormap(4,:));
%         f1.MarkerFaceAlpha = 0.5;
%         f1 = scatter([1:1:8],arr{ii,jj}(5,:),'MarkerFaceColor',colormap(4,:),'MarkerEdgeColor',colormap(4,:));
%         f1.MarkerFaceAlpha = 0.5;
%         f1 = scatter([1:1:8],arr{ii,jj}(6,:),'MarkerFaceColor',colormap(4,:),'MarkerEdgeColor',colormap(4,:));        
%         f1.MarkerFaceAlpha = 0.5;
%         
%         err = std(arr{ii,jj}(7:end,:))/sqrt(3);
%         errorbar(nanmean(arr{ii,jj}(7:end,:)),err,'LineWidth',2,'Color',colormap(6,:));
%         xlabel('Trial Blocks')
%         ylabel('Performance')
%         hold on
%         f1 = scatter([1:1:8],arr{ii,jj}(7,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));
%         f1.MarkerFaceAlpha = 0.5;
%         f1 = scatter([1:1:8],arr{ii,jj}(8,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));
%         f1.MarkerFaceAlpha = 0.5;
%         if ii==3 || jj ==2
%             f1 = scatter([1:1:8],arr{ii,jj}(9,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));        
%             f1.MarkerFaceAlpha = 0.5;
%         end
%         ylim([0.5 1.1])
%         xlim([0 9])
%         title(strcat('Region: ',reg{ii},' Zone:',zone{jj}));
%     end
% end

figure
for ii = 1:size(arr,1)
    for jj = 1:size(arr,2)
        if ii <3
            
            subplot(3,2,2*(ii-1)+jj)
            err = std(arr{ii,jj}(1:3,:))/sqrt(3);
            errorbar(nanmean(arr{ii,jj}(1:3,:)),err,'LineWidth',2,'Color',colormap(3,:));
            xlabel('Trial Blocks')
            ylabel('Performance')
            hold on
            f1 = scatter([1:1:2],arr{ii,jj}(1,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));
            f1.MarkerFaceAlpha = 0.5;
            f1 = scatter([1:1:2],arr{ii,jj}(2,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));
            f1.MarkerFaceAlpha = 0.5;
            f1 = scatter([1:1:2],arr{ii,jj}(3,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));        
            f1.MarkerFaceAlpha = 0.5;

    %         err = std(arr{ii,jj}(4:6,:))/sqrt(3);
    %         errorbar(nanmean(arr{ii,jj}(4:6,:)),err,'LineWidth',2,'Color',colormap(4,:));
    %         xlabel('Trial Blocks')
    %         ylabel('Performance')
    %         hold on
    %         f1 = scatter([1:1:2],arr{ii,jj}(4,:),'MarkerFaceColor',colormap(4,:),'MarkerEdgeColor',colormap(4,:));
    %         f1.MarkerFaceAlpha = 0.5;
    %         f1 = scatter([1:1:2],arr{ii,jj}(5,:),'MarkerFaceColor',colormap(4,:),'MarkerEdgeColor',colormap(4,:));
    %         f1.MarkerFaceAlpha = 0.5;
    %         f1 = scatter([1:1:2],arr{ii,jj}(6,:),'MarkerFaceColor',colormap(4,:),'MarkerEdgeColor',colormap(4,:));        
    %         f1.MarkerFaceAlpha = 0.5;
    %
            err = std(arr{ii,jj}(4:end,:))/sqrt(size(arr{ii,jj},1)-3);
            errorbar(nanmean(arr{ii,jj}(4:end,:)),err,'LineWidth',2,'Color',colormap(6,:));
            xlabel('Trial Blocks')
            ylabel('Performance')
            hold on
            f1 = scatter([1:1:2],arr{ii,jj}(4,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));
            f1.MarkerFaceAlpha = 0.5;
            if (size(arr{ii,jj},1)-3) == 2
                f1 = scatter([1:1:2],arr{ii,jj}(5,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));
                f1.MarkerFaceAlpha = 0.5;
            end
            if (size(arr{ii,jj},1)-3) == 3
                f1 = scatter([1:1:2],arr{ii,jj}(6,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));        
                f1.MarkerFaceAlpha = 0.5;
            end
    %
    %
    %         err = std(arr{ii,jj}(7:end,:))/sqrt(3);
    %         errorbar(nanmean(arr{ii,jj}(7:end,:)),err,'LineWidth',2,'Color',colormap(6,:));
    %         xlabel('Trial Blocks')
    %         ylabel('Performance')
    %         hold on
    %         f1 = scatter([1:1:2],arr{ii,jj}(7,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));
    %         f1.MarkerFaceAlpha = 0.5;
    %         f1 = scatter([1:1:2],arr{ii,jj}(8,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));
    %         f1.MarkerFaceAlpha = 0.5;
    %         if ii==3 || jj ==2
    %             f1 = scatter([1:1:2],arr{ii,jj}(9,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));        
    %             f1.MarkerFaceAlpha = 0.5;
    %         end
            ylim([0.5 1.1])
            xlim([0 3])
            title(strcat('Region: ',reg{ii},' Zone:',zone{jj}));
            
        else
            
            subplot(3,2,2*(ii-1)+jj)
            err = std(arr{ii,jj}(1:2,:))/sqrt(2);
            errorbar(nanmean(arr{ii,jj}(1:3,:)),err,'LineWidth',2,'Color',colormap(3,:));
            xlabel('Trial Blocks')
            ylabel('Performance')
            hold on
            f1 = scatter([1:1:2],arr{ii,jj}(1,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));
            f1.MarkerFaceAlpha = 0.5;
            f1 = scatter([1:1:2],arr{ii,jj}(2,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));
            f1.MarkerFaceAlpha = 0.5;
%             f1 = scatter([1:1:2],arr{ii,jj}(3,:),'MarkerFaceColor',colormap(3,:),'MarkerEdgeColor',colormap(3,:));        
%             f1.MarkerFaceAlpha = 0.5;

            err = std(arr{ii,jj}(3:end,:))/sqrt(size(arr{ii,jj},1)-2);
            errorbar(nanmean(arr{ii,jj}(3:end,:)),err,'LineWidth',2,'Color',colormap(6,:));
            xlabel('Trial Blocks')
            ylabel('Performance')
            hold on
            f1 = scatter([1:1:2],arr{ii,jj}(3,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));
            f1.MarkerFaceAlpha = 0.5;
            if (size(arr{ii,jj},1)-3) == 2
                f1 = scatter([1:1:2],arr{ii,jj}(4,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));
                f1.MarkerFaceAlpha = 0.5;
            end
            if (size(arr{ii,jj},1)-3) == 3
                f1 = scatter([1:1:2],arr{ii,jj}(5,:),'MarkerFaceColor',colormap(6,:),'MarkerEdgeColor',colormap(6,:));        
                f1.MarkerFaceAlpha = 0.5;
            end

            ylim([0.5 1.1])
            xlim([0 3])
            title(strcat('Region: ',reg{ii},' Zone:',zone{jj}));
        end
    end
end
