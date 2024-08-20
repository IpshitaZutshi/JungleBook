function FigS17_UMAPdeliberation

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[600 32 1147 905]);

numrow = 5;
numcol = 4;

sess = {'IZ39_220714_sess18','IZ43_220828_sess4','IZ44_220830_sess7','IZ47_230707_sess24','IZ48_230628_sess17'};
decodingSess = {'IZ39\Final\IZ39_220714_sess18','IZ43\Final\IZ43_220828_sess4','IZ44\Final\IZ44_220830_sess7','IZ47\Final\IZ47_230707_sess24','IZ48\Final\IZ48_230628_sess17'};

A = [-181.1478,-182,363.1123,-74,179];
E = [90,-2,6.4,90,90];

framelag = [-30 0 15];


for ss = 1:length(sess)

    sessName  = strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\',decodingSess{ss});
    cd(sessName)
    file = dir('*.Tracking.Behavior.mat');
    load(file.name)
    file = dir('*.TrialBehavior.Behavior.mat');
    load(file.name)
    Dec = findDecelerationPoints('plotfig',false);

    umap_path = strcat('Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\',sess{ss},'\manifold');
    cd(umap_path)     
    umap_name = 'Umap_behavior_speed_1_smooth_5_bin_0.1.csv';
    file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
    behav_file = file.name;
    
    acc = gradient(tracking.position.v)/0.033;
    for ii = 1:4
        if ii < 4
            decTS = Dec.ts(Dec.decType==ii);
        else
            % If 4, randomly find timepoints when the mouse is near a port
            numsample = sum(Dec.decType==2);
            toneTrialidx = find(behavTrials.linTrial==0);            
            idxLocs = find(((tracking.position.y>28 & tracking.position.y<34) | (tracking.position.y>53 & tracking.position.y<59) | ...
            (tracking.position.y>76 & tracking.position.y <83) | (tracking.position.y >100 & tracking.position.y <108)) & ...
            (tracking.timestamps>behavTrials.timestamps(toneTrialidx(1),1) & tracking.timestamps<behavTrials.timestamps(toneTrialidx(end),2)) & ...
            tracking.position.vy>25 & acc>5);
            subLocs = randsample(idxLocs,numsample);
            decTS = tracking.timestamps(subLocs);
        end
        plotwithinTrialUMAP(umap_name,behav_file,A(ss), E(ss), decTS,tracking,fig2,numrow,numcol,4*(ss-1)+ii,framelag)
    end

end

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure17_UMAPDeliberation.png'));
saveas(gcf,strcat(expPath,'SupFigure17_UMAPDeliberation.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure17_UMAPDeliberation.fig'));


end

function plotwithinTrialUMAP(umap_name,behav_file,A, E, decTS,tracking,fighandle,numrows,numcol,plotloc,framelag)

TRIAL_TYPE = 0:5;
% load Umap result
Umap_results = readtable(umap_name);
Umap_results = table2array(Umap_results);
    
% load position direction and other information
load(behav_file);

plot_ind = [];
    
for tt = 1:length(TRIAL_TYPE)
    plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & probe_ds==0 & correct_ds==1 & position_y_all>2 & speed_all'>2)];   
end

% Plot manifold in gray
ax1 = subplot(numrows,numcol,plotloc,'Parent',fighandle);
scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),3,[0.9 0.9 0.9],'filled','MarkerFaceAlpha',1);
view(A,E)
grid off;
axis off;
axis tight
hold on   

% Now plot the timepoints of interest
for dt = 1:length(decTS)
    [~,frameidx] = min(abs(tracking.timestamps-decTS(dt)));
    startTime = tracking.timestamps(frameidx+framelag(1));
    endTime = tracking.timestamps(frameidx+framelag(3));
    
    [~,startidx] = min(abs(timestamp_beh-startTime));
    [~,endidx] = min(abs(timestamp_beh-endTime));
    plot_ind_pos = zeros(1,length(position_y_all));
    plot_ind_pos(startidx:1:endidx) = 1;
    plot_ind_final = plot_ind_pos & position_y_all>1 & speed_all' >1;
    tsAxis = linspace(1,sum(plot_ind_final),sum(plot_ind_final));
    scatter3(Umap_results(plot_ind_final,1),Umap_results(plot_ind_final,2),Umap_results(plot_ind_final,3),3,tsAxis','filled');    
    colormap(ax1,'magma');    
end
%colorbar

end
