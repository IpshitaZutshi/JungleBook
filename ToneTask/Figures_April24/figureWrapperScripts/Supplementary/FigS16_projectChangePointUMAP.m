function FigS16_projectChangePointUMAP

fig2  = figure;
set(fig2,'Renderer','painters')
set(fig2,'Color','w')
set(fig2,'Position',[600 58 671 905]);

numrow = 5;
numcol = 2;

sess = {'IZ39_220714_sess18','IZ43_220828_sess4','IZ44_220830_sess7','IZ47_230707_sess24','IZ48_230628_sess17'};
decodingSess = {'IZ39\Final\IZ39_220714_sess18','IZ43\Final\IZ43_220828_sess4','IZ44\Final\IZ44_220830_sess7','IZ47\Final\IZ47_230707_sess24','IZ48\Final\IZ48_230628_sess17'};

A = [-181.1478,-182,363.1123,-74,179];
E = [90,-2,6.4,90,90];

col = [83/255 0/255 0/255;...
    184/255 15/255 10/255;...
    241/255 114/255 42/255;...
    249/255 197/255 81/255;...
    143/255 189/255 107/255;...
    87/255 116/255 144/255];

for ss = 1:length(sess)
    umap_path = strcat('Z:\Buzsakilabspace\LabShare\WinnieYang\Ipshita\NatureRevisions\',sess{ss},'\manifold');
    cd(umap_path)     
    umap_name = 'Umap_behavior_speed_1_smooth_5_bin_0.1.csv';
    file = dir('*.position_behavior_speed_1_smooth_5_bin_0.1.mat');
    behav_file = file.name;
    
    % Load behavior file
    load(strcat('Z:\Homes\zutshi01\Recordings\Auditory_Task\',decodingSess{ss},'\',sess{ss},'.TrialBehavior.Behavior.mat'))
    if ~isfield(behavTrials,'probe')
        behavTrials.probe(1:length(behavTrials.linTrial))=0;
    end

    % Load decoding
    decodingPath = strcat('Z:\Homes\zz737\ipshita_data\Auditory_Task\',decodingSess{ss},'\py_data\theta_decoding_lickLoc_y\up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].nc');
    changePointName = strcat('Z:\Homes\zz737\ipshita_data\Auditory_Task\',decodingSess{ss},'\py_data\theta_decoding_lickLoc_y\change_point_posterior_up_samp_binsize[0.01]movement_var[25]sticky_p[0.999].mat');

    posterior_goal =  ncread(decodingPath,'x_position') ;
    posterior_pos =  ncread(decodingPath,'y_position') ;
    post_time =  ncread(decodingPath,'time') ;
    post_pos =  ncread(decodingPath,'y_position_value') ;
    post_goal =  ncread(decodingPath,'x_position_value') ;

    load(changePointName);
    
    ts_dec = [];
    trial_dec = [];

    %% Within each trial, find the change point time that the goal is decoded 
     for tt = 1:length(behavTrials.lickLoc)
        if behavTrials.linTrial(tt)==0 && behavTrials.probe(tt) == 0 && behavTrials.correct(tt)==1
    
            [~,idxstart] = min(abs(post_time-behavTrials.timestamps(tt,1)));
            if post_time(idxstart)<behavTrials.timestamps(tt,1) %Take the next index
                idxstart = idxstart+1;
            end        
            [~,idxend] = min(abs(post_time-behavTrials.timestamps(tt,2)));
            if post_time(idxend)>behavTrials.timestamps(tt,2) %Take the previous index
                idxend = idxend-1;
            end   
    
            [~,decGoal] = max(posterior_goal(:,idxstart:idxend));
    
            %% Get last change point for that trial
            if sum(trial==(tt-1))>0
                curChanges = change_point{trial==(tt-1)};
    
                idxGoal = curChanges(end);
                trialDecGoal = mode(decGoal(curChanges(end)+1:end));
    
                if trialDecGoal==(behavTrials.lickLoc(tt)+1)             
                    ts_dec = [ts_dec post_time(idxGoal+idxstart)];
                    trial_dec = [trial_dec tt];
                end
    
            end
        end
     end

     %% load Umap result
    Umap_results = readtable(umap_name);
    Umap_results = table2array(Umap_results);

    %% load position direction and other information
    load(behav_file);
    TRIAL_TYPE = 0:5;
    plot_ind = [];
    
    for tt = 1:length(TRIAL_TYPE)
        plot_ind =  [plot_ind,find(lick_loc_ds==TRIAL_TYPE(tt) & probe_ds==0 & correct_ds ==1 & position_y_all>5 & speed_all'>1)];   
    end    
       
    lick_plot = lick_loc_ds; 
    lick_plot(isnan(lick_plot))=0; % deal with nan
    lick_plot = lick_plot(plot_ind);

    ax2 = subplot(numrow,numcol,numcol*(ss-1)+1);
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,lick_plot,'filled','MarkerFaceAlpha',1);
    colormap(ax2,col);
    view(A(ss),E(ss))
    grid off;
    axis off;
    axis tight
    
    subplot(numrow,numcol,numcol*(ss-1)+2);
    scatter3(Umap_results(plot_ind,1),Umap_results(plot_ind,2),Umap_results(plot_ind,3),5,[0.9 0.9 0.9],'filled','MarkerFaceAlpha',1);
    view(A(ss),E(ss))
    grid off;
    axis off;
    axis tight
    hold on   
    
    for ii = 1:length(ts_dec)
        [~,idxUMAP] = min(abs(timestamp_beh-ts_dec(ii)));
        colid = behavTrials.lickLoc(trial_dec(ii))+1;
        subplot(numrow,numcol,numcol*(ss-1)+2)
        scatter3(Umap_results(idxUMAP-1,1),Umap_results(idxUMAP-1,2),Umap_results(idxUMAP-1,3),35,col(colid,:),'filled');
    end

end

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\Figures_April2024\SuppFigures\';
saveas(gcf,strcat(expPath,'SupFigure16_projectChangePoint.png'));
saveas(gcf,strcat(expPath,'SupFigure15_projectChangePoint.eps'),'epsc');
saveas(gcf,strcat(expPath,'SupFigure15_projectChangePoint.fig'));


end

