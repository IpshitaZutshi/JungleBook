function ripple_HSE_score_super = decode_replay_w_score_wrapper_iz(options,basepath,saving_folder,ripple_HSE_dec,template_beh,behavior,intervals,direction)
%% 1.1 get transition  matrix
[trans_mtrx, ~] = behavior_trans_mtrx(behavior,'interval',intervals,'position_bins',size(template_beh,2),'dim_names', {'lin'});

%  figure; imagesc(trans_mtrx)
%  colorbar;

%% 1.2 get state distance matrix
dist_states_matrix_super = get_states_distance_matrix(trans_mtrx,options);
% figure; imagesc(dist_states_matrix_super);
% colorbar;

%caxis([0 10]);
%% 1.3 calclate z score for both distance correlolation and weighted correlation



ripple_HSE_score_super =  caculateCorrandZ_super(ripple_HSE_dec.spikecounts,dist_states_matrix_super,template_beh,options);


%% 1.4 calculate Monte Carlo p-value 

for e = 1:length(ripple_HSE_score_super)
    [MonteP_super_Rw(1),~] = SignificanceTest(abs(ripple_HSE_score_super(e).Rw), abs(ripple_HSE_score_super(e).Rw_null{1}));
    [MonteP_super_Rw(2),~] = SignificanceTest(abs(ripple_HSE_score_super(e).Rw), abs(ripple_HSE_score_super(e).Rw_null{2}));
    [MonteP_super_Rw(3),~] = SignificanceTest(abs(ripple_HSE_score_super(e).Rw), abs(ripple_HSE_score_super(e).Rw_null{3}));
    %[p_super_Rw(2,1),~] = SignificanceTest(abs(example_results_super(2).Rw), abs(example_results_super(2).Rw_null{1}));
    ripple_HSE_score_super(e).MonteP_Rw = MonteP_super_Rw;
    
    
    [MonteP_super_Rwd(1),~] = SignificanceTest(ripple_HSE_score_super(e).Rwd, ripple_HSE_score_super(e).Rwd_null{1});
    [MonteP_super_Rwd(2),~] = SignificanceTest(ripple_HSE_score_super(e).Rwd, ripple_HSE_score_super(e).Rwd_null{2});
    [MonteP_super_Rwd(3),~] = SignificanceTest(ripple_HSE_score_super(e).Rwd, ripple_HSE_score_super(e).Rwd_null{3});


    ripple_HSE_score_super(e).MonteP_Rwd = MonteP_super_Rwd;

    %[p_super_Rwd(2,1),~] = SignificanceTest(example_results_super(2).Rwd, example_results_super(2).Rwd_null{1});
end

% [p_unsuper_Rwd(1,1),~] = SignificanceTest(example_results_unsuper(1).Rwd, example_results_unsuper(1).Rwd_null{1});
% [p_unsuper_Rwd(2,1),~] = SignificanceTest(example_results_unsuper(2).Rwd, example_results_unsuper(2).Rwd_null{1});
% T = table(p_super_Rw,p_super_Rwd,'Rownames',{'Example1';'Example2'})
% T = table(p_super_Rw,p_super_Rwd,'Rownames',{'Example1';'Example2'});

%% 1.5 save the result


basename = bz_BasenameFromBasepath(basepath);
save([basepath,saving_folder,basename,'.ripple_HSE_score_super_',direction,'.mat'],'ripple_HSE_score_super')


end
