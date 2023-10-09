function [results_dec, template_all] = bz_decodeBehavior_wy(basepath,spikes,behavior,varargin)
% Build decoder and test via 
% TODO: Make this generalizable by adding custom field selection for
% behavior
%
% INPUTS
%   'spikes'    buzcode compatible 'spikes.cellinfo' struct
%   'behavior'  buzcode compatible '.behavior' struct
%
%   (optional)
%	'algorithm'	Currently supported: 'bayes', 'PVcorr', default 'bayes'
%   'dim_names' Names of the behavioral dimensions to be decoded. These
%               have to be subfields of behavior.position. By default, this
%               will use every non-empty subfield of position.
%   'dim_bins'  Can be ndims x 1 cell with the edges for dimension bins or
%               an ndims x 1 vector with the desired number of bins for
%               each dimension, which will be automatically spread from
%               the min- to max-value.
%   'intervals' Subset of behavior ts to be used for decoding
%   'wdw'       Decoding time bin width in seconds, default 0.05s
%   'frequency' Frequency at which to sample and decode position for
%               testing the decoder, default 1Hz
%   'Ncrossval' Number of cross-validations, default 10-fold
%   'UIDs'      Sub-select specific units from the dataset for decoding.
%               Use all units by default.
%
% OUTPUT
%   Outvar:   	Description
%
% EXAMPLE CALLS
% [] = DefaultTemplateBuzcode(pwd)
%
% Thomas Hainmueller, 2020, Buzsakilab
% Winnie Yang, 2021, Buzsakilab 
% TODO: make better plots for decoding result 
%% Input handling
p = inputParser;
addParameter(p,'algorithm','bayes',@ischar);
addParameter(p,'intervals',[0 Inf],@isnumeric);
addParameter(p,'dim_names',[],@iscell);
addParameter(p,'dim_bins',[]);
addParameter(p,'frequency',1,@isnumeric);
addParameter(p,'Ncrossval',10,@isnumeric);
addParameter(p,'UIDs',[],@isnumeric);
addParameter(p,'minspeed',2,@isnumeric); % minimum speed / second (default 5 m/s)
addParameter(p,'templ_filters',{},@iscell); % nfilters x 2 cell with filter - dimension pairs of 1D filters.
addParameter(p,'time_smooth',0.5,@isnumeric); % Temporal smoothing / binning (seconds);
addParameter(p,'min_occup',0,@isnumeric);
addParameter(p,'templ_offset',.01,@isnumeric); % Addition factor to template for avoiding 0 cancellation.
addParameter(p,'dim_dtypes',{},@iscell); % For error calculation
addParameter(p,'plot_results',true,@islogical);
addParameter(p,'save_plots',true,@islogical);

parse(p,varargin{:})

algorithm = p.Results.algorithm;
intervals = p.Results.intervals;
dim_names = p.Results.dim_names;
dim_bins = p.Results.dim_bins;
frequency = p.Results.frequency;
Ncrossval = p.Results.Ncrossval;
UIDs = p.Results.UIDs;
minspeed = p.Results.minspeed;
templ_filters = p.Results.templ_filters;
time_smooth = p.Results.time_smooth;
min_occup = p.Results.min_occup;
templ_offset = p.Results.templ_offset;
dim_dtypes = p.Results.dim_dtypes;
plot_results = p.Results.plot_results;
save_plots = p.Results.save_plots;

%% Checking inputs
if ~isempty(dim_bins) && length(dim_bins) ~= length(dim_names)
    warning('Bins must be given for all or none of the dimensions.')
    return
end

if ~isempty(dim_dtypes) && length(dim_dtypes) ~= length(dim_names)
    warning('Data Types must be specified for all or none of the dimensions.');
    return
end

if ~iscell(dim_bins) && ~isnumeric(dim_bins)
    warning('Dim_bins must be either numeric or cell')
    return
end

%% Generating defaults
tidcs = find(InIntervals(behavior.timestamps,intervals));
wdw = 1/behavior.samplingRate;

if isempty(dim_names)
    fields = fieldnames(behavior.position);
    for f = 1:length(fields)
        if ~isempty(behavior.position.(fields{f}))
            dim_names{end+1} = fields{f};
        end
    end
end
ndims = length(dim_names);

if isempty(dim_bins)
    % Default: 20 bins per dimension.
    dim_bins(1:ndims) = 20;
end

if isempty(dim_dtypes)
    dim_dtypes(1:ndims) = {'lin'};
end

% Convert dim_bins to cell with edges.
if isnumeric(dim_bins)
    for d = 1:ndims
        % use the 2.5 percentile of distribution to determine lower bound
        tempsort = sort(behavior.position.(dim_names{d})(tidcs));
        lbound = tempsort(round(0.025*length(tempsort)));
        ubound = max(behavior.position.(dim_names{d})(tidcs)) +...
             .001*abs(max(behavior.position.(dim_names{d})(tidcs)));
        % Set slightly higher to include edges, particular for categorical! 
        new_dim_bins{d} = linspace(lbound,ubound,dim_bins(d)+1);
    end
    dim_bins = new_dim_bins;
end

% Select all units if no sub-selection is given.
if isempty(UIDs)
    if isfield(spikes,'ids')
        UIDs = 1:length(spikes.ids);
    elseif isfield(spikes,'UID')
        UIDs = 1:length(spikes.UID);
    else
        warning('Spikes has no recognized unit index field, provide UIDs manually')
        return
    end
end

%% Get firing rates in the relevant bins

% Get spike counts for the decoding time bins using uniform, nonoverlapping 
% bins at behavior acquisition frequency.
% Alternative with arbitrary bins: results.spikecounts = get_spikecounts(spikes,behavior.timestamps(smp_idcs),wdw);
spk_edges = [behavior.timestamps(tidcs)' behavior.timestamps(tidcs(end))+1/behavior.samplingRate];
for n = length(UIDs):-1:1
    spikecounts(n,:) = histcounts(spikes.times{UIDs(n)},spk_edges);
end

% Normalize counts to accomodate uneven bins, dropped frames and recording gaps.
spikecounts = spikecounts ./diff(spk_edges) / behavior.samplingRate;

% Smooth spikecounts to 'simulate' wider bins
spikecounts = movmean(spikecounts,behavior.samplingRate*time_smooth,2);

% Find non-moving periods for exclusion
for d = length(dim_names):-1:1
    v(d,:) = abs([0 diff(behavior.position.(dim_names{d})(tidcs))'])*behavior.samplingRate; % TODO: rho is still explicit!
    v(d,:) = conv(v(d,:),fspecial('gauss',[1 round(3*behavior.samplingRate)],behavior.samplingRate),'same'); % Filter at 1Hz standard
end
v = sqrt(nansum(v.^2,1)); % Euclidean distance
moving = v>minspeed & (behavior.position.(dim_names{d})(tidcs)>lbound)';

% Sub-sample timestamps for increased performance
% smp_idcs = round(tidcs(1) : behavior.samplingRate / frequency : tidcs(end));
% smp_idcs = intersect(smp_idcs,tidcs(moving));
% [~,spk_idcs] = intersect(tidcs,smp_idcs); % Spikecounts has len tidcs
smp_idcs = tidcs(moving);
[~,spk_idcs] = intersect(tidcs,smp_idcs); % Spikecounts has len tidcs

%% Initialize results struct
% Copy subsampled timestamps and position to results.
results_dec.frequency = frequency;
results_dec.UIDs = UIDs;
results_dec.timestamps = behavior.timestamps(smp_idcs);
results_dec.spikecounts = spikecounts(:,spk_idcs);%smp_idcs); % Changed 201010

% Copy position vectors
for d = 1:ndims
    results_dec.position.(dim_names{d}) = behavior.position.(dim_names{d})(smp_idcs);
    dim_size(d) = length(dim_bins{d})-1;
    dim_idcs(d) = {1:dim_size(d)};
end

% Initialize posterior-probability matrices
% template.count = zeros([dim_size size(spikecounts,1)]);

results_dec.pPost = zeros([dim_size length(smp_idcs)]);
results_dec.pPostM = zeros(length(dim_size),length(smp_idcs));
results_dec.pPost_linear = zeros(prod(dim_size), length(smp_idcs));
results_dec.pPostM_linear = zeros(1,length(smp_idcs));

%% Iterative cross-validation

for tp = 1:Ncrossval
    % Split timestamps into template and decoding using logical array
    test_idcs = false(length(smp_idcs),1);
    test_idcs(round((tp-1)/Ncrossval*length(smp_idcs))+1:...
        round(tp/Ncrossval*length(smp_idcs))) = true;
    templ_idcs = ~test_idcs;
    res_idcs = find(test_idcs);
    
    % Generating template
    for d = 1:ndims
        position_v(:,d) = behavior.position.(dim_names{d})(smp_idcs(templ_idcs));
    end
    
    template = buildTemplate(behavior.timestamps(smp_idcs(templ_idcs)),...
        results_dec.spikecounts(:,templ_idcs), position_v,'position_bins',...
        dim_bins, 'fsample', behavior.samplingRate, 'dim_names', dim_names,...
        'min_occup', min_occup, 'filters', templ_filters, 'offset', templ_offset);

    results_dec.templates{tp} = template;
    
%     % Filter template
%     for f = 1:size(templ_filters,1)
%         template = filter_template(template,templ_filters{f,1},...
%             templ_filters{f,2},'preserve_zeros',false);
%     end
%     
    % Linearize and post-process template
    templ_dims = size(template.rate);
    templ_dims = templ_dims(1:end-1); % Last dimension is cell number
    templ_v = template.linear.rate;
    templ_v = templ_v./max(templ_v,[],1);
    [~,cidcs] = max(templ_v,[],1);
    [~,cidcs] = sort(cidcs);
    results_dec.templates{tp}.linear.rate  = templ_v(:,cidcs)';
    
%     templ_v = templ_v + templ_offset;
    
    % Decoding
    if strcmp(algorithm,'bayes')
        [pPostM, pPost] = bayes_decode(templ_v, results_dec.spikecounts(:,res_idcs), wdw);
    elseif strcmp(algorithm,'PVcorr')
        [pPostM, pPost] = PVcorr_decode(templ_v, results_dec.spikecounts(:,res_idcs));
    end
    
    % Append decoding results_dec to struct
    results_dec.pPost_linear(:,res_idcs) = pPost;
    results_dec.pPostM_linear(1,res_idcs) = pPostM;
    
    % Reshape to original format
    [pPostM2{1:length(templ_dims)}] = ind2sub(templ_dims,pPostM);
    pPostM2 = cat(1,pPostM2{:});
    results_dec.pPost(dim_idcs{:},res_idcs) = reshape(pPost,[templ_dims length(res_idcs)]);
    results_dec.pPostM(:,res_idcs) = pPostM2;

    clear pPostM2 position_v
end
%% generate template using all data, used as template for decoding events later
position_v = [];
for d = 1:ndims
    position_v(:,d) = behavior.position.(dim_names{d})(smp_idcs);
end

template = buildTemplate(behavior.timestamps(smp_idcs(:)),...
    results_dec.spikecounts(:,:), position_v,'position_bins',...
    dim_bins, 'fsample', behavior.samplingRate, 'dim_names', dim_names,...
    'min_occup', min_occup, 'filters', templ_filters, 'offset', templ_offset);
templ_v = template.linear.rate;
% templ_v = templ_v./max(templ_v,[],1);
template_all = templ_v';
% [~,cidcs] = max(templ_v,[],1);
% [~,cidcs] = sort(cidcs);
% template_all = templ_v(:,cidcs)';

%% Calculate errors
for d = 1:size(results_dec.pPostM,1)
    dec_pos = dim_bins{d}(results_dec.pPostM(d,:));
    true_pos = results_dec.position.(dim_names{d})';
    rand_pos = randsample(true_pos,length(true_pos));
    results_dec.decodedPosition.(dim_names{d}) = dec_pos;
    
    switch dim_dtypes{d}
        case 'lin'
            results_dec.errors(d,:) = abs(dec_pos - true_pos);
            results_dec.errors_rand(d,:) = abs(dec_pos - rand_pos);
%              results_dec.errors(d,:) = abs(results_dec.decodedPosition.(dim_names{d})-...
%                 results_dec.position.(dim_names{d})');
        case 'cat'
            results_dec.errors(d,:) = double(dec_pos - true_pos ~= 0);
            results_dec.errors_rand(d,:) = double(dec_pos - rand_pos ~= 0);
%              results_dec.errors(d,:) = double(results_dec.decodedPosition.(dim_names{d})-...
%                 results_dec.position.(dim_names{d})' ~= 0);
        case 'circ'
            results_dec.errors(d,:) = min(cat(1, abs(dec_pos - true_pos), ...
                mod(dec_pos + true_pos, dim_bins{d}(end))),[],1);
            results_dec.errors_rand(d,:) = min(cat(1, abs(dec_pos - rand_pos), ...
                mod(dec_pos + rand_pos, dim_bins{d}(end))),[],1);
%             results_dec.errors(d,:) = min(cat(1,...
%                 abs(results_dec.decodedPosition.(dim_names{d})-...
%                 results_dec.position.(dim_names{d})'), ...
%                 mod(results_dec.decodedPosition.(dim_names{d})+...
%                 results_dec.position.(dim_names{d})', dim_bins{d}(end))),[],1);
    end
end     
%%
if ~exist(strcat(basepath,'\BayesDecoding'), 'dir')
   mkdir('BayesDecoding')
end
%% Plot results_dec
if plot_results
    
    % Get regular timestamps and linear behavior
    ts_plot = (0:length(results_dec.timestamps)-1)/results_dec.frequency;
    
    for d = 1:ndims
        position_res{d} = discretize(results_dec.position.(dim_names{d}),dim_bins{d});
    end
    results_dec.position.linear = sub2ind(size(template.time),position_res{:});
    
    % Display last template
    h{1} = figure('color','white');
    imagesc(template_all);
    xlabel('Spatial bin #');
    ylabel('Cell ID');
    saveas(gcf,[basepath,'\BayesDecoding\template.pdf'])
    
    % Display linear position
    h{2} = figure('color','white');
    imagesc(ts_plot,1:prod(dim_size),results_dec.pPost_linear);
    % y-axis label with the different dimensions would be great here.
    caxis([0 prctile(results_dec.pPost_linear,99,'all')]); colormap(flipud(gray));
    hold on
    scatter(ts_plot,results_dec.position.linear,'r','.');
    scatter(ts_plot,results_dec.pPostM_linear,'g','.');
    xlabel('Time (s)');
    ylabel('Spatial bin #');
    saveas(gcf,[basepath,'\BayesDecoding\linearPosition.pdf'])
    
    % Display results_dec in separate dimensions
    h{3} = figure('color','white');
    h{4} = figure('color','white');
    for d = 1:ndims
        cat_dims = 1:ndims;
        cat_dims(d) = [];
        
        % Plot results_dec separate for each dimension
        figure(h{3})
        ax{d} = subplot(ndims,1,d); hold on;
        if ~isempty(cat_dims)
            im_plot = squeeze(sum(results_dec.pPost,cat_dims));
        else % Special case 1D
            im_plot = results_dec.pPost;
        end
        imagesc(ts_plot, dim_bins{d}(1:end-1), im_plot);
        colormap(flipud(gray)); caxis([0 prctile(im_plot,95,'all')]);
        plot(ts_plot,results_dec.position.(dim_names{d}),'r');
        scatter(ts_plot,results_dec.decodedPosition.(dim_names{d}),'g');
        plot(ts_plot,results_dec.errors(d,:));
        ylim([min(dim_bins{d}(1),0) dim_bins{d}(end)]);
        xlim([0 max(ts_plot)]);
        ylabel(dim_names{d});
        xlabel('Time (s)');
        legend('Actual Positon','Decoded Position','Decoding error')
        saveas(gcf,[basepath,'\BayesDecoding\template.pdf'])
        
        % Plot error distributions
        figure(h{4})
        subplot(ndims,1,d); hold on;
        pl(:,1) = cdfplot(results_dec.errors(d,:));
        pl(:,2) = cdfplot(results_dec.errors_rand(d,:));
        set(pl(:,1),'Color','r');
        set(pl(:,2),'Color','k');
        title(sprintf('%s - decoding errors',dim_names{d}))
        xlabel(sprintf('%s error (cm)',dim_names{d}));
        ylabel('CDF');
        saveas(gcf,[basepath,'\BayesDecoding\cdf_behavior.pdf'])
        
        % TODO: Confusion matrices
    end
    
    linkaxes([ax{:}],'x');
end

end