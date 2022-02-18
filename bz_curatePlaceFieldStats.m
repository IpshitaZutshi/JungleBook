function [placeFieldStats_curated] = bz_curatePlaceFieldStats(placeFieldStats, positions, spikes, varargin)

% USAGE
% [placeFieldStats_curated] = bz_curatePlaceFieldStats(placeFieldStats, varargin)
% Takes already computed placeFieldStats and curates them. It takes out
% incomplete place fields (too near the maze sides) and place fields with 
% more than one peak. Thought to be used when computing phase precession
%
% INPUTS
%
%   placeFieldStats         Output from bz_findPlaceFields1D
%   positions               (1 x #conditions) cell array, where each cell is a
%                           Nx2 array with timestamps (sec) and linear positions
%   <options>               optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties        Values
%    -------------------------------------------------------------------------
%     'input_mapstats'      placeFieldStats.mapStats can be input instead
%                           of placeFieldStats, if this variable is set to
%                           true. If true, it will return mapStats
%                           structure.
%                           field is going to be discarded
%     'min_slope'           Threshold (% of maximum slope) below which place 
%                           field is going to be discarded
%     'height_boundaries'   Threshold (% of maximum firing rate) above which
%                           two local maxima will be treated as one
%     'do_plot'             Plot firing maps of all units and all conditions in
%                           one figure, both old (blue) and new (red) PFs 
%    =========================================================================
%
%   OUTPUTS
%
%   placeFieldStats cellinfo structure with the following fields
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%		.UID                unit ids
%		.sessionName        name of session
%		.params             parameters:
%         .sizeMaze
%         .threshold
%         .minSize
%         .maxSize              
%         .sepEdge
%         .minPeak
%         .minPeak2nd
%         .verbose
%         .saveMat
%       .mapStats           Statistics of the Firing Map
%         .x
%         .field                (updated)
%         .size                 (updated)
%         .peak
%         .mean
%         .fieldX               (updated)
%         .specificity
%         .m
%         .r
%         .mode
%         .k
%         .y
%         .fieldY
%
%    =========================================================================

% Andrea Navas-Olive, 2021

    p=inputParser;
    addParameter(p,'input_mapstats',false,@islogical);
    addParameter(p,'min_slope',0.3,@isnumeric);
    addParameter(p,'height_boundaries',0.8,@isnumeric);
    addParameter(p,'do_plot',false,@islogical);

    parse(p,varargin{:});
    input_mapstats = p.Results.input_mapstats;
    min_slope = p.Results.min_slope;
    height_boundaries = p.Results.height_boundaries;
    do_plot = p.Results.do_plot;


    % Define new variables
    % - if first input is placeFieldStats.mapStats
    if input_mapstats
        boundaries = placeFieldStats;
    % - if first input is placeFieldStats
    else
        % Check it has the mapStats field
        if ~isfield(placeFieldStats,'mapStats')
            warning('First input is not "placeFieldStats". Make sure "input_mapstats=true" if you have input "placeFieldStats.mapStats" instead of "placeFieldStats"')
        end
        boundaries = placeFieldStats.mapStats;
    end
    n_units = length(boundaries);
    n_conditions = length(positions);

    % Boundaries to be changed
    boundaries_curated = boundaries;

    % Take number of bins used for this session
    nBins_all = [];
    for unit = 1:n_units
        for c = 1:n_conditions   
            nBins_all = [nBins_all size(boundaries{unit}{c}.field,1)];
        end
    end
    nBins = mode(nBins_all(nBins_all>0));

    % Plot
    if do_plot
        figure('pos', [10,10,1500,1200])
        subx = ceil(sqrt(n_units));
        suby = ceil(n_units/subx);
    end

    % Go through all firing maps
    for unit = 1:n_units        
        if do_plot
            subplot(subx, suby, unit)
            hold on
        end
        for c = 1:n_conditions        
            for bd = 1:size(boundaries{unit}{c}.fieldX, 1)

                % Boudary place fields (normalized between 0 and 1)
                boundaryData = boundaries{unit}{c}.fieldX(bd, :) / nBins;
                % Firing curve
                curve = FiringCurve(positions{c}, spikes.times{unit}, 'nbins',nBins, 'smooth',2, 'type','ll');
                curvex_norm = curve.x / max(curve.x);
                curvey_norm = curve.rate / max(curve.rate);

                % FIRST CURATION: if place field is too near the maze sides
                % so it is incomplete, we remove it
                is_not_placefield = isnan(boundaryData(1));
                if boundaries{unit}{c}.fieldX(bd,1)>0
                    is_in_border1 = (boundaryData(1) < 0.03) && (curve.rate(boundaries{unit}{c}.fieldX(bd,1)) > max(curve.rate)/2);
                    is_in_border2 = (boundaryData(2) > 0.97) && (curve.rate(boundaries{unit}{c}.fieldX(bd,2)) > max(curve.rate)/2);
                else
                    is_in_border1  = 1;
                end
                if is_not_placefield || is_in_border1 || is_in_border2
                    if bd == 1
                        % Plot firing map even if there is no place field
                        boundaries_curated{unit}{c}.fieldX = [nan nan];
                        boundaries_curated{unit}{c}.field = [];
                        boundaries_curated{unit}{c}.size = 0;
                    end
                    if do_plot
                        plot(curve.x/max(curve.x), c + curve.rate/max(curve.rate), 'k')
                    end
                    continue
                end

                % Indexes of place field
                idxsPF = find(boundaries{unit}{c}.field(:,bd)==1);
                idxsPF_curated = idxsPF;
                field_curated = boundaries{unit}{c}.field(:,bd);
                boundaryData_curated = boundaries{unit}{c}.fieldX(bd,:);
                center_curated = round(mean(boundaryData_curated));
                size_curated = round(diff(boundaryData_curated));

                % SECOND CURATION: If more than two peaks, take maximum
                % and the two nearest local minima
                if length(findpeaks_Andrea(curvey_norm(idxsPF)).loc ) > 1
                    % Maximum
                    [~, idxmax] = max(curvey_norm(idxsPF));
                    idxmax = idxsPF(idxmax);
                    % Take local minima nearest to peak
                    % - right
                    local_minima_before = find(diff(diff(curvey_norm(1:idxmax))>0)==1);
                    if isempty(local_minima_before)
                        idxmin1 =  boundaries{unit}{c}.fieldX(bd,1);
                    else
                        idxmin1 = 1 + local_minima_before(end);
                    end
                    % - left
                    local_minima_after = find(diff(diff(curvey_norm(idxmax:end))>0)==1);
                    if isempty(local_minima_after)
                        idxmin2 =  boundaries{unit}{c}.fieldX(bd,end);
                    else
                        idxmin2 = idxmax-1 + local_minima_after(1);
                    end
                    % Check than minimums low enough to be considered two
                    % different peaks and not noise
                    if (curvey_norm(idxmin1) < height_boundaries) && (curvey_norm(idxmin2) < height_boundaries)
                        % New boundaries
                        idxsPF_curated = [idxmin1:idxmin2];
                        boundaryData_curated = [idxmin1 idxmin2];
                        field_curated = zeros(nBins,1);
                        field_curated(idxsPF_curated) = 1;
                        center_curated = round(mean(boundaryData_curated));
                        size_curated = round(diff(boundaryData_curated));
                    end
                end
                
                % THIRD CURATION: if boundaries are too flat, remove it
                dydx = [abs(diff(curvey_norm)) 0];
                max_dydx = max(abs(dydx));
                % - Side 1
                idxs_PFside1 = round([boundaryData_curated(1) : center_curated-size_curated/10]);
                if any( dydx(idxs_PFside1) < max_dydx*min_slope )
                    % Define new boundary, as the first index with big slope
                    idxup1 = idxs_PFside1(find(dydx(idxs_PFside1) >= max_dydx*min_slope, 1, 'first'));
                    % If all slopes are less than required, then take the
                    % most conservative choice
                    if isempty(idxup1)
                        idxup1 = idxs_PFside1(end);
                    end
                    % New boundaries
                    idxsPF_curated = [idxup1:boundaryData_curated(2)];
                    boundaryData_curated(1) = idxup1;
                    field_curated = zeros(nBins,1);
                    field_curated(idxsPF_curated) = 1;
                    center_curated = round(mean(boundaryData_curated));
                    size_curated = round(diff(boundaryData_curated));
                end
                % - Side 2
                idxs_PFside2 = round([center_curated+size_curated/10 : boundaryData_curated(end)]);
                if any( dydx(idxs_PFside2) < max_dydx*min_slope )
                    % Define new boundary, as the last index with big slope
                    idxup2 = idxs_PFside2(find(dydx(idxs_PFside2) >= max_dydx*min_slope, 1, 'last'));
                    % If all slopes are less than required, then take the
                    % most conservative choice
                    if isempty(idxup2)
                        idxup2 = idxs_PFside2(1);
                    end
                    % New boundaries
                    idxsPF_curated = [boundaryData_curated(1):idxup2];
                    boundaryData_curated(2) = idxup2;
                    field_curated = zeros(nBins,1);
                    field_curated(idxsPF_curated) = 1;
                    center_curated = round(mean(boundaryData_curated));
                    size_curated = round(diff(boundaryData_curated));
                end
                
                % Plot
                if do_plot            
                    % Original Firing map
                    if bd == 1
                        plot(curvex_norm, c+curvey_norm, 'k')
                    end
                    plot(curvex_norm(idxsPF), c+curvey_norm(idxsPF), '-b', 'linewidth', 2)
                    % New firing map
                    plot(curvex_norm(idxsPF_curated), c+curvey_norm(idxsPF_curated), '-r')
                    % Axis
                    ylim([1 1+n_conditions])
                    title(sprintf('Unit %d', unit))
                end

                % Save 
                % - ini and end
                boundaries_curated{unit}{c}.fieldX(bd, :) = boundaryData_curated;
                % - binary vector
                boundaries_curated{unit}{c}.field(:,bd) = field_curated;
                % - mean and size
                boundaries_curated{unit}{c}.x(bd) = center_curated;
                boundaries_curated{unit}{c}.size(bd) = size_curated;
            end
        end
    end
    
    if input_mapstats
        placeFieldStats_curated = boundaries_curated;
    else
        placeFieldStats_curated = placeFieldStats;
        placeFieldStats_curated.mapStats = boundaries_curated;
    end
    
    
end 