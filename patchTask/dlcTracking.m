
function [DLC_tracking] = dlcTracking(varargin)

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'plotfig',true,@islogical);
addParameter(p,'forceRun',true,@islogical);
%addParameter(p,'updatedIntan',true,@islogical);

parse(p,varargin{:});
saveMat = p.Results.saveMat;
plotfig = p.Results.plotfig;
forceRun = p.Results.forceRun;
%updatedIntan = p.Results.updatedIntan;

basepath = pwd;
[~, currentFolderName] = fileparts(basepath);


% dlc_file = dir('*T22'); 
% tracking_data = csvread(dlc_file, 2, 1);
tracking_data = csvread('\\research-cifs.nyumc.org\research\buzsakilab\Buzsakilabspace\LabShare\ZutshiI\Cue_task-Lucy-2025-02-24\T22_41022_150631 - T22_41022_150631.csv', 2, 0);

for i = 1:size(tracking_data, 1)
    if tracking_data(i, 4) >= 0.98 && tracking_data(i, 7) >= 0.98 && tracking_data(i, 10) >= 0.98 ...
            && tracking_data(i, 13) >= 0.98 && tracking_data(i, 16) >= 0.98
        cutoff = i;
        break
    else
        continue
    end
end

tracking_data = tracking_data(cutoff:end, :);

tracking_data_og=tracking_data;

sz = size(tracking_data, 1);
cut_rows = [];

% delete uncertainty
j = 1;
while j <= sz
    if tracking_data(j, 13) < 0.98
        tracking_data(j, :) = [];
        sz = sz - 1;
    else
        j = j + 1;
    end
end

tracking_data_og_2=tracking_data;

j=2;
while j < size(tracking_data, 1)
    if tracking_data(j, 11) >= tracking_data(j-1, 11)+10 || tracking_data(j, 12) >= tracking_data(j-1, 12)+10 ...
            tracking_data(j, 11) <= tracking_data(j-1, 11)-10 || tracking_data(j, 12) <= tracking_data(j-1, 12)-10;
        start_cut = j;
        cut_pt = tracking_data(j, 11);
        
        k=j+1;
        % while k < size(tracking_data, 1)
        %     if tracking_data(k, 11) >= tracking_data(k-1, 11)+10 || tracking_data(k, 12) >= tracking_data(k-1, 12)+10 || ...
        %             tracking_data(k, 11) <= tracking_data(k-1, 11)-10 || tracking_data(k, 12) <= tracking_data(k-1, 12)-10
        %         break;  % End of bad block (jump back)
        %     else
        %         k = k + 1;
        %     end
        % end

        while k < size(tracking_data, 1)
            if tracking_data(k, 12) <= 850 && tracking_data(k, 11) > 680 ... 
                    && tracking_data(k, 11) < 730 % these are rough numbers -> and only for this session
                break;  % End of bad block (jump back)
            elseif tracking_data(k, 12) > 850 && tracking_data(k, 12) < 890 
                break;  % End of bad block (jump back)
            else
                k = k + 1;
            end
        end

        end_cut = k - 1;

        % Delete the bad block
        tracking_data(start_cut:end_cut, :) = [];
        cut_rows = [cut_rows; (start_cut:end_cut)'];

        % Adjust j after deletion
        j = start_cut;
        j = j + 1;
    else
        j = j + 1;
    end
end

%{
figure
colormap(jet)
plot(tracking_data(:, 10), tracking_data(:, 11), tracking_data(:, 10))

x_pos = tracking_data(:, 10);
y_pos = tracking_data(:, 11);
c = linspace(1, 10, length(x_pos));

figure
h = surface([x_pos; x_pos], [y_pos: y_pos], [zeros(size(x_pos)); zeros(size(x_pos))], [c; c], ...
    'FaceColor', 'no', 'EdgeColor', 'interp', 'LineWidth', 2);
colormap(jet)
%}


figure
colormap(jet)
scatter(tracking_data(:,11), tracking_data(:,12), 20, tracking_data(:,11), 'filled') % 20 is point size
colorbar

end



















