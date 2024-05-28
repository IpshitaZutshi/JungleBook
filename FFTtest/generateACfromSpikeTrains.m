%% Plot the autocorrelation. The spike train must be a column vector or else
%% it gives an error

function generateACfromSpikeTrains(tSp)

t_bin = 0.005; %% Time resolution for the autocorrelation

[cor, lag] = CrossCorr(tSp, tSp);
[cor, lag, smooth] = SmoothCor(cor, lag, t_bin); % solve for smoothed signal, and chop off center peak and negative half of autocorrelation
[S, f] = getChronuxSpectrum(cor, t_bin);

% Sall = S(f>6.2&f<14);
% fall = f(f>6.2&f<14);
% [~,freqP] = max(Sall);
% freqPeak(currf) = fall(freqP);

figure
bar(lag, cor, 'FaceColor', [.2 .2 1], 'EdgeColor', [ 0 0 .75]), hold on;
line(lag, smooth, 'Color', 'k', 'linewidth', 1.5);
        
if lag(end)>0
    xlim([0 lag(end)]);
end

if max(smooth)>0
    ylim([0 max(smooth)*1.1]);
end

end

function [cor, lag, smooth] = SmoothCor(cor, lag, t_bin)
    
    std_smooth_kernel = .005;
    
    kernel = pdf('normal', -std_smooth_kernel*10/t_bin:std_smooth_kernel*10/t_bin, 0, std_smooth_kernel/t_bin); % convolve with gaussian

    if isempty(cor), smooth = []; return; % if no signals, then send it back
    
    else
        
        smooth = zeros(size(cor));
        
        for i = 1:size(cor,2)  

            smooth(:,i) = ndnanfilter(cor(:,i), kernel(:), []);
        
        end
        
        lag = lag(round(end/2):end,:);
        
        cor = cor(round(end/2):end,:);
        
        smooth = smooth(round(end/2):end, :);

    end 
    
end

function [S, f] = getChronuxSpectrum(cor, t_bin)

    params.tapers = [1.25, 1];
    params.fpass = [];
    params.Fs = 1/t_bin;
    params.pad = 6;  

    [S,f]=mtspectrumc(cor,params); % compute spectrum of cor

end