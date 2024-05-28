 function [cor, lag, smooth] = SmoothCor(cor, lag, t_bin)
    
    std_smooth_kernel = .01;
    
    kernel = pdf('normal', -std_smooth_kernel*10/t_bin:std_smooth_kernel*10/t_bin, 0, std_smooth_kernel/t_bin); % convolve with gaussian

    if isempty(cor), smooth = []; return; % if no signals, then send it back
    
    else
        
        smooth = zeros(size(cor));
        
        for i = 1:size(cor,2)  

            smooth(:,i) = ndnanfilter(cor(:,i), kernel(:), []);
        
        end
        
        % lag = lag(round(end/2):end,:);
        % 
        % cor = cor(round(end/2):end,:);
        % 
        % smooth = smooth(round(end/2):end, :);

    end 
    
end