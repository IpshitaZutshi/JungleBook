function [trajectory_numbers_direction, trajectory_nums_no_direction] = getTrajectories(lickedPorts)

ports = 1:7;  

%% IGNORE DIRECTION

combinations = nchoosek(ports, 2);
% lickedPorts = behavTrials.port;
trajectory_nums_no_direction = zeros(length(lickedPorts), 1); 

for k = 2:length(lickedPorts) 
    previous_lick = lickedPorts(k-1);
    current_lick = lickedPorts(k);
    
    sorted_pair = sort([previous_lick, current_lick]);
    for n = 1:size(combinations, 1)
        if combinations(n, 1) == sorted_pair(1) && ...
           combinations(n, 2) == sorted_pair(2)
           trajectory_nums_no_direction(k) = n; 
           break;
        end
    end
end


%% TAKE INTO ACCOUNT DIRECTION

permutations = reshape(combinations(:,perms(1:2)),[],2);
trajectory_numbers_direction = zeros(length(lickedPorts), 1);

for k = 2:length(lickedPorts)
    previous_lick = lickedPorts(k-1);
    current_lick = lickedPorts(k);
    
    for n = 1:length(permutations)
        if permutations(n, 1) == previous_lick && ...
           permutations(n, 2) == current_lick
            trajectory_numbers_direction(k) = n;
            break;
        end
    end
end
end