file = dir('*.ArmChoice.Events.mat');
load(file.name);

numtrials = length(armChoice.arm);

for ii = 1:fix(numtrials/10)
    
        if ii < fix(numtrials/10)
            armChoice.successBlocks(ii,1) = sum(armChoice.armAlternation(((10*(ii-1))+2):((10*ii)+1)))/10;
        elseif ii == fix(numtrials/10)
            if mod(numtrials,10) == 0
                armChoice.successBlocks(ii,1) = sum(armChoice.armAlternation(((10*(ii-1))+2):end))/10;
            else
                armChoice.successBlocks(ii,1) = sum(armChoice.armAlternation(((10*(ii-1))+2):((10*ii)+1)))/10;
            end
        end
    
end

save(file.name, 'armChoice');