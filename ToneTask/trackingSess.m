function framesDropped = trackingSess

warning('off','all')
allSess = dir('*_sess*');
framesDropped(1:size(allSess,1)) = nan;

for ii = 7:(size(allSess,1))
    cd(strcat(allSess(ii).folder,'\',allSess(ii).name));
    tracking = getToneTracking('forceReload',true);
    framesDropped(ii) = tracking.framesDropped{1}; 
end

end