function [xtileInds,xtileBounds,countXtiles,meanXtile,stdXtile] = Xtiles(values,xtile,space)
%% Sort ripple duration and amplitude by desired xtile 
% Input:
% vector of values
% xtile - desired percentile number of devisions (e.g. quartiles = 4) 
% space - 'log' or 'lin' 

% Output: 
% xtileInds - xtile ind for each ripple, from small to large 

% Later: write option to do in log or lin space, for now just doing in
% logbase 10 space!!!
% Write more general version of this for any vector 

%% Determine xtile xtileBounds

z=[0:100/xtile:100];
for ii = 1:length(z)-1
    ints(ii,1) = z(ii);
    ints(ii,2) = z(ii+1)-.01;
end

%% Split into xtiles

if strcmp(space,'log')
    % Do in log base 10 space 
    durationlog=log10(values); %10.^values = actual durations 
    % Find xtileBounds
    for dd = 1:xtile
        pdurlog(dd,1) = prctile(durationlog,ints(dd,1)); 
        pdurlog(dd,2) = prctile(durationlog,ints(dd,2));
    end
    % Values of each timebin
    xtileBounds = reshape(10.^(reshape(pdurlog,numel(pdurlog),1)),xtile,2);
    
else strcmp(space,'lin')
    % Find xtileBounds
    for dd = 1:xtile
        xtileBounds(dd,1) = prctile(values,ints(dd,1)); 
        xtileBounds(dd,2) = prctile(values,ints(dd,2));
    end
end

% Find indices of ripples within each of these values - in ms
xtileInds = nan(numel(values),1);
for i=1:xtile
    r=find(values >= xtileBounds(i,1) & values <= xtileBounds(i,2));
    xtileInds(r) = i;
    countXtiles(i) = length(r); 
    meanXtile(i) = mean(values(r)); 
    stdXtile(i) = std(values(r)); 
end

% Plot 


end
