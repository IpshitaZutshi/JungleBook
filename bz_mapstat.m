% Shannon information, sparseness, and selectivity
function [information,sparsity,selectivity] = bz_mapstat(map,posPDF);
n = size(map,1);
meanrate = nansum(nansum( map .* posPDF ));
meansquarerate = nansum(nansum( (map.^2) .* posPDF ));
if meansquarerate == 0
    sparsity = NaN;
else
    sparsity = meanrate^2 / meansquarerate;
end
maxrate = max(max(map));
if meanrate == 0;
    selectivity = NaN;
else
    selectivity = maxrate/meanrate;
end
[i1, i2] = find( (map>0) & (posPDF>0) );  % the limit of x*log(x) as x->0 is 0
if ~isempty(i1)
    akksum = 0;
    for i = 1:length(i1);
        ii1 = i1(i);
        ii2 = i2(i);
        akksum = akksum + posPDF(ii1,ii2) * (map(ii1,ii2)/meanrate) * log2( map(ii1,ii2) / meanrate );
    end
    information = akksum;
else
    information = NaN;
end