% Finds the position to the spikes
function [xSp,ySp,tSp,iSp] = bz_spikePos2(tSp0,x,y,t)
nSp = length(tSp0);
count = 0;
xSp = []; ySp = []; iSp = []; tSp = [];
for sp = 1:nSp
    tdiff = (t-tSp0(sp)).^2;
    [m,ind] = min(tdiff);
    % Check if spike occurs within a short time interval from a valid
    % position point
    if m < 0.04 %0.2 *0.2 sec
        count = count + 1;
        xSp = cat(1,xSp,x(ind(1)));
        ySp = cat(1,ySp,y(ind(1)));
        tSp = cat(1,tSp,tSp0(sp));
        iSp = cat(1,iSp,ind(1));
    end
end