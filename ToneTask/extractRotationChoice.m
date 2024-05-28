file = dir('*.rotation.mat');

for ii = 1:length(file)

    load(file(ii).name)
    locs = trackRotation.lickLoc(trackRotation.linTrial==0);
    rot = trackRotation.direction(trackRotation.linTrial==0);

    for rr = 1:6
        clock(ii,rr) = sum((rot(locs==rr))==1)./length(rot(locs==rr));
        counterclock(ii,rr) = sum((rot(locs==rr))==-1)./length(rot(locs==rr));        
    end
end