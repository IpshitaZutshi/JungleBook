function map = bz_calRateMap(spikeMap,timeMap,minOccupancy)
if nargin<3 || isempty(minOccupancy)
    minOccupancy = 0.15;
end
box = [0.0025 0.0125 0.0200 0.0125 0.0025;...
    0.0125 0.0625 0.1000 0.0625 0.0125;...
    0.0200 0.1000 0.1600 0.1000 0.0200;...
    0.0125 0.0625 0.1000 0.0625 0.0125;...
    0.0025 0.0125 0.0200 0.0125 0.0025];
timeMap(timeMap < minOccupancy) = NaN; %eliminate pixels with very low occupancy, Nov. 30, 2010; Increased from 0.01 to 0.150, Jan 10, 2014
rawMap = spikeMap./timeMap;
rawMap = padarray(rawMap,[2 2],NaN);
for i = 3:size(rawMap,1)-2
    for j = 3:size(rawMap,2)-2
        current_sum = 0;
        current_box = 0;
        box_i = 0;
        for ii = i-2:i+2
            box_i = box_i+1;
            box_j = 0;
            for jj = j-2:j+2
                box_j = box_j+1;
                if ~isnan(rawMap(ii,jj))
                    current_sum = current_sum+rawMap(ii,jj)*box(box_i,box_j);
                    current_box = current_box+box(box_i,box_j);
                end
            end
        end
        map(i,j) = current_sum/current_box;
    end
end
map = map(3:end,3:end);

