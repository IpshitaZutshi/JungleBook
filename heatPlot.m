function map = heatPlot(x,y,xRange,yRange,mapSize,mapSizeSmooth,kernal)
if nargin<7 || isempty(kernal)
    kernal =  [0.0025 0.0125 0.0200 0.0125 0.0025;...
        0.0125 0.0625 0.1000 0.0625 0.0125;...
        0.0200 0.1000 0.1600 0.1000 0.0200;...
        0.0125 0.0625 0.1000 0.0625 0.0125;...
        0.0025 0.0125 0.0200 0.0125 0.0025];
end
if length(mapSize)==1
    mapSize = [mapSize mapSize];
end
if length(mapSizeSmooth)==1
    mapSizeSmooth = [mapSizeSmooth mapSizeSmooth];
end
xBinSize = diff(xRange)/mapSize(1);
yBinSize = diff(yRange)/mapSize(2);
xEdges = xRange(1):xBinSize:xRange(2);
yEdges = yRange(1):yBinSize:yRange(2);
xBinSizeSmooth = diff(xRange)/mapSizeSmooth(1);
yBinSizeSmooth = diff(yRange)/mapSizeSmooth(2);
xEdgesSmooth = xRange(1):xBinSizeSmooth:xRange(2);
yEdgesSmooth = yRange(1):yBinSizeSmooth:yRange(2);

[~,xBin] = histc(x,xEdges);
[~,yBin] = histc(y,yEdges);
map = zeros(length(yEdges),length(xEdges));
nY = length(yEdges);
nX = length(xEdges);
for iy = 1:nY
    for ix = 1:nX
        map(iy,ix) = length(find(xBin==ix & yBin==iy));
    end
end
map = conv2(map,kernal,'same');
[xM,yM] = meshgrid(xEdgesSmooth,yEdgesSmooth);
map = interp2(xEdges(1:end),yEdges(1:end),map,xM,yM);
imagesc(map,'XData',xEdgesSmooth,'YData',yEdgesSmooth)
set(gca,'YDir','Normal')