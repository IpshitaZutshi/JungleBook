function sessRippleCortexDS

expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
pathToSessionsAll = {'IZ12\IZ12_288um_200121_sess2','IZ12\IZ12_288um_200122_sess3','IZ12\IZ12_288um_200127_sess4','IZ12\IZ12_432um_200129_sess5','IZ12\IZ12_576um_200131_sess6',...
    'IZ13\IZ13_216um_200304_sess3','IZ13\IZ13_360um_200306_sess4','IZ13\IZ13_504um_200309_sess6','IZ16\IZ16_144um_200616_sess1','IZ18\IZ18_0um_200707_sess1',...
    'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3'};%,'IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3',...
%     'IZ27\IZ27_0um_201009_sess1','IZ27\IZ27_0um_201015_sess2','IZ33\IZ33_0um_210222_sess1','IZ34\IZ34_0um_210222_sess1'};

win = [-0.75 0.75];
rast_x = []; rast_y = [];
lengthDS = [];
numDS = 1;
for ii=1:size(pathToSessionsAll,2)
    
    cd(strcat(expPath,pathToSessionsAll{ii}));
    file = dir('*.SlowWaves.events.mat');
    load(file.name);
    file = dir('*.ripples.events.mat');
    load(file.name);    
    
    st = SlowWaves.ints.DOWN(:,1);
    % extract all 
    
    for kk = 1:length(st)
        temp_rast = ripples.peaks - st(kk);
        temp_rast = temp_rast(temp_rast>win(1) & temp_rast<win(2));
        tempTS = zeros(1,2*(win(2)*1000)+1);
        if ~isempty(temp_rast)
            tempRastTS = round(temp_rast*1000)+(win(2)*1000)+1;
            tempTS(tempRastTS) = 1;
            rast_x = [rast_x; tempTS];
            rast_y = [rast_y; numDS*ones(1,length(tempTS))];
            lengthDS = [lengthDS (SlowWaves.ints.DOWN(kk,2)-SlowWaves.ints.DOWN(kk,1))];
            numDS = numDS+1;
        end
    end
    
    %[stccg, t] = CCG({ripples.peak st},[],'binSize',0.01,'duration',2);    
end
figure
[~,idx] =  sort(lengthDS);
plot(rast_x, rast_y,'.','MarkerSize',2)
end