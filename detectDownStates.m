
expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
pathToSessionsAll = {'IZ12\IZ12_288um_200121_sess2','IZ12\IZ12_288um_200122_sess3','IZ12\IZ12_288um_200127_sess4','IZ12\IZ12_432um_200129_sess5','IZ12\IZ12_576um_200131_sess6',...
    'IZ13\IZ13_216um_200304_sess3','IZ13\IZ13_360um_200306_sess4','IZ13\IZ13_504um_200309_sess6','IZ16\IZ16_144um_200616_sess1','IZ18\IZ18_0um_200707_sess1',...
    'IZ24\IZ24_0um_200903_sess1','IZ24\IZ24_288um_200915_sess3','IZ26\IZ26_0um_201003_sess1','IZ26\IZ26_0um_201006_sess2','IZ26\IZ26_0um_201015_sess3',...
    'IZ27\IZ27_0um_201015_sess2','IZ33\IZ33_0um_210222_sess1','IZ34\IZ34_0um_210222_sess1'};
ripCh = [4 4 4 62 18 2 30 3 4 62 37 2 101 101 118 37 125 119];
ctxCh = {{19 54},{19 54},{19 54},{27 49},{8 56},{0 54},{31 57},{27 40},{30 49},{5 53},{4 55},{19 58},{4 37 104},{4 37 104},{0 50 99},{40},{63 104},{48 99}};

for i=1:size(pathToSessionsAll,2)
    tic
    disp(num2str(i))
    cd(strcat(expPath,pathToSessionsAll{i}))
    basePath = cd;
    baseName = bz_BasenameFromBasepath(basePath);
    load([baseName '.sessionInfo.mat'])
    
    % Identify channels
    CTXchans = [];
    for jj = 1:(size(sessionInfo.AnatGrps,2)-1)
        idx = find(sessionInfo.AnatGrps(jj).Channels==ctxCh{i}{jj});
        CTXchans = [CTXchans sessionInfo.AnatGrps(jj).Channels(1:idx)];
    end
    
%    badChannels = [sessionInfo.AnatGrps(size(sessionInfo.AnatGrps,2)).Channels];
%    SleepScoreMaster(basePath,'rejectChannels',badChannels); 
    
   %  [SlowWaves] = DetectSlowWaves(basePath,'CTXChans',CTXchans,'noPrompts',true);
     [Ripples] = rippleMasterDetectorIZ('rippleChannel',ripCh(i));
    %[SlowWaves] = DetectSlowWavesMUA(basePath,SWChan,'smoothwin',smoothwin,'startbins',startbins,'refineDipEstimate',refineDipEstimate);
    toc
end