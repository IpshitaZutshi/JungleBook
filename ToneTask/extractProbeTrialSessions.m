function Summary = extractProbeTrialSessions

sessloc = {'IZ47\Pre-implantation\IZ47_230513_161152','IZ47\Pre-implantation\IZ47_230515_183010','IZ47\Pre-implantation\IZ47_230516_113908',...
    'IZ48\Pre-implantation\IZ48_230504_125124','IZ48\Pre-implantation\IZ48_230508_190336','IZ48\Pre-implantation\IZ48_230516_190310',...
    'IZ47\Pre-implantation\IZ47_230511_180039','IZ47\Pre-implantation\IZ47_230514_151540','IZ47\Pre-implantation\IZ47_230515_160841',...
    'IZ48\Pre-implantation\IZ48_230502_110530','IZ48\Pre-implantation\IZ48_230503_113520','IZ48\Pre-implantation\IZ48_230505_170717',...
    'IZ47\Pre-implantation\IZ47_230510_190225','IZ47\Pre-implantation\IZ47_230511_121543','IZ47\Pre-implantation\IZ47_230512_182734',...
    'IZ48\Pre-implantation\IZ48_230504_192514','IZ48\Pre-implantation\IZ48_230505_101009','IZ48\Pre-implantation\IZ48_230508_141932'};

freq = [4 4 4 4 4 4 8 8 8 8 8 8 16 16 16 16 16 16];

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

for ii = 1:2
    Summary.PortDiffDist4{ii} = [];
    Summary.PortDiffDist8{ii} = [];
    Summary.PortDiffDist16{ii} = [];
end

for ii = 1:length(sessloc)
    cd(strcat(expPath,sessloc{ii}))
    [behavTrials] = getToneBehavior7ports;
    
    Summary.freq(ii) = freq(ii);
    
    %% Average fraction of no-probe trials (ports 3-6)
    Summary.fractProbe(ii) = sum(behavTrials.toneGain>1 &behavTrials.linTrial==0&behavTrials.probe==1)./sum(behavTrials.toneGain>1 &behavTrials.linTrial==0);
    
    %% Extract performance for probe and no probe trials
    Summary.perf(ii,1) = sum(behavTrials.toneGain>1 & behavTrials.linTrial==0 & behavTrials.probe==0 & behavTrials.correct==1)./...
        sum(behavTrials.toneGain>1 & behavTrials.linTrial==0 & behavTrials.probe==0);
    Summary.perf(ii,2) = sum(behavTrials.toneGain>1 & behavTrials.linTrial==0 & behavTrials.probe==1 & behavTrials.correct==1)./...
        sum(behavTrials.toneGain>1 & behavTrials.linTrial==0 & behavTrials.probe==1);
    
    %% Extract difference between target and lick
    idxNoProbe = behavTrials.toneGain>1 & behavTrials.linTrial==0 & behavTrials.probe==0;% & behavTrials.correct==0;
    idxProbe = behavTrials.toneGain>1 & behavTrials.linTrial==0 & behavTrials.probe==1;% & behavTrials.correct==0;
    if freq(ii) == 4
        Summary.PortDiffDist4{1} = [Summary.PortDiffDist4{1}; behavTrials.toneGain(idxNoProbe)-behavTrials.lickLoc(idxNoProbe)];
        Summary.PortDiffDist4{2} = [Summary.PortDiffDist4{2}; behavTrials.toneGain(idxProbe)-behavTrials.lickLoc(idxProbe)];
    elseif freq(ii) == 8
        Summary.PortDiffDist8{1} = [Summary.PortDiffDist8{1}; behavTrials.toneGain(idxNoProbe )-behavTrials.lickLoc(idxNoProbe)];
        Summary.PortDiffDist8{2} = [Summary.PortDiffDist8{2}; behavTrials.toneGain(idxProbe)-behavTrials.lickLoc(idxProbe)];        
    elseif freq(ii) == 16
        Summary.PortDiffDist16{1} = [Summary.PortDiffDist16{1}; behavTrials.toneGain(idxNoProbe )-behavTrials.lickLoc(idxNoProbe)];
        Summary.PortDiffDist16{2} = [Summary.PortDiffDist16{2}; behavTrials.toneGain(idxProbe)-behavTrials.lickLoc(idxProbe)];        
    end
end
    
end
