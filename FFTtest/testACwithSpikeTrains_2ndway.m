%% This function generates a rhythmic spike train based on provided
% frequencies, taking into account bursting properties of a cell. The peak
% frequency of the FFT of the autocorrelation of the cell is then
% determined to compare it to the frequency the spike train is modeled on.
% December 2017

function testACwithSpikeTrains_2ndway

%%Input variables:

%Burst: Maximum number of spikes possible within a burst. 
        %e.g.: if Burst = 3, each cycle will randomly be assigned 1, 2 or 3
        %spikes. 

%Frequency: Vector with the list of frequencies used to generate spike
        %trains. The spiketrains with be modeled with interspike intervals
        %corresponding to these frequencies.

%Cyclenum: Number of 'cycles'. Assuming a cell is not bursting, i.e.,
        %number of spikes = Cyclenum. Currently set at 1500 cycles. 

%timediff: Interspike interval within a burst. Currently set at 20 ms. 

%startSp: Arbitrary start timestamp for the spiketrain. Set at 1500, but
        %does not affect the analysis. 

%%Initialize variables
Burst = 3;
Frequency = [9:0.005:14];
Cyclenum = 1500;
timediff = 0.020;
startSp = 1500; 

% Generate a random distribution to assign 1, 2 or 3 spikes per burst per
% cycle.
r2 = randi(Burst,Cyclenum,1);
r2(1) = 1; %% Assume the first spike in the spike train is not a burst. 

% Generate another random array to account for cases in which there are 2
% spikes in a burst. If the corresponding value of randBurst is 1, the
% 'second' spike occurs before the spike occuring at the timelag
% determined by the frequency. If randBurst is 2, it occurs after. 

randBurst = randi(2,Cyclenum,1);

t_bin = 0.005; %% Time resolution for the autocorrelation

for currf = 1:length(Frequency)
   
    %%Intitalize
    tSp = [];
    n = 1;
    
    %Determine timelag for the current input frequency
    timelag = (1/Frequency(currf));

    %%generateSpikeTrain
    
    for t = 1:Cyclenum
        
        switch r2(t)
            case 1 %Burst with 1 spike
                if n ==1 
                    tSp(n) = startSp;
                else
                    tSp(n) = startSp+(t*timelag);
                end
                n = n+1;
            case 2 %Burst with 2 spikes
                 if randBurst(t) == 1
                        tSp(n+1) = startSp+(t*timelag);
                        tSp(n) = tSp(n+1)-timediff;
                 elseif randBurst(t) == 2
                        tSp(n) = startSp+(t*timelag);
                        tSp(n+1) = tSp(n)+timediff;
                end 
                 n = n+2;
            case 3 %Burst with 3 spikes
                tSp(n+1) = startSp+(t*timelag);
                tSp(n) = tSp(n+1)-timediff;
                tSp(n+2) = tSp(n+1)+timediff;
                n = n+3;
        end

    end
    
    %%Now spikeTrain is generated. Find the FFT of the spiketrain. 

%      [cor, lag] = CrossCorr(tSp', tSp');
%      [cor, lag, smooth] = SmoothCor(cor, lag, t_bin); % solve for smoothed signal, and chop off center peak and negative half of autocorrelation
%      [S, f] = getChronuxSpectrum(cor, t_bin);
    Fs = 1000;
    p.Fs = Fs;
    p.tapers = [1,1];
    p.pad = 1;
    p.fPass = [0 125];
    [ACValues ACBins] = getTemporalACValues(tSp,Fs,500);
    if max(ACValues) > 0
        [S,f] = mtspectrumpb(ACValues-mean(ACValues),p);
    end
     Sall = S(f>6.2&f<14);
     fall = f(f>6.2&f<14);
     [~,freqP] = max(Sall);
     freqPeak(currf) = fall(freqP);
    
end

%% Plot the results across all frequencies
figure
scatter(Frequency,freqPeak,'ko','filled')
xlabel('Input Frequency')
ylabel('Frequency from FFT')
refline(1)
end

function [cor, lag, smooth] = SmoothCor(cor, lag, t_bin)
    
    std_smooth_kernel = .005;
    
    kernel = pdf('normal', -std_smooth_kernel*10/t_bin:std_smooth_kernel*10/t_bin, 0, std_smooth_kernel/t_bin); % convolve with gaussian

    if isempty(cor), smooth = []; return; % if no signals, then send it back
    
    else
        
        smooth = zeros(size(cor));
        
        for i = 1:size(cor,2)  

            smooth(:,i) = ndnanfilter(cor(:,i), kernel(:), []);
        
        end
        
        lag = lag(round(end/2):end,:);
        
        cor = cor(round(end/2):end,:);
        
        smooth = smooth(round(end/2):end, :);

    end 
    
end

function [S, f] = getChronuxSpectrum(cor, t_bin)

    params.tapers = [1.25, 1];
    params.fpass = [];
    params.Fs = 1/t_bin;
    params.pad = 6;  

    [S,f]=mtspectrumc(cor,params); % compute spectrum of cor

end