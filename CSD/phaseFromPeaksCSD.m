function phase = phaseFromPeaksCSD(signal,ts)
phase = nan(length(signal),1);
[~,iPk] = findpeaks(signal);
[~,iTrf] = findpeaks(-signal);
phase(iPk) = pi;
phase(iTrf) = 0;

if isequal(iPk(1),min(iPk(1),iTrf(1)))
    peak = true;
else
    peak = false;
end

iCyc = sort([iPk iTrf]);

for i = 1:length(iCyc)-1
   tCyc = ts(iCyc(i):iCyc(i+1));
   tCyc = tCyc-min(tCyc);
   tCyc = tCyc/max(tCyc);
   if peak
       phase(iCyc(i):iCyc(i+1)) = -(1-tCyc)*pi;
   else
       phase(iCyc(i):iCyc(i+1)) = tCyc*pi;
   end
   peak = ~peak;
end

