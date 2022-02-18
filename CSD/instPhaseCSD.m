function [signal_phase signal_amp] = instPhaseCSD(signal)
% signal_phase = LFP.InstPhase(signal);
%
% Returns instantaneous phase using the Hilbert transform for 'signal'
%
% andrew bogaard 27 sept 2010

H = hilbert(signal);
signal_phase = angle(H);
signal_amp = abs(H);

end