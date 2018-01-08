function [signal,phaseatstart,newfreq] = phasereset_jz (time, fsample, wavefr, minfr, maxfr, startpos, resetpos, numones)

%
% Function generates signal composed of a single sinusoid whose phase is
% being reset. The frequency of the sinusoid is chosen randomly from a
% specified range [minfr, maxfr]. The initial phase of the sinusoid is chosen randomly.
%
% Required inputs:
%  time  - 
%  minfr - minimum frequency of the sinusoid which is being reset
%  maxfr - maximum frequency of the sinusoid which is being reset
% Optional inputs:
%  position - position of the reset [in frames]; default: frames/2 => in the middle
%  tjitter - stdev of time jitter of the reset [in frames]; default: 0 => no jitter
% Output:
%  signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
% Implemented by: Rafal Bogacz and Nick Yeung, September 2006


startpos=round(startpos);
resetpos=round(resetpos);

%generating the wave
signal = nan (size(time));
if isempty(wavefr)
  wavefr = rand(1) * (maxfr-minfr) + minfr; % wave frequency
end
initphase = 2*pi*rand(1)-pi;

phaseatstart=wrapToPi(wavefr*2*pi*time(startpos) + initphase);
phaseatreset=acos(mean([cos(phaseatstart) ones(1,numones)]));  % aiming for phase of 0degrees 

% possible new frequencies that would bring phase to desired by resetpos
newfreqposs=((phaseatreset+[-2*pi 0 2*pi 4*pi] -phaseatstart)/[(resetpos-startpos)/fsample])/(2*pi);
newfreq = newfreqposs(newfreqposs> minfr & newfreqposs<maxfr);
if isempty(newfreq) % if none in range, then select closest to period of interval over which to shift
  newfreq = newfreqposs(dsearchn(newfreqposs', fsample/(resetpos-startpos)));
end

phaseatstartdiff=(wavefr-newfreq)*2*pi*time(startpos);

finalfreq=mean([newfreq wavefr]);

phaseatresetdiff=(newfreq-finalfreq)*2*pi*time(resetpos)+phaseatstartdiff;

signal(1:startpos-1)=cos(wavefr*2*pi*time(1:startpos-1) + initphase);
signal(startpos:resetpos)=cos(newfreq*2*pi*time(startpos:resetpos) + initphase + phaseatstartdiff) ;
signal(resetpos:length(time))=cos(finalfreq*2*pi*time(resetpos:length(time)) + initphase + phaseatresetdiff);

