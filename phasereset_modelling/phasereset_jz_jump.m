function [signal,phaseatstart,phaseatreset] = phasereset_jz_jump (time, wavefr, startpos, resetpos, resetwin)

%
% Function generates signal composed of a single sinusoid whose phase is
% being reset. The initial phase of the sinusoid is chosen randomly.
%
% Required inputs:
%  time  - vector (s)
%  wavefr - frequency of the sinusoid which is being reset
%  startpos - position of the reset start [sample of the time input]
%  resetpos - position of the final/actual reset [sample of the time input]
%  numones - number of 'ones' to go into mean of amplitude to approach one
%            at time of resetpos; The higher the number, the tighter the
%            range of values across trials that amplitude at resetpos will be
% Output:
%  signal - simulated EEG signal; vector: 1 by frames*epochs containing concatenated trials
%  phaseatstart - ongoing random phase at startpos
%
% Implemented by: Rafal Bogacz and Nick Yeung, September 2006
% Modified largely by: Johanna Zumer, September 2015


startpos=round(startpos);
resetpos=round(resetpos + rand*resetwin);

%generating the wave
signal = nan (size(time));
initphase = 2*pi*rand(1)-pi;

phaseatstart=wrapToPi(wavefr*2*pi*time(startpos) + initphase);
phaseatreset=acos(mean([cos(phaseatstart) 0.5]));  % aiming for phase that corresponds to number in second position

ampatstart=cos(phaseatstart);
ampatreset=cos(phaseatreset);

% fit line to start/reset points
mslope=(ampatreset-ampatstart)/(time(resetpos)-time(startpos));
bintercept=ampatreset-mslope*time(resetpos);
if mslope>0 % uprise of cosine means phase between 180->360deg
  phaseatreset=-phaseatreset;
end

% possible new frequencies that would bring phase to desired by resetpos
% newfreqposs=((phaseatreset+[-2*pi 0 2*pi 4*pi] -phaseatstart)/[(resetpos-startpos)/fsample])/(2*pi);
% newfreq = newfreqposs(newfreqposs> 1 & newfreqposs<12);
% if isempty(newfreq) % if none in range, then select closest to period of interval over which to shift
%   newfreq = newfreqposs(dsearchn(newfreqposs', fsample/(resetpos-startpos)));
% end
% phaseatstartdiff=(wavefr-newfreq)*2*pi*time(startpos);
% finalfreq=mean([newfreq wavefr]);
% phaseatresetdiff=(newfreq-finalfreq)*2*pi*time(resetpos)+phaseatstartdiff;

resetphasediff=wrapToPi(wavefr*2*pi*time(resetpos))-phaseatreset;

signal(1:startpos-1)=cos(wavefr*2*pi*time(1:startpos-1) + initphase);
signal(startpos:resetpos)=mslope*time(startpos:resetpos)+bintercept;
signal(resetpos:length(time))=cos(wavefr*2*pi*time(resetpos:length(time)) -resetphasediff);
% signal(startpos:resetpos)=cos(newfreq*2*pi*time(startpos:resetpos) + initphase + phaseatstartdiff) ;
% signal(resetpos:length(time))=cos(finalfreq*2*pi*time(resetpos:length(time)) + initphase + phaseatresetdiff);

