function [signal,phaseatstart,phaseatreset] = phasereset_jz_linphase (time, wavefr, startpos, resetpos, resetwin)

%
% Function generates signal composed of a single sinusoid whose phase is
% being reset. The initial phase of the sinusoid is chosen randomly.
%
% Required inputs:
%  time  - vector (s)
%  wavefr - frequency of the sinusoid which is being reset
%  startpos - position of the reset start [sample of the time input]
%  resetpos - position of the final/actual reset [sample of the time input]
%  resetwin - duration of flexible window after input resetpos that actual
%             resetpos will occur

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

phaseatstart=wrapTo2Pi(wavefr*2*pi*time(startpos) + initphase);
%   phaseatreset=acos(mean([cos(phaseatstart) 0.5]));  % aiming for phase that corresponds to number in second position
phaseatreset=2*pi;  % cos(0) = 1; thus, peak of curve

ampatstart=cos(phaseatstart);
ampatreset=cos(phaseatreset);


freqneed=[(1000/(2*pi))*(phaseatreset-phaseatstart-2*pi)/(resetpos-startpos)-wavefr (1000/(2*pi))*(phaseatreset-phaseatstart)/(resetpos-startpos)-wavefr (1000/(2*pi))*(phaseatreset-phaseatstart+2*pi)/(resetpos-startpos)-wavefr]

[mn,mnind]=min(abs(freqneed));
mnind
if phaseatstart<(2*pi) && phaseatstart>5.4
  mnind=3
end
if mnind==1 && (phaseatstart+2*pi)>(2*pi)
  mnind=2
end
[phaseatstart+2*pi phaseatstart  phaseatstart-2*pi]
if mnind==1
  phaseatstart=phaseatstart+2*pi;
elseif mnind==2
elseif mnind==3
  phaseatstart=phaseatstart-2*pi;
end
  


postwavefr=wavefr+0.5*randn(1);
resetphasediff=wrapToPi(postwavefr*2*pi*time(resetpos))-phaseatreset;

signal(1:startpos-1)=cos(wavefr*2*pi*time(1:startpos-1) + initphase);
signal(startpos:resetpos)=cos(linspace(phaseatstart,phaseatreset,numel(startpos:resetpos)));
signal(resetpos:length(time))=cos(postwavefr*2*pi*time(resetpos:length(time)) -resetphasediff);

% figure;plot(time(startpos-200:resetpos+200),signal(startpos-200:resetpos+200))

% signal(startpos:resetpos)=cos(newfreq*2*pi*time(startpos:resetpos) + initphase + phaseatstartdiff) ;
% signal(resetpos:length(time))=cos(finalfreq*2*pi*time(resetpos:length(time)) + initphase + phaseatresetdiff);

