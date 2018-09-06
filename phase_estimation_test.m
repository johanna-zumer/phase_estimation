% Testing / simulating examples for estimating phase
%
% Question: What is best between bandpass filter + Hilbert versus FFT/taper
% versus wavelet?
%
% Overview
% 1) Signal has one primary oscillatory component with mild noise added
% 2) Signal has 1, 2, 4, 8 Hz components; how does phase of one bleed on to other?
% 3) Signal has amplitude peak at either 4,5,6, or 7 Hz
% 4) Create new data, with ERP just after time point of interest, with additive model
% 5) Create new data, with ERP just after time point of interest, with phase-resetting model
% 6) Create new data, with Kc evoked at time point of interest, with additive model (nearly same as 4)
% 7) Add alpha ERD (and theta ERS?) on to ERP model
% 8) Signal has frequency randomly assigned between 8-12 Hz (constant though for all 100 trials); Use Hilbert to determine this frequency then FFT at that frequency
% 8a) like 4, ERP additive
% 8b) like 5, ERP phase-reset

clear all

% as of 31 aug, 2015: rerun 2-6 for fix of FIR cfg.fouse

run1=0;
run2=0;
run3=0;
run4=0;
run5=0;
run6=0;
run7=0;
run8=0;
run9=1;

ft_defaults;

fsample=1000;
time=0:1/fsample:8; % 8s long trial
wdr='D:\phase_estimation\';
addpath('D:\phase_estimation\phasereset_modelling')
wavwidth=[4 5 6 7 8];
wavgwidth=[3 4];

% FIXME / To-Do list:
% 1) Double check that correction/calculation for 1/2 period shift is computed correctly!
% 2) Discuss with Simon, Ali M, Ole, (also Rik Henson, Olivier David) and literature searches
% 3) Keep backround oscillation frequency constant but vary phase and noise
% added
% 4) Vary amplitude of ERP and Kc added (is it just amplitude effect or
% frequency effect, or ..?);  phase-reset model can be varied via tightness
% of window of desired phase
% 5) Report separately the bias and the variance (not just RMS)
% 6) plot Hanning, Hamming, Gaussian to see how they differ
% 7) complete simulation 7 (alpha ERD as well as ERP additive or ERP
% phase-reset)
% 8) Ensure filter order for FIR computed correctly (not 2*fsample/4)...set
% it explicity perhaps?  Set it equal to that of FIRWS?
% 9) See David's email!  (Huang-Hilbert and shaping bandpass of FFT)


%% 0) Test properties of FFT and Hilbert and Wavelet in terms of what angle '0' means

spectrum=hilbert(cos(2*pi*time));
angle(spectrum(2001)) % 0 degrees
spectrum=hilbert(sin(2*pi*time));
angle(spectrum(2001)) % -90 degrees
% Conclusion: use cos and then Hilbert gives angle=0deg

[spectrum,~,freqoi,timeoi]=ft_specest_mtmconvol(cos(2*pi*time),time,'timeoi',2,'taper','hanning','timwin',3*ones(1,3),'freqoi',1:3);
angle(spectrum(1)) % 0 degrees
[spectrum,~,freqoi,timeoi]=ft_specest_mtmconvol(sin(2*pi*time),time,'timeoi',2,'taper','hanning','timwin',3*ones(1,3),'freqoi',1:3);
angle(spectrum(1)) % -90 degrees
% Conclusion: use cos and then FFT gives angle=0deg

[spectrum,freqoi,timeoi]=ft_specest_wavelet(cos(2*pi*time),time,'timeoi',2,'freqoi',1:3,'width',3);
angle(spectrum(1)) % 0 degrees
[spectrum,freqoi,timeoi]=ft_specest_wavelet(sin(2*pi*time),time,'timeoi',2,'freqoi',1:3,'width',3);
angle(spectrum(1)) % -90 degrees
% Conclusion: use cos and then Wavelet gives angle=0deg

%% 0.1) How to report results? (added Jan 2018)

phasedeg(:,1)=20*randn(100,1); % centre 0, 20 vairance
phasedeg(:,2)=50*randn(100,1); % centre 0, 50 vairance
phasedeg(:,3)=20*randn(100,1)+20; % centre 20, 20 vairance
phasedeg(:,4)=50*randn(100,1)+20; % centre 20, 50 vairance
figure;subplot(2,2,1);rose(deg2rad(phasedeg(:,1)),100);
subplot(2,2,2);rose(deg2rad(phasedeg(:,2)),100);
subplot(2,2,3);rose(deg2rad(phasedeg(:,3)),100);
subplot(2,2,4);rose(deg2rad(phasedeg(:,4)),100);

rmsp=rms(phasedeg,1);
circvar=1-abs(sum(exp(i*deg2rad(phasedeg)),1)/100);
meanang=rad2deg(angle(sum(exp(i*deg2rad(phasedeg)))/100));


%% 1) Signal has one primary oscillatory component with mild noise added

if run1
  close all;
  clearvars -except run* time wdr fsample wav*width
  cd(wdr);
  
  halftimeind=round(length(time)/2);
  numtrials=150;
  raw.label{1}='test';
  raw.dimord='chan_time';
  rawpn=raw;
  for ss=1:numtrials
    cd('D:\fieldtrip_svn\utilities\private');
    state{ss}=randomseed(13+ss);
    cd(wdr);
    
    
    phaseshift(ss)=2*pi*rand(1)-pi; % random phase added, but fixed per frequency over all tests
    if ss<51
      frequse(ss)=3*rand(1)+1; % choose a single frequency between 1-4 Hz
    elseif ss<101
      frequse(ss)=3.5*rand(1)+4.5; % choose a single frequency between 4.5-8 Hz
    elseif ss<151
      frequse(ss)=3.5*rand(1)+8.5; % choose a single frequency between 8.5-12 Hz
    end
    trueang_raw(ss,:)=rad2deg(wrapToPi(frequse(ss)*2*pi*time+phaseshift(ss)*ones(size(time))));
    
    raw.time{ss}=time;
    rawpn.time{ss}=time;
    
    raw.trial{ss}=cos(frequse(ss)*2*pi*time+phaseshift(ss)*ones(size(time))); % add random phase
    rawpn.trial{ss}=raw.trial{ss}+10*pinknoise(length(time));
    
    
  end
  
  cfg=[];
  tlock=ft_timelockanalysis(cfg,raw);
  tlockpn=ft_timelockanalysis(cfg,rawpn);
  
  if 0
    cfg=[];
    ft_databrowser(cfg,raw);
    ft_databrowser(cfg,rawpn);
  end
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1a) Bandpass filter + Hilbert
  
  % test effect of time window length
  timwin=[4 2 1 .5 .25]; % duration in seconds; must be ordered from longest to shortest
  freqwin=[1 4; 4.5 8; 8.5 12];
  timeinduse=halftimeind-fsample*min(timwin)/2:halftimeind+fsample*min(timwin)/2-1;
  timeuse=time(timeinduse);
  for tt=1:length(timwin)
    
    cfg=[];
    cfg.latency=time([halftimeind-fsample*timwin(tt)/2 halftimeind+fsample*timwin(tt)/2]);
    rawuse=ft_selectdata(cfg,raw);
    rawpnuse=ft_selectdata(cfg,rawpn);
    
    for ff=1:size(freqwin,1)
      
      
      % Bandpass filter: Butterworth
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=freqwin(ff,:);
      cfg.bpfilttype='but';
      cfg.plotfiltresp='yes';
      cfg.fouse=2:4;
      cfg.figind=0;
      cfg.plotflag=0;
      raw_bpbutr=filter4phase_estim8(cfg,rawuse);
      rawpn_bpbutr=filter4phase_estim8(cfg,rawpnuse);
      cfg.hilbert='complex';
      raw_bpbut=filter4phase_estim8(cfg,rawuse);
      rawpn_bpbut=filter4phase_estim8(cfg,rawpnuse);
      
      % Bandpass filter: FIR (Matlab 'fir1')
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=freqwin(ff,:);
      cfg.bpfilttype='fir';
      cfg.plotfiltresp='yes';
      %       cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4]; % 3* is default
      cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
      cfg.figind=10;
      cfg.plotflag=0;
      raw_bpfirr=filter4phase_estim8(cfg,rawuse);
      rawpn_bpfirr=filter4phase_estim8(cfg,rawpnuse);
      cfg.hilbert='complex';
      raw_bpfir=filter4phase_estim8(cfg,rawuse);
      rawpn_bpfir=filter4phase_estim8(cfg,rawpnuse);
      
      % Bandpass filter: FIRWS (Matlab 'fir1')
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=freqwin(ff,:);
      cfg.bpfilttype='firws';
      cfg.plotfiltresp='yes';
      cfg.fouse=2; % this is actually cfg.bpfiltdf the transition width in Hz; 2Hz is effectively default chosen for low frequencies.
      cfg.figind=20;
      cfg.plotflag=0;
      raw_bpfirwsr=filter4phase_estim8(cfg,rawuse);
      rawpn_bpfirwsr=filter4phase_estim8(cfg,rawpnuse);
      cfg.hilbert='complex';
      raw_bpfirws=filter4phase_estim8(cfg,rawuse);
      rawpn_bpfirws=filter4phase_estim8(cfg,rawpnuse);
      if length(cfg.fouse)<3
        raw_bpfirws{3,2}=[];
        rawpn_bpfirws{3,2}=[];
      end
      
      if 0
        
        cfg=[];
        ft_databrowser(cfg,raw_bpbut{4,2});
        ft_databrowser(cfg,rawpn_bpbut{4,2});
        ft_databrowser(cfg,raw_bpfir{4,2});
        ft_databrowser(cfg,rawpn_bpfir{4,2});
      end
      
      if tt==1 % use same timephaseind for all tt, as chosen by largest time window during tt=1
        
        timephase(:,ff)=nan(numtrials,1);
        timephaseind(:,ff)=nan(numtrials,1);
        
        for ss=1:numtrials
          % assess difference across filter order setting as measure of variability
          try
            for fo=1:3
              angdiff(fo,1,:)=wrapToPi(angle(raw_bpbut{fo,2}.trial{ss})  -angle(raw_bpbut{fo+1,2}.trial{ss}));
              %               angdiff(fo,2,:)=wrapToPi(angle(raw_bpfir{fo,2}.trial{ss})  -angle(raw_bpfir{fo+1,2}.trial{ss}));
              %               angdiff(fo,3,:)=wrapToPi(angle(raw_bpfirws{fo,2}.trial{ss})-angle(raw_bpfirws{fo+1,2}.trial{ss}));
            end
          catch
            try
              for fo=1:2
                angdiff(fo,1,:)=wrapToPi(angle(raw_bpbut{fo,2}.trial{ss})-angle(raw_bpbut{fo+1,2}.trial{ss}));
                %               angdiff(fo,2,:)=wrapToPi(angle(raw_bpfir{fo,2}.trial{ss})-angle(raw_bpfir{fo+1,2}.trial{ss}));
                %               angdiff(fo,3,:)=wrapToPi(angle(raw_bpfirws{fo,2}.trial{ss})-angle(raw_bpfirws{fo+1,2}.trial{ss}));
              end
            catch
              for fo=1
                angdiff(fo,1,:)=wrapToPi(angle(raw_bpbut{fo,2}.trial{ss})-angle(raw_bpbut{fo+1,2}.trial{ss}));
                %               angdiff(fo,2,:)=wrapToPi(angle(raw_bpfir{fo,2}.trial{ss})-angle(raw_bpfir{fo+1,2}.trial{ss}));
                %               angdiff(fo,3,:)=wrapToPi(angle(raw_bpfirws{fo,2}.trial{ss})-angle(raw_bpfirws{fo+1,2}.trial{ss}));
              end
            end
          end
          if ss<26 || [ss>50 && ss<76] || [ss>100 && ss<126]  % first 50 is delta band; next 50 is theta band; last 50 is alpha band
            [mval,mind]=min(rms(reshape(angdiff(:,:,dsearchn(rawuse.time{1}',timeuse')),[size(angdiff,1)*size(angdiff,2) length(timeuse)]))); % avoid edge artifacts and so that later can have at least a 500ms window around it
          else
            [mval,mind]=max(rms(reshape(angdiff(:,:,dsearchn(rawuse.time{1}',timeuse')),[size(angdiff,1)*size(angdiff,2) length(timeuse)])));
          end
          timephaseindorig(ss,ff)=mind+timeinduse(1)-1;
          timephase(ss,ff)=time(timephaseindorig(ss,ff));
          trueang(ss,ff)=trueang_raw(ss,timephaseindorig(ss,ff));
        end %ss
        clear angdiff
        %       else
        %         timephaseind(:,ff)=dsearchn(rawuse.time{1}',timephase(:,ff));
      end
      timephaseind(:,ff)=dsearchn(rawuse.time{1}',timephase(:,ff));        % specific to tt/ff
      
      % reset these for every ff and tt
      angdiffbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angdifffir=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angdifffirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angdiffpnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angdiffpnfir=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angdiffpnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeepbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeepfir=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeepfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeeppnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeeppnfir=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeeppnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      for ss=1:numtrials
        for fo=1:size(raw_bpbut,1)
          for fd=1:size(raw_bpbut,2)
            if ~isempty(raw_bpbut{fo,fd})
              angkeepbut(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpbut{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
              angkeeppnbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpbut{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
            end
            if ~isempty(raw_bpfir{fo,fd})
              angkeepfir(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfir{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
              angkeeppnfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfir{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
            end
            if ~isempty(raw_bpfirws{fo,fd})
              angkeepfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfirws{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
              angkeeppnfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfirws{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
            end
          end
        end
        % don't need to index these angdiff* over ff as we don't save them
        angdiffbut(:,:,ss)=anglediff(angkeepbut(:,:,ss),trueang(ss,ff),1);
        angdifffir(:,:,ss)=anglediff(angkeepfir(:,:,ss),trueang(ss,ff),1);
        angdifffirws(:,:,ss)=anglediff(angkeepfirws(:,:,ss),trueang(ss,ff),1);
        angdiffpnbut(:,:,ss)=anglediff(angkeeppnbut(:,:,ss),trueang(ss,ff),1);
        angdiffpnfir(:,:,ss)=anglediff(angkeeppnfir(:,:,ss),trueang(ss,ff),1);
        angdiffpnfirws(:,:,ss)=anglediff(angkeeppnfirws(:,:,ss),trueang(ss,ff),1);
        %         angdiffbut(:,:,ss)=angkeepbut(:,:,ss)-trueang(ss,ff);
        %         angdifffir(:,:,ss)=angkeepfir(:,:,ss)-trueang(ss,ff);
        %         angdiffpnbut(:,:,ss)=angkeeppnbut(:,:,ss)-trueang(ss,ff);
        %         angdiffpnfir(:,:,ss)=angkeeppnfir(:,:,ss)-trueang(ss,ff);
      end
      
      % Only these summary RMS values keep/save out:
      %     angbutte(:,:,tt,ff)=rms(angdiffbut(:,:,1:25),3);
      %     angbutth(:,:,tt,ff)=rms(angdiffbut(:,:,26:50),3);
      %     angbutae(:,:,tt,ff)=rms(angdiffbut(:,:,51:75),3);
      %     angbutah(:,:,tt,ff)=rms(angdiffbut(:,:,76:100),3);
      angbutda(:,:,tt,ff)=rms(angdiffbut(:,:,1:50),3);
      angbutta(:,:,tt,ff)=rms(angdiffbut(:,:,51:100),3);
      angbutaa(:,:,tt,ff)=rms(angdiffbut(:,:,101:150),3);
      angbuta(:,:,tt,ff)=rms(angdiffbut(:,:,:),3);
      
      %     angfirte(:,:,tt,ff)=rms(angdifffir(:,:,1:25),3);
      %     angfirth(:,:,tt,ff)=rms(angdifffir(:,:,26:50),3);
      %     angfirae(:,:,tt,ff)=rms(angdifffir(:,:,51:75),3);
      %     angfirah(:,:,tt,ff)=rms(angdifffir(:,:,76:100),3);
      angfirda(:,:,tt,ff)=rms(angdifffir(:,:,1:50),3);
      angfirta(:,:,tt,ff)=rms(angdifffir(:,:,51:100),3);
      angfiraa(:,:,tt,ff)=rms(angdifffir(:,:,101:150),3);
      angfira(:,:,tt,ff)=rms(angdifffir(:,:,:),3);
      
      angfirwsda(:,:,tt,ff)=rms(angdifffirws(:,:,1:50),3);
      angfirwsta(:,:,tt,ff)=rms(angdifffirws(:,:,51:100),3);
      angfirwsaa(:,:,tt,ff)=rms(angdifffirws(:,:,101:150),3);
      angfirwsa(:,:,tt,ff)=rms(angdifffirws(:,:,:),3);
      %     angpnbutte(:,:,tt,ff)=rms(angdiffpnbut(:,:,1:25),3);
      %     angpnbutth(:,:,tt,ff)=rms(angdiffpnbut(:,:,26:50),3);
      %     angpnbutae(:,:,tt,ff)=rms(angdiffpnbut(:,:,51:75),3);
      %     angpnbutah(:,:,tt,ff)=rms(angdiffpnbut(:,:,76:100),3);
      angpnbutda(:,:,tt,ff)=rms(angdiffpnbut(:,:,1:50),3);
      angpnbutta(:,:,tt,ff)=rms(angdiffpnbut(:,:,51:100),3);
      angpnbutaa(:,:,tt,ff)=rms(angdiffpnbut(:,:,101:150),3);
      angpnbuta(:,:,tt,ff)=rms(angdiffpnbut(:,:,:),3);
      
      %     angpnfirte(:,:,tt,ff)=rms(angdiffpnfir(:,:,1:25),3);
      %     angpnfirth(:,:,tt,ff)=rms(angdiffpnfir(:,:,26:50),3);
      %     angpnfirae(:,:,tt,ff)=rms(angdiffpnfir(:,:,51:75),3);
      %     angpnfirah(:,:,tt,ff)=rms(angdiffpnfir(:,:,76:100),3);
      angpnfirda(:,:,tt,ff)=rms(angdiffpnfir(:,:,1:50),3);
      angpnfirta(:,:,tt,ff)=rms(angdiffpnfir(:,:,51:100),3);
      angpnfiraa(:,:,tt,ff)=rms(angdiffpnfir(:,:,101:150),3);
      angpnfira(:,:,tt,ff)=rms(angdiffpnfir(:,:,:),3);
      
      angpnfirwsda(:,:,tt,ff)=rms(angdiffpnfirws(:,:,1:50),3);
      angpnfirwsta(:,:,tt,ff)=rms(angdiffpnfirws(:,:,51:100),3);
      angpnfirwsaa(:,:,tt,ff)=rms(angdiffpnfirws(:,:,101:150),3);
      angpnfirwsa(:,:,tt,ff)=rms(angdiffpnfirws(:,:,:),3);
    end % ff
    
  end % tt
  
  save(['D:\phase_estimation\anglermsrun1.mat'],'ang*a','timephase','-append')
  
  figure(1);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);    bar(squeeze(angbutda(:,2,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Delta only trials');end
    subplot(4,length(timwin),tt+5);  bar(squeeze(angbutta(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Theta only trials');end
    subplot(4,length(timwin),tt+10); bar(squeeze(angbutaa(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Alpha only trials');end
    subplot(4,length(timwin),tt+15); bar(squeeze(angbuta(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('1-12Hz trials');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(1,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  figure(2);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);    bar(squeeze(angfirda(:,1,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Delta only trials');end
    subplot(4,length(timwin),tt+5);  bar(squeeze(angfirta(:,1,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Theta only trials');end
    subplot(4,length(timwin),tt+10); bar(squeeze(angfiraa(:,1,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Alpha only trials');end
    subplot(4,length(timwin),tt+15); bar(squeeze(angfira(:,1,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('1-12Hz trials');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(2,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  figure(3);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);    bar(squeeze(angfirwsda(:,1,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Delta only trials');end
    subplot(4,length(timwin),tt+5);  bar(squeeze(angfirwsta(:,1,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Theta only trials');end
    subplot(4,length(timwin),tt+10); bar(squeeze(angfirwsaa(:,1,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Alpha only trials');end
    subplot(4,length(timwin),tt+15); bar(squeeze(angfirwsa(:,1,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('1-12Hz trials');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(3,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  figure(4);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);    bar(squeeze(angpnbutda(:,2,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Delta only trials');end
    subplot(4,length(timwin),tt+5);  bar(squeeze(angpnbutta(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Theta only trials');end
    subplot(4,length(timwin),tt+10); bar(squeeze(angpnbutaa(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Alpha only trials');end
    subplot(4,length(timwin),tt+15); bar(squeeze(angpnbuta(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('1-12Hz trials');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(4,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  figure(5);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);    bar(squeeze(angpnfirda(:,2,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Delta only trials');end
    subplot(4,length(timwin),tt+5);  bar(squeeze(angpnfirta(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Theta only trials');end
    subplot(4,length(timwin),tt+10); bar(squeeze(angpnfiraa(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Alpha only trials');end
    subplot(4,length(timwin),tt+15); bar(squeeze(angpnfira(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('1-12Hz trials');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(5,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  figure(6);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);    bar(squeeze(angpnfirwsda(:,2,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Delta only trials');end
    subplot(4,length(timwin),tt+5);  bar(squeeze(angpnfirwsta(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Theta only trials');end
    subplot(4,length(timwin),tt+10); bar(squeeze(angpnfirwsaa(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Alpha only trials');end
    subplot(4,length(timwin),tt+15); bar(squeeze(angpnfirwsa(:,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('1-12Hz trials');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(6,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1b) FFT + taper
  
  cfg=[];
  cfg.method='mtmconvol';
  cfg.output='fourier';
  cfg.taper='hanning';
  cfg.foi=1:12;
  cfg.keeptrials='yes';
  
  
  timephase_fft(1:50)=timephase(1:50,1); %optimal for delta
  timephase_fft(51:100)=timephase(51:100,2); % optimal for theta
  timephase_fft(101:150)=timephase(101:150,3); % optimal for alpha
  trueang_fft(1:50)=trueang(1:50,1); %optimal for delta
  trueang_fft(51:100)=trueang(51:100,2); % optimal for theta
  trueang_fft(101:150)=trueang(101:150,3); % optimal for alpha
  
  timephaseind_fft=dsearchn(time',timephase_fft');
  
  %  timwin=[4 2 1 .5 .25]; % duration in seconds; must be ordered from longest to shortest
  t_ftimwin{1}=timwin(1)*ones(size(cfg.foi)); % full length of data
  t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 2 periods 4 Hz
  t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{6}=4./cfg.foi; %
  t_ftimwin{7}=3./cfg.foi; %
  t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
  t_ftimwin{9}=1./cfg.foi; % one period for each frequency
  t_ftimwin{10}=0.5./cfg.foi; %
  
  %   angfreq=nan(numtrials,length(cfg.foi),length(toi),length(pad),length(t_ftimwin));
  angfreq=nan(numtrials,length(cfg.foi),2,1,length(t_ftimwin));
  angfreqdiff=nan(size(angfreq));
  
  for ss=1:numtrials
    cfg.trials=ss;
    
    % first, centre on time of interest
    toi{1}=timephase_fft(ss);
    % second, centre half period before time of interest and add pi
    toi{2}=timephase_fft(ss)-0.5./cfg.foi;
    
    pad{1}=time(end)+1;
    %     pad{2}=2*time(end);
    %     pad{3}=2^13/1000; % 8.192
    
    
    for tt=1:length(toi)
      cfg.toi=toi{tt};
      for pp=1:length(pad)
        cfg.pad=pad{pp};
        for tf=1:length(t_ftimwin)
          cfg.t_ftimwin=t_ftimwin{tf};
          freq{tt,pp,tf}  =ft_freqanalysis(cfg, raw);
          freqpn{tt,pp,tf}=ft_freqanalysis(cfg, rawpn);
          if tt==1
            %           angfreq(:,:,tt,pp,tf)=angle(squeeze(freq{tt,pp,tf}.fourierspctrm))/(2*pi)*360;
            angfreq(ss,:,tt,pp,tf)=rad2deg(wrapToPi(angle(squeeze(freq{tt,pp,tf}.fourierspctrm))));
          elseif tt==2
            angfreq(ss,:,tt,pp,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freq{tt,pp,tf}.fourierspctrm(1,1,:,:)))+pi  )));
            %           for ss=1:numtrials
            %                angfreq(ss,:,tt,pp,tf)=diag(angle(squeeze(freq{tt,pp,tf}.fourierspctrm(ss,1,:,:)))/(2*pi)*360+180);
            %           end
          end
        end
      end
    end
  end
  
  for ss=1:numtrials
    angfreqdiff(ss,:,:,:,:)=anglediff(angfreq(ss,:,:,:,:),trueang_raw(ss,timephaseind_fft(ss)),1);
  end
  angfrms=squeeze(rms(angfreqdiff,1));
  angfrmsd=squeeze(rms(angfreqdiff(1:50,:,:,:,:),1));
  angfrmst=squeeze(rms(angfreqdiff(51:100,:,:,:,:),1));
  angfrmsa=squeeze(rms(angfreqdiff(101:150,:,:,:,:),1));
  
  %   figure(100);
  %   for tf=1:length(t_ftimwin),subplot(3,length(t_ftimwin),tf);bar(angfrms(:,:,1,tf));end
  %   for tf=1:length(t_ftimwin),subplot(3,length(t_ftimwin),tf+length(t_ftimwin));bar(angfrms(:,:,2,tf));end
  %   for tf=1:length(t_ftimwin),subplot(3,length(t_ftimwin),tf+2*length(t_ftimwin));bar(angfrms(:,:,3,tf));end
  
  % pad{1} is fine; over subplots for timwin length, then each has frequency and time-centre-point
  figure(100);
  for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(angfrms(:,:,tf));axis([-inf inf 0 160]);end
  figure(101);
  for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(angfrmsd(:,:,tf));axis([-inf inf 0 160]);end
  figure(102);
  for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(angfrmst(:,:,tf));axis([-inf inf 0 160]);end
  figure(103);
  for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(angfrmsa(:,:,tf));axis([-inf inf 0 160]);end
  
  
  save(['D:\phase_estimation\anglermsrun1.mat'],'ang*rms*','ang*a','timephase_fft','-append')
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 1c) Wavelet
  
  
  error('do this over diff gwidth; see run2');
  
  load(['D:\phase_estimation\anglermsrun1.mat'],'timephase_fft')
  
  foi=1:12;
  angwave=nan(numtrials,length(foi),2,length(wavwidth),3);
  
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  cfg.foi=foi;
  
  for ss=1:numtrials
    cfg.trials=ss;
    % first, centre on time of interest
    toi{1}=timephase_fft(ss);
    % second, centre half period before time of interest and add pi
    toi{2}=timephase_fft(ss)-0.5./cfg.foi;
    for tt=1:length(toi)
      cfg.toi=toi{tt};
      for ww=1:length(wavwidth)
        cfg.width=wavwidth(ww);
        wavefoi{tt,ww}=ft_freqanalysis(cfg,raw);
        if tt==1
          angwave(ss,:,tt,ww,1)= rad2deg(wrapToPi(angle(squeeze(wavefoi{tt,ww}.fourierspctrm))));
        elseif tt==2
          angwave(ss,:,tt,ww,1)= rad2deg(diag(wrapToPi(   angle(squeeze(wavefoi{tt,ww}.fourierspctrm(1,1,:,:)))+pi   )));
        end
      end
    end
  end
  
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  
  for ff=1:size(freqwin,1)
    
    cfg.foilim=freqwin(ff,:);
    
    for ss=1:numtrials
      cfg.trials=ss;
      % first, centre on time of interest
      toi{1}=timephase_fft(ss);
      % second, centre half period before time of interest and add pi
      toi{2}=timephase_fft(ss)-0.5./foi;
      
      for tt=1:length(toi)
        if tt==1
          cfg.toi=toi{tt};
        elseif tt==2
          cfg.toi=toi{tt}(floor(cfg.foilim(1)):ceil(cfg.foilim(2)));
        end
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoilim{tt,ww}=ft_freqanalysis(cfg,raw);
          % this produces .freq of length25
          
          % one option: take the middle freq within the range and use that angle.
          freqwave2=round(wavefoilim{tt,ww}.freq(ceil(length(wavefoilim{tt,ww}.freq)/2)));
          
          if tt==1
            %             angwave(ss,freqwave2,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(wavefoilim{tt,ww}.fourierspctrm(:,:,ceil(length(wavefoilim{tt,ww}.freq)/2))))));
            angwave(ss,freqwave2,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(wavefoilim{tt,ww}.fourierspctrm(:,:,dsearchn(wavefoilim{tt,ww}.freq',freqwave2)  )))));
          elseif tt==2
            angwave(ss,freqwave2,tt,ww,2)=rad2deg(wrapToPi( (angle(squeeze(wavefoilim{tt,ww}.fourierspctrm(1,1,dsearchn(wavefoilim{tt,ww}.freq',freqwave2) ,dsearchn([floor(cfg.foilim(1)):ceil(cfg.foilim(2))]',freqwave2)))))+pi ));
            %             freqwave2=round(wavefoilim{tt,ww}.freq(dsearchn(wavefoilim{tt,ww}.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'))); % rounded frequencies
            %             angwave(ss,min(freqwave2):max(freqwave2),tt,ww,2)=rad2deg(wrapToPi( diag(angle(squeeze(wavefoilim{tt,ww}.fourierspctrm(1,1,dsearchn(wavefoilim{tt,ww}.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi ));
            %             freqwave2=round(wavefoilim{tt,ww}.freq(dsearchn(wavefoilim{tt,ww}.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'))); % rounded frequencies
            %             angwave(ss,min(freqwave2):max(freqwave2),tt,ww,2)=rad2deg(wrapToPi( diag(angle(squeeze(wavefoilim{tt,ww}.fourierspctrm(1,1,ceil(length(wavefoilim{tt,ww}.freq)/2),:))))+pi ));
          end
          % another option: take average over
          if tt==1
            freqwave2=round(wavefoilim{tt,ww}.freq(ceil(length(wavefoilim{tt,ww}.freq)/2)));
            angwave(ss,freqwave2,tt,ww,3)=mean(rad2deg(wrapToPi(angle(squeeze(wavefoilim{tt,ww}.fourierspctrm)))));
          elseif tt==2
            freqwave2=round(wavefoilim{tt,ww}.freq(ceil(length(wavefoilim{tt,ww}.freq)/2)));
            angwave(ss,freqwave2,tt,ww,3)=mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavefoilim{tt,ww}.fourierspctrm(1,1,dsearchn(wavefoilim{tt,ww}.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
          end
          
        end % ww
      end % tt
    end  % ss
  end  % ff
  
  
  for ss=1:numtrials
    angwavediff(ss,:,:,:,:)=anglediff(angwave(ss,:,:,:,:),trueang_raw(ss,timephaseind_fft(ss)),1);
  end
  
  angwrms=squeeze(rms(angwavediff,1));
  angwrmsd=squeeze(rms(angwavediff(1:50,:,:,:,:),1));
  angwrmst=squeeze(rms(angwavediff(51:100,:,:,:,:),1));
  angwrmsa=squeeze(rms(angwavediff(101:150,:,:,:,:),1));
  
  figure(200);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrms(:,:,ww,1));axis([-inf inf 0 160]);end
  figure(201);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmsd(:,:,ww,1));axis([-inf inf 0 160]);end
  figure(202);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmst(:,:,ww,1));axis([-inf inf 0 160]);end
  figure(203);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmsa(:,:,ww,1));axis([-inf inf 0 160]);end
  
  figure(204);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrms(:,:,ww,2));axis([-inf inf 0 160]);end
  figure(205);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmsd(:,:,ww,2));axis([-inf inf 0 160]);end
  figure(206);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmst(:,:,ww,2));axis([-inf inf 0 160]);end
  figure(207);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmsa(:,:,ww,2));axis([-inf inf 0 160]);end
  
  figure(208);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrms(:,:,ww,3));axis([-inf inf 0 160]);end
  figure(209);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmsd(:,:,ww,3));axis([-inf inf 0 160]);end
  figure(210);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmst(:,:,ww,3));axis([-inf inf 0 160]);end
  figure(211);
  for ww=1:length(wavwidth),subplot(1,5,ww);bar(angwrmsa(:,:,ww,3));axis([-inf inf 0 160]);end
  
  
  
  
  
  
  
  
  
  
  
  save(['D:\phase_estimation\anglermsrun1.mat'],'ang*rms*','ang*a','-append')
  
  
end
%% 2) Signal has 1, 2, 4, 8 Hz components; how does phase of one bleed on to other?
% Each trial progressively has more of phase difference between frequency components.

if run2
  
  close all;
  clearvars -except run* time wdr fsample wav*width
  cd(wdr);
  
  halftimeind=round(length(time)/2);
  numtrials=90;
  raw.label{1}='test';
  raw.dimord='chan_time';
  rawpn=raw;
  frequse=[1 2 4 8];
  
  for ss=1:numtrials
    %     cd('D:\fieldtrip_svn\utilities\private');
    %     state{ss}=randomseed(13+ss);
    %     cd(wdr);
    
    % with each incremental ss, phase difference between frequencies is larger
    phaseshift(ss,1)=0;
    for pp=2:length(frequse)
      phaseshift(ss,pp)=wrapToPi(phaseshift(ss,1)+(ss/90)*[(pp-1)/length(frequse)]*2*pi);
    end
    
    raw.trial{ss}=0;
    for pp=1:length(frequse)
      trueang_raw(ss,pp,:)=frequse(pp)*2*pi*time+phaseshift(ss,pp);
      raw.trial{ss}=raw.trial{ss}+cos(squeeze(trueang_raw(ss,pp,:)))';
      esttrueang_raw(ss,:)=angle(hilbert(raw.trial{ss}));
    end
    
    raw.time{ss}=time;
    rawpn.time{ss}=time;
    
    rawpn.trial{ss}=raw.trial{ss}+10*pinknoise(length(time));
  end
  
  cfg=[];
  tlock=ft_timelockanalysis(cfg,raw);
  tlockpn=ft_timelockanalysis(cfg,rawpn);
  
  if 0
    cfg=[];
    ft_databrowser(cfg,raw);
    ft_databrowser(cfg,rawpn);
  end
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 2a) Bandpass filter + Hilbert
  
  
  % test effect of time window length
  timwin=[4 2 1 .5 .25]; % duration in seconds; must be ordered from longest to shortest
  freqwin=[1 4; 4.5 8; 8.5 12];
  timeinduse=halftimeind-fsample*min(timwin)/2:halftimeind+fsample*min(timwin)/2-1;
  timeuse=time(timeinduse);
  for tt=1:length(timwin)
    
    cfg=[];
    cfg.latency=time([halftimeind-fsample*timwin(tt)/2 halftimeind+fsample*timwin(tt)/2]);
    rawuse=ft_selectdata(cfg,raw);
    rawpnuse=ft_selectdata(cfg,rawpn);
    
    for ff=1:size(freqwin,1)
      
      
      % Bandpass filter: Butterworth
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=freqwin(ff,:);
      cfg.bpfilttype='but';
      cfg.plotfiltresp='yes';
      cfg.fouse=1:4;
      cfg.figind=0;
      cfg.plotflag=0;
      raw_bpbutr=filter4phase_estim8(cfg,rawuse);
      rawpn_bpbutr=filter4phase_estim8(cfg,rawpnuse);
      cfg.hilbert='complex';
      raw_bpbut=filter4phase_estim8(cfg,rawuse);
      rawpn_bpbut=filter4phase_estim8(cfg,rawpnuse);
      
      % Bandpass filter: FIR (Matlab 'fir1')
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=freqwin(ff,:);
      cfg.bpfilttype='fir';
      cfg.plotfiltresp='yes';
      %       cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4 5*fsample/4]; % 3* is default
      cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
      cfg.figind=10;
      cfg.plotflag=0;
      raw_bpfirr=filter4phase_estim8(cfg,rawuse);
      rawpn_bpfirr=filter4phase_estim8(cfg,rawpnuse);
      cfg.hilbert='complex';
      raw_bpfir=filter4phase_estim8(cfg,rawuse);
      rawpn_bpfir=filter4phase_estim8(cfg,rawpnuse);
      
      % Bandpass filter: FIRWS
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=freqwin(ff,:);
      cfg.bpfilttype='firws';
      cfg.plotfiltresp='yes';
      cfg.fouse=2;
      cfg.figind=20;
      cfg.plotflag=0;
      raw_bpfirwsr=filter4phase_estim8(cfg,rawuse);
      rawpn_bpfirwsr=filter4phase_estim8(cfg,rawpnuse);
      cfg.hilbert='complex';
      raw_bpfirws=filter4phase_estim8(cfg,rawuse);
      rawpn_bpfirws=filter4phase_estim8(cfg,rawpnuse);
      
      
      if 0
        
        cfg=[];
        ft_databrowser(cfg,raw_bpbut{4,2});
        ft_databrowser(cfg,rawpn_bpbut{4,2});
        ft_databrowser(cfg,raw_bpfir{4,2});
        ft_databrowser(cfg,rawpn_bpfir{4,2});
      end
      
      if tt==1 % use same timephaseind for all tt, as chosen by largest time window during tt=1
        
        timephase(:,ff)=nan(numtrials,1);
        timephaseind(:,ff)=nan(numtrials,1);
        
        for ss=1:numtrials
          % assess difference across filter order setting as measure of variability
          try
            for fo=1:3
              angdiff(fo,1,:)=wrapToPi(angle(raw_bpbut{fo,2}.trial{ss})-angle(raw_bpbut{fo+1,2}.trial{ss}));
              %               angdiff(fo,2,:)=wrapToPi(angle(raw_bpfir{fo,2}.trial{ss})-angle(raw_bpfir{fo+1,2}.trial{ss}));
              %               angdiff(fo,3,:)=wrapToPi(angle(raw_bpfirws{fo,2}.trial{ss})-angle(raw_bpfirws{fo+1,2}.trial{ss}));
            end
          catch
            for fo=1:2
              angdiff(fo,1,:)=wrapToPi(angle(raw_bpbut{fo,2}.trial{ss})-angle(raw_bpbut{fo+1,2}.trial{ss}));
              %               angdiff(fo,2,:)=wrapToPi(angle(raw_bpfir{fo,2}.trial{ss})-angle(raw_bpfir{fo+1,2}.trial{ss}));
              %               angdiff(fo,3,:)=wrapToPi(angle(raw_bpfirws{fo,2}.trial{ss})-angle(raw_bpfirws{fo+1,2}.trial{ss}));
            end
          end
          if ss<26 || [ss>50 && ss<76] || [ss>100 && ss<126]  % first 50 is delta band; next 50 is theta band; last 50 is alpha band
            [mval,mind]=min(rms(reshape(angdiff(:,:,dsearchn(rawuse.time{1}',timeuse')),[size(angdiff,1)*size(angdiff,2) length(timeuse)]))); % avoid edge artifacts and so that later can have at least a 500ms window around it
          else
            [mval,mind]=max(rms(reshape(angdiff(:,:,dsearchn(rawuse.time{1}',timeuse')),[size(angdiff,1)*size(angdiff,2) length(timeuse)])));
          end
          timephaseindorig(ss,ff)=mind+timeinduse(1)-1;
          timephase(ss,ff)=time(timephaseindorig(ss,ff));
          trueang(ss,:,ff)=trueang_raw(ss,:,timephaseindorig(ss,ff));
          esttrueang(ss,ff)=esttrueang_raw(ss,timephaseindorig(ss,ff));
        end %ss
        clear angdiff
        %       else
        %         timephaseind(:,ff)=dsearchn(rawuse.time{1}',timephase(:,ff));
      end
      timephaseind(:,ff)=dsearchn(rawuse.time{1}',timephase(:,ff));        % specific to tt/ff
      
      % reset these for every ff and tt
      angdiffbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),size(trueang_raw,2)+1,numtrials);
      angdifffir=nan(size(raw_bpbut,1),size(raw_bpbut,2),size(trueang_raw,2)+1,numtrials);
      angdifffirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),size(trueang_raw,2)+1,numtrials);
      angdiffpnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),size(trueang_raw,2)+1,numtrials);
      angdiffpnfir=nan(size(raw_bpbut,1),size(raw_bpbut,2),size(trueang_raw,2)+1,numtrials);
      angdiffpnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),size(trueang_raw,2)+1,numtrials);
      
      angkeepbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeepfir=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeepfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeeppnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeeppnfir=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      angkeeppnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
      %       trueang_raw(ss,pp,:)=frequse(pp)*2*pi*time+phaseshift(ss,pp);
      %       raw.trial{ss}=raw.trial{ss}+cos(squeeze(trueang_raw(ss,pp,:)))';
      %       esttrueang_raw(ss,:)=angle(hilbert(raw.trial{ss}));
      
      for ss=1:numtrials
        for fo=1:size(raw_bpbut,1)
          for fd=1:size(raw_bpbut,2)
            if ~isempty(raw_bpbut{fo,fd})
              angkeepbut(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpbut{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
              angkeeppnbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpbut{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
            end
            if ~isempty(raw_bpfir{fo,fd})
              angkeepfir(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfir{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
              angkeeppnfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfir{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
            end
            if ~isempty(raw_bpfirws{fo,fd})
              angkeepfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfirws{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
              angkeeppnfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfirws{fo,fd}.trial{ss}(timephaseind(ss,ff)))));
            end
          end
        end
        % don't need to index these angdiff* over ff as we don't save them
        for pp=1:size(trueang,2)
          angdiffbut(:,:,pp,ss)=anglediff(angkeepbut(:,:,ss),trueang(ss,pp,ff),1);
          angdifffir(:,:,pp,ss)=anglediff(angkeepfir(:,:,ss),trueang(ss,pp,ff),1);
          angdifffirws(:,:,pp,ss)=anglediff(angkeepfirws(:,:,ss),trueang(ss,pp,ff),1);
          angdiffpnbut(:,:,pp,ss)=anglediff(angkeeppnbut(:,:,ss),trueang(ss,pp,ff),1);
          angdiffpnfir(:,:,pp,ss)=anglediff(angkeeppnfir(:,:,ss),trueang(ss,pp,ff),1);
          angdiffpnfirws(:,:,pp,ss)=anglediff(angkeeppnfirws(:,:,ss),trueang(ss,pp,ff),1);
          angdiffbut(:,:,pp,ss)=anglediff(angkeepbut(:,:,ss),trueang(ss,pp,ff),1);
          angdifffir(:,:,pp,ss)=anglediff(angkeepfir(:,:,ss),trueang(ss,pp,ff),1);
          angdifffirws(:,:,pp,ss)=anglediff(angkeepfirws(:,:,ss),trueang(ss,pp,ff),1);
          angdiffpnbut(:,:,pp,ss)=anglediff(angkeeppnbut(:,:,ss),trueang(ss,pp,ff),1);
          angdiffpnfir(:,:,pp,ss)=anglediff(angkeeppnfir(:,:,ss),trueang(ss,pp,ff),1);
          angdiffpnfirws(:,:,pp,ss)=anglediff(angkeeppnfirws(:,:,ss),trueang(ss,pp,ff),1);
        end
        
        angdiffbut(:,:,length(frequse)+1,ss)=anglediff(angkeepbut(:,:,ss),esttrueang(ss,ff),1);
        angdifffir(:,:,length(frequse)+1,ss)=anglediff(angkeepfir(:,:,ss),esttrueang(ss,ff),1);
        angdifffirws(:,:,length(frequse)+1,ss)=anglediff(angkeepfirws(:,:,ss),esttrueang(ss,ff),1);
        angdiffpnbut(:,:,length(frequse)+1,ss)=anglediff(angkeeppnbut(:,:,ss),esttrueang(ss,ff),1);
        angdiffpnfir(:,:,length(frequse)+1,ss)=anglediff(angkeeppnfir(:,:,ss),esttrueang(ss,ff),1);
        angdiffpnfirws(:,:,length(frequse)+1,ss)=anglediff(angkeeppnfirws(:,:,ss),esttrueang(ss,ff),1);
        angdiffbut(:,:,length(frequse)+1,ss)=anglediff(angkeepbut(:,:,ss),esttrueang(ss,ff),1);
        angdifffir(:,:,length(frequse)+1,ss)=anglediff(angkeepfir(:,:,ss),esttrueang(ss,ff),1);
        angdifffirws(:,:,length(frequse)+1,ss)=anglediff(angkeepfirws(:,:,ss),esttrueang(ss,ff),1);
        angdiffpnbut(:,:,length(frequse)+1,ss)=anglediff(angkeeppnbut(:,:,ss),esttrueang(ss,ff),1);
        angdiffpnfir(:,:,length(frequse)+1,ss)=anglediff(angkeeppnfir(:,:,ss),esttrueang(ss,ff),1);
        angdiffpnfirws(:,:,length(frequse)+1,ss)=anglediff(angkeeppnfirws(:,:,ss),esttrueang(ss,ff),1);
      end
      
      % Only these summary RMS values keep/save out:
      %       angbuta(:,:,:,tt,ff)=rms(angdiffbut,4);
      %       angfira(:,:,:,tt,ff)=rms(angdifffir,4);
      %       angpnbuta(:,:,:,tt,ff)=rms(angdiffpnbut,4);
      %       angpnfira(:,:,:,tt,ff)=rms(angdiffpnfir,4);
      
      
      angbuta(:,:,:,tt,ff)=rms(angdiffbut(:,:,:,:),4);
      angbutae(:,:,:,tt,ff)=rms(angdiffbut(:,:,:,1:25),4);
      angbutam(:,:,:,tt,ff)=rms(angdiffbut(:,:,:,26:50),4);
      angbutal(:,:,:,tt,ff)=rms(angdiffbut(:,:,:,51:90),4);
      
      angfira(:,:,:,tt,ff)=rms(angdifffir(:,:,:,:),4);
      angfirae(:,:,:,tt,ff)=rms(angdifffir(:,:,:,1:25),4);
      angfiram(:,:,:,tt,ff)=rms(angdifffir(:,:,:,26:50),4);
      angfiral(:,:,:,tt,ff)=rms(angdifffir(:,:,:,51:90),4);
      
      angfirwsa(:,:,:,tt,ff)=rms(angdifffirws(:,:,:,:),4);
      angfirwsae(:,:,:,tt,ff)=rms(angdifffirws(:,:,:,1:25),4);
      angfirwsam(:,:,:,tt,ff)=rms(angdifffirws(:,:,:,26:50),4);
      angfirwsal(:,:,:,tt,ff)=rms(angdifffirws(:,:,:,51:90),4);
      
    end % ff
    
  end % tt
  
  save(['D:\phase_estimation\anglermsrun2.mat'],'angbut*','angfir*','timephase*','-append')
  
  figind=1;figure(figind);
  for tt=1:length(timwin),
    subplot(5,length(timwin),tt);    bar(squeeze(angbuta(:,2,1,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Angle 1Hz');end
    subplot(5,length(timwin),tt+5);  bar(squeeze(angbuta(:,2,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 2Hz');end
    subplot(5,length(timwin),tt+10); bar(squeeze(angbuta(:,2,3,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 4Hz');end
    subplot(5,length(timwin),tt+15); bar(squeeze(angbuta(:,2,4,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 8Hz');end
    subplot(5,length(timwin),tt+20); bar(squeeze(angbuta(:,2,5,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle Overall');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  figind=11;figure(figind);
  for tt=1:length(timwin),
    subplot(5,length(timwin),tt);    bar(squeeze(angbutae(:,2,1,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Angle 1Hz');end
    subplot(5,length(timwin),tt+5);  bar(squeeze(angbutae(:,2,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 2Hz');end
    subplot(5,length(timwin),tt+10); bar(squeeze(angbutae(:,2,3,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 4Hz');end
    subplot(5,length(timwin),tt+15); bar(squeeze(angbutae(:,2,4,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 8Hz');end
    subplot(5,length(timwin),tt+20); bar(squeeze(angbutae(:,2,5,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle Overall');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  figind=21;figure(figind);
  for tt=1:length(timwin),
    subplot(5,length(timwin),tt);    bar(squeeze(angbutam(:,2,1,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Angle 1Hz');end
    subplot(5,length(timwin),tt+5);  bar(squeeze(angbutam(:,2,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 2Hz');end
    subplot(5,length(timwin),tt+10); bar(squeeze(angbutam(:,2,3,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 4Hz');end
    subplot(5,length(timwin),tt+15); bar(squeeze(angbutam(:,2,4,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 8Hz');end
    subplot(5,length(timwin),tt+20); bar(squeeze(angbutam(:,2,5,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle Overall');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  figind=31;figure(figind);
  for tt=1:length(timwin),
    subplot(5,length(timwin),tt);    bar(squeeze(angbutal(:,2,1,tt,:)));axis([-inf inf 0 200])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('Angle 1Hz');end
    subplot(5,length(timwin),tt+5);  bar(squeeze(angbutal(:,2,2,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 2Hz');end
    subplot(5,length(timwin),tt+10); bar(squeeze(angbutal(:,2,3,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 4Hz');end
    subplot(5,length(timwin),tt+15); bar(squeeze(angbutal(:,2,4,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle 8Hz');end
    subplot(5,length(timwin),tt+20); bar(squeeze(angbutal(:,2,5,tt,:)));axis([-inf inf 0 200])
    if tt==1,ylabel('Angle Overall');end
  end
  legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  for fd=1:2
    figind=2+fd-1;figure(figind);
    for tt=1:length(timwin),
      subplot(5,length(timwin),tt);    bar(squeeze(angfira(:,fd,1,tt,:)));axis([-inf inf 0 200])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Angle 1Hz');end
      subplot(5,length(timwin),tt+5);  bar(squeeze(angfira(:,fd,2,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 2Hz');end
      subplot(5,length(timwin),tt+10); bar(squeeze(angfira(:,fd,3,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 4Hz');end
      subplot(5,length(timwin),tt+15); bar(squeeze(angfira(:,fd,4,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 8Hz');end
      subplot(5,length(timwin),tt+20); bar(squeeze(angfira(:,fd,5,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle Overall');end
    end
    legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
    set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    figind=12+fd-1;figure(figind);
    for tt=1:length(timwin),
      subplot(5,length(timwin),tt);    bar(squeeze(angfirae(:,fd,1,tt,:)));axis([-inf inf 0 200])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Angle 1Hz');end
      subplot(5,length(timwin),tt+5);  bar(squeeze(angfirae(:,fd,2,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 2Hz');end
      subplot(5,length(timwin),tt+10); bar(squeeze(angfirae(:,fd,3,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 4Hz');end
      subplot(5,length(timwin),tt+15); bar(squeeze(angfirae(:,fd,4,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 8Hz');end
      subplot(5,length(timwin),tt+20); bar(squeeze(angfirae(:,fd,5,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle Overall');end
    end
    legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
    set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    figind=22+fd-1;figure(figind);
    for tt=1:length(timwin),
      subplot(5,length(timwin),tt);    bar(squeeze(angfiram(:,fd,1,tt,:)));axis([-inf inf 0 200])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Angle 1Hz');end
      subplot(5,length(timwin),tt+5);  bar(squeeze(angfiram(:,fd,2,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 2Hz');end
      subplot(5,length(timwin),tt+10); bar(squeeze(angfiram(:,fd,3,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 4Hz');end
      subplot(5,length(timwin),tt+15); bar(squeeze(angfiram(:,fd,4,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 8Hz');end
      subplot(5,length(timwin),tt+20); bar(squeeze(angfiram(:,fd,5,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle Overall');end
    end
    legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
    set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    figind=32+fd-1;figure(figind);
    for tt=1:length(timwin),
      subplot(5,length(timwin),tt);    bar(squeeze(angfiral(:,fd,1,tt,:)));axis([-inf inf 0 200])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Angle 1Hz');end
      subplot(5,length(timwin),tt+5);  bar(squeeze(angfiral(:,fd,2,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 2Hz');end
      subplot(5,length(timwin),tt+10); bar(squeeze(angfiral(:,fd,3,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 4Hz');end
      subplot(5,length(timwin),tt+15); bar(squeeze(angfiral(:,fd,4,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 8Hz');end
      subplot(5,length(timwin),tt+20); bar(squeeze(angfiral(:,fd,5,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle Overall');end
    end
    legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
    set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    
    figind=2+fd+1;figure(figind);
    for tt=1:length(timwin),
      subplot(5,length(timwin),tt);    bar(squeeze(angfirwsa(:,fd,1,tt,:)));axis([-inf inf 0 200])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Angle 1Hz');end
      subplot(5,length(timwin),tt+5);  bar(squeeze(angfirwsa(:,fd,2,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 2Hz');end
      subplot(5,length(timwin),tt+10); bar(squeeze(angfirwsa(:,fd,3,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 4Hz');end
      subplot(5,length(timwin),tt+15); bar(squeeze(angfirwsa(:,fd,4,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 8Hz');end
      subplot(5,length(timwin),tt+20); bar(squeeze(angfirwsa(:,fd,5,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle Overall');end
    end
    legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
    set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    figind=12+fd+1;figure(figind);
    for tt=1:length(timwin),
      subplot(5,length(timwin),tt);    bar(squeeze(angfirwsae(:,fd,1,tt,:)));axis([-inf inf 0 200])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Angle 1Hz');end
      subplot(5,length(timwin),tt+5);  bar(squeeze(angfirwsae(:,fd,2,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 2Hz');end
      subplot(5,length(timwin),tt+10); bar(squeeze(angfirwsae(:,fd,3,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 4Hz');end
      subplot(5,length(timwin),tt+15); bar(squeeze(angfirwsae(:,fd,4,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 8Hz');end
      subplot(5,length(timwin),tt+20); bar(squeeze(angfirwsae(:,fd,5,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle Overall');end
    end
    legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
    set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    figind=22+fd+1;figure(figind);
    for tt=1:length(timwin),
      subplot(5,length(timwin),tt);    bar(squeeze(angfirwsam(:,fd,1,tt,:)));axis([-inf inf 0 200])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Angle 1Hz');end
      subplot(5,length(timwin),tt+5);  bar(squeeze(angfirwsam(:,fd,2,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 2Hz');end
      subplot(5,length(timwin),tt+10); bar(squeeze(angfirwsam(:,fd,3,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 4Hz');end
      subplot(5,length(timwin),tt+15); bar(squeeze(angfirwsam(:,fd,4,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 8Hz');end
      subplot(5,length(timwin),tt+20); bar(squeeze(angfirwsam(:,fd,5,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle Overall');end
    end
    legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
    set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    figind=32+fd+1;figure(figind);
    for tt=1:length(timwin),
      subplot(5,length(timwin),tt);    bar(squeeze(angfirwsal(:,fd,1,tt,:)));axis([-inf inf 0 200])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Angle 1Hz');end
      subplot(5,length(timwin),tt+5);  bar(squeeze(angfirwsal(:,fd,2,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 2Hz');end
      subplot(5,length(timwin),tt+10); bar(squeeze(angfirwsal(:,fd,3,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 4Hz');end
      subplot(5,length(timwin),tt+15); bar(squeeze(angfirwsal(:,fd,4,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle 8Hz');end
      subplot(5,length(timwin),tt+20); bar(squeeze(angfirwsal(:,fd,5,tt,:)));axis([-inf inf 0 200])
      if tt==1,ylabel('Angle Overall');end
    end
    legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
    set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
  end
  
  %   figind=3;figure(figind);
  %   for tt=1:length(timwin),
  %     subplot(5,length(timwin),tt);    bar(squeeze(angpnbuta(:,2,1,tt,:)));axis([-inf inf 0 200])
  %     title(['Time Window length ' num2str(timwin(tt))])
  %     if tt==1,ylabel('Angle 1Hz');end
  %     subplot(5,length(timwin),tt+5);  bar(squeeze(angpnbuta(:,2,2,tt,:)));axis([-inf inf 0 200])
  %     if tt==1,ylabel('Angle 2Hz');end
  %     subplot(5,length(timwin),tt+10); bar(squeeze(angpnbuta(:,2,3,tt,:)));axis([-inf inf 0 200])
  %     if tt==1,ylabel('Angle 4Hz');end
  %     subplot(5,length(timwin),tt+15); bar(squeeze(angpnbuta(:,2,4,tt,:)));axis([-inf inf 0 200])
  %     if tt==1,ylabel('Angle 8Hz');end
  %     subplot(5,length(timwin),tt+20); bar(squeeze(angpnbuta(:,2,5,tt,:)));axis([-inf inf 0 200])
  %     if tt==1,ylabel('Angle Overall');end
  %   end
  %   legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  %   set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  %
  %   figind=4;figure(figind);
  %   for tt=1:length(timwin),
  %     subplot(5,length(timwin),tt);    bar(squeeze(angpnfira(:,2,1,tt,:)));axis([-inf inf 0 200])
  %     title(['Time Window length ' num2str(timwin(tt))])
  %     if tt==1,ylabel('Angle 1Hz');end
  %     subplot(5,length(timwin),tt+5);  bar(squeeze(angpnfira(:,2,2,tt,:)));axis([-inf inf 0 200])
  %     if tt==1,ylabel('Angle 2Hz');end
  %     subplot(5,length(timwin),tt+10); bar(squeeze(angpnfira(:,2,3,tt,:)));axis([-inf inf 0 200])
  %     if tt==1,ylabel('Angle 4Hz');end
  %     subplot(5,length(timwin),tt+15); bar(squeeze(angpnfira(:,2,4,tt,:)));axis([-inf inf 0 200])
  %     if tt==1,ylabel('Angle 8Hz');end
  %     subplot(5,length(timwin),tt+20); bar(squeeze(angpnfira(:,2,5,tt,:)));axis([-inf inf 0 200])
  %     if tt==1,ylabel('Angle Overall');end
  %   end
  %   legend({'delta-filtered', 'theta-filtered', 'alpha-filtered'})
  %   set(get(figind,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 2b) FFT + taper
  
  
  cfg=[];
  cfg.method='mtmconvol';
  cfg.output='fourier';
  cfg.taper='hanning';
  cfg.keeptrials='yes';
  
  % timephase is 90x3 (numtrials x freqwin)
  % trueang is 90x4x3 (numtrials x frequse x freqwin)
  
  %   timephase_fft(1:50)=timephase(1:50,1); %optimal for delta
  %   timephase_fft(51:100)=timephase(51:100,2); % optimal for theta
  %   timephase_fft(101:150)=timephase(101:150,3); % optimal for alpha
  %   trueang_fft(1:50)=trueang(1:50,1); %optimal for delta
  %   trueang_fft(51:100)=trueang(51:100,2); % optimal for theta
  %   trueang_fft(101:150)=trueang(101:150,3); % optimal for alpha
  
  
  %   angfreq=nan(numtrials,length(cfg.foi),length(toi),length(pad),length(t_ftimwin));
  angfreq=nan(numtrials,12,2,1,10);
  angfreqdiff=nan(size(angfreq));
  
  for ff=1:size(freqwin,1)
    cfg.foi=round(freqwin(ff,1):freqwin(ff,2));
    
    %  timwin=[4 2 1 .5 .25]; % duration in seconds; must be ordered from longest to shortest
    t_ftimwin{1}=timwin(1)*ones(size(cfg.foi)); % full length of data
    t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 2 periods 4 Hz
    t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
    t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
    t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
    t_ftimwin{6}=4./cfg.foi; %
    t_ftimwin{7}=3./cfg.foi; %
    t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
    t_ftimwin{9}=1./cfg.foi; % one period for each frequency
    t_ftimwin{10}=0.5./cfg.foi; %
    
    for ss=1:numtrials
      cfg.trials=ss;
      
      timephaseind_fft(ss,ff)=dsearchn(time',timephase(ss,ff)');
      
      % first, centre on time of interest
      toi{1}=timephase(ss,ff);
      % second, centre half period before time of interest and add pi
      toi{2}=timephase(ss,ff)-0.5./cfg.foi;
      
      pad{1}=time(end)+1;
      %     pad{2}=2*time(end);
      %     pad{3}=2^13/1000; % 8.192
      
      
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for pp=1:length(pad)
          cfg.pad=pad{pp};
          for tf=1:length(t_ftimwin)
            cfg.t_ftimwin=t_ftimwin{tf};
            freq  =ft_freqanalysis(cfg, raw);
            freqpn=ft_freqanalysis(cfg, rawpn);
            if tt==1
              angfreq(ss,cfg.foi,tt,pp,tf)=rad2deg(wrapToPi(angle(squeeze(freq.fourierspctrm))));
            elseif tt==2
              angfreq(ss,cfg.foi,tt,pp,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freq.fourierspctrm(1,1,:,:)))+pi  )));
            end
          end
        end
      end
    end
  end
  
  for ff=1:size(freqwin,1)
    foi=round(freqwin(ff,1):freqwin(ff,2));
    for ss=1:numtrials
      for pf=1:length(frequse)
        angfreqdiff(ss,foi,:,:,:,pf)=anglediff(angfreq(ss,foi,:,:,:,:),trueang_raw(ss,pf,timephaseind_fft(ss,ff)),1);
      end
      angfreqdiff(ss,foi,:,:,:,pf+1)=anglediff(angfreq(ss,foi,:,:,:,:),esttrueang_raw(ss,timephaseind_fft(ss,ff)) ,1);
    end
  end
  angfrms=squeeze(rms(angfreqdiff,1));
  angfrmse=squeeze(rms(angfreqdiff(1:25,:,:,:,:,:),1));
  angfrmsm=squeeze(rms(angfreqdiff(26:50,:,:,:,:,:),1));
  angfrmsl=squeeze(rms(angfreqdiff(51:90,:,:,:,:,:),1));
  
  % angfrms is foi x toistart x timwin x angfreqcompare
  % assuming only 1 pad value, gets squeezed out
  % 0,1,2,3 is all, early, middle, late
  % 10,20,30,40,50 is versus 1hz, 2hz, 4hz, 8hz, or combined
  for pf=1:size(angfrms,4)
    figure(100+pf*10);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(angfrms(:,:,tf,pf));axis([-inf inf 0 160]);end
    figure(101+pf*10);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(angfrmse(:,:,tf,pf));axis([-inf inf 0 160]);end
    figure(102+pf*10);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(angfrmsm(:,:,tf,pf));axis([-inf inf 0 160]);end
    figure(103+pf*10);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(angfrmsl(:,:,tf,pf));axis([-inf inf 0 160]);end
  end
  
  
  save(['D:\phase_estimation\anglermsrun2.mat'],'ang*rms*','ang*a','angbut*','angfir*','timephase*','-append')
  
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 2c) Wavelet
  
  load(['D:\phase_estimation\anglermsrun2.mat'],'timephase*')
  
  foi=1:12;
  % numtrials x foi x toitime x width x 3 compute_options
  angwave=nan(numtrials,length(foi),2,length(wavwidth),length(wavgwidth),3);
  
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  
  
  for ff=1:size(freqwin,1)  % don't need to keep angwave as 'ff' indexed as that's wrapped into cfg.foi
    cfg.foi=round(freqwin(ff,1):freqwin(ff,2));
    for ss=1:numtrials
      cfg.trials=ss;
      % first, centre on time of interest
      toi{1}=timephase(ss,ff);
      % second, centre half period before time of interest and add pi
      toi{2}=timephase(ss,ff)-0.5./cfg.foi;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          for gw=1:length(wavgwidth)
            cfg.gwidth=wavgwidth(gw);
            wavefoi=ft_freqanalysis(cfg,raw);
            if tt==1
              angwave(ss,cfg.foi,tt,ww,gw,1)= rad2deg(wrapToPi(angle(squeeze(wavefoi.fourierspctrm))));
            elseif tt==2
              angwave(ss,cfg.foi,tt,ww,gw,1)= rad2deg(diag(wrapToPi(   angle(squeeze(wavefoi.fourierspctrm(1,1,:,:)))+pi   )));
            end
          end
        end
      end
    end
  end
  
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  
  for ff=1:size(freqwin,1)
    
    cfg.foilim=freqwin(ff,:);
    
    for ss=1:numtrials
      cfg.trials=ss;
      % first, centre on time of interest
      toi{1}=timephase(ss,ff);
      % second, centre half period before time of interest and add pi
      toi{2}=timephase(ss,ff)-0.5./foi;
      
      for tt=1:length(toi)
        if tt==1
          cfg.toi=toi{tt};
        elseif tt==2
          cfg.toi=toi{tt}(floor(cfg.foilim(1)):ceil(cfg.foilim(2)));
        end
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          for gw=1:length(wavgwidth)
            cfg.gwidth=wavgwidth(gw);
            wavefoilim=ft_freqanalysis(cfg,raw);
            % this produces .freq of length25
            
            % one option: take the middle freq within the range and use that angle.
            freqwave2=round(wavefoilim.freq(ceil(length(wavefoilim.freq)/2)));
            
            if tt==1
              angwave(ss,freqwave2,tt,ww,gw,2)=rad2deg(wrapToPi(angle(squeeze(wavefoilim.fourierspctrm(:,:,dsearchn(wavefoilim.freq',freqwave2)  )))));
            elseif tt==2
              angwave(ss,freqwave2,tt,ww,gw,2)=rad2deg(wrapToPi( (angle(squeeze(wavefoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq',freqwave2) ,dsearchn([floor(cfg.foilim(1)):ceil(cfg.foilim(2))]',freqwave2)))))+pi ));
            end
            % another option: take average over
            if tt==1
              freqwave2=round(wavefoilim.freq(ceil(length(wavefoilim.freq)/2)));
              angwave(ss,freqwave2,tt,ww,gw,3)=mean(rad2deg(wrapToPi(angle(squeeze(wavefoilim.fourierspctrm)))));
            elseif tt==2
              freqwave2=round(wavefoilim.freq(ceil(length(wavefoilim.freq)/2)));
              angwave(ss,freqwave2,tt,ww,gw,3)=mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavefoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            end
            
          end % gw
        end % ww
      end % tt
    end  % ss
  end  % ff
  
  for ff=1:size(freqwin,1)
    foi=round(freqwin(ff,1):freqwin(ff,2));
    for ss=1:numtrials
      for pf=1:length(frequse)
        angwavediff(ss,foi,:,:,:,:,pf)=anglediff(angwave(ss,foi,:,:,:,:),trueang_raw(ss,pf,timephaseind_fft(ss,ff)),1);
      end
      angwavediff(ss,foi,:,:,:,:,pf+1)=anglediff(angwave(ss,foi,:,:,:,:),esttrueang_raw(ss,timephaseind_fft(ss,ff)),1);
    end
  end
  
  angwrms=squeeze(rms(angwavediff,1));
  angwrmse=squeeze(rms(angwavediff(1:25,:,:,:,:,:,:),1));
  angwrmsm=squeeze(rms(angwavediff(26:50,:,:,:,:,:,:),1));
  angwrmsl=squeeze(rms(angwavediff(51:90,:,:,:,:,:,:),1));
  
  % angwrms is:  foi x toistart x wavwidth x wavgwidth x freqmethod x angfreqcompare
  % 0,1,2,3 is all, early, middle, late
  % 10,20,30,40,50 is versus 1hz, 2hz, 4hz, 8hz, or combined
  %   for pf=1:size(angwrms,6)
  %     figure(200+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrms(:,:,ww,gw,1,pf));axis([-inf inf 0 160]);end;end
  %     figure(201+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmse(:,:,ww,gw,1,pf));axis([-inf inf 0 160]);end;end
  %     figure(202+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmsm(:,:,ww,gw,1,pf));axis([-inf inf 0 160]);end;end
  %     figure(203+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmsl(:,:,ww,gw,1,pf));axis([-inf inf 0 160]);end;end
  %   end
  %
  %   for pf=1:size(angwrms,6)
  %     figure(204+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrms(:,:,ww,gw,2,pf));axis([-inf inf 0 160]);end;end
  %     figure(205+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmse(:,:,ww,gw,2,pf));axis([-inf inf 0 160]);end;end
  %     figure(206+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmsm(:,:,ww,gw,2,pf));axis([-inf inf 0 160]);end;end
  %     figure(207+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmsl(:,:,ww,gw,2,pf));axis([-inf inf 0 160]);end;end
  %   end
  %
  %   for pf=1:size(angwrms,6)
  %     figure(254+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrms(:,:,ww,gw,3,pf));axis([-inf inf 0 160]);end;end
  %     figure(255+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmse(:,:,ww,gw,3,pf));axis([-inf inf 0 160]);end;end
  %     figure(256+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmsm(:,:,ww,gw,3,pf));axis([-inf inf 0 160]);end;end
  %     figure(257+pf*10);
  %     for ww=1:length(wavwidth),for gw=1:length(wavgwidth),subplot(length(wavgwidth),length(wavwidth),(gw-1)*length(wavwidth)+ww);bar(angwrmsl(:,:,ww,gw,3,pf));axis([-inf inf 0 160]);end;end
  %   end
  
  save(['D:\phase_estimation\anglermsrun2.mat'],'ang*rms*','ang*a','timephase*','-append')
  
  % Can ignore wavwidth and wavgwidth: results are identical (see above commented out plots)
  % angwrms is:  foi x toistart x wavwidth x wavgwidth x freqmethod x angfreqcompare
  ww=4;
  gw=1;
  figure(400);
  for pf=1:size(angwrms,6),for fm=1:size(angwrms,5),subplot(size(angwrms,6),size(angwrms,5),(pf-1)*size(angwrms,5)+fm);bar(angwrms(:,:,ww,gw,fm,pf));axis([-inf inf 0 160]);end;end
  figure(401);
  for pf=1:size(angwrms,6),for fm=1:size(angwrms,5),subplot(size(angwrms,6),size(angwrms,5),(pf-1)*size(angwrms,5)+fm);bar(angwrmse(:,:,ww,gw,fm,pf));axis([-inf inf 0 160]);end;end
  figure(402);
  for pf=1:size(angwrms,6),for fm=1:size(angwrms,5),subplot(size(angwrms,6),size(angwrms,5),(pf-1)*size(angwrms,5)+fm);bar(angwrmsm(:,:,ww,gw,fm,pf));axis([-inf inf 0 160]);end;end
  figure(403);
  for pf=1:size(angwrms,6),for fm=1:size(angwrms,5),subplot(size(angwrms,6),size(angwrms,5),(pf-1)*size(angwrms,5)+fm);bar(angwrmsl(:,:,ww,gw,fm,pf));axis([-inf inf 0 160]);end;end
  
  
  
  
end

%% 3) Signal has amplitude peak at either 4,5,6, or 7 Hz

if run3
  % Create data
  
  
  % params for whole simulation
  clearvars -except run* time wdr fsample wav*width
  close all;
  cd(wdr);
  timeask=0;
  halftimeind=round(length(time)/2);
  % Effect of time window length
  timwin=[round(length(time)/fsample) 2 1 .5 .25]; % 4, 2, 1 period of 4Hz
  
  angrms=nan(4,2,3,length(timwin),20);     % FO, [But FIR], timwin (4s, 1s, 500ms, 250ms), 20 simulations
  angtrms=nan(4,2,3,length(timwin),20);    % FO, [But FIR], timwin (4s, 1s, 500ms, 250ms), 20 simulations
  
  angfrms=nan(4,2,3,10,20);  % freq, toi, pad, t_ftimwin, 20 simulations
  angftrms=nan(4,2,3,10,20);
  angwrms=nan(4,2,length(wavwidth),2,20);  % freq, toi, wavwidth, foi/foilim, 20 simulations
  angwtrms=nan(4,2,length(wavwidth),2,20);
  
  for ss=1:20 % 20 difference simulations: 10 where easy phase is picked, 10 where difficult phase is picked.
    cd('D:\fieldtrip_svn\utilities\private');
    state{ss}=randomseed(13+ss);
    cd(wdr);
    
    freqtest=[1:12];
    clear raw
    phaseshift=2*pi*rand(1,length(freqtest))-pi; % random phase added, but fixed per frequency over all tests
    for rr=1:[length(freqtest)+1]  % rr+1 has no dominant frequency (sum of all equally)
      raw{rr}.label{1}='test';
      raw{rr}.dimord='chan_time';
      rawpn{rr}=raw{rr};
      for tr=1:length(freqtest) % index of frequency used
        raw{rr}.time{tr}=time;
        rawpn{rr}.time{tr}=time;
        % Always use 'cos' to generate signals, as that's what Hilbert finds.
        raw{rr}.trial{tr}=cos(freqtest(tr)*2*pi*time+phaseshift(tr)*ones(size(time))); % add random phase
        if rr==tr  % triple amplitude
          raw{rr}.trial{tr}=2*raw{rr}.trial{tr}; % weight this one as double
        end
        rawpn{rr}.trial{tr}=raw{rr}.trial{tr}+5*pinknoise(length(time));
      end
    end
    for rr=1:[length(freqtest)+1]  % rr+1 has no dominant frequency (sum of all equally)
      cfg=[];
      tlock{rr}=ft_timelockanalysis(cfg,raw{rr});
      cfg.trials=4:7;
      tlock_theta{rr}=ft_timelockanalysis(cfg,raw{rr}); % only computed for purpose of getting 'true' phase
      tlockpn{rr}=ft_timelockanalysis(cfg,rawpn{rr});
      if 0
        figure;plot(tlock{rr}.time,tlock{rr}.avg);
        figure;plot(tlock{rr}.time,tlockpn{rr}.avg);
      end
    end
    
    % % computation of true phase
    % Theoretically
    for tr=1:length(freqtest)+1
      if tr<=length(freqtest)
        trueang_raw(tr,:)=rad2deg(wrapToPi(freqtest(tr)*2*pi*time+phaseshift(tr)*ones(size(time))));
      end
      trueang_tlocktheta(tr,:)=rad2deg(wrapToPi(ft_preproc_hilbert(tlock_theta{tr}.avg,'angle')));
      trueang_tlock(tr,:)=rad2deg(wrapToPi(ft_preproc_hilbert(tlock{tr}.avg,'angle')));
    end
    
    if 0
      figure;
      subplot(5,4,1);plot(time,raw{4}.trial{4});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,2);plot(time,raw{4}.trial{5});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,3);plot(time,raw{4}.trial{6});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,4);plot(time,raw{4}.trial{7});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,5);plot(time,raw{5}.trial{4});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,6);plot(time,raw{5}.trial{5});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,7);plot(time,raw{5}.trial{6});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,8);plot(time,raw{5}.trial{7});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,9);plot(time,raw{6}.trial{4});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,10);plot(time,raw{6}.trial{5});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,11);plot(time,raw{6}.trial{6});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,12);plot(time,raw{6}.trial{7});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,13);plot(time,raw{7}.trial{4});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,14);plot(time,raw{7}.trial{5});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,15);plot(time,raw{7}.trial{6});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,16);plot(time,raw{7}.trial{7});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,17);plot(time,raw{13}.trial{4});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,18);plot(time,raw{13}.trial{5});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,19);plot(time,raw{13}.trial{6});axis([-.5 4.5 -2.5 2.5])
      subplot(5,4,20);plot(time,raw{13}.trial{7});axis([-.5 4.5 -2.5 2.5])
      
      figure;
      subplot(5,1,1);plot(time,tlock_theta{4}.avg);axis([-.5 4.5 -2.5 2.5])
      subplot(5,1,2);plot(time,tlock_theta{5}.avg);axis([-.5 4.5 -2.5 2.5])
      subplot(5,1,3);plot(time,tlock_theta{6}.avg);axis([-.5 4.5 -2.5 2.5])
      subplot(5,1,4);plot(time,tlock_theta{7}.avg);axis([-.5 4.5 -2.5 2.5])
      subplot(5,1,5);plot(time,tlock_theta{13}.avg);axis([-.5 4.5 -2.5 2.5])
    end
    
    
    close all
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3a) Bandpass filter + Hilbert
    
    %
    rrtest=[4:7 13];
    
    
    for tf=1:length(timwin)
      
      if tf==1
        timephase(ss)=time(halftimeind); % this will get updated for tf=1 and be something else for subsequent tf
      end
      
      cfg=[];
      cfg.latency=[timephase(ss)-.5*timwin(tf) timephase(ss)+.5*timwin(tf)];  % 4-period of 4Hz time window centred at time of interest
      for rr=rrtest
        tlocktw{rr}=ft_selectdata(cfg,tlock{rr});
      end
      
      % Bandpass filter: Butterworth
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=[4 7];
      cfg.bpfilttype='but';
      cfg.plotfiltresp='yes';
      cfg.hilbert='complex';
      cfg.fouse=1:4;
      cfg.figind=0;
      cfg.plotflag=0;
      tlock_bpbut=filter4phase_estimate(cfg,tlocktw);
      
      % Bandpass filter: FIR (Matlab 'fir1')
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=[4 7];
      cfg.bpfilttype='fir';
      cfg.plotfiltresp='yes';
      cfg.hilbert='complex';
      %       cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4 5*fsample/4]; % 3* is default
      cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
      cfg.figind=10;
      cfg.plotflag=0;
      tlock_bpfir=filter4phase_estimate(cfg,tlocktw);
      
      % Bandpass filter: FIRWS
      cfg=[];
      cfg.bpfilter='yes';
      cfg.bpfreq=[4 7];
      cfg.bpfilttype='firws';
      cfg.plotfiltresp='yes';
      cfg.hilbert='complex';
      cfg.fouse=2;
      cfg.figind=20;
      cfg.plotflag=0;
      tlock_bpfirws=filter4phase_estimate(cfg,tlocktw);
      
      if tf==1
        % quantitatively assess variance of phase over different methods and use
        % to objectively choose time points where phase variance is most and least
        if timeask
          if ss<11
            timephase(ss)=input('Which time point to use as possibly EASY to determine phase?');
          else
            timephase(ss)=input('Which time point to use as possibly HARD to determine phase?');
          end
          timephaseind=dsearchn(time',timephase(ss));
        else
          for rr=rrtest
            for fo=1:3
              angdiff(fo,rr,:)=anglediff(rad2deg(wrapToPi(angle(tlock_bpfir{rr,fo,2}.avg))),rad2deg(wrapToPi(angle(tlock_bpfir{rr,fo+1,2}.avg))),1);
            end
          end
          if ss<11  % choose middle 2 seconds to avoid edge artifacts and so that later can have at least a 500ms window around it
            [mval,mind]=min(rms(reshape(angdiff(:,rrtest,halftimeind-1000:halftimeind+999),[15 2000])));
          else
            [mval,mind]=max(rms(reshape(angdiff(:,rrtest,halftimeind-1000:halftimeind+999),[15 2000])));
          end
          timephaseind=mind+halftimeind-1000;
          timephase(ss)=time(timephaseind);
        end
        clear angdiff
      else
        timephaseind=dsearchn(tlock_bpfir{4,1,1}.time',timephase(ss));
      end
      
      for fo=[1 2 3 4] % can't go to 5 with this datalength
        for fd=1:size(tlock_bpbut,3)
          for rr=rrtest
            angkeepbut(fo,fd,rr)=rad2deg(wrapToPi(angle(tlock_bpbut{rr,fo,fd}.avg(timephaseind))));
            angkeepfir(fo,fd,rr)=rad2deg(wrapToPi(angle(tlock_bpfir{rr,fo,fd}.avg(timephaseind))));
            try
              angkeepfirws(fo,fd,rr)=rad2deg(wrapToPi(angle(tlock_bpfirws{rr,fo,fd}.avg(timephaseind))));
            catch
              angkeepfirws(fo,fd,rr)=nan;
            end
          end
        end
      end
      
      if 0
        figure;imagesc(anglediff(squeeze(angkeepbut(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,dsearchn(time',timephase(ss)))',[4 1]),1));colorbar;title('filter order by dominant freq; diff from truth');caxis([-180 180])
        figure;imagesc(anglediff(squeeze(angkeepfir(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,dsearchn(time',timephase(ss)))',[4 1]),1));colorbar;title('filter order by dominant freq; diff from truth');caxis([-180 180])
      end
      
      % more relevant for frequency-specific FFT
      angrms(:,:,1,tf,ss)=rms(anglediff(squeeze(angkeepbut(:,:,4:7)),  permute(repmat(trueang_raw(4:7,dsearchn(time',timephase(ss))),[1 4 2]),[2 3 1]),1),3);
      angrms(:,:,2,tf,ss)=rms(anglediff(squeeze(angkeepfir(:,:,4:7)),  permute(repmat(trueang_raw(4:7,dsearchn(time',timephase(ss))),[1 4 2]),[2 3 1]),1),3);
      angrms(:,:,3,tf,ss)=rms(anglediff(squeeze(angkeepfirws(:,:,4:7)),permute(repmat(trueang_raw(4:7,dsearchn(time',timephase(ss))),[1 4 2]),[2 3 1]),1),3);
      
      % more relevant for time-domain filter & Hilbert
      angtrms(:,:,1,tf,ss)=rms(anglediff(squeeze(angkeepbut(:,:,rrtest)),  permute(repmat(trueang_tlocktheta(rrtest,dsearchn(time',timephase(ss))),[1 4 2]),[2 3 1]),1),3);
      angtrms(:,:,2,tf,ss)=rms(anglediff(squeeze(angkeepfir(:,:,rrtest)),  permute(repmat(trueang_tlocktheta(rrtest,dsearchn(time',timephase(ss))),[1 4 2]),[2 3 1]),1),3);
      angtrms(:,:,3,tf,ss)=rms(anglediff(squeeze(angkeepfirws(:,:,rrtest)),permute(repmat(trueang_tlocktheta(rrtest,dsearchn(time',timephase(ss))),[1 4 2]),[2 3 1]),1),3);
      
    end % tf
    
    
    
    
    
    %     cfg=[];
    %     cfg.latency=[timephase(ss)-.5*t_ftimwin(tf) timephase(ss)+.5*t_ftimwin(tf)];  % 4-period of 4Hz time window centred at time of interest
    %     for rr=1:length(tlock)
    %       tlocktw{rr,tf}=ft_selectdata(cfg,tlock{rr});
    %     end
    %
    %     % Bandpass filter: Butterworth
    %     cfg=[];
    %     cfg.bpfilter='yes';
    %     cfg.bpfreq=[4 7];
    %     cfg.bpfilttype='but';
    %     cfg.plotfiltresp='yes';
    %     cfg.hilbert='complex';
    %     cfg.fouse=1:4;
    %     cfg.figind=0;
    %     cfg.plotflag=0;
    %     tlocktw_bpbut=filter4phase_estimate(cfg,tlocktw);
    %
    %     % Bandpass filter: FIR (Matlab 'fir1')
    %     cfg=[];
    %     cfg.bpfilter='yes';
    %     cfg.bpfreq=[4 7];
    %     cfg.bpfilttype='fir';
    %     cfg.plotfiltresp='yes';
    %     cfg.hilbert='complex';
    %     cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4 5*fsample/4]; % 3* is default
    %       cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
    %     cfg.figind=10;
    %     cfg.plotflag=0;
    %     tlocktw4_bpfir=filter4phase_estimate(cfg,tlocktw4);
    %
    %     timephaseindtw4=dsearchn(tlocktw4_bpbut{4,fo,fd}.time',timephase(ss));
    %     for fo=[1 2 3 4] % can't go to 5 with this datalength
    %       for fd=1:size(tlocktw4_bpbut,3)
    %         for rr=rrtest
    %           angkeeptwbut4(fo,fd,rr)=angle(tlocktw4_bpbut{rr,fo,fd}.avg(timephaseindtw4))/(2*pi)*360;
    %           angkeeptwfir4(fo,fd,rr)=angle(tlocktw4_bpfir{rr,fo,fd}.avg(timephaseindtw4))/(2*pi)*360;
    %         end
    %       end
    %     end
    %     figure;imagesc(anglediff(squeeze(angkeeptwbut4(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1));colorbar;title('filter order by dominant freq; diff from truth');caxis([-180 180])
    %     figure;imagesc(anglediff(squeeze(angkeeptwfir4(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1));colorbar;title('filter order by dominant freq; diff from truth');caxis([-180 180])
    %
    %     % more relevant for frequency-specific FFT
    %     angrms(:,1,2,ss)=rms(anglediff(squeeze(angkeeptwbut4(:,2,4:7)),repmat(trueang_raw(4:7,timephaseind)',[4 1]),1),2);
    %     angrms(:,2,2,ss)=rms(anglediff(squeeze(angkeeptwfir4(:,2,4:7)),repmat(trueang_raw(4:7,timephaseind)',[4 1]),1),2);
    %
    %     % more relevant for time-domain filter & Hilbert
    %     angtrms(:,1,2,ss)=rms(anglediff(squeeze(angkeeptwbut4(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1),2);
    %     angtrms(:,2,2,ss)=rms(anglediff(squeeze(angkeeptwfir4(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1),2);
    %
    %
    %     % %
    %     cfg=[];
    %     cfg.latency=[timephase(ss)-.25 timephase(ss)+.25];   % 2-period of 4Hz time window centred at time of interest
    %     for rr=1:length(tlock)
    %       tlocktw2{rr}=ft_selectdata(cfg,tlock{rr});
    %     end
    %
    %     % Bandpass filter: Butterworth
    %     cfg=[];
    %     cfg.bpfilter='yes';
    %     cfg.bpfreq=[4 7];
    %     cfg.bpfilttype='but';
    %     cfg.plotfiltresp='yes';
    %     cfg.hilbert='complex';
    %     cfg.fouse=1:4;
    %     cfg.figind=0;
    %     cfg.plotflag=0;
    %     tlocktw2_bpbut=filter4phase_estimate(cfg,tlocktw2);
    %
    %     % Bandpass filter: FIR (Matlab 'fir1')
    %     cfg=[];
    %     cfg.bpfilter='yes';
    %     cfg.bpfreq=[4 7];
    %     cfg.bpfilttype='fir';
    %     cfg.plotfiltresp='yes';
    %     cfg.hilbert='complex';
    %     cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4 5*fsample/4]; % 3* is default
    %       cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
    %     cfg.figind=10;
    %     cfg.plotflag=0;
    %     tlocktw2_bpfir=filter4phase_estimate(cfg,tlocktw2);
    %
    %     timephaseindtw2=dsearchn(tlocktw2_bpbut{4,fo,fd}.time',timephase(ss));
    %     for fo=[1 2 3 4] % can't go to 5 with this datalength
    %       for fd=1:size(tlocktw2_bpbut,3)
    %         for rr=rrtest
    %           angkeeptwbut2(fo,fd,rr)=angle(tlocktw2_bpbut{rr,fo,fd}.avg(timephaseindtw2))/(2*pi)*360;
    %           angkeeptwfir2(fo,fd,rr)=angle(tlocktw2_bpfir{rr,fo,fd}.avg(timephaseindtw2))/(2*pi)*360;
    %         end
    %       end
    %     end
    %
    %     figure;imagesc(anglediff(squeeze(angkeeptwbut2(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1));colorbar;title('filter order by dominant freq; diff from truth');caxis([-180 180])
    %     figure;imagesc(anglediff(squeeze(angkeeptwfir2(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1));colorbar;title('filter order by dominant freq; diff from truth');caxis([-180 180])
    %
    %     % more relevant for frequency-specific FFT
    %     angrms(:,1,3,ss)=rms(anglediff(squeeze(angkeeptwbut2(:,2,4:7)),repmat(trueang_raw(4:7,timephaseind)',[4 1]),1),2);
    %     angrms(:,2,3,ss)=rms(anglediff(squeeze(angkeeptwfir2(:,2,4:7)),repmat(trueang_raw(4:7,timephaseind)',[4 1]),1),2);
    %
    %     % more relevant for time-domain filter & Hilbert
    %     angtrms(:,1,3,ss)=rms(anglediff(squeeze(angkeeptwbut2(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1),2);
    %     angtrms(:,2,3,ss)=rms(anglediff(squeeze(angkeeptwfir2(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1),2);
    %
    %
    %     % %
    %     cfg=[];
    %     cfg.latency=[timephase(ss)-.125 timephase(ss)+.125];  % 1-period of 4Hz time window centred at time of interest
    %     for rr=1:length(tlock)
    %       tlocktw1{rr}=ft_selectdata(cfg,tlock{rr});
    %     end
    %
    %     % Bandpass filter: Butterworth
    %     cfg=[];
    %     cfg.bpfilter='yes';
    %     cfg.bpfreq=[4 7];
    %     cfg.bpfilttype='but';
    %     cfg.plotfiltresp='yes';
    %     cfg.hilbert='complex';
    %     cfg.fouse=1:4;
    %     cfg.figind=0;
    %     cfg.plotflag=0;
    %     tlocktw1_bpbut=filter4phase_estimate(cfg,tlocktw1);
    %
    %     % Bandpass filter: FIR (Matlab 'fir1')
    %     cfg=[];
    %     cfg.bpfilter='yes';
    %     cfg.bpfreq=[4 7];
    %     cfg.bpfilttype='fir';
    %     cfg.plotfiltresp='yes';
    %     cfg.hilbert='complex';
    %     cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4 5*fsample/4]; % 3* is default
    %       cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
    %     cfg.figind=10;
    %     cfg.plotflag=0;
    %     tlocktw1_bpfir=filter4phase_estimate(cfg,tlocktw1);
    %
    %     timephaseindtw1=dsearchn(tlocktw1_bpbut{4,fo,fd}.time',timephase(ss));
    %     for fo=[1 2 3 4] % can't go to 5 with this datalength
    %       for fd=1:size(tlocktw1_bpbut,3)
    %         for rr=rrtest
    %           angkeeptwbut1(fo,fd,rr)=angle(tlocktw1_bpbut{rr,fo,fd}.avg(timephaseindtw1))/(2*pi)*360;
    %           angkeeptwfir1(fo,fd,rr)=angle(tlocktw1_bpfir{rr,fo,fd}.avg(timephaseindtw1))/(2*pi)*360;
    %         end
    %       end
    %     end
    %     figure;imagesc(anglediff(squeeze(angkeeptwbut1(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1));colorbar;title('filter order by dominant freq; diff from truth');caxis([-180 180])
    %     figure;imagesc(anglediff(squeeze(angkeeptwfir1(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1));colorbar;title('filter order by dominant freq; diff from truth');caxis([-180 180])
    %
    %     % more relevant for frequency-specific FFT
    %     angrms(:,1,4,ss)=rms(anglediff(squeeze(angkeeptwbut1(:,2,4:7)),repmat(trueang_raw(4:7,timephaseind)',[4 1]),1),2);
    %     angrms(:,2,4,ss)=rms(anglediff(squeeze(angkeeptwfir1(:,2,4:7)),repmat(trueang_raw(4:7,timephaseind)',[4 1]),1),2);
    %
    %     % more relevant for time-domain filter & Hilbert
    %     angtrms(:,1,4,ss)=rms(anglediff(squeeze(angkeeptwbut1(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1),2);
    %     angtrms(:,2,4,ss)=rms(anglediff(squeeze(angkeeptwfir1(:,2,rrtest)),repmat(trueang_tlocktheta(rrtest,timephaseind)',[4 1]),1),2);
    
    
    
    
    if 0
      figure(1);for ee=1:size(angtrms,3),subplot(1,size(angtrms,3),ee);bar(angtrms(:,:,ee,ss));axis([0 5 0 180]);
        if ee==1,ylabel('RMS degrees'),  end;
        if ee==2,title('tlocktheta:  For 4s, 1s, 500ms, 250ms sized windows');
          legend({'Butterworth' 'FIR1' 'FIRWS'});xlabel('Filter order'); end; end
      
      figure(2);for ee=1:size(angtrms,3),subplot(1,size(angtrms,3),ee);bar(angrms(:,:,ee,ss));axis([0 5 0 180]);
        if ee==1,ylabel('RMS degrees'),  end;
        if ee==2,title('tlock:  For 4s, 1s, 500ms, 250ms sized windows');
          legend({'Butterworth' 'FIR1' 'FIRWS'});xlabel('Filter order'); end; end
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3b) FFT + taper
    
    
    cfg=[];
    cfg.method='mtmconvol';
    cfg.output='fourier';
    cfg.taper='hanning';
    cfg.foi=1:13;
    
    % first, centre on time of interest
    toi{1}=timephase(ss);
    % second, centre half period before time of interest and add pi
    toi{2}=timephase(ss)-0.5./cfg.foi;
    
    pad{1}=time(end)+1;
    pad{2}=2*time(end);
    pad{3}=2^13/1000; % 8.192
    
    t_ftimwin{1}=timwin(1)/2*ones(size(cfg.foi)); % to match Hilbert calculations  % 16 periods 4 Hz
    t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 8 periods 4 Hz
    t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 4 periods 4 Hz
    t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 2 periods 4 Hz
    t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
    t_ftimwin{6}=4./cfg.foi; %
    t_ftimwin{7}=3./cfg.foi; %
    t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
    t_ftimwin{9}=1./cfg.foi; % one period for each frequency
    t_ftimwin{10}=0.5./cfg.foi; %
    
    angfreq=nan(4,max(rrtest),length(toi),length(pad),length(t_ftimwin));
    for rr=rrtest % over differnt datasets with each of these frequencies as more dominant
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for pp=1:length(pad)
          cfg.pad=pad{pp};
          for tf=1:length(t_ftimwin)
            cfg.t_ftimwin=t_ftimwin{tf};
            freq{rr,tt,pp,tf}=ft_freqanalysis(cfg,tlock{rr});
            if tt==1
              angfreq(:,rr,tt,pp,tf)=rad2deg(wrapToPi(angle(squeeze(freq{rr,tt,pp,tf}.fourierspctrm(1,1,4:7,1)))));
            elseif tt==2
              angfreq(:,rr,tt,pp,tf)=rad2deg(wrapToPi(   diag(angle(squeeze(freq{rr,tt,pp,tf}.fourierspctrm(1,1,4:7,4:7)))+pi)  ));
            end
          end
        end
      end
    end
    
    % more relevant for frequency-specific FFT
    angfrms(:,:,:,:,ss)=squeeze(rms(anglediff(squeeze(angfreq(:,4:7,:,:,:)),repmat(trueang_raw(4:7,dsearchn(time',timephase(ss)))',[4 1 2 3 10]),1),2));
    
    % more relevant for time-domain filter & Hilbert
    angftrms(:,:,:,:,ss)=squeeze(rms(anglediff(squeeze(angfreq(:,rrtest,:,:,:)),repmat(trueang_tlocktheta(rrtest,dsearchn(time',timephase(ss)))',[4 1 2 3 10]),1),2));
    
    if 0
      figure(100);for tf=1:10,subplot(2,5,tf);bar(angfrms(:,:,3,tf,ss));axis([0 5 0 180]);
        set(get(100,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
        if tf==1,ylabel('RMS degrees'),  end;
        if tf==2,title('raw:  For 8 diff sized windows');
          legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
      
      figure(101);for tf=1:10,subplot(2,5,tf);bar(angftrms(:,:,3,tf,ss));axis([0 5 0 180]);
        set(get(101,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
        if tf==1,ylabel('RMS degrees'),  end;
        if tf==2,title('tlocktheta:  For 8 diff sized windows');
          legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
    end
    
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3c) Wavelet
    
    cfg=[];
    cfg.method='wavelet';
    cfg.output='fourier';
    cfg.foi=1:13;
    
    % first, centre on time of interest
    toi{1}=timephase(ss);
    % second, centre half period before time of interest and add pi
    toi{2}=timephase(ss)-0.5./cfg.foi;
    
    for rr=rrtest % over differnt datasets with each of these frequencies as more dominant
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoi{rr,tt,ww}=ft_freqanalysis(cfg,tlock{rr});
          if tt==1
            angwave(:,rr,tt,ww,1)=rad2deg(wrapToPi(angle(squeeze(wavefoi{rr,tt,ww}.fourierspctrm(1,1,4:7)))));
          elseif tt==2
            angwave(:,rr,tt,ww,1)=rad2deg(wrapToPi(diag(angle(squeeze(wavefoi{rr,tt,ww}.fourierspctrm(1,1,4:7,4:7))))+pi));
          end
        end
      end
    end
    
    cfg=[];
    cfg.method='wavelet';
    cfg.output='fourier';
    cfg.foilim=[4 7];
    
    for rr=rrtest % over differnt datasets with each of these frequencies as more dominant
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoilim{rr,tt,ww}=ft_freqanalysis(cfg,tlock{rr});
          if tt==1
            angwave(:,rr,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(wavefoilim{rr,tt,ww}.fourierspctrm(1,1,4:7)))));
          elseif tt==2
            angwave(:,rr,tt,ww,2)=rad2deg(wrapToPi(diag(angle(squeeze(wavefoilim{rr,tt,ww}.fourierspctrm(1,1,4:7,4:7))))+pi  ));
          end
        end
      end
    end
    
    % more relevant for frequency-specific FFT
    angwrms(:,:,:,:,ss)=squeeze(rms(anglediff(squeeze(angwave(:,4:7,:,:,:)),repmat(trueang_raw(4:7,dsearchn(time',timephase(ss)))',[4 1 2 length(wavwidth) 2]),1),2));
    
    % more relevant for time-domain filter & Hilbert
    angwtrms(:,:,:,:,ss)=squeeze(rms(anglediff(squeeze(angwave(:,rrtest,:,:,:)),repmat(trueang_tlocktheta(rrtest,dsearchn(time',timephase(ss)))',[4 1 2 length(wavwidth) 2]),1),2));
    
    if 0
      figure(200);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwrms(:,:,ww,1,ss));axis([0 5 0 180]);
        set(get(200,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
        if tf==1,ylabel('RMS degrees'),  end;
        if tf==2,title('raw:  For 8 diff sized windows');
          legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
      
      figure(201);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwtrms(:,:,ww,1,ss));axis([0 5 0 180]);
        set(get(201,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
        if tf==1,ylabel('RMS degrees'),  end;
        if tf==2,title('tlocktheta:  For 8 diff sized windows');
          legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
      
      figure(202);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwrms(:,:,ww,2,ss));axis([0 5 0 180]);
        set(get(202,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
        if tf==1,ylabel('RMS degrees'),  end;
        if tf==2,title('raw:  For 8 diff sized windows');
          legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
      
      figure(203);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwtrms(:,:,ww,2,ss));axis([0 5 0 180]);
        set(get(203,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
        if tf==1,ylabel('RMS degrees'),  end;
        if tf==2,title('tlocktheta:  For 8 diff sized windows');
          legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
    end
    
    
  end % ss
  % save timephase, ang*rms, state
  
  save(['D:\phase_estimation\anglermsrun3.mat'],'ang*rms')
  
  angrmse=mean(angrms(:,:,:,:,1:10),5);
  angrmsh=mean(angrms(:,:,:,:,11:20),5);
  angrmses=std(angrms(:,:,:,:,1:10),[],5);
  angrmshs=std(angrms(:,:,:,:,11:20),[],5);
  angrmsa=mean(angrms(:,:,:,:,1:20),5);
  angrmsas=std(angrms(:,:,:,:,1:20),[],5);
  
  angtrmse=mean(angtrms(:,:,:,:,1:10),5);
  angtrmsh=mean(angtrms(:,:,:,:,11:20),5);
  angtrmses=std(angtrms(:,:,:,:,1:10),[],5);
  angtrmshs=std(angtrms(:,:,:,:,11:20),[],5);
  angtrmsa=mean(angtrms(:,:,:,:,1:20),5);
  angtrmsas=std(angtrms(:,:,:,:,1:20),[],5);
  
  angfrmse=mean(angfrms(:,:,:,:,1:10),5);
  angfrmsh=mean(angfrms(:,:,:,:,11:20),5);
  angfrmses=std(angfrms(:,:,:,:,1:10),[],5);
  angfrmshs=std(angfrms(:,:,:,:,11:20),[],5);
  angfrmsa=mean(angfrms(:,:,:,:,1:20),5);
  angfrmsas=std(angfrms(:,:,:,:,1:20),[],5);
  
  angftrmse=mean(angftrms(:,:,:,:,1:10),5);
  angftrmsh=mean(angftrms(:,:,:,:,11:20),5);
  angftrmses=std(angftrms(:,:,:,:,1:10),[],5);
  angftrmshs=std(angftrms(:,:,:,:,11:20),[],5);
  angftrmsa=mean(angftrms(:,:,:,:,1:20),5);
  angftrmsas=std(angftrms(:,:,:,:,1:20),[],5);
  
  angwrmse=mean(angwrms(:,:,:,:,1:10),5);
  angwrmsh=mean(angwrms(:,:,:,:,11:20),5);
  angwrmses=std(angwrms(:,:,:,:,1:10),[],5);
  angwrmshs=std(angwrms(:,:,:,:,11:20),[],5);
  angwrmsa=mean(angwrms(:,:,:,:,1:20),5);
  angwrmsas=std(angwrms(:,:,:,:,1:20),[],5);
  
  angwtrmse=mean(angwtrms(:,:,:,:,1:10),5);
  angwtrmsh=mean(angwtrms(:,:,:,:,11:20),5);
  angwtrmses=std(angwtrms(:,:,:,:,1:10),[],5);
  angwtrmshs=std(angwtrms(:,:,:,:,11:20),[],5);
  angwtrmsa=mean(angwtrms(:,:,:,:,1:20),5);
  angwtrmsas=std(angwtrms(:,:,:,:,1:20),[],5);
  
  
  for fd=1:2
    figure(1+10*fd);for ee=1:size(angrmse,4),subplot(1,size(angrmse,4),ee);bar(squeeze(angrmse(:,fd,:,ee)));axis([0 5 0 180]);
      if ee==1,ylabel('RMS degrees'),  end;
      if ee==2,title('raw:  For 8s, 2s, 1s, 500ms, 250ms sized windows');
        legend({'Butterworth' 'FIR1' 'FIRWS'});xlabel('Filter order'); end; end
    
    figure(2+10*fd);for ee=1:size(angrmse,4),subplot(1,size(angrmse,4),ee);bar(squeeze(angrmsh(:,fd,:,ee)));axis([0 5 0 180]);
      if ee==1,ylabel('RMS degrees'),  end;
      if ee==2,title('raw:  For 8s, 2s, 1s, 500ms, 250ms sized windows');
        legend({'Butterworth' 'FIR1' 'FIRWS'});xlabel('Filter order'); end; end
    
    figure(3+10*fd);for ee=1:size(angrmse,4),subplot(1,size(angrmse,4),ee);bar(squeeze(angrmsa(:,fd,:,ee)));axis([0 5 0 180]);
      if ee==1,ylabel('RMS degrees'),  end;
      if ee==2,title('raw:  For 8s, 2s, 1s, 500ms, 250ms sized windows');
        legend({'Butterworth' 'FIR1' 'FIRWS'});xlabel('Filter order'); end; end
    
    figure(4+10*fd);for ee=1:size(angrmse,4),subplot(1,size(angrmse,4),ee);bar(squeeze(angtrmse(:,fd,:,ee)));axis([0 5 0 180]);
      if ee==1,ylabel('RMS degrees'),  end;
      if ee==2,title('tlocktheta:  For 8s, 2s, 1s, 500ms, 250ms sized windows');
        legend({'Butterworth' 'FIR1' 'FIRWS'});xlabel('Filter order'); end; end
    
    figure(5+10*fd);for ee=1:size(angrmse,4),subplot(1,size(angrmse,4),ee);bar(squeeze(angtrmsh(:,fd,:,ee)));axis([0 5 0 180]);
      if ee==1,ylabel('RMS degrees'),  end;
      if ee==2,title('tlocktheta:  For 8s, 2s, 1s, 500ms, 250ms sized windows');
        legend({'Butterworth' 'FIR1' 'FIRWS'});xlabel('Filter order'); end; end
    
    figure(6+10*fd);for ee=1:size(angrmse,4),subplot(1,size(angrmse,4),ee);bar(squeeze(angtrmsa(:,fd,:,ee)));axis([0 5 0 180]);
      if ee==1,ylabel('RMS degrees'),  end;
      if ee==2,title('tlocktheta:  For 8s, 2s, 1s, 500ms, 250ms sized windows');
        legend({'Butterworth' 'FIR1' 'FIRWS'});xlabel('Filter order'); end; end
  end
  
  figure(100);for tf=[1 2 3 8 7 6 5 4],subplot(1,8,tf);bar(angfrmse(:,:,3,tf));axis([0 5 0 180]);
    set(get(100,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 8 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(101);for tf=[1 2 3 8 7 6 5 4],subplot(1,8,tf);bar(angfrmsh(:,:,3,tf));axis([0 5 0 180]);
    set(get(101,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 8 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(102);for tf=[1 2 3 8 7 6 5 4],subplot(1,8,tf);bar(angfrmsa(:,:,3,tf));axis([0 5 0 180]);
    set(get(102,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 8 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(103);for tf=[1 2 3 8 7 6 5 4],subplot(1,8,tf);bar(angftrmse(:,:,3,tf));axis([0 5 0 180]);
    set(get(103,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 8 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(104);for tf=[1 2 3 8 7 6 5 4],subplot(1,8,tf);bar(angftrmsh(:,:,3,tf));axis([0 5 0 180]);
    set(get(104,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 8 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(105);for tf=[1 2 3 8 7 6 5 4],subplot(1,8,tf);bar(angftrmsa(:,:,3,tf));axis([0 5 0 180]);
    set(get(105,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 8 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  
  figure(200);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwrmse(:,:,ww,1));axis([0 5 0 180]);
    set(get(200,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(201);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwrmsh(:,:,ww,1));axis([0 5 0 180]);
    set(get(201,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(202);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwrmsa(:,:,ww,1));axis([0 5 0 180]);
    set(get(202,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(203);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwtrmse(:,:,ww,1));axis([0 5 0 180]);
    set(get(203,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(204);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwtrmsh(:,:,ww,1));axis([0 5 0 180]);
    set(get(204,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(205);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwtrmsa(:,:,ww,1));axis([0 5 0 180]);
    set(get(205,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(220);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwrmse(:,:,ww,2));axis([0 5 0 180]);
    set(get(220,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(221);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwrmsh(:,:,ww,2));axis([0 5 0 180]);
    set(get(221,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(222);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwrmsa(:,:,ww,2));axis([0 5 0 180]);
    set(get(222,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('raw:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(223);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwtrmse(:,:,ww,2));axis([0 5 0 180]);
    set(get(223,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(224);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwtrmsh(:,:,ww,2));axis([0 5 0 180]);
    set(get(224,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  figure(225);for ww=1:length(wavwidth),subplot(1,length(wavwidth),ww);bar(angwtrmsa(:,:,ww,2));axis([0 5 0 180]);
    set(get(225,'Children'),'xTickLabel',{'4', '5', '6', '7Hz'})
    if tf==1,ylabel('RMS degrees'),  end;
    if tf==2,title('tlocktheta:  For 7 diff sized windows');
      legend({'Centre on 0' 'Centre -pi'});xlabel('Frequency (4,5,6,7) Hz'); end; end
  
  
end

%% 4) Create new data, with ERP just after time point of interest, with additive model
% One raw structure with 100 trials: each at slight diff phase and freq
% within alpha band; then additive ERP

if run4
  
  % params for whole simulation
  clearvars -except run* time wdr fsample wav*width
  close all;
  numtrials=100;
  halftimeind=round(length(time)/2);
  plotflag=1;
  foilim=[8 12];
  
  cd('D:\fieldtrip_svn\utilities\private');
  state=randomseed(13);
  cd(wdr);
  
  
  
  % New way:
  % base frequency 8,9,10,11,or 12 Hz
  % vary realisation of noise added (x100)
  % vary amplitude of ERP added (3 levels?)
  
  if 0 % old way
    clear raw
    raw.label{1}='test';
    raw.dimord='chan_time';
    [raw.time{1:numtrials}]=deal(time);
    rawpn=raw;
    rawerp=raw;
    rawpnerp=raw;
    %     phaseshift=2*pi*rand(1,length(freqtest))-pi; % random phase added, but fixed per frequency over all tests
    for tr=1:numtrials % simulate 100 trials
      % Always use 'cos' to generate signals
      % Each trial has a base frequency, pink noise, and ERP added
      %       frequse=ceil(12*rand);
      frequse(tr) = diff(foilim)*rand(1)+foilim(1);
      phaseshift(tr) = wrapToPi(2*pi*rand(1));
      trueang(tr,:) = frequse(tr)*2*pi*time+phaseshift(tr)*ones(size(time)); % random freq and phase
      raw.trial{tr}=cos(trueang(tr,:));
      trueang_prenoise(tr,:)=rad2deg(wrapToPi(trueang(tr,:)));
      rawpn.trial{tr}=raw.trial{tr}+10*pinknoise(length(time));
      trueang_postnoise(tr,:)=rad2deg(wrapToPi(angle(hilbert(rawpn.trial{tr})))); % after noise added.
      % add sinusoid (starting at amplitude zero) from time 0-100ms, at 10.1Hz frequency
      rawpnerp.trial{tr}=rawpn.trial{tr}+[zeros(1,halftimeind-1), gausswin(301,2)'.*sin(10.1*2*pi*time(halftimeind:halftimeind+300)-10.1*2*pi*time(halftimeind)), zeros(1,halftimeind-301)];
      rawerp.trial{tr}  =raw.trial{tr}  +[zeros(1,halftimeind-1), gausswin(301,2)'.*sin(10.1*2*pi*time(halftimeind:halftimeind+300)-10.1*2*pi*time(halftimeind)), zeros(1,halftimeind-301)];
    end
    cfg=[];
    tlock=ft_timelockanalysis(cfg,raw);
    tlockpn=ft_timelockanalysis(cfg,rawpn);
    tlockerp=ft_timelockanalysis(cfg,rawerp);
    tlockpnerp=ft_timelockanalysis(cfg,rawpnerp);
    figure;plot(tlock.time,tlock.avg)
    hold on;plot(tlockpn.time,tlockpn.avg,'g');axis([-inf inf -2 2])
    hold on;plot(tlockerp.time,tlockerp.avg,'r');axis([-inf inf -2 2])
    hold on;plot(tlockerp.time,tlockpnerp.avg,'k');axis([-inf inf -2 2])
    if 0
      cfg=[];
      ft_databrowser(cfg,raw);
      ft_databrowser(cfg,rawerp);
    end
  end
  
  % New way:
  clear raw
  erpamp=[0.5 1 2];
  for frequse=8:12
    for erpampind=1:length(erpamp)
      raw{frequse,erpampind}.label{1}='test';
      raw{frequse,erpampind}.dimord='chan_time';
      [raw{frequse,erpampind}.time{1:numtrials}]=deal(time);
      rawpn{frequse,erpampind}=raw{frequse,erpampind};
      rawerp{frequse,erpampind}=raw{frequse,erpampind};
      rawpnerp{frequse,erpampind}=raw{frequse,erpampind};
      for tr=1:numtrials % simulate 100 trials
        % Always use 'cos' to generate signals
        % Each trial has a base frequency, pink noise, and ERP added
        phaseshift(frequse,erpampind,tr) = wrapToPi(2*pi*rand(1));
        trueang(:,frequse,erpampind,tr) = frequse*2*pi*time+phaseshift(frequse,erpampind,tr)*ones(size(time)); % random freq and phase
        raw{frequse,erpampind}.trial{tr}=cos(trueang(:,frequse,erpampind,tr))';
        trueang_prenoise(:,frequse,erpampind,tr)=rad2deg(wrapToPi(trueang(:,frequse,erpampind,tr)));
        rawpn{frequse,erpampind}.trial{tr}=raw{frequse,erpampind}.trial{tr}+10*pinknoise(length(time));
        trueang_postnoise(:,frequse,erpampind,tr)=rad2deg(wrapToPi(angle(hilbert(rawpn{frequse,erpampind}.trial{tr})))); % after noise added.
        % add sinusoid (starting at amplitude zero) from time 0-100ms, at 10.1Hz frequency
        rawpnerp{frequse,erpampind}.trial{tr}=rawpn{frequse,erpampind}.trial{tr}+[zeros(1,halftimeind-1), erpamp(erpampind)*gausswin(301,2)'.*sin(10.1*2*pi*time(halftimeind:halftimeind+300)-10.1*2*pi*time(halftimeind)), zeros(1,halftimeind-301)];
        rawerp{frequse,erpampind}.trial{tr}  =raw{frequse,erpampind}.trial{tr}  +[zeros(1,halftimeind-1), erpamp(erpampind)*gausswin(301,2)'.*sin(10.1*2*pi*time(halftimeind:halftimeind+300)-10.1*2*pi*time(halftimeind)), zeros(1,halftimeind-301)];
      end
      if plotflag
        cfg=[];
        tlock=ft_timelockanalysis(cfg,raw{frequse,erpampind});
        tlockpn=ft_timelockanalysis(cfg,rawpn{frequse,erpampind});
        tlockerp=ft_timelockanalysis(cfg,rawerp{frequse,erpampind});
        tlockpnerp=ft_timelockanalysis(cfg,rawpnerp{frequse,erpampind});
        figure;plot(tlock.time,tlock.avg)
        hold on;plot(tlockpn.time,tlockpn.avg,'g');axis([-inf inf -2 2])
        hold on;plot(tlockerp.time,tlockerp.avg,'r');axis([-inf inf -2 2])
        hold on;plot(tlockerp.time,tlockpnerp.avg,'k');axis([-inf inf -2 2])
      end
    end
  end
  
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 4a) Bandpass filter + Hilbert
  
  close all
  %     freqtest=[1 4; 4 7; 8 12];
  freqtest=[8 12];
  
  timwin=[4 2 1 .5 .25]; % duration in seconds; must be ordered from longest to shortest
  
  %   pn=0; % use rawpnuse (else use rawuse)
  
  
  for tt=1:length(timwin)
    
    for frequse=8:12
      for erpampind=1:length(erpamp)
        
        cfg=[];
        cfg.latency=time([halftimeind-fsample*timwin(tt)/2 halftimeind+fsample*timwin(tt)/2]);
        rawuse=ft_selectdata(cfg,raw{frequse,erpampind});
        rawpnuse=ft_selectdata(cfg,rawpn{frequse,erpampind});
        rawerpuse=ft_selectdata(cfg,rawerp{frequse,erpampind});
        rawpnerpuse=ft_selectdata(cfg,rawpnerp{frequse,erpampind});
        halftimeinduse(tt)=dsearchn(rawuse.time{1}',time(halftimeind));
        
        for ff=1:size(freqtest,1)
          
          if 0
            % Bandpass filter: Butterworth
            rawpn_bpbut{3,2}=[];
            raw_bpbut{3,2}=[];
            rawerp_bpbut{3,2}=[];
            rawpnerp_bpbut{3,2}=[];
            cfg=[];
            cfg.bpfilter='yes';
            cfg.bpfreq=freqtest(ff,:);
            cfg.bpfilttype='but';
            cfg.plotfiltresp='yes';
            cfg.fouse=2:4;
            cfg.figind=0;
            cfg.plotflag=0;
            rawpn_bpbutr=filter4phase_estim8(cfg,rawpnuse);
            raw_bpbutr=filter4phase_estim8(cfg,rawuse);
            rawerp_bpbutr=filter4phase_estim8(cfg,rawerpuse);
            rawpnerp_bpbutr=filter4phase_estim8(cfg,rawpnerpuse);
            cfg.hilbert='complex';
            rawpn_bpbut=filter4phase_estim8(cfg,rawpnuse);
            raw_bpbut=filter4phase_estim8(cfg,rawuse);
            rawerp_bpbut=filter4phase_estim8(cfg,rawerpuse);
            rawpnerp_bpbut=filter4phase_estim8(cfg,rawpnerpuse);
          end
          
          % Bandpass filter: FIR (Matlab 'fir1')
          rawpn_bpfir{3,2}=[];
          raw_bpfir{3,2}=[];
          rawerp_bpfir{3,2}=[];
          rawpnerp_bpfir{3,2}=[];
          cfg=[];
          cfg.bpfilter='yes';
          cfg.bpfreq=freqtest(ff,:);
          cfg.bpfilttype='fir';
          cfg.plotfiltresp='yes';
          %       cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4 ]; % 3* is default
          %       cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
          cfg.fouse=[3*fsample/cfg.bpfreq(1)]; % 3* is default
          cfg.figind=10;
          cfg.plotflag=0;
          rawpn_bpfirr=filter4phase_estim8(cfg,rawpnuse);
          raw_bpfirr=filter4phase_estim8(cfg,rawuse);
          rawerp_bpfirr=filter4phase_estim8(cfg,rawerpuse);
          rawpnerp_bpfirr=filter4phase_estim8(cfg,rawpnerpuse);
          cfg.hilbert='complex';
          rawpn_bpfir=filter4phase_estim8(cfg,rawpnuse);
          raw_bpfir=filter4phase_estim8(cfg,rawuse);
          rawerp_bpfir=filter4phase_estim8(cfg,rawerpuse);
          rawpnerp_bpfir=filter4phase_estim8(cfg,rawpnerpuse);
          
          if 0
            % Bandpass filter: FIRWS
            cfg=[];
            cfg.bpfilter='yes';
            cfg.bpfreq=freqtest(ff,:);
            cfg.bpfilttype='firws';
            cfg.plotfiltresp='yes';
            cfg.fouse=2;
            cfg.figind=20;
            cfg.plotflag=0;
            rawpn_bpfirwsr=filter4phase_estim8(cfg,rawpnuse);
            raw_bpfirwsr=filter4phase_estim8(cfg,rawuse);
            rawerp_bpfirwsr=filter4phase_estim8(cfg,rawerpuse);
            rawpnerp_bpfirwsr=filter4phase_estim8(cfg,rawpnerpuse);
            cfg.hilbert='complex';
            rawpn_bpfirws=filter4phase_estim8(cfg,rawpnuse);
            raw_bpfirws=filter4phase_estim8(cfg,rawuse);
            rawerp_bpfirws=filter4phase_estim8(cfg,rawerpuse);
            rawpnerp_bpfirws=filter4phase_estim8(cfg,rawpnerpuse);
            if length(cfg.fouse)<3
              rawpn_bpfirws{3,2}=[];
              raw_bpfirws{3,2}=[];
              rawerp_bpfirws{3,2}=[];
              rawpnerp_bpfirws{3,2}=[];
            end
          end
          
          % reset these for every ff and tt
          if 0
            angdiffbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdifferpbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnerpbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            
            angdifffirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdifferpfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnerpfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeepbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeperpbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnerpbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeepfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeperpfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnerpfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
          end
          
          angdifffir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdiffpnfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdifferpfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdiffpnerpfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeepfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeeppnfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeeperpfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeeppnerpfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          
          
          for ss=1:numtrials
            for fo=1:size(raw_bpfir,1)
              for fd=1:size(raw_bpfir,2)
                if 0
                  if ~isempty(raw_bpbut{fo,fd})
                    angkeepbut(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeperpbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawerp_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnerpbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnerp_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  end
                  if ~isempty(raw_bpfirws{fo,fd})
                    angkeepfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeperpfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawerp_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnerpfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnerp_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  end
                end
                if ~isempty(raw_bpfir{fo,fd})
                  angkeepfir(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeeppnfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeeperpfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawerp_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeeppnerpfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnerp_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                end
              end
            end % fo
            % don't need to index these angdiff* over ff as we don't save them
            if 0
              angdiffbutpre(:,:,ss)=anglediff(angkeepbut(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdiffpnbutpre(:,:,ss)=anglediff(angkeeppnbut(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdifferpbutpre(:,:,ss)=anglediff(angkeeperpbut(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdiffpnerpbutpre(:,:,ss)=anglediff(angkeeppnerpbut(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdifffirwspre(:,:,ss)=anglediff(angkeepfirws(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdiffpnfirwspre(:,:,ss)=anglediff(angkeeppnfirws(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdifferpfirwspre(:,:,ss)=anglediff(angkeeperpfirws(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdiffpnerpfirwspre(:,:,ss)=anglediff(angkeeppnerpfirws(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              
              angdiffbutpost(:,:,ss)=anglediff(angkeepbut(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdiffpnbutpost(:,:,ss)=anglediff(angkeeppnbut(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdifferpbutpost(:,:,ss)=anglediff(angkeeperpbut(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdiffpnerpbutpost(:,:,ss)=anglediff(angkeeppnerpbut(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdifffirwspost(:,:,ss)=anglediff(angkeepfirws(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdiffpnfirwspost(:,:,ss)=anglediff(angkeeppnfirws(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdifferpfirwspost(:,:,ss)=anglediff(angkeeperpfirws(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdiffpnerpfirwspost(:,:,ss)=anglediff(angkeeppnerpfirws(:,:,ss),trueang_postnoise(ss,halftimeind),1);
            end
            
            angdifffirpre(:,:,ss)=anglediff(angkeepfir(:,:,ss),trueang_prenoise(halftimeind,frequse,erpampind,ss),1);
            angdiffpnfirpre(:,:,ss)=anglediff(angkeeppnfir(:,:,ss),trueang_prenoise(halftimeind,frequse,erpampind,ss),1);
            angdifferpfirpre(:,:,ss)=anglediff(angkeeperpfir(:,:,ss),trueang_prenoise(halftimeind,frequse,erpampind,ss),1);
            angdiffpnerpfirpre(:,:,ss)=anglediff(angkeeppnerpfir(:,:,ss),trueang_prenoise(halftimeind,frequse,erpampind,ss),1);
            
            angdifffirpost(:,:,ss)=anglediff(angkeepfir(:,:,ss),trueang_postnoise(halftimeind,frequse,erpampind,ss),1);
            angdiffpnfirpost(:,:,ss)=anglediff(angkeeppnfir(:,:,ss),trueang_postnoise(halftimeind,frequse,erpampind,ss),1);
            angdifferpfirpost(:,:,ss)=anglediff(angkeeperpfir(:,:,ss),trueang_postnoise(halftimeind,frequse,erpampind,ss),1);
            angdiffpnerpfirpost(:,:,ss)=anglediff(angkeeppnerpfir(:,:,ss),trueang_postnoise(halftimeind,frequse,erpampind,ss),1);
            
          end % ss
          
          if 0
            angbutpre(:,:,tt,ff)=rms(angdiffbutpre,3);
            angfirwspre(:,:,tt,ff)=rms(angdifffirwspre,3);
            angpnbutpre(:,:,tt,ff)=rms(angdiffpnbutpre,3);
            angpnfirwspre(:,:,tt,ff)=rms(angdiffpnfirwspre,3);
            angerpbutpre(:,:,tt,ff)=rms(angdifferpbutpre,3);
            angerpfirwspre(:,:,tt,ff)=rms(angdifferpfirwspre,3);
            angpnerpbutpre(:,:,tt,ff)=rms(angdiffpnerpbutpre,3);
            angpnerpfirwspre(:,:,tt,ff)=rms(angdiffpnerpfirwspre,3);
            
            angbutpost(:,:,tt,ff)=rms(angdiffbutpost,3);
            angfirwspost(:,:,tt,ff)=rms(angdifffirwspost,3);
            angpnbutpost(:,:,tt,ff)=rms(angdiffpnbutpost,3);
            angpnfirwspost(:,:,tt,ff)=rms(angdiffpnfirwspost,3);
            angerpbutpost(:,:,tt,ff)=rms(angdifferpbutpost,3);
            angerpfirwspost(:,:,tt,ff)=rms(angdifferpfirwspost,3);
            angpnerpbutpost(:,:,tt,ff)=rms(angdiffpnerpbutpost,3);
            angpnerpfirwspost(:,:,tt,ff)=rms(angdiffpnerpfirwspost,3);
            
            angfirpre(:,:,tt,ff)=rms(angdifffirpre,3);
            angpnfirpre(:,:,tt,ff)=rms(angdiffpnfirpre,3);
            angerpfirpre(:,:,tt,ff)=rms(angdifferpfirpre,3);
            angpnerpfirpre(:,:,tt,ff)=rms(angdiffpnerpfirpre,3);
            angfirpost(:,:,tt,ff)=rms(angdifffirpost,3);
            angpnfirpost(:,:,tt,ff)=rms(angdiffpnfirpost,3);
            angerpfirpost(:,:,tt,ff)=rms(angdifferpfirpost,3);
            angpnerpfirpost(:,:,tt,ff)=rms(angdiffpnerpfirpost,3);
          end
          
          % Note: not saving over 'ff' (but this is just alpha band anyway)
          % figure;rose(deg2rad(squeeze(angdiffpnerpfirpre(1,1,:))))
          angfirprem(tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdifffirpre(1,1,:))),3) /100);
          angfirprep(tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdifffirpre(1,1,:))),3) /100));
          angpnfirprem(tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdiffpnfirpre(1,1,:))),3) /100);
          angpnfirprep(tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnfirpre(1,1,:))),3) /100));
          angerpfirprem(tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdifferpfirpre(1,1,:))),3) /100);
          angerpfirprep(tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdifferpfirpre(1,1,:))),3) /100));
          angpnerpfirprem(tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdiffpnerpfirpre(1,1,:))),3) /100);
          angpnerpfirprep(tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnerpfirpre(1,1,:))),3) /100));
          
          
          %       angfirprem(tt,frequse,erpampind)=mean(abs(squeeze(angdifffirpre(1,1,:))));
          %       angfirprev(tt,frequse,erpampind)=std(abs(squeeze(angdifffirpre(1,1,:))));
          %       angpnfirprem(tt,frequse,erpampind)=mean(abs(squeeze(angdiffpnfirpre(1,1,:))));
          %       angpnfirprev(tt,frequse,erpampind)=std(abs(squeeze(angdiffpnfirpre(1,1,:))));
          %       angerpfirprem(tt,frequse,erpampind)=mean(abs(squeeze(angdifferpfirpre(1,1,:))));
          %       angerpfirprev(tt,frequse,erpampind)=std(abs(squeeze(angdifferpfirpre(1,1,:))));
          %       angpnerpfirprem(tt,frequse,erpampind)=mean(abs(squeeze(angdiffpnerpfirpre(1,1,:))));
          %       angpnerpfirprev(tt,frequse,erpampind)=std(abs(squeeze(angdiffpnerpfirpre(1,1,:))));
          %
          %       angfirpostm(tt,frequse,erpampind)=mean(abs(squeeze(angdifffirpost(1,1,:))));
          %       angfirpostv(tt,frequse,erpampind)=std(abs(squeeze(angdifffirpost(1,1,:))));
          %       angpnfirpostm(tt,frequse,erpampind)=mean(abs(squeeze(angdiffpnfirpost(1,1,:))));
          %       angpnfirpostv(tt,frequse,erpampind)=std(abs(squeeze(angdiffpnfirpost(1,1,:))));
          %       angerpfirpostm(tt,frequse,erpampind)=mean(abs(squeeze(angdifferpfirpost(1,1,:))));
          %       angerpfirpostv(tt,frequse,erpampind)=std(abs(squeeze(angdifferpfirpost(1,1,:))));
          %       angpnerpfirpostm(tt,frequse,erpampind)=mean(abs(squeeze(angdiffpnerpfirpost(1,1,:))));
          %       angpnerpfirpostv(tt,frequse,erpampind)=std(abs(squeeze(angdiffpnerpfirpost(1,1,:))));
          
        end % ff
      end % erpampind
    end % frequse
  end % tt
  clear angdiff*
  
  if 0
    save(['D:\phase_estimation\anglermsrun4.mat'],'ang*but*','ang*fir*','frequse','-append')
  else
    try
      save(['D:\phase_estimation\anglermsrun4b.mat'],'ang*fir*','-append')
    catch
      save(['D:\phase_estimation\anglermsrun4b.mat'],'ang*fir*')
    end
  end
  
  load(['D:\phase_estimation\anglermsrun4b.mat'])
  
  % "1 - circular-mean" is the same as circular-variance
  figure(1);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);                     bar(1-squeeze(angfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('No PN, No ERP');end
    subplot(4,length(timwin),tt+length(timwin));      bar(1-squeeze(angpnfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    if tt==1,ylabel('PN, No ERP');end
    subplot(4,length(timwin),tt+2*length(timwin));    bar(1-squeeze(angerpfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    if tt==1,ylabel('No PN, ERP');end
    subplot(4,length(timwin),tt+3*length(timwin));    bar(1-squeeze(angpnerpfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    if tt==1,ylabel('PN, ERP');end
  end
  legend({'8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  set(get(1,'Children'),'xTickLabel',{'ERPamp 1' 'ERPamp 2' 'ERPamp 3'})
  
  figure(2);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);                     bar(squeeze(angfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('No PN, No ERP');end
    subplot(4,length(timwin),tt+length(timwin));      bar(squeeze(angpnfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('PN, No ERP');end
    subplot(4,length(timwin),tt+2*length(timwin));    bar(squeeze(angerpfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('No PN, ERP');end
    subplot(4,length(timwin),tt+3*length(timwin));    bar(squeeze(angpnerpfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('PN, ERP');end
  end
  legend({'8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  set(get(1,'Children'),'xTickLabel',{'ERPamp 1' 'ERPamp 2' 'ERPamp 3'})
  
  figure(3);
  for tt=1:length(timwin),
    subplot(2,length(timwin),tt+0*length(timwin));    bar(1-squeeze(angpnerpfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('Circ-Var');end
    subplot(2,length(timwin),tt+1*length(timwin));    bar(squeeze(angpnerpfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('Mean-Angle');end
  end
  legend({'8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  set(get(3,'Children'),'xTickLabel',{'ERPamp 1' 'ERPamp 2' 'ERPamp 3'})
  
  if 0
    %'Compared to trueang_prenoise')
    figure(1);
    for tt=1:length(timwin),
      subplot(3,length(timwin),tt);    bar([squeeze(angbutpre(:,2,tt,ff)), squeeze(angbutpre(:,2,tt,ff)), squeeze(angerpbutpre(:,2,tt,ff)), squeeze(angpnerpbutpre(:,2,tt,ff))]);axis([-inf inf 0 160])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Butterworth');end
      subplot(3,length(timwin),tt+5);    bar([squeeze(angfirpre(:,2,tt,ff)), squeeze(angfirpre(:,2,tt,ff)), squeeze(angerpfirpre(:,2,tt,ff)), squeeze(angpnerpfirpre(:,2,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR');end
      subplot(3,length(timwin),tt+10);    bar([squeeze(angfirwspre(:,2,tt,ff)), squeeze(angfirwspre(:,2,tt,ff)), squeeze(angerpfirwspre(:,2,tt,ff)), squeeze(angpnerpfirwspre(:,2,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIRWS');end
    end
    legend({'Orig signal' 'Noise, but no ERP added', 'ERP added', 'Noise and ERP added'})
    set(get(1,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    %     %'Compared to trueang_postnoise'
    %   figure(2);
    %   for tt=1:length(timwin),
    %     subplot(2,length(timwin),tt);  bar([squeeze(angbutpost(:,2,tt,ff)), squeeze(angerpbutpost(:,2,tt,ff))]);axis([-inf inf 0 160])
    %     title(['Time Window length ' num2str(timwin(tt))])
    %     if tt==1,ylabel('Butterworth');end
    %     subplot(2,length(timwin),tt+5);  bar([squeeze(angfirpost(:,2,tt,ff)), squeeze(angerpfirpost(:,2,tt,ff))]);axis([-inf inf 0 160])
    %     if tt==1,ylabel('FIR');end
    %   end
    %   legend({'No ERP added', 'ERP added'})
    %   set(get(2,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    %'Compared to trueang_prenoise')
    figure(3);
    for tt=1:length(timwin),
      subplot(3,length(timwin),tt);    bar([squeeze(angbutpre(:,1,tt,ff)), squeeze(angpnbutpre(:,1,tt,ff)), squeeze(angerpbutpre(:,1,tt,ff)), squeeze(angpnerpbutpre(:,1,tt,ff))]);axis([-inf inf 0 160])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Butterworth');end
      subplot(3,length(timwin),tt+5);    bar([squeeze(angfirpre(:,1,tt,ff)), squeeze(angpnfirpre(:,1,tt,ff)), squeeze(angerpfirpre(:,1,tt,ff)), squeeze(angpnerpfirpre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR');end
      subplot(3,length(timwin),tt+10);    bar([squeeze(angfirwspre(:,1,tt,ff)), squeeze(angpnfirwspre(:,1,tt,ff)), squeeze(angerpfirwspre(:,1,tt,ff)), squeeze(angpnerpfirwspre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR');end
    end
    legend({'Orig signal' 'Noise, but no ERP added', 'ERP added', 'Noise and ERP added'})
    set(get(3,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    figure(4);
    for tt=1:length(timwin),
      subplot(3,length(timwin),tt);    bar([squeeze(angbutpre(:,2,tt,ff)), squeeze(angpnbutpre(:,2,tt,ff)), squeeze(angerpbutpre(:,2,tt,ff)), squeeze(angpnerpbutpre(:,2,tt,ff))]);axis([-inf inf 0 160])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Butterworth 2pass');end
      subplot(3,length(timwin),tt+5);    bar([squeeze(angfirpre(:,1,tt,ff)), squeeze(angpnfirpre(:,1,tt,ff)), squeeze(angerpfirpre(:,1,tt,ff)), squeeze(angpnerpfirpre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR 1pass-0phase');end
      subplot(3,length(timwin),tt+10);    bar([squeeze(angfirwspre(:,1,tt,ff)), squeeze(angpnfirwspre(:,1,tt,ff)), squeeze(angerpfirwspre(:,1,tt,ff)), squeeze(angpnerpfirwspre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIRWS 1pass-0phase');end
    end
    legend({'Orig signal' 'Noise, but no ERP added', 'ERP added', 'Noise and ERP added'})
    set(get(3,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  end
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 4b) FFT + taper
  
  cfg=[];
  cfg.method='mtmconvol';
  cfg.output='fourier';
  cfg.taper='hanning';
  cfg.foi=foilim(1):foilim(end);
  cfg.keeptrials='yes';
  
  t_ftimwin{1}=timwin(1)*ones(size(cfg.foi)); % full length of data
  t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 2 periods 4 Hz
  t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{6}=4./cfg.foi; %
  t_ftimwin{7}=3./cfg.foi; %
  t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
  t_ftimwin{9}=1./cfg.foi; % one period for each frequency
  t_ftimwin{10}=0.5./cfg.foi; %
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angfreq=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqerp=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpn=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpnerp=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  
  angfreqspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqerpspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpnspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpnerpspec=nan(numtrials,length(toi),length(t_ftimwin));
  
  angfreq10=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqerp10=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpn10=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpnerp10=nan(numtrials,length(toi),length(t_ftimwin));
  
  for frequse=8:12
    for erpampind=1:length(erpamp)
      for ss=1:numtrials
        cfg.trials=ss;
        for tt=1:length(toi)
          cfg.toi=toi{tt};
          for tf=1:length(t_ftimwin)
            cfg.t_ftimwin=t_ftimwin{tf};
            freq  =ft_freqanalysis(cfg, raw{frequse,erpampind});
            freqerp=ft_freqanalysis(cfg, rawerp{frequse,erpampind});
            freqpn  =ft_freqanalysis(cfg, rawpn{frequse,erpampind});
            freqpnerp=ft_freqanalysis(cfg, rawpnerp{frequse,erpampind});
            if tt==1
              %           angfreq(:,:,tt,pp,tf)=angle(squeeze(freq{tt,pp,tf}.fourierspctrm))/(2*pi)*360;
              angfreq(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freq.fourierspctrm))));
              angfreqerp(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqerp.fourierspctrm))));
              angfreqpn(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freqpn.fourierspctrm))));
              angfreqpnerp(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqpnerp.fourierspctrm))));
            elseif tt==2
              angfreq(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freq.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqerp(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqerp.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpn(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freqpn.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpnerp(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqpnerp.fourierspctrm(1,1,:,:)))+pi  )));
            end
          end % tf
        end % tt
        if 0
          % given that we know exactly which freq was used for a given trial
          angfreqspec(ss,:,:)   =angfreq(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqerpspec(ss,:,:)=angfreqerp(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqpnspec(ss,:,:)   =angfreqpn(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqpnerpspec(ss,:,:)=angfreqpnerp(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          % use the freq that was closest to the ERP (10.1Hz)
          angfreq10(ss,:,:)   =angfreq(ss,3,:,:);
          angfreqerp10(ss,:,:)=angfreqerp(ss,3,:,:);
          angfreqpn10(ss,:,:)   =angfreqpn(ss,3,:,:);
          angfreqpnerp10(ss,:,:)=angfreqpnerp(ss,3,:,:);
          
          angfreqspecdiff(ss,:,:,:)   =anglediff(angfreqspec(ss,:,:,:),   trueang_prenoise(ss,halftimeind),1);
          angfreqerpspecdiff(ss,:,:,:)=anglediff(angfreqerpspec(ss,:,:,:),trueang_prenoise(ss,halftimeind),1);
          angfreqpnspecdiff(ss,:,:,:)   =anglediff(angfreqpnspec(ss,:,:,:),   trueang_prenoise(ss,halftimeind),1);
          angfreqpnerpspecdiff(ss,:,:,:)=anglediff(angfreqpnerpspec(ss,:,:,:),trueang_prenoise(ss,halftimeind),1);
          
          angfreq10diff(ss,:,:,:)     =anglediff(angfreq10(ss,:,:,:),     trueang_prenoise(ss,halftimeind),1);
          angfreqerp10diff(ss,:,:,:)  =anglediff(angfreqerp10(ss,:,:,:),  trueang_prenoise(ss,halftimeind),1);
          angfreqpn10diff(ss,:,:,:)     =anglediff(angfreqpn10(ss,:,:,:),     trueang_prenoise(ss,halftimeind),1);
          angfreqpnerp10diff(ss,:,:,:)  =anglediff(angfreqpnerp10(ss,:,:,:),  trueang_prenoise(ss,halftimeind),1);
        end
        
        % compute differences
        angfreqdiff(ss,:,:,:)       =anglediff(angfreq(ss,:,:,:),       trueang_prenoise(halftimeind,frequse,erpampind,ss),1);
        angfreqerpdiff(ss,:,:,:)    =anglediff(angfreqerp(ss,:,:,:),    trueang_prenoise(halftimeind,frequse,erpampind,ss),1);
        angfreqpndiff(ss,:,:,:)       =anglediff(angfreqpn(ss,:,:,:),       trueang_prenoise(halftimeind,frequse,erpampind,ss),1);
        angfreqpnerpdiff(ss,:,:,:)    =anglediff(angfreqpnerp(ss,:,:,:),    trueang_prenoise(halftimeind,frequse,erpampind,ss),1);
        
      end % ss
      if 0
        angfrms  =squeeze(rms(angfreqdiff,1));
        angferms =squeeze(rms(angfreqerpdiff,1));
        angfprms  =squeeze(rms(angfreqpndiff,1));
        angfperms =squeeze(rms(angfreqpnerpdiff,1));
        
        angfsrms =squeeze(rms(angfreqspecdiff,1));
        angfesrms=squeeze(rms(angfreqerpspecdiff,1));
        angfpsrms =squeeze(rms(angfreqpnspecdiff,1));
        angfpesrms=squeeze(rms(angfreqpnerpspecdiff,1));
        
        angftrms =squeeze(rms(angfreq10diff,1));
        angfetrms=squeeze(rms(angfreqerp10diff,1));
        angfptrms =squeeze(rms(angfreqpn10diff,1));
        angfpetrms=squeeze(rms(angfreqpnerp10diff,1));
      end
      
      %   angfmean(:,:,:,frequse,erpampind)=squeeze(mean(angfreqdiff,1));
      %   angfemean(:,:,:,frequse,erpampind)=squeeze(mean(angfreqerpdiff,1));
      %   angfpmean(:,:,:,frequse,erpampind)=squeeze(mean(angfreqpndiff,1));
      %   angfpemean(:,:,:,frequse,erpampind)=squeeze(mean(angfreqpnerpdiff,1));
      
      angfmean(:,:,:,frequse,erpampind)=  abs(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials);
      angfemean(:,:,:,frequse,erpampind)=abs(sum(exp(i*deg2rad(angfreqerpdiff)),1)/numtrials);
      angfpmean(:,:,:,frequse,erpampind)=abs(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials);
      angfpemean(:,:,:,frequse,erpampind)=abs(sum(exp(i*deg2rad(angfreqpnerpdiff)),1)/numtrials);
      angfang(:,:,:,frequse,erpampind)=  rad2deg(angle(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials));
      angfeang(:,:,:,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqerpdiff)),1)/numtrials));
      angfpang(:,:,:,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials));
      angfpeang(:,:,:,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpnerpdiff)),1)/numtrials));
      
    end % erpampind
  end % frequse
  
  if 0
    save(['D:\phase_estimation\anglermsrun4.mat'],'ang*rms*','frequse','-append')
  else
    save(['D:\phase_estimation\anglermsrun4b.mat'],'ang*mean','ang*ang','-append')
  end
  
  load(['D:\phase_estimation\anglermsrun4b.mat'])
  
  for erpampind=1:3
    for tt=1:2
      figure(100+tt+(erpampind-1)*2);
      for tf=1:length(t_ftimwin),
        subplot(4,length(t_ftimwin),tf);                     bar(1-squeeze(angfmean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['TW ' num2str(10*t_ftimwin{tf}(3)) ' per.'])
        end
        if tf==1,ylabel('No PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+length(t_ftimwin));      bar(1-squeeze(angfpmean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+2*length(t_ftimwin));    bar(1-squeeze(angfemean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('No PN, ERP');end
        subplot(4,length(t_ftimwin),tf+3*length(t_ftimwin));    bar(1-squeeze(angfpemean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('PN, ERP');end
      end
      legend({'cfg.foi 8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
      set(get(100+tt+(erpampind-1)*2,'Children'),'xTickLabel',{'Sim 8' '9' '10' '11' '12'})
    end
  end
  
  for erpampind=1:3
    for tt=1:2
      figure(110+tt+(erpampind-1)*2);
      for tf=1:length(t_ftimwin),
        subplot(4,length(t_ftimwin),tf);                     bar(squeeze(angfang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['TW ' num2str(10*t_ftimwin{tf}(3)) ' per.'])
        end
        if tf==1,ylabel('No PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+length(t_ftimwin));      bar(squeeze(angfpang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+2*length(t_ftimwin));    bar(squeeze(angfeang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('No PN, ERP');end
        subplot(4,length(t_ftimwin),tf+3*length(t_ftimwin));    bar(squeeze(angfpeang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('PN, ERP');end
      end
      legend({'cfg.foi 8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
      set(get(110+tt+(erpampind-1)*2,'Children'),'xTickLabel',{'Sim 8' '9' '10' '11' '12'})
    end
  end
  
  for erpampind=1:3
    for tt=1:2
      figure(120+tt+(erpampind-1)*2);
      for tf=1:length(t_ftimwin),
        subplot(4,5,tf);       bar(1-squeeze(angfpemean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['TW ' num2str(10*t_ftimwin{tf}(3)) ' per.'])
        end
        if tf==1,ylabel('Circ-Var');end
        subplot(4,5,tf+10);    bar(squeeze(angfpeang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('Mean Angle');end
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['TW ' num2str(10*t_ftimwin{tf}(3)) ' per.'])
        end
      end
      %       set(get(120+tt+(erpampind-1)*2,'Children'),'xTickLabel',{'Sim 8' '9' '10' '11' '12'})
      set(get(120+tt+(erpampind-1)*2,'Children'),'XTickLabel',{'Sim 8' '9' '10' '11' '12'})
      legend({'cfg.foi 8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
    end
  end
  
  if 0
    figure(100);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfrms(:,:,tf));axis([-inf inf 0 160]);end
    figure(101);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angferms(:,:,tf));axis([-inf inf 0 160]);end
    figure(102);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfprms(:,:,tf));axis([-inf inf 0 160]);end
    figure(103);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfperms(:,:,tf));axis([-inf inf 0 160]);end
    figure(104);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([angfsrms(:,tf),angfesrms(:,tf),angftrms(:,tf),angfetrms(:,tf), angfpsrms(:,tf),angfpesrms(:,tf),angfptrms(:,tf),angfpetrms(:,tf)]');axis([-inf inf 0 160]);end
    set(get(104,'Children'),'xTickLabel',{'spec', 'ERP spec', '10hz', 'ERP 10hz', 'PN spec', 'PN ERP spec', 'PN 10hz', 'PN ERP 10hz'})
    
    % no point in using '10hz' as it's already in 8-12.
    figure(110);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 13 ],[angfrms(:,:,tf); angfsrms(:,tf)'; ]);axis([-inf inf 0 160]);end
    set(get(110,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(111);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 13 ],[angferms(:,:,tf); angfesrms(:,tf)'; ]);axis([-inf inf 0 160]);end
    set(get(111,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(112);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 13],[angfprms(:,:,tf); angfpsrms(:,tf)';]);axis([-inf inf 0 160]);end
    set(get(112,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(113);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 13 ],[angfperms(:,:,tf); angfpesrms(:,tf)'; ]);axis([-inf inf 0 160]);end
    set(get(113,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
  end
  
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 4c) Wavelet
  
  foi=foilim(1):foilim(end);
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  cfg.foi=foi;
  cfg.pad=time(end);
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angwave   =nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwaveerp=nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwavepn   =nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwavepnerp=nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  
  angwavediff    = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwaveerpdiff = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwavepndiff    = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwavepnerpdiff = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  
  angwavespecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwaveerpspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepnspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepnerpspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  
  angwave10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwaveerp10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepn10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepnerp10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  
  
  for tf=1:length(timwin) % more important than gwidth?
    scfg=[];
    scfg.latency=time([halftimeind-fsample*timwin(tf)/2 halftimeind+fsample*timwin(tf)/2]);
    rawuse=ft_selectdata(scfg,raw);
    rawerpuse=ft_selectdata(scfg,rawerp);
    rawpnuse=ft_selectdata(scfg,rawpn);
    rawpnerpuse=ft_selectdata(scfg,rawpnerp);
    halftimeinduse(tf)=dsearchn(rawuse.time{1}',time(halftimeind));
    
    for ss=1:numtrials
      cfg.trials=ss;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoi=ft_freqanalysis(cfg,rawuse);
          waveerpfoi=ft_freqanalysis(cfg,rawerpuse);
          wavepnfoi=ft_freqanalysis(cfg,rawpnuse);
          wavepnerpfoi=ft_freqanalysis(cfg,rawpnerpuse);
          indfoi1=unique(dsearchn(foi',wavefoi.freq'));
          if tt==1
            angwave(ss,indfoi1,tt,ww,1)   = rad2deg(wrapToPi(angle(squeeze(wavefoi.fourierspctrm))));
            angwaveerp(ss,indfoi1,tt,ww,1)= rad2deg(wrapToPi(angle(squeeze(waveerpfoi.fourierspctrm))));
            angwavepn(ss,indfoi1,tt,ww,1)   = rad2deg(wrapToPi(angle(squeeze(wavepnfoi.fourierspctrm))));
            angwavepnerp(ss,indfoi1,tt,ww,1)= rad2deg(wrapToPi(angle(squeeze(wavepnerpfoi.fourierspctrm))));
          elseif tt==2
            angwave(ss,indfoi1,tt,ww,1)   = rad2deg(diag(wrapToPi(   angle(squeeze(wavefoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwaveerp(ss,indfoi1,tt,ww,1)= rad2deg(diag(wrapToPi(   angle(squeeze(waveerpfoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwavepn(ss,indfoi1,tt,ww,1)   = rad2deg(diag(wrapToPi(   angle(squeeze(wavepnfoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwavepnerp(ss,indfoi1,tt,ww,1)= rad2deg(diag(wrapToPi(   angle(squeeze(wavepnerpfoi.fourierspctrm(1,1,:,:)))+pi   )));
          end
        end % ww
      end % tt
      % given that we know exactly which freq was used for a given trial
      angwavespec(ss,:,:)   =angwave(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwaveerpspec(ss,:,:)=angwaveerp(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwavepnspec(ss,:,:)   =angwavepn(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwavepnerpspec(ss,:,:)=angwavepnerp(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      % use the freq that was closest to the ERP (10.1Hz)
      angwave10(ss,:,:)   =angwave(ss,3,:,:,1);
      angwaveerp10(ss,:,:)=angwaveerp(ss,3,:,:,1);
      angwavepn10(ss,:,:)   =angwavepn(ss,3,:,:,1);
      angwavepnerp10(ss,:,:)=angwavepnerp(ss,3,:,:,1);
      
      % compute differences
      angwavediff(ss,:,:,:,tf,1)       =anglediff(angwave(ss,:,:,:,1),       trueang_prenoise(ss,halftimeind),1);
      angwaveerpdiff(ss,:,:,:,tf,1)    =anglediff(angwaveerp(ss,:,:,:,1),    trueang_prenoise(ss,halftimeind),1);
      angwavepndiff(ss,:,:,:,tf,1)       =anglediff(angwavepn(ss,:,:,:,1),       trueang_prenoise(ss,halftimeind),1);
      angwavepnerpdiff(ss,:,:,:,tf,1)    =anglediff(angwavepnerp(ss,:,:,:,1),    trueang_prenoise(ss,halftimeind),1);
      
      angwavespecdiff(ss,:,:,tf,1)   =anglediff(angwavespec(ss,:,:),   trueang_prenoise(ss,halftimeind),1);
      angwaveerpspecdiff(ss,:,:,tf,1)=anglediff(angwaveerpspec(ss,:,:),trueang_prenoise(ss,halftimeind),1);
      angwavepnspecdiff(ss,:,:,tf,1)   =anglediff(angwavepnspec(ss,:,:),   trueang_prenoise(ss,halftimeind),1);
      angwavepnerpspecdiff(ss,:,:,tf,1)=anglediff(angwavepnerpspec(ss,:,:),trueang_prenoise(ss,halftimeind),1);
      
      angwave10diff(ss,:,:,tf,1)     =anglediff(angwave10(ss,:,:),     trueang_prenoise(ss,halftimeind),1);
      angwaveerp10diff(ss,:,:,tf,1)  =anglediff(angwaveerp10(ss,:,:),  trueang_prenoise(ss,halftimeind),1);
      angwavepn10diff(ss,:,:,tf,1)     =anglediff(angwavepn10(ss,:,:),     trueang_prenoise(ss,halftimeind),1);
      angwavepnerp10diff(ss,:,:,tf,1)  =anglediff(angwavepnerp10(ss,:,:),  trueang_prenoise(ss,halftimeind),1);
    end % ss
  end % tf
  
  
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  cfg.foilim=foilim;
  cfg.pad=time(end);
  
  for tf=1:length(timwin) % more important than gwidth
    scfg=[];
    scfg.latency=time([halftimeind-fsample*timwin(tf)/2 halftimeind+fsample*timwin(tf)/2]);
    rawuse=ft_selectdata(scfg,raw);
    rawerpuse=ft_selectdata(scfg,rawerp);
    rawpnuse=ft_selectdata(scfg,rawpn);
    rawpnerpuse=ft_selectdata(scfg,rawpnerp);
    halftimeinduse(tf)=dsearchn(rawuse.time{1}',time(halftimeind));
    
    for ss=1:numtrials
      cfg.trials=ss;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoilim=ft_freqanalysis(cfg,rawuse);
          waveerpfoilim=ft_freqanalysis(cfg,rawerpuse);
          wavepnfoilim=ft_freqanalysis(cfg,rawpnuse);
          wavepnerpfoilim=ft_freqanalysis(cfg,rawpnerpuse);
          % this produces .freq of length varying with timwin (spacing is 1/timwin)
          
          % one option: take the middle freq within the range and use that angle.
          indfoi2=dsearchn(foi', round(wavefoilim.freq(ceil(length(wavefoilim.freq)/2))));
          if tt==1
            angwave(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi(angle(squeeze(wavefoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwaveerp(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(waveerpfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwavepn(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi(angle(squeeze(wavepnfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwavepnerp(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(wavepnerpfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
          elseif tt==2
            angwave(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi( angle(squeeze(wavefoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwaveerp(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi( angle(squeeze(waveerpfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwavepn(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi( angle(squeeze(wavepnfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwavepnerp(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi( angle(squeeze(wavepnerpfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
          end
          % another option: take average over
          if tt==1
            angwave(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(angle(squeeze(wavefoilim.fourierspctrm)))));
            angwaveerp(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(angle(squeeze(waveerpfoilim.fourierspctrm)))));
            angwavepn(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(angle(squeeze(wavepnfoilim.fourierspctrm)))));
            angwavepnerp(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(angle(squeeze(wavepnerpfoilim.fourierspctrm)))));
          elseif tt==2
            angwave(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavefoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwaveerp(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(  diag(angle(squeeze(waveerpfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwavepn(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavepnfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwavepnerp(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavepnerpfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
          end
          
        end % ww
      end % tt
      
      % compute differences
      for dd=2:3
        angwavediff(ss,indfoi2,:,:,tf,dd)       =anglediff(angwave(ss,indfoi2,:,:,dd),       trueang_prenoise(ss,halftimeind),1);
        angwaveerpdiff(ss,indfoi2,:,:,tf,dd)    =anglediff(angwaveerp(ss,indfoi2,:,:,dd),    trueang_prenoise(ss,halftimeind),1);
        angwavepndiff(ss,indfoi2,:,:,tf,dd)       =anglediff(angwavepn(ss,indfoi2,:,:,dd),       trueang_prenoise(ss,halftimeind),1);
        angwavepnerpdiff(ss,indfoi2,:,:,tf,dd)    =anglediff(angwavepnerp(ss,indfoi2,:,:,dd),    trueang_prenoise(ss,halftimeind),1);
      end
    end  % ss
  end % tf
  
  
  
  angwrms  =squeeze(rms(angwavediff,1));
  angwerms =squeeze(rms(angwaveerpdiff,1));
  angwprms  =squeeze(rms(angwavepndiff,1));
  angwperms =squeeze(rms(angwavepnerpdiff,1));
  
  angwsrms =squeeze(rms(angwavespecdiff,1));
  angwesrms=squeeze(rms(angwaveerpspecdiff,1));
  angwpsrms =squeeze(rms(angwavepnspecdiff,1));
  angwpesrms=squeeze(rms(angwavepnerpspecdiff,1));
  
  angwtrms =squeeze(rms(angwave10diff,1));
  angwetrms=squeeze(rms(angwaveerp10diff,1));
  angwptrms =squeeze(rms(angwavepn10diff,1));
  angwpetrms=squeeze(rms(angwavepnerp10diff,1));
  
  save(['D:\phase_estimation\anglermsrun4.mat'],'ang*rms*','frequse','-append')
  
  % 3rd dim index 4 means wavwidth=7
  for dd=1:3 % 3 diff ways of computing
    figure(200+dd*10);
    for tf=1:length(timwin),subplot(2,5,tf);bar(foi,angwrms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(201+dd*10);
    for tf=1:length(timwin),subplot(2,5,tf);bar(foi,angwerms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(202+dd*10);
    for tf=1:length(timwin),subplot(2,5,tf);bar(foi,angwprms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(203+dd*10);
    for tf=1:length(timwin),subplot(2,5,tf);bar(foi,angwperms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
  end
  figure(204);
  for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([angwsrms(:,tf),angwesrms(:,tf),angwtrms(:,tf),angwetrms(:,tf), angwpsrms(:,tf),angwpesrms(:,tf),angwptrms(:,tf),angwpetrms(:,tf)]');axis([-inf inf 0 160]);end
  set(get(204,'Children'),'xTickLabel',{'spec', 'ERP spec', '10hz', 'ERP 10hz', 'PN spec', 'PN ERP spec', 'PN 10hz', 'PN ERP 10hz'})
  
  % 3rd dim index 4 means wavwidth=7
  for dd=1:3 % 3 diff ways of computing
    figure(250+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar([foi 13],[angwrms(:,:,4,tf,dd); angwsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(250+dd*10,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(251+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar([foi 13],[angwerms(:,:,4,tf,dd); angwesrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(251+dd*10,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(252+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar([foi 13],[angwprms(:,:,4,tf,dd); angwpsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(252+dd*10,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(253+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar([foi 13],[angwperms(:,:,4,tf,dd); angwpesrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(253+dd*10,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
  end
  
end

%% 5) Create new data, with ERP just after time point of interest, with phase-resetting model

if run5
  
  clearvars -except run* time wdr fsample wav*width
  addpath('D:\phase_estimation\phasereset_modelling')
  halftimeind=round(length(time)/2);
  numtrials=100;
  
  foilim=[8 12]; % should be same as foilim in run4
  
  %   for ss=1:20
  % inspired by
  % phasereset (frames, epochs, srate, minfr, maxfr, position, tjitter)
  % from % Implemented by: Rafal Bogacz and Nick Yeung, September 2006
  
  numones=[1 4 7];
  resetwin=[40 20 0];
  
  % Model 1:  frequency shift to meet phase at later time; return to new frequency
  % Model 12: Same effectively as 1, except modelled as linear phase procession (at new frequency not determinted explicitly); return to new frequency
  % Model 2:  frequency shift to meet phase at later time; return to original frequency
  % Model 22: Same effectively as 2, except modelled as linear phase procession (at new frequency not determinted explicitly); return to original frequency
  % Model 3:  short time to meet new phase via line;
  % Model 4:  short time to meet new phase via line; new phase determined by 'phase-response-curve'
  modeluse=12;
  
  for frequse=foilim(1):foilim(2)
    for erpampind=1:length(resetwin)
      cd('D:\fieldtrip_svn\utilities\private');
      state=randomseed(13);
      cd(wdr);
      
      % control with no phase reset
      raw{frequse,erpampind}.label{1}='test';
      raw{frequse,erpampind}.dimord='chan_time';
      [raw{frequse,erpampind}.time{1:numtrials}]=deal(time);
      rawpn{frequse,erpampind}=raw{frequse,erpampind};
      rawpr{frequse,erpampind}=raw{frequse,erpampind};
      rawpnpr{frequse,erpampind}=raw{frequse,erpampind};
      for tr=1:numtrials
        switch modeluse
          case 1
            [raw{frequse,erpampind}.trial{tr}] = phasereset_jz (time, fsample, frequse, foilim(1), foilim(2), length(time)-1, length(time), 1);
          case 12
            [raw{frequse,erpampind}.trial{tr}] = phasereset_jz_linphase_postjitter (time, frequse, length(time)-1, length(time), 0 );
          case 3
            %       [raw{frequse,erpampind}.trial{tr}] = phasereset_jz_jump (time, frequse, length(time)-1, length(time), numones(erpampind) );
            [raw{frequse,erpampind}.trial{tr}] = phasereset_jz_jump (time, frequse, length(time)-1, length(time), 0 );
          case 2
            [raw{frequse,erpampind}.trial{tr}] = phasereset_jz_steadybasefreq (time, fsample, frequse, foilim(1), foilim(2), length(time)-1, length(time), 1);
          case 22
            [raw{frequse,erpampind}.trial{tr}] = phasereset_jz_linphase (time, frequse, length(time)-1, length(time), 0 );
          case 4 % phase response curve
            [raw{frequse,erpampind}.trial{tr}] = phasereset_jz_jump_prc (time, frequse, length(time)-1, length(time), 0 );
        end
        rawpn{frequse,erpampind}.trial{tr}=raw{frequse,erpampind}.trial{tr}+10*pinknoise(length(time));
      end
      
      % set back again to exact same state, thus same phase at time points until phase reset
      cd('D:\fieldtrip_svn\utilities\private');
      state=randomseed(13);
      cd(wdr);
      
      % do proper reset at halfway point
      for tr=1:numtrials
        switch modeluse
          case 1
            [rawpr{frequse,erpampind}.trial{tr},trueang_prenoise(tr,frequse,erpampind),newfreq(tr,frequse,erpampind)] = phasereset_jz (time, fsample, frequse, foilim(1), foilim(2), halftimeind, halftimeind+.1*fsample, numones(erpampind));
          case 12
            [rawpr{frequse,erpampind}.trial{tr},trueang_prenoise(tr,frequse,erpampind),resetang(tr,frequse,erpampind)] = phasereset_jz_linphase_postjitter (time, frequse, halftimeind, halftimeind+20, resetwin(erpampind) );
          case 3
            %       [rawpr{frequse,erpampind}.trial{tr},trueang_prenoise(tr,frequse,erpampind)] = phasereset_jz_jump (time, frequse, halftimeind, halftimeind+20, numones(erpampind) );
            [rawpr{frequse,erpampind}.trial{tr},trueang_prenoise(tr,frequse,erpampind),resetang(tr,frequse,erpampind)] = phasereset_jz_jump (time, frequse, halftimeind, halftimeind+20, resetwin(erpampind) );
          case 2
            [rawpr{frequse,erpampind}.trial{tr},trueang_prenoise(tr,frequse,erpampind)] = phasereset_jz_steadybasefreq (time, fsample, frequse, foilim(1), foilim(2), halftimeind, halftimeind+.1*fsample, numones(erpampind));
          case 22
            [rawpr{frequse,erpampind}.trial{tr},trueang_prenoise(tr,frequse,erpampind),resetang(tr,frequse,erpampind)] = phasereset_jz_linphase (time, frequse, halftimeind, halftimeind+20, resetwin(erpampind) );
          case 4
        end
        rawpnpr{frequse,erpampind}.trial{tr}=rawpr{frequse,erpampind}.trial{tr}+10*pinknoise(length(time));
      end
      
      prc(:,frequse,erpampind)=resetang(:,frequse,erpampind)-trueang_prenoise(:,frequse,erpampind);
      trueang_prenoise(:,frequse,erpampind)=rad2deg(wrapToPi(trueang_prenoise(:,frequse,erpampind)));
      
      cfg=[];
      tlock=ft_timelockanalysis(cfg,raw{frequse,erpampind});
      tlockpr=ft_timelockanalysis(cfg,rawpr{frequse,erpampind});
      tlockpn=ft_timelockanalysis(cfg,rawpn{frequse,erpampind});
      tlockpnpr=ft_timelockanalysis(cfg,rawpnpr{frequse,erpampind});
      
      if 0
        figure;plot(time,[raw{frequse,erpampind}.trial{tr}; rawpr{frequse,erpampind}.trial{tr}]);
        figure;plot(time,[raw{frequse,erpampind}.trial{1}; rawpr{frequse,erpampind}.trial{1}]);
      end
      if 1
        figure;plot(time,[tlock.avg; tlockpr.avg; tlockpn.avg; tlockpnpr.avg]);axis([-inf inf -1 1])
      end
      
    end % erpampind
  end % frequse
  
  
  cfg=[];
  for erpampind=1:3
    tlockpnpr=ft_timelockanalysis(cfg,rawpnpr{10,erpampind});
    figure;plot(time,[tlockpnpr.avg]);axis([-inf inf -1 1])
  end
  
  
  startang=resetang-prc;
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 5a) Bandpass filter + Hilbert
  
  close all
  freqtest=[8 12];
  
  timwin=[4 2 1 .5 .25]; % duration in seconds; must be ordered from longest to shortest
  
  for tt=1:length(timwin)
    
    for frequse=8:12
      for erpampind=1:length(resetwin)
        
        cfg=[];
        cfg.latency=time([halftimeind-fsample*timwin(tt)/2 halftimeind+fsample*timwin(tt)/2]);
        rawuse=ft_selectdata(cfg,raw{frequse,erpampind});
        rawpnuse=ft_selectdata(cfg,rawpn{frequse,erpampind});
        rawpruse=ft_selectdata(cfg,rawpr{frequse,erpampind});
        rawpnpruse=ft_selectdata(cfg,rawpnpr{frequse,erpampind});
        halftimeinduse(tt)=dsearchn(rawuse.time{1}',time(halftimeind));
        
        for ff=1:size(freqtest,1)
          
          if 0
            % Bandpass filter: Butterworth
            rawpn_bpbut{3,2}=[];
            raw_bpbut{3,2}=[];
            rawpr_bpbut{3,2}=[];
            rawpnpr_bpbut{3,2}=[];
            cfg=[];
            cfg.bpfilter='yes';
            cfg.bpfreq=freqtest(ff,:);
            cfg.bpfilttype='but';
            cfg.plotfiltresp='yes';
            cfg.fouse=2:4;
            cfg.figind=0;
            cfg.plotflag=0;
            rawpn_bpbutr=filter4phase_estim8(cfg,rawpnuse);
            raw_bpbutr=filter4phase_estim8(cfg,rawuse);
            rawpr_bpbutr=filter4phase_estim8(cfg,rawpruse);
            rawpnpr_bpbutr=filter4phase_estim8(cfg,rawpnpruse);
            cfg.hilbert='complex';
            rawpn_bpbut=filter4phase_estim8(cfg,rawpnuse);
            raw_bpbut=filter4phase_estim8(cfg,rawuse);
            rawpr_bpbut=filter4phase_estim8(cfg,rawpruse);
            rawpnpr_bpbut=filter4phase_estim8(cfg,rawpnpruse);
          end
          
          % FIR
          rawpn_bpfir{3,2}=[];
          raw_bpfir{3,2}=[];
          rawpr_bpfir{3,2}=[];
          rawpnpr_bpfir{3,2}=[];
          % Bandpass filter: FIR (Matlab 'fir1')
          cfg=[];
          cfg.bpfilter='yes';
          cfg.bpfreq=freqtest(ff,:);
          cfg.bpfilttype='fir';
          cfg.plotfiltresp='yes';
          %       cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4]; % 3* is default
          %      cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
          %           cfg.fouse=[3*fsample/cfg.bpfreq(1)]; % 3* is default
          cfg.fouse=round([1 .8 .6].*length(rawpnuse.time{1})/3);  % see ft_preproc_bandpassfilter 'fir':   if N > floor( (size(dat,2) - 1) / 3);   N=floor(size(dat,2)/3) - 1;  end
          cfg.figind=10;
          cfg.plotflag=0;
          rawpn_bpfirr=filter4phase_estim8(cfg,rawpnuse);
          raw_bpfirr=filter4phase_estim8(cfg,rawuse);
          rawpr_bpfirr=filter4phase_estim8(cfg,rawpruse);
          rawpnpr_bpfirr=filter4phase_estim8(cfg,rawpnpruse);
          cfg.hilbert='complex';
          rawpn_bpfir=filter4phase_estim8(cfg,rawpnuse);
          raw_bpfir=filter4phase_estim8(cfg,rawuse);
          rawpr_bpfir=filter4phase_estim8(cfg,rawpruse);
          rawpnpr_bpfir=filter4phase_estim8(cfg,rawpnpruse);
          
          if 0
            % Bandpass filter: FIRWS
            cfg=[];
            cfg.bpfilter='yes';
            cfg.bpfreq=freqtest(ff,:);
            cfg.bpfilttype='firws';
            cfg.plotfiltresp='yes';
            cfg.fouse=2;
            cfg.figind=20;
            cfg.plotflag=0;
            rawpn_bpfirwsr=filter4phase_estim8(cfg,rawpnuse);
            raw_bpfirwsr=filter4phase_estim8(cfg,rawuse);
            rawpr_bpfirwsr=filter4phase_estim8(cfg,rawpruse);
            rawpnpr_bpfirwsr=filter4phase_estim8(cfg,rawpnpruse);
            cfg.hilbert='complex';
            rawpn_bpfirws=filter4phase_estim8(cfg,rawpnuse);
            raw_bpfirws=filter4phase_estim8(cfg,rawuse);
            rawpr_bpfirws=filter4phase_estim8(cfg,rawpruse);
            rawpnpr_bpfirws=filter4phase_estim8(cfg,rawpnpruse);
            if length(cfg.fouse)<3
              rawpn_bpfirws{3,2}=[];
              raw_bpfirws{3,2}=[];
              rawpr_bpfirws{3,2}=[];
              rawpnpr_bpfirws{3,2}=[];
            end
          end
          
          
          % reset these for every ff and tt
          if 0
            angdiffbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffprbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnprbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdifffirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffprfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnprfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            
            angkeepbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeepfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeepprbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeepprfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnprbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnprfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
          end
          
          angdifffir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdiffpnfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdiffprfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdiffpnprfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          
          angkeepfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeeppnfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeepprfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeeppnprfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          
          for ss=1:numtrials
            for fo=1:size(raw_bpfir,1)
              for fd=1:size(raw_bpfir,2)
                if 0
                  if ~isempty(raw_bpbut{fo,fd})
                    angkeepbut(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeepprbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpr_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnprbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnpr_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  end
                  if ~isempty(raw_bpfirws{fo,fd})
                    angkeepfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeepprfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpr_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnprfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnpr_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  end
                end
                if ~isempty(raw_bpfir{fo,fd})  % FIR
                  angkeepfir(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeeppnfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeepprfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpr_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeeppnprfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnpr_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                end
              end
            end % fo
            % don't need to index these angdiff* over ff as we don't save them
            if 0
              angdiffbutpre(:,:,ss)=anglediff(angkeepbut(:,:,ss),trueang_prenoise(ss),1);
              angdifffirwspre(:,:,ss)=anglediff(angkeepfirws(:,:,ss),trueang_prenoise(ss),1);
              angdiffpnbutpre(:,:,ss)=anglediff(angkeeppnbut(:,:,ss),trueang_prenoise(ss),1);
              angdiffpnfirwspre(:,:,ss)=anglediff(angkeeppnfirws(:,:,ss),trueang_prenoise(ss),1);
              angdiffprbutpre(:,:,ss)=anglediff(angkeepprbut(:,:,ss),trueang_prenoise(ss),1);
              angdiffprfirwspre(:,:,ss)=anglediff(angkeepprfirws(:,:,ss),trueang_prenoise(ss),1);
              angdiffpnprbutpre(:,:,ss)=anglediff(angkeeppnprbut(:,:,ss),trueang_prenoise(ss),1);
              angdiffpnprfirwspre(:,:,ss)=anglediff(angkeeppnprfirws(:,:,ss),trueang_prenoise(ss),1);
            end
            
            angdifffirpre(:,:,ss)=anglediff(angkeepfir(:,:,ss),trueang_prenoise(ss,frequse,erpampind),1);
            angdiffpnfirpre(:,:,ss)=anglediff(angkeeppnfir(:,:,ss),trueang_prenoise(ss,frequse,erpampind),1);
            angdiffprfirpre(:,:,ss)=anglediff(angkeepprfir(:,:,ss),trueang_prenoise(ss,frequse,erpampind),1);
            angdiffpnprfirpre(:,:,ss)=anglediff(angkeeppnprfir(:,:,ss),trueang_prenoise(ss,frequse,erpampind),1);
          end % ss
          
          if 0
            angbutpre(:,:,tt,ff)=rms(angdiffbutpre,3);
            angfirpre(:,:,tt,ff)=rms(angdifffirpre,3);
            angfirwspre(:,:,tt,ff)=rms(angdifffirwspre,3);
            angpnbutpre(:,:,tt,ff)=rms(angdiffpnbutpre,3);
            angpnfirpre(:,:,tt,ff)=rms(angdiffpnfirpre,3);
            angpnfirwspre(:,:,tt,ff)=rms(angdiffpnfirwspre,3);
            angprbutpre(:,:,tt,ff)=rms(angdiffprbutpre,3);
            angprfirpre(:,:,tt,ff)=rms(angdiffprfirpre,3);
            angprfirwspre(:,:,tt,ff)=rms(angdiffprfirwspre,3);
            angpnprbutpre(:,:,tt,ff)=rms(angdiffpnprbutpre,3);
            angpnprfirpre(:,:,tt,ff)=rms(angdiffpnprfirpre,3);
            angpnprfirwspre(:,:,tt,ff)=rms(angdiffpnprfirwspre,3);
          end
          
          if 0
            angfirprem(tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdifffirpre(1,1,:))),3) /100);
            angfirprep(tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdifffirpre(1,1,:))),3) /100));
            angpnfirprem(tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdiffpnfirpre(1,1,:))),3) /100);
            angpnfirprep(tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnfirpre(1,1,:))),3) /100));
            angerpfirprem(tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdiffprfirpre(1,1,:))),3) /100);
            angerpfirprep(tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffprfirpre(1,1,:))),3) /100));
            angpnerpfirprem(tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdiffpnprfirpre(1,1,:))),3) /100);
            angpnerpfirprep(tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnprfirpre(1,1,:))),3) /100));
          else
            angfirprem(:,tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdifffirpre(:,1,:))),3) /100);
            angfirprep(:,tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdifffirpre(:,1,:))),3) /100));
            angpnfirprem(:,tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdiffpnfirpre(:,1,:))),3) /100);
            angpnfirprep(:,tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnfirpre(:,1,:))),3) /100));
            angerpfirprem(:,tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdiffprfirpre(:,1,:))),3) /100);
            angerpfirprep(:,tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffprfirpre(:,1,:))),3) /100));
            angpnerpfirprem(:,tt,frequse,erpampind)=abs(sum(exp(i*deg2rad(angdiffpnprfirpre(:,1,:))),3) /100);
            angpnerpfirprep(:,tt,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnprfirpre(:,1,:))),3) /100));
          end
          
        end % ff
      end % erpampind
    end % frequse
  end % tt
  clear angdiff*
  
  if 0
    save(['D:\phase_estimation\anglermsrun5.mat'],'angbut*pre','angp*but*pre','angfir*pre','angp*fir*pre','frequse','trueang_prenoise')
  else
    try
      save(['D:\phase_estimation\anglermsrun5b_model' num2str(modeluse) '.mat'],'ang*fir*','-append')
    catch
      save(['D:\phase_estimation\anglermsrun5b_model' num2str(modeluse) '.mat'],'ang*fir*')
    end
  end
  
  load(['D:\phase_estimation\anglermsrun5b_model' num2str(modeluse) '.mat'])
  
  % "1 - circular-mean" is the same as circular-variance
  figure(1);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);                     bar(1-squeeze(angfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('No PN, No ERP');end
    subplot(4,length(timwin),tt+length(timwin));      bar(1-squeeze(angpnfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    if tt==1,ylabel('PN, No ERP');end
    subplot(4,length(timwin),tt+2*length(timwin));    bar(1-squeeze(angerpfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    if tt==1,ylabel('No PN, ERP');end
    subplot(4,length(timwin),tt+3*length(timwin));    bar(1-squeeze(angpnerpfirprem(tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    if tt==1,ylabel('PN, ERP');end
  end
  legend({'8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  set(get(1,'Children'),'xTickLabel',{'ERPamp 1' 'ERPamp 2' 'ERPamp 3'})
  
  figure(2);
  for tt=1:length(timwin),
    subplot(4,length(timwin),tt);                     bar(squeeze(angfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    title(['Time Window length ' num2str(timwin(tt))])
    if tt==1,ylabel('No PN, No ERP');end
    subplot(4,length(timwin),tt+length(timwin));      bar(squeeze(angpnfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('PN, No ERP');end
    subplot(4,length(timwin),tt+2*length(timwin));    bar(squeeze(angerpfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('No PN, ERP');end
    subplot(4,length(timwin),tt+3*length(timwin));    bar(squeeze(angpnerpfirprep(tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('PN, ERP');end
  end
  legend({'8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  set(get(1,'Children'),'xTickLabel',{'ERPamp 1' 'ERPamp 2' 'ERPamp 3'})
  
  figure(3);
  for tt=1:length(timwin),
    subplot(2,length(timwin),tt+0*length(timwin));    bar(1-squeeze(angpnerpfirprem(1,tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('Circ-Var');end
    subplot(2,length(timwin),tt+1*length(timwin));    bar(1-squeeze(angpnerpfirprep(1,tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('Mean-Angle');end
  end
  legend({'8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  %   set(get(3,'Children'),'XTickLabel',{'ERPamp 1' 'ERPamp 2' 'ERPamp 3'})
  
  figure(4);
  for tt=1:length(timwin),
    subplot(2,length(timwin),tt+0*length(timwin));    bar(1-squeeze(angpnerpfirprem(2,tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('Circ-Var');end
    subplot(2,length(timwin),tt+1*length(timwin));    bar(1-squeeze(angpnerpfirprep(2,tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('Mean-Angle');end
  end
  legend({'8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  %   set(get(4,'Children'),'XTickLabel',{'ERPamp 1' 'ERPamp 2' 'ERPamp 3'})
  
  figure(5);
  for tt=1:length(timwin),
    subplot(2,length(timwin),tt+0*length(timwin));    bar(1-squeeze(angpnerpfirprem(3,tt,8:12,:))')    ;axis([-inf inf 0 0.5])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('Circ-Var');end
    subplot(2,length(timwin),tt+1*length(timwin));    bar(1-squeeze(angpnerpfirprep(3,tt,8:12,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('Mean-Angle');end
  end
  legend({'8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  %   set(get(5,'Children'),'XTickLabel',{'ERPamp 1' 'ERPamp 2' 'ERPamp 3'})
  
  if 0
    %'Compared to trueang_prenoise')
    figure(1);
    for tt=1:length(timwin),
      subplot(3,length(timwin),tt);    bar([squeeze(angbutpre(:,2,tt,ff)), squeeze(angbutpre(:,2,tt,ff)), squeeze(angprbutpre(:,2,tt,ff)), squeeze(angpnprbutpre(:,2,tt,ff))]);axis([-inf inf 0 160])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Butterworth');end
      subplot(3,length(timwin),tt+5);    bar([squeeze(angfirpre(:,2,tt,ff)), squeeze(angfirpre(:,2,tt,ff)), squeeze(angprfirpre(:,2,tt,ff)), squeeze(angpnprfirpre(:,2,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR');end
      subplot(3,length(timwin),tt+10);    bar([squeeze(angfirwspre(:,2,tt,ff)), squeeze(angfirwspre(:,2,tt,ff)), squeeze(angprfirwspre(:,2,tt,ff)), squeeze(angpnprfirwspre(:,2,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIRWS');end
    end
    legend({'Orig signal' 'Noise, but no pr added', 'PR added', 'Noise and pr added'})
    set(get(1,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    %     %'Compared to trueang_postnoise'
    %   figure(2);
    %   for tt=1:length(timwin),
    %     subplot(2,length(timwin),tt);  bar([squeeze(angbutpost(:,2,tt,ff)), squeeze(angprbutpost(:,2,tt,ff))]);axis([-inf inf 0 160])
    %     title(['Time Window length ' num2str(timwin(tt))])
    %     if tt==1,ylabel('Butterworth');end
    %     subplot(2,length(timwin),tt+5);  bar([squeeze(angfirpost(:,2,tt,ff)), squeeze(angprfirpost(:,2,tt,ff))]);axis([-inf inf 0 160])
    %     if tt==1,ylabel('FIR');end
    %   end
    %   legend({'No pr added', 'pr added'})
    %   set(get(2,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    %'Compared to trueang_prenoise')
    figure(3);
    for tt=1:length(timwin),
      subplot(3,length(timwin),tt);    bar([squeeze(angbutpre(:,1,tt,ff)), squeeze(angpnbutpre(:,1,tt,ff)), squeeze(angprbutpre(:,1,tt,ff)), squeeze(angpnprbutpre(:,1,tt,ff))]);axis([-inf inf 0 160])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Butterworth');end
      subplot(3,length(timwin),tt+5);    bar([squeeze(angfirpre(:,1,tt,ff)), squeeze(angpnfirpre(:,1,tt,ff)), squeeze(angprfirpre(:,1,tt,ff)), squeeze(angpnprfirpre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR');end
      subplot(3,length(timwin),tt+10);    bar([squeeze(angfirwspre(:,1,tt,ff)), squeeze(angpnfirwspre(:,1,tt,ff)), squeeze(angprfirwspre(:,1,tt,ff)), squeeze(angpnprfirwspre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIRWS');end
    end
    legend({'Orig signal' 'Noise, but no PR added', 'PR added', 'Noise and PR added'})
    set(get(3,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    %'Compared to trueang_prenoise')
    figure(4);
    for tt=1:length(timwin),
      subplot(3,length(timwin),tt);    bar([squeeze(angbutpre(:,1,tt,ff)), squeeze(angpnbutpre(:,2,tt,ff)), squeeze(angprbutpre(:,2,tt,ff)), squeeze(angpnprbutpre(:,2,tt,ff))]);axis([-inf inf 0 160])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Butterworth 2pass');end
      subplot(3,length(timwin),tt+5);    bar([squeeze(angfirpre(:,1,tt,ff)), squeeze(angpnfirpre(:,1,tt,ff)), squeeze(angprfirpre(:,1,tt,ff)), squeeze(angpnprfirpre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR 1pass-0phase');end
      subplot(3,length(timwin),tt+10);    bar([squeeze(angfirwspre(:,1,tt,ff)), squeeze(angpnfirwspre(:,1,tt,ff)), squeeze(angprfirwspre(:,1,tt,ff)), squeeze(angpnprfirwspre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIRWS 1pass-0phase');end
    end
    legend({'Orig signal' 'Noise, but no PR added', 'PR added', 'Noise and PR added'})
    set(get(4,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  end
  
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 5b) FFT + taper
  
  cfg=[];
  cfg.method='mtmconvol';
  cfg.output='fourier';
  cfg.taper='hanning';
  cfg.foi=foilim(1):foilim(end);
  cfg.keeptrials='yes';
  
  t_ftimwin{1}=timwin(1)*ones(size(cfg.foi)); % full length of data
  t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 2 periods 4 Hz
  t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{6}=4./cfg.foi; %
  t_ftimwin{7}=3./cfg.foi; %
  t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
  t_ftimwin{9}=1./cfg.foi; % one period for each frequency
  t_ftimwin{10}=0.5./cfg.foi; %
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angfreq=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpr=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpn=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpnpr=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  
  angfreqspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqprspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpnspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpnprspec=nan(numtrials,length(toi),length(t_ftimwin));
  
  %   angfreq10=nan(numtrials,length(toi),length(t_ftimwin));
  %   angfreqpr10=nan(numtrials,length(toi),length(t_ftimwin));
  %   angfreqpn10=nan(numtrials,length(toi),length(t_ftimwin));
  %   angfreqpnpr10=nan(numtrials,length(toi),length(t_ftimwin));
  
  for frequse=8:12
    for erpampind=1:length(resetwin)
      for ss=1:numtrials
        cfg.trials=ss;
        for tt=1:length(toi)
          cfg.toi=toi{tt};
          for tf=1:length(t_ftimwin)
            cfg.t_ftimwin=t_ftimwin{tf};
            freq  =ft_freqanalysis(cfg, raw{frequse,erpampind});
            freqpr=ft_freqanalysis(cfg, rawpr{frequse,erpampind});
            freqpn  =ft_freqanalysis(cfg, rawpn{frequse,erpampind});
            freqpnpr=ft_freqanalysis(cfg, rawpnpr{frequse,erpampind});
            if tt==1
              %           angfreq(:,:,tt,pp,tf)=angle(squeeze(freq{tt,pp,tf}.fourierspctrm))/(2*pi)*360;
              angfreq(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freq.fourierspctrm))));
              angfreqpr(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqpr.fourierspctrm))));
              angfreqpn(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freqpn.fourierspctrm))));
              angfreqpnpr(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqpnpr.fourierspctrm))));
            elseif tt==2
              angfreq(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freq.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpr(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqpr.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpn(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freqpn.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpnpr(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqpnpr.fourierspctrm(1,1,:,:)))+pi  )));
            end
          end % tf
        end % tt
        if 0
          % given that we know exactly which freq was used for a given trial
          angfreqspec(ss,:,:)   =angfreq(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqprspec(ss,:,:)=angfreqpr(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqpnspec(ss,:,:)   =angfreqpn(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqpnprspec(ss,:,:)=angfreqpnpr(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          %     % use the freq that was closest to the pr (10.1Hz)
          %     angfreq10(ss,:,:)   =angfreq(ss,3,:,:);
          %     angfreqpr10(ss,:,:)=angfreqpr(ss,3,:,:);
          %     angfreqpn10(ss,:,:)   =angfreqpn(ss,3,:,:);
          %     angfreqpnpr10(ss,:,:)=angfreqpnpr(ss,3,:,:);
          angfreqspecdiff(ss,:,:,:)   =anglediff(angfreqspec(ss,:,:,:),   trueang_prenoise(ss),1);
          angfreqprspecdiff(ss,:,:,:)=anglediff(angfreqprspec(ss,:,:,:),trueang_prenoise(ss),1);
          angfreqpnspecdiff(ss,:,:,:)   =anglediff(angfreqpnspec(ss,:,:,:),   trueang_prenoise(ss),1);
          angfreqpnprspecdiff(ss,:,:,:)=anglediff(angfreqpnprspec(ss,:,:,:),trueang_prenoise(ss),1);
        end
        
        % compute differences
        angfreqdiff(ss,:,:,:)       =anglediff(angfreq(ss,:,:,:),       trueang_prenoise(ss,frequse,erpampind),1);
        angfreqprdiff(ss,:,:,:)    =anglediff(angfreqpr(ss,:,:,:),    trueang_prenoise(ss,frequse,erpampind),1);
        angfreqpndiff(ss,:,:,:)       =anglediff(angfreqpn(ss,:,:,:),       trueang_prenoise(ss,frequse,erpampind),1);
        angfreqpnprdiff(ss,:,:,:)    =anglediff(angfreqpnpr(ss,:,:,:),    trueang_prenoise(ss,frequse,erpampind),1);
        
        
        %     angfreq10diff(ss,:,:,:)     =anglediff(angfreq10(ss,:,:,:),     trueang_prenoise(ss,halftimeind),1);
        %     angfreqpr10diff(ss,:,:,:)  =anglediff(angfreqpr10(ss,:,:,:),  trueang_prenoise(ss,halftimeind),1);
        %     angfreqpn10diff(ss,:,:,:)     =anglediff(angfreqpn10(ss,:,:,:),     trueang_prenoise(ss,halftimeind),1);
        %     angfreqpnpr10diff(ss,:,:,:)  =anglediff(angfreqpnpr10(ss,:,:,:),  trueang_prenoise(ss,halftimeind),1);
      end
      if 0
        angfrms  =squeeze(rms(angfreqdiff,1));
        angfrrms =squeeze(rms(angfreqprdiff,1));
        angfprms  =squeeze(rms(angfreqpndiff,1));
        angfprrms =squeeze(rms(angfreqpnprdiff,1));
        
        angfsrms =squeeze(rms(angfreqspecdiff,1));
        angfrsrms=squeeze(rms(angfreqprspecdiff,1));
        angfpsrms =squeeze(rms(angfreqpnspecdiff,1));
        angfprsrms=squeeze(rms(angfreqpnprspecdiff,1));
        
        %   angftrms =squeeze(rms(angfreq10diff,1));
        %   angfetrms=squeeze(rms(angfreqerp10diff,1));
        %   angfptrms =squeeze(rms(angfreqpn10diff,1));
        %   angfpetrms=squeeze(rms(angfreqpnerp10diff,1));
      end
      
      angfmean(:,:,:,frequse,erpampind)=  abs(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials);
      angfemean(:,:,:,frequse,erpampind)=abs(sum(exp(i*deg2rad(angfreqprdiff)),1)/numtrials);
      angfpmean(:,:,:,frequse,erpampind)=abs(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials);
      angfpemean(:,:,:,frequse,erpampind)=abs(sum(exp(i*deg2rad(angfreqpnprdiff)),1)/numtrials);
      angfang(:,:,:,frequse,erpampind)=  rad2deg(angle(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials));
      angfeang(:,:,:,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqprdiff)),1)/numtrials));
      angfpang(:,:,:,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials));
      angfpeang(:,:,:,frequse,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpnprdiff)),1)/numtrials));
      
    end % erpampind
  end % frequse
  
  if 0
    save(['D:\phase_estimation\anglermsrun5.mat'],'ang*rms*','frequse','-append')
  else
    save(['D:\phase_estimation\anglermsrun5b_model' num2str(modeluse) '.mat'],'ang*mean','ang*ang','-append')
  end
  
  load(['D:\phase_estimation\anglermsrun5b_model' num2str(modeluse) '.mat'])
  
  for erpampind=1:3
    for tt=1:2
      figure(100+tt+(erpampind-1)*2);
      for tf=1:length(t_ftimwin),
        subplot(4,length(t_ftimwin),tf);                     bar(1-squeeze(angfmean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        title(['Time Window length ' num2str(t_ftimwin{tf}(3))])
        if tf==1,ylabel('No PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+length(t_ftimwin));      bar(1-squeeze(angfpmean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+2*length(t_ftimwin));    bar(1-squeeze(angfemean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('No PN, ERP');end
        subplot(4,length(t_ftimwin),tf+3*length(t_ftimwin));    bar(1-squeeze(angfpemean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('PN, ERP');end
      end
      legend({'cfg.foi 8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
      set(get(100+tt+(erpampind-1)*2,'Children'),'xTickLabel',{'Sim 8' '9' '10' '11' '12'})
    end
  end
  
  for erpampind=1:3
    for tt=1:2
      figure(110+tt+(erpampind-1)*2);
      for tf=1:length(t_ftimwin),
        subplot(4,length(t_ftimwin),tf);                     bar(squeeze(angfang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        title(['Time Window length ' num2str(t_ftimwin{tf}(1))])
        if tf==1,ylabel('No PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+length(t_ftimwin));      bar(squeeze(angfpang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+2*length(t_ftimwin));    bar(squeeze(angfeang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('No PN, ERP');end
        subplot(4,length(t_ftimwin),tf+3*length(t_ftimwin));    bar(squeeze(angfpeang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('PN, ERP');end
      end
      legend({'cfg.foi 8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
      set(get(110+tt+(erpampind-1)*2,'Children'),'xTickLabel',{'Sim 8' '9' '10' '11' '12'})
    end
  end
  
  
  
  for erpampind=1:3
    for tt=1:2
      figure(120+tt+(erpampind-1)*2);
      for tf=1:length(t_ftimwin),
        subplot(4,5,tf);       bar(1-squeeze(angfpemean(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf 0 0.5])
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['TW ' num2str(10*t_ftimwin{tf}(3)) ' per.'])
        end
        if tf==1,ylabel('Circ-Var');end
        subplot(4,5,tf+10);    bar(squeeze(angfpeang(:,tt,tf,8:12,erpampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('Mean Angle');end
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['TW ' num2str(10*t_ftimwin{tf}(3)) ' per.'])
        end
      end
      %       set(get(120+tt+(erpampind-1)*2,'Children'),'xTickLabel',{'Sim 8' '9' '10' '11' '12'})
      set(get(120+tt+(erpampind-1)*2,'Children'),'XTickLabel',{'Sim 8' '9' '10' '11' '12'})
      legend({'cfg.foi 8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
    end
  end
  
  
  if 0
    %   figure(100);
    %   for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfrms(:,:,tf));axis([-inf inf 0 160]);end
    %   figure(101);
    %   for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfrrms(:,:,tf));axis([-inf inf 0 160]);end
    %   figure(102);
    %   for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfprms(:,:,tf));axis([-inf inf 0 160]);end
    %   figure(103);
    %   for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfprrms(:,:,tf));axis([-inf inf 0 160]);end
    %
    %   figure(104);
    %   %   for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([angfsrms(:,tf),angfrsrms(:,tf),angftrms(:,tf),angfrtrms(:,tf), angfpsrms(:,tf),angfprsrms(:,tf),angfptrms(:,tf),angfprtrms(:,tf)]');axis([-inf inf 0 160]);end
    %   %   set(get(104,'Children'),'xTickLabel',{'spec', 'ERP spec', '10hz', 'ERP 10hz', 'PN spec', 'PN ERP spec', 'PN 10hz', 'PN ERP 10hz'})
    %   for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([angfsrms(:,tf),angfrsrms(:,tf), angfpsrms(:,tf),angfprsrms(:,tf)]');axis([-inf inf 0 160]);end
    %   set(get(104,'Children'),'xTickLabel',{'spec', 'ERP spec', 'PN spec', 'PN ERP spec'})
    
    figure(110);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 13],[angfrms(:,:,tf); angfsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(110,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(111);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 13],[angfrrms(:,:,tf); angfrsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(111,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(112);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 13],[angfprms(:,:,tf); angfpsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(112,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(113);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 13],[angfprrms(:,:,tf); angfprsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(113,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
  end
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 5c) Wavelet
  
  foi=foilim(1):foilim(end);
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  cfg.foi=foi;
  cfg.pad=time(end);
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angwave   =nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwavepr=nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwavepn   =nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwavepnpr=nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  
  angwavediff    = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwaveprdiff = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwavepndiff    = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwavepnprdiff = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  
  angwavespecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwaveprspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepnspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepnprspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  
  %   angwave10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  %   angwavepr10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  %   angwavepn10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  %   angwavepnpr10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  
  
  for tf=1:length(timwin) % more important than gwidth?
    scfg=[];
    scfg.latency=time([halftimeind-fsample*timwin(tf)/2 halftimeind+fsample*timwin(tf)/2]);
    rawuse=ft_selectdata(scfg,raw);
    rawpruse=ft_selectdata(scfg,rawpr);
    rawpnuse=ft_selectdata(scfg,rawpn);
    rawpnpruse=ft_selectdata(scfg,rawpnpr);
    halftimeinduse(tf)=dsearchn(rawuse.time{1}',time(halftimeind));
    
    for ss=1:numtrials
      cfg.trials=ss;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoi=ft_freqanalysis(cfg,rawuse);
          waveprfoi=ft_freqanalysis(cfg,rawpruse);
          wavepnfoi=ft_freqanalysis(cfg,rawpnuse);
          wavepnprfoi=ft_freqanalysis(cfg,rawpnpruse);
          indfoi1=unique(dsearchn(foi',wavefoi.freq'));
          if tt==1
            angwave(ss,indfoi1,tt,ww,1)   = rad2deg(wrapToPi(angle(squeeze(wavefoi.fourierspctrm))));
            angwavepr(ss,indfoi1,tt,ww,1)= rad2deg(wrapToPi(angle(squeeze(waveprfoi.fourierspctrm))));
            angwavepn(ss,indfoi1,tt,ww,1)   = rad2deg(wrapToPi(angle(squeeze(wavepnfoi.fourierspctrm))));
            angwavepnpr(ss,indfoi1,tt,ww,1)= rad2deg(wrapToPi(angle(squeeze(wavepnprfoi.fourierspctrm))));
          elseif tt==2
            angwave(ss,indfoi1,tt,ww,1)   = rad2deg(diag(wrapToPi(   angle(squeeze(wavefoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwavepr(ss,indfoi1,tt,ww,1)= rad2deg(diag(wrapToPi(   angle(squeeze(waveprfoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwavepn(ss,indfoi1,tt,ww,1)   = rad2deg(diag(wrapToPi(   angle(squeeze(wavepnfoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwavepnpr(ss,indfoi1,tt,ww,1)= rad2deg(diag(wrapToPi(   angle(squeeze(wavepnprfoi.fourierspctrm(1,1,:,:)))+pi   )));
          end
        end % ww
      end % tt
      % given that we know exactly which freq was used for a given trial
      angwavespec(ss,:,:)   =angwave(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwaveprspec(ss,:,:)=angwavepr(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwavepnspec(ss,:,:)   =angwavepn(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwavepnprspec(ss,:,:)=angwavepnpr(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      %       % use the freq that was closest to the pr (10.1Hz)
      %       angwave10(ss,:,:)   =angwave(ss,3,:,:,1);
      %       angwavepr10(ss,:,:)=angwavepr(ss,3,:,:,1);
      %       angwavepn10(ss,:,:)   =angwavepn(ss,3,:,:,1);
      %       angwavepnpr10(ss,:,:)=angwavepnpr(ss,3,:,:,1);
      
      % compute differences
      angwavediff(ss,:,:,:,tf,1)       =anglediff(angwave(ss,:,:,:,1),       trueang_prenoise(ss),1);
      angwaveprdiff(ss,:,:,:,tf,1)    =anglediff(angwavepr(ss,:,:,:,1),    trueang_prenoise(ss),1);
      angwavepndiff(ss,:,:,:,tf,1)       =anglediff(angwavepn(ss,:,:,:,1),       trueang_prenoise(ss),1);
      angwavepnprdiff(ss,:,:,:,tf,1)    =anglediff(angwavepnpr(ss,:,:,:,1),    trueang_prenoise(ss),1);
      
      angwavespecdiff(ss,:,:,tf,1)   =anglediff(angwavespec(ss,:,:),   trueang_prenoise(ss),1);
      angwaveprspecdiff(ss,:,:,tf,1)=anglediff(angwaveprspec(ss,:,:),trueang_prenoise(ss),1);
      angwavepnspecdiff(ss,:,:,tf,1)   =anglediff(angwavepnspec(ss,:,:),   trueang_prenoise(ss),1);
      angwavepnprspecdiff(ss,:,:,tf,1)=anglediff(angwavepnprspec(ss,:,:),trueang_prenoise(ss),1);
      
      %       angwave10diff(ss,:,:,tf,1)     =anglediff(angwave10(ss,:,:),     trueang_prenoise(ss,halftimeind),1);
      %       angwavepr10diff(ss,:,:,tf,1)  =anglediff(angwavepr10(ss,:,:),  trueang_prenoise(ss,halftimeind),1);
      %       angwavepn10diff(ss,:,:,tf,1)     =anglediff(angwavepn10(ss,:,:),     trueang_prenoise(ss,halftimeind),1);
      %       angwavepnpr10diff(ss,:,:,tf,1)  =anglediff(angwavepnpr10(ss,:,:),  trueang_prenoise(ss,halftimeind),1);
    end % ss
  end % tf
  
  
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  cfg.foilim=foilim;
  cfg.pad=time(end);
  
  for tf=1:length(timwin) % more important than gwidth
    scfg=[];
    scfg.latency=time([halftimeind-fsample*timwin(tf)/2 halftimeind+fsample*timwin(tf)/2]);
    rawuse=ft_selectdata(scfg,raw);
    rawpruse=ft_selectdata(scfg,rawpr);
    rawpnuse=ft_selectdata(scfg,rawpn);
    rawpnpruse=ft_selectdata(scfg,rawpnpr);
    halftimeinduse(tf)=dsearchn(rawuse.time{1}',time(halftimeind));
    
    for ss=1:numtrials
      cfg.trials=ss;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoilim=ft_freqanalysis(cfg,rawuse);
          waveprfoilim=ft_freqanalysis(cfg,rawpruse);
          wavepnfoilim=ft_freqanalysis(cfg,rawpnuse);
          wavepnprfoilim=ft_freqanalysis(cfg,rawpnpruse);
          % this produces .freq of length varying with timwin (spacing is 1/timwin)
          
          % one option: take the middle freq within the range and use that angle.
          indfoi2=dsearchn(foi', round(wavefoilim.freq(ceil(length(wavefoilim.freq)/2))));
          if tt==1
            angwave(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi(angle(squeeze(wavefoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwavepr(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(waveprfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwavepn(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi(angle(squeeze(wavepnfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwavepnpr(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(wavepnprfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
          elseif tt==2
            angwave(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi( angle(squeeze(wavefoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwavepr(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi( angle(squeeze(waveprfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwavepn(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi( angle(squeeze(wavepnfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwavepnpr(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi( angle(squeeze(wavepnprfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
          end
          % another option: take average over
          if tt==1
            angwave(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(angle(squeeze(wavefoilim.fourierspctrm)))));
            angwavepr(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(angle(squeeze(waveprfoilim.fourierspctrm)))));
            angwavepn(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(angle(squeeze(wavepnfoilim.fourierspctrm)))));
            angwavepnpr(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(angle(squeeze(wavepnprfoilim.fourierspctrm)))));
          elseif tt==2
            angwave(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavefoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwavepr(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(  diag(angle(squeeze(waveprfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwavepn(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavepnfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwavepnpr(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavepnprfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
          end
          
        end % ww
      end % tt
      
      % compute differences
      for dd=2:3
        angwavediff(ss,indfoi2,:,:,tf,dd)       =anglediff(angwave(ss,indfoi2,:,:,dd),       trueang_prenoise(ss),1);
        angwaveprdiff(ss,indfoi2,:,:,tf,dd)    =anglediff(angwavepr(ss,indfoi2,:,:,dd),    trueang_prenoise(ss),1);
        angwavepndiff(ss,indfoi2,:,:,tf,dd)       =anglediff(angwavepn(ss,indfoi2,:,:,dd),       trueang_prenoise(ss),1);
        angwavepnprdiff(ss,indfoi2,:,:,tf,dd)    =anglediff(angwavepnpr(ss,indfoi2,:,:,dd),    trueang_prenoise(ss),1);
      end
    end  % ss
  end % tf
  
  
  
  angwrms  =squeeze(rms(angwavediff,1));
  angwrrms =squeeze(rms(angwaveprdiff,1));
  angwprms  =squeeze(rms(angwavepndiff,1));
  angwprrms =squeeze(rms(angwavepnprdiff,1));
  
  angwsrms =squeeze(rms(angwavespecdiff,1));
  angwrsrms=squeeze(rms(angwaveprspecdiff,1));
  angwpsrms =squeeze(rms(angwavepnspecdiff,1));
  angwprsrms=squeeze(rms(angwavepnprspecdiff,1));
  
  %   angwtrms =squeeze(rms(angwave10diff,1));
  %   angwetrms=squeeze(rms(angwavepr10diff,1));
  %   angwptrms =squeeze(rms(angwavepn10diff,1));
  %   angwpetrms=squeeze(rms(angwavepnpr10diff,1));
  
  save(['D:\phase_estimation\anglermsrun5.mat'],'ang*rms*','frequse','-append')
  
  % 3rd dim index 4 means wavwidth=7
  for dd=1:3 % 3 diff ways of computing
    figure(200+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar(foi,angwrms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(201+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar(foi,angwrrms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(202+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar(foi,angwprms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(203+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar(foi,angwprrms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
  end
  figure(204);
  %   for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([angwsrms(:,tf),angwrsrms(:,tf),angwtrms(:,tf),angwrtrms(:,tf), angwpsrms(:,tf),angwprsrms(:,tf),angwptrms(:,tf),angwprtrms(:,tf)]');axis([-inf inf 0 160]);end
  %   set(get(204,'Children'),'xTickLabel',{'spec', 'ERP spec', '10hz', 'ERP 10hz', 'PN spec', 'PN ERP spec', 'PN 10hz', 'PN ERP 10hz'})
  for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([angwsrms(:,tf),angwrsrms(:,tf), angwpsrms(:,tf),angwprsrms(:,tf)]');axis([-inf inf 0 160]);end
  set(get(204,'Children'),'xTickLabel',{'spec', 'PR spec', 'PN spec', 'PN PR spec'})
  
  
  % 3rd dim index 4 means wavwidth=7
  for dd=1:3 % 3 diff ways of computing
    figure(250+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar([foi 13],[angwrms(:,:,4,tf,dd); angwsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(250+dd*10,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(251+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar([foi 13],[angwrrms(:,:,4,tf,dd); angwrsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(251+dd*10,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(252+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar([foi 13],[angwprms(:,:,4,tf,dd); angwpsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(252+dd*10,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
    figure(253+dd*10);
    for tf=1:length(timwin),subplot(1,5,tf);bar([foi 13],[angwprrms(:,:,4,tf,dd); angwprsrms(:,tf)']);axis([-inf inf 0 160]);end
    set(get(253+dd*10,'Children'),'xTickLabel',{'8', '9', '10', '11', '12', 'custom'})
  end
  
  
  
  
  %   end % ss
end

%% 6) Create new data, with Kc evoked at time point of interest, with additive model
% nearly same as (4)

if run6
  
  % params for whole simulation
  clearvars -except run* time wdr fsample wav*width
  close all;
  cd('D:\fieldtrip_svn\utilities\private');
  state=randomseed(13);
  cd(wdr);
  
  %   angrms=nan(4,2,4,20);     % FO, [But FIR], timwin (4s, 1s, 500ms, 250ms), 20 simulations
  %   angtrms=nan(4,2,4,20);    % FO, [But FIR], timwin (4s, 1s, 500ms, 250ms), 20 simulations
  %   angfrms=nan(4,2,3,10,20);  % freq, toi, pad, timwin, 20 simulations
  %   angftrms=nan(4,2,3,10,20);
  %   angwrms=nan(4,2,length(wavwidth),2,20);  % freq, toi, wavwidth, foi/foilim, 20 simulations
  %   angwtrms=nan(4,2,length(wavwidth),2,20);
  
  if 0
    foilim=[1 8]; % main frequencies during N2
    raw.label{1}='test';
    raw.dimord='chan_time';
    [raw.time{1:numtrials}]=deal(time);
    rawpn=raw;
    rawerp=raw;
    rawpnerp=raw;
  end
  
  clear raw
  numtrials=100;
  halftimeind=round(length(time)/2);
  kcamp=[3 6 9];
  
  
  for frequse=1:7
    for kcampind=1:length(kcamp)
      raw{frequse,kcampind}.label{1}='test';
      raw{frequse,kcampind}.dimord='chan_time';
      [raw{frequse,kcampind}.time{1:numtrials}]=deal(time);
      rawpn{frequse,kcampind}=raw{frequse,kcampind};
      rawerp{frequse,kcampind}=raw{frequse,kcampind};
      rawpnerp{frequse,kcampind}=raw{frequse,kcampind};
      
      for tr=1:numtrials % simulate 100 trials
        % Always use 'cos' to generate signals
        % Each trial has a base frequency, pink noise, and ERP added
        if 0
          frequse(tr) = diff(foilim)*rand(1)+foilim(1);
        end
        phaseshift(tr,frequse,kcampind) = wrapToPi(2*pi*rand(1));
        if 0
          trueang(tr,:) = frequse(tr)*2*pi*time+phaseshift(tr)*ones(size(time)); % random freq and phase
        else
          trueang(:,tr,frequse,kcampind) = frequse*2*pi*time+phaseshift(tr,frequse,kcampind)*ones(size(time)); % random phase
        end
        raw{frequse,kcampind}.trial{tr}=cos(trueang(:,tr,frequse,kcampind)');
        trueang_prenoise(:,tr,frequse,kcampind)=rad2deg(wrapToPi(trueang(:,tr,frequse,kcampind)));
        rawpn{frequse,kcampind}.trial{tr}=raw{frequse,kcampind}.trial{tr}+10*pinknoise(length(time));
        trueang_postnoise(:,tr,frequse,kcampind)=rad2deg(wrapToPi(angle(hilbert(rawpn{frequse,kcampind}.trial{tr})))); % after noise added.
        % add sinusoid (starting at amplitude zero) from time 0-100ms, at 1.25Hz frequency
        rawpnerp{frequse,kcampind}.trial{tr}=rawpn{frequse,kcampind}.trial{tr}+[zeros(1,halftimeind-1), kcamp(kcampind)*gausswin(801,1)'.*sin(1.25*2*pi*time(halftimeind:halftimeind+800)-1.25*2*pi*time(halftimeind)), zeros(1,halftimeind-801)];
        rawerp{frequse,kcampind}.trial{tr}  =raw{frequse,kcampind}.trial{tr}  +[zeros(1,halftimeind-1), kcamp(kcampind)*gausswin(801,1)'.*sin(1.25*2*pi*time(halftimeind:halftimeind+800)-1.25*2*pi*time(halftimeind)), zeros(1,halftimeind-801)];
      end
      cfg=[];
      tlock=ft_timelockanalysis(cfg,raw{frequse,kcampind});
      tlockpn=ft_timelockanalysis(cfg,rawpn{frequse,kcampind});
      tlockerp=ft_timelockanalysis(cfg,rawerp{frequse,kcampind});
      tlockpnerp=ft_timelockanalysis(cfg,rawpnerp{frequse,kcampind});
      
      if 0
        figure;plot(tlock.time,tlock.avg)
        hold on;plot(tlockpn.time,tlockpn.avg,'g');axis([-inf inf -3 3])
        hold on;plot(tlockerp.time,tlockerp.avg,'r');axis([-inf inf -3 3])
        hold on;plot(tlockerp.time,tlockpnerp.avg,'k');axis([-inf inf -3 3])
      end
      
      if 1
        figure;plot(tlock.time,raw{frequse,kcampind}.trial{1})
        hold on;plot(tlockpn.time,rawpn{frequse,kcampind}.trial{1},'g');axis([-inf inf -10 10])
        hold on;plot(tlockerp.time,rawerp{frequse,kcampind}.trial{1},'r');axis([-inf inf -10 10])
        hold on;plot(tlockerp.time,rawpnerp{frequse,kcampind}.trial{1},'k');axis([-inf inf -10 10])
      end
      
      if 0
        cfg=[];
        ft_databrowser(cfg,raw);
        ft_databrowser(cfg,rawerp);
      end
    end
  end
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 6a) Bandpass filter + Hilbert
  
  close all
  %     freqtest=[1 4; 4 7; 8 12];
  freqtest=[1 4];
  
  timwin=[4 2 1 .5]; % duration in seconds; must be ordered from longest to shortest
  
  for tt=1:length(timwin)
    
    for frequse=1:7
      for kcampind=1:length(kcamp)
        
        cfg=[];
        cfg.latency=time([halftimeind-fsample*timwin(tt)/2 halftimeind+fsample*timwin(tt)/2]);
        rawuse=ft_selectdata(cfg,raw{frequse,kcampind});
        rawpnuse=ft_selectdata(cfg,rawpn{frequse,kcampind});
        rawerpuse=ft_selectdata(cfg,rawerp{frequse,kcampind});
        rawpnerpuse=ft_selectdata(cfg,rawpnerp{frequse,kcampind});
        halftimeinduse(tt)=dsearchn(rawuse.time{1}',time(halftimeind));
        
        for ff=1:size(freqtest,1)
          
          if 0
            % Bandpass filter: Butterworth
            rawpn_bpbut{3,2}=[];
            raw_bpbut{3,2}=[];
            rawerp_bpbut{3,2}=[];
            rawpnerp_bpbut{3,2}=[];
            cfg=[];
            cfg.bpfilter='yes';
            cfg.bpfreq=freqtest(ff,:);
            cfg.bpfilttype='but';
            cfg.plotfiltresp='yes';
            cfg.fouse=2:4;
            cfg.figind=0;
            cfg.plotflag=0;
            rawpn_bpbutr=filter4phase_estim8(cfg,rawpnuse);
            raw_bpbutr=filter4phase_estim8(cfg,rawuse);
            rawerp_bpbutr=filter4phase_estim8(cfg,rawerpuse);
            rawpnerp_bpbutr=filter4phase_estim8(cfg,rawpnerpuse);
            cfg.hilbert='complex';
            rawpn_bpbut=filter4phase_estim8(cfg,rawpnuse);
            raw_bpbut=filter4phase_estim8(cfg,rawuse);
            rawerp_bpbut=filter4phase_estim8(cfg,rawerpuse);
            rawpnerp_bpbut=filter4phase_estim8(cfg,rawpnerpuse);
          end
          
          % Bandpass filter: FIR (Matlab 'fir1')
          rawpn_bpfir{3,2}=[];
          raw_bpfir{3,2}=[];
          rawerp_bpfir{3,2}=[];
          rawpnerp_bpfir{3,2}=[];
          cfg=[];
          cfg.bpfilter='yes';
          cfg.bpfreq=freqtest(ff,:);
          cfg.bpfilttype='fir';
          cfg.plotfiltresp='yes';
          %       cfg.fouse=[2*fsample/4 3*fsample/4 4*fsample/4]; % 3* is default
          %           cfg.fouse=[2*fsample/cfg.bpfreq(1) 3*fsample/cfg.bpfreq(1) 4*fsample/cfg.bpfreq(1)]; % 3* is default
          %           cfg.fouse=[3*fsample/cfg.bpfreq(1)]; % 3* is default
          cfg.fouse=round([1 .8 .6].*length(rawpnuse.time{1})/3);  % see ft_preproc_bandpassfilter 'fir':   if N > floor( (size(dat,2) - 1) / 3);   N=floor(size(dat,2)/3) - 1;  end
          cfg.figind=10;
          cfg.plotflag=0;
          rawpn_bpfirr=filter4phase_estim8(cfg,rawpnuse);
          raw_bpfirr=filter4phase_estim8(cfg,rawuse);
          rawerp_bpfirr=filter4phase_estim8(cfg,rawerpuse);
          rawpnerp_bpfirr=filter4phase_estim8(cfg,rawpnerpuse);
          cfg.hilbert='complex';
          rawpn_bpfir=filter4phase_estim8(cfg,rawpnuse);
          raw_bpfir=filter4phase_estim8(cfg,rawuse);
          rawerp_bpfir=filter4phase_estim8(cfg,rawerpuse);
          rawpnerp_bpfir=filter4phase_estim8(cfg,rawpnerpuse);
          
          if 0
            % Bandpass filter: FIRWS
            cfg=[];
            cfg.bpfilter='yes';
            cfg.bpfreq=freqtest(ff,:);
            cfg.bpfilttype='firws';
            cfg.plotfiltresp='yes';
            cfg.fouse=2;
            cfg.figind=20;
            cfg.plotflag=0;
            rawpn_bpfirwsr=filter4phase_estim8(cfg,rawpnuse);
            raw_bpfirwsr=filter4phase_estim8(cfg,rawuse);
            rawerp_bpfirwsr=filter4phase_estim8(cfg,rawerpuse);
            rawpnerp_bpfirwsr=filter4phase_estim8(cfg,rawpnerpuse);
            cfg.hilbert='complex';
            rawpn_bpfirws=filter4phase_estim8(cfg,rawpnuse);
            raw_bpfirws=filter4phase_estim8(cfg,rawuse);
            rawerp_bpfirws=filter4phase_estim8(cfg,rawerpuse);
            rawpnerp_bpfirws=filter4phase_estim8(cfg,rawpnerpuse);
            if length(cfg.fouse)<3
              rawpn_bpfirws{3,2}=[];
              raw_bpfirws{3,2}=[];
              rawerp_bpfirws{3,2}=[];
              rawpnerp_bpfirws{3,2}=[];
            end
          end
          
          % reset these for every ff and tt
          if 0
            angdiffbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdifferpbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnerpbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            
            angdifffirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdifferpfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angdiffpnerpfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeepbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeperpbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnerpbut=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeepfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeperpfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
            angkeeppnerpfirws=nan(size(raw_bpbut,1),size(raw_bpbut,2),numtrials);
          end
          
          angdifffir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdiffpnfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdifferpfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angdiffpnerpfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeepfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeeppnfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeeperpfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          angkeeppnerpfir=nan(size(raw_bpfir,1),size(raw_bpfir,2),numtrials);
          
          for ss=1:numtrials
            for fo=1:size(raw_bpfir,1)
              for fd=1:size(raw_bpfir,2)
                if 0
                  if ~isempty(raw_bpbut{fo,fd})
                    angkeepbut(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeperpbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawerp_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnerpbut(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnerp_bpbut{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  end
                  if ~isempty(raw_bpfirws{fo,fd})
                    angkeepfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeperpfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawerp_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                    angkeeppnerpfirws(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnerp_bpfirws{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  end
                end
                if ~isempty(raw_bpfir{fo,fd})
                  angkeepfir(fo,fd,ss)=rad2deg(wrapToPi(angle(raw_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeeppnfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpn_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeeperpfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawerp_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                  angkeeppnerpfir(fo,fd,ss)=rad2deg(wrapToPi(angle(rawpnerp_bpfir{fo,fd}.trial{ss}(halftimeinduse(tt)))));
                end
              end
            end % fo
            % don't need to index these angdiff* over ff as we don't save them
            if 0
              angdiffbutpre(:,:,ss)=anglediff(angkeepbut(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdiffpnbutpre(:,:,ss)=anglediff(angkeeppnbut(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdifferpbutpre(:,:,ss)=anglediff(angkeeperpbut(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdiffpnerpbutpre(:,:,ss)=anglediff(angkeeppnerpbut(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdifffirwspre(:,:,ss)=anglediff(angkeepfirws(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdiffpnfirwspre(:,:,ss)=anglediff(angkeeppnfirws(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdifferpfirwspre(:,:,ss)=anglediff(angkeeperpfirws(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              angdiffpnerpfirwspre(:,:,ss)=anglediff(angkeeppnerpfirws(:,:,ss),trueang_prenoise(ss,halftimeind),1);
              
              angdiffbutpost(:,:,ss)=anglediff(angkeepbut(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdiffpnbutpost(:,:,ss)=anglediff(angkeeppnbut(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdifferpbutpost(:,:,ss)=anglediff(angkeeperpbut(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdiffpnerpbutpost(:,:,ss)=anglediff(angkeeppnerpbut(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdifffirwspost(:,:,ss)=anglediff(angkeepfirws(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdiffpnfirwspost(:,:,ss)=anglediff(angkeeppnfirws(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdifferpfirwspost(:,:,ss)=anglediff(angkeeperpfirws(:,:,ss),trueang_postnoise(ss,halftimeind),1);
              angdiffpnerpfirwspost(:,:,ss)=anglediff(angkeeppnerpfirws(:,:,ss),trueang_postnoise(ss,halftimeind),1);
            end
            
            angdifffirpre(:,:,ss)=anglediff(angkeepfir(:,:,ss),trueang_prenoise(halftimeind,ss,frequse,kcampind),1);
            angdiffpnfirpre(:,:,ss)=anglediff(angkeeppnfir(:,:,ss),trueang_prenoise(halftimeind,ss,frequse,kcampind),1);
            angdifferpfirpre(:,:,ss)=anglediff(angkeeperpfir(:,:,ss),trueang_prenoise(halftimeind,ss,frequse,kcampind),1);
            angdiffpnerpfirpre(:,:,ss)=anglediff(angkeeppnerpfir(:,:,ss),trueang_prenoise(halftimeind,ss,frequse,kcampind),1);
            
            angdifffirpost(:,:,ss)=anglediff(angkeepfir(:,:,ss),trueang_postnoise(halftimeind,ss,frequse,kcampind),1);
            angdiffpnfirpost(:,:,ss)=anglediff(angkeeppnfir(:,:,ss),trueang_postnoise(halftimeind,ss,frequse,kcampind),1);
            angdifferpfirpost(:,:,ss)=anglediff(angkeeperpfir(:,:,ss),trueang_postnoise(halftimeind,ss,frequse,kcampind),1);
            angdiffpnerpfirpost(:,:,ss)=anglediff(angkeeppnerpfir(:,:,ss),trueang_postnoise(halftimeind,ss,frequse,kcampind),1);
          end % ss
          if 0
            angbutpre(:,:,tt,ff)=rms(angdiffbutpre,3);
            angfirpre(:,:,tt,ff)=rms(angdifffirpre,3);
            angfirwspre(:,:,tt,ff)=rms(angdifffirwspre,3);
            angpnbutpre(:,:,tt,ff)=rms(angdiffpnbutpre,3);
            angpnfirpre(:,:,tt,ff)=rms(angdiffpnfirpre,3);
            angpnfirwspre(:,:,tt,ff)=rms(angdiffpnfirwspre,3);
            angerpbutpre(:,:,tt,ff)=rms(angdifferpbutpre,3);
            angerpfirpre(:,:,tt,ff)=rms(angdifferpfirpre,3);
            angerpfirwspre(:,:,tt,ff)=rms(angdifferpfirwspre,3);
            angpnerpbutpre(:,:,tt,ff)=rms(angdiffpnerpbutpre,3);
            angpnerpfirpre(:,:,tt,ff)=rms(angdiffpnerpfirpre,3);
            angpnerpfirwspre(:,:,tt,ff)=rms(angdiffpnerpfirwspre,3);
            
            angbutpost(:,:,tt,ff)=rms(angdiffbutpost,3);
            angfirpost(:,:,tt,ff)=rms(angdifffirpost,3);
            angfirwspost(:,:,tt,ff)=rms(angdifffirwspost,3);
            angpnbutpost(:,:,tt,ff)=rms(angdiffpnbutpost,3);
            angpnfirpost(:,:,tt,ff)=rms(angdiffpnfirpost,3);
            angpnfirwspost(:,:,tt,ff)=rms(angdiffpnfirwspost,3);
            angerpbutpost(:,:,tt,ff)=rms(angdifferpbutpost,3);
            angerpfirpost(:,:,tt,ff)=rms(angdifferpfirpost,3);
            angerpfirwspost(:,:,tt,ff)=rms(angdifferpfirwspost,3);
            angpnerpbutpost(:,:,tt,ff)=rms(angdiffpnerpbutpost,3);
            angpnerpfirpost(:,:,tt,ff)=rms(angdiffpnerpfirpost,3);
            angpnerpfirwspost(:,:,tt,ff)=rms(angdiffpnerpfirwspost,3);
          end
          
          %           % Note: not saving over 'ff' (but this is just alpha band anyway)
          %           angfirprem(tt,frequse,kcampind)=abs(sum(exp(i*deg2rad(angdifffirpre(1,1,:))),3) /100);
          %           angfirprep(tt,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angdifffirpre(1,1,:))),3) /100));
          %           angpnfirprem(tt,frequse,kcampind)=abs(sum(exp(i*deg2rad(angdiffpnfirpre(1,1,:))),3) /100);
          %           angpnfirprep(tt,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnfirpre(1,1,:))),3) /100));
          %           angerpfirprem(tt,frequse,kcampind)=abs(sum(exp(i*deg2rad(angdifferpfirpre(1,1,:))),3) /100);
          %           angerpfirprep(tt,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angdifferpfirpre(1,1,:))),3) /100));
          %           angpnerpfirprem(tt,frequse,kcampind)=abs(sum(exp(i*deg2rad(angdiffpnerpfirpre(1,1,:))),3) /100);
          %           angpnerpfirprep(tt,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnerpfirpre(1,1,:))),3) /100));
          
          % Note: not saving over 'ff' ;  % here we are saving all FO options
          angfirprem(:,tt,frequse,kcampind)=abs(sum(exp(i*deg2rad(angdifffirpre(:,1,:))),3) /100);
          angfirprep(:,tt,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angdifffirpre(:,1,:))),3) /100));
          angpnfirprem(:,tt,frequse,kcampind)=abs(sum(exp(i*deg2rad(angdiffpnfirpre(:,1,:))),3) /100);
          angpnfirprep(:,tt,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnfirpre(:,1,:))),3) /100));
          angerpfirprem(:,tt,frequse,kcampind)=abs(sum(exp(i*deg2rad(angdifferpfirpre(:,1,:))),3) /100);
          angerpfirprep(:,tt,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angdifferpfirpre(:,1,:))),3) /100));
          angpnerpfirprem(:,tt,frequse,kcampind)=abs(sum(exp(i*deg2rad(angdiffpnerpfirpre(:,1,:))),3) /100);
          angpnerpfirprep(:,tt,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angdiffpnerpfirpre(:,1,:))),3) /100));
          
        end % ff
      end % kcampind
    end % frequse
  end % tt
  clear angdiff*
  
  if 0
    save(['D:\phase_estimation\anglermsrun6.mat'],'ang*but*pre','ang*fir*pre','frequse','trueang_prenoise','-append');
  else
    try
      save(['D:\phase_estimation\anglermsrun6b.mat'],'ang*fir*','trueang_prenoise','-append');
    catch
      save(['D:\phase_estimation\anglermsrun6b.mat'],'ang*fir*','trueang_prenoise');
    end
  end
  
  load(['D:\phase_estimation\anglermsrun6b.mat']);
  
  
  %   figure(1);
  %   for tt=1:length(timwin),
  %     subplot(4,length(timwin),tt);                     bar(1-squeeze(angfirprem(tt,:,:))')    ;axis([-inf inf 0 0.5])
  %     title(['Time Window length ' num2str(timwin(tt))])
  %     if tt==1,ylabel('No PN, No Kc');end
  %     subplot(4,length(timwin),tt+length(timwin));      bar(1-squeeze(angpnfirprem(tt,:,:))')    ;axis([-inf inf 0 0.5])
  %     if tt==1,ylabel('PN, No Kc');end
  %     subplot(4,length(timwin),tt+2*length(timwin));    bar(1-squeeze(angerpfirprem(tt,:,:))')    ;axis([-inf inf 0 0.5])
  %     if tt==1,ylabel('No PN, Kc');end
  %     subplot(4,length(timwin),tt+3*length(timwin));    bar(1-squeeze(angpnerpfirprem(tt,:,:))')    ;axis([-inf inf 0 0.5])
  %     if tt==1,ylabel('PN, Kc');end
  %   end
  %   legend({'1 Hz' '2 Hz' '3 Hz' '4 Hz', '5 Hz', '6 Hz' '7 Hz'})
  %   set(get(1,'Children'),'xTickLabel',{'KCamp 1' 'KCamp 2' 'KCamp 3'})
  %
  %
  %   figure(2);
  %   for tt=1:length(timwin),
  %     subplot(4,length(timwin),tt);                     bar(squeeze(angfirprep(tt,:,:))')    ;axis([-inf inf -20 20])
  %     title(['Time Window length ' num2str(timwin(tt))])
  %     if tt==1,ylabel('No PN, No Kc');end
  %     subplot(4,length(timwin),tt+length(timwin));      bar(squeeze(angpnfirprep(tt,:,:))')    ;axis([-inf inf -20 20])
  %     if tt==1,ylabel('PN, No Kc');end
  %     subplot(4,length(timwin),tt+2*length(timwin));    bar(squeeze(angerpfirprep(tt,:,:))')    ;axis([-inf inf -20 20])
  %     if tt==1,ylabel('No PN, Kc');end
  %     subplot(4,length(timwin),tt+3*length(timwin));    bar(squeeze(angpnerpfirprep(tt,:,:))')    ;axis([-inf inf -20 20])
  %     if tt==1,ylabel('PN, Kc');end
  %   end
  %   legend({'1 Hz' '2 Hz' '3 Hz' '4 Hz', '5 Hz', '6 Hz' '7 Hz'})
  %   set(get(2,'Children'),'xTickLabel',{'KCamp 1' 'KCamp 2' 'KCamp 3'})
  
  
  figure(3);
  for tt=1:length(timwin),
    subplot(2,length(timwin),tt+0*length(timwin));    bar(1-squeeze(angpnerpfirprem(1,tt,1:4,:))')    ;axis([-inf inf 0 0.5])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('Circ-Var');end
    subplot(2,length(timwin),tt+1*length(timwin));    bar(squeeze(angpnerpfirprep(1,tt,1:4,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('Mean-Angle');end
  end
  legend({'1 Hz' '2 Hz' '3 Hz' '4 Hz'})
  %   set(get(3,'Children'),'xTickLabel',{'KCamp 1' 'KCamp 2' 'KCamp 3'})
  
  figure(4);
  for tt=1:length(timwin),
    subplot(2,length(timwin),tt+0*length(timwin));    bar(1-squeeze(angpnerpfirprem(2,tt,1:4,:))')    ;axis([-inf inf 0 0.5])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('Circ-Var');end
    subplot(2,length(timwin),tt+1*length(timwin));    bar(squeeze(angpnerpfirprep(2,tt,1:4,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('Mean-Angle');end
  end
  legend({'1 Hz' '2 Hz' '3 Hz' '4 Hz'})
  %   set(get(4,'Children'),'xTickLabel',{'KCamp 1' 'KCamp 2' 'KCamp 3'})
  
  figure(5);
  for tt=1:length(timwin),
    subplot(2,length(timwin),tt+0*length(timwin));    bar(1-squeeze(angpnerpfirprem(3,tt,1:4,:))')    ;axis([-inf inf 0 0.5])
    title(['TimWin ' num2str(timwin(tt)) ' s'])
    if tt==1,ylabel('Circ-Var');end
    subplot(2,length(timwin),tt+1*length(timwin));    bar(squeeze(angpnerpfirprep(3,tt,1:4,:))')    ;axis([-inf inf -20 20])
    if tt==1,ylabel('Mean-Angle');end
  end
  legend({'1 Hz' '2 Hz' '3 Hz' '4 Hz'})
  %   set(get(5,'Children'),'xTickLabel',{'KCamp 1' 'KCamp 2' 'KCamp 3'})
  
  
  
  if 0
    %'Compared to trueang_prenoise')
    figure(1);
    for tt=1:length(timwin),
      subplot(3,length(timwin),tt);    bar([squeeze(angbutpre(:,2,tt,ff)), squeeze(angbutpre(:,2,tt,ff)), squeeze(angerpbutpre(:,2,tt,ff)), squeeze(angpnerpbutpre(:,2,tt,ff))]);axis([-inf inf 0 160])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Butterworth');end
      subplot(3,length(timwin),tt+length(timwin));    bar([squeeze(angfirpre(:,2,tt,ff)), squeeze(angfirpre(:,2,tt,ff)), squeeze(angerpfirpre(:,2,tt,ff)), squeeze(angpnerpfirpre(:,2,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR');end
      subplot(3,length(timwin),tt+2*length(timwin));    bar([squeeze(angfirwspre(:,2,tt,ff)), squeeze(angfirwspre(:,2,tt,ff)), squeeze(angerpfirwspre(:,2,tt,ff)), squeeze(angpnerpfirwspre(:,2,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIRWS');end
    end
    legend({'Orig signal' 'Noise, but no Kc added', 'Kc added', 'Noise and Kc added'})
    set(get(1,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    %     %'Compared to trueang_postnoise'
    %   figure(2);
    %   for tt=1:length(timwin),
    %     subplot(2,length(timwin),tt);  bar([squeeze(angbutpost(:,2,tt,ff)), squeeze(angerpbutpost(:,2,tt,ff))]);axis([-inf inf 0 160])
    %     title(['Time Window length ' num2str(timwin(tt))])
    %     if tt==1,ylabel('Butterworth');end
    %     subplot(2,length(timwin),tt+5);  bar([squeeze(angfirpost(:,2,tt,ff)), squeeze(angerpfirpost(:,2,tt,ff))]);axis([-inf inf 0 160])
    %     if tt==1,ylabel('FIR');end
    %   end
    %   legend({'No ERP added', 'ERP added'})
    %   set(get(2,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
    
    %'Compared to trueang_prenoise')
    figure(3);
    for tt=1:length(timwin),
      subplot(3,length(timwin),tt);    bar([squeeze(angbutpre(:,1,tt,ff)), squeeze(angpnbutpre(:,1,tt,ff)), squeeze(angerpbutpre(:,1,tt,ff)), squeeze(angpnerpbutpre(:,1,tt,ff))]);axis([-inf inf 0 160])
      title(['Time Window length ' num2str(timwin(tt))])
      if tt==1,ylabel('Butterworth');end
      subplot(3,length(timwin),tt+length(timwin));    bar([squeeze(angfirpre(:,1,tt,ff)), squeeze(angpnfirpre(:,1,tt,ff)), squeeze(angerpfirpre(:,1,tt,ff)), squeeze(angpnerpfirpre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR');end
      subplot(3,length(timwin),tt+2*length(timwin));    bar([squeeze(angfirwspre(:,1,tt,ff)), squeeze(angpnfirwspre(:,1,tt,ff)), squeeze(angerpfirwspre(:,1,tt,ff)), squeeze(angpnerpfirwspre(:,1,tt,ff))]);axis([-inf inf 0 160])
      if tt==1,ylabel('FIR');end
    end
    legend({'Orig signal' 'Noise, but no Kc added', 'Kc added', 'Noise and Kc added'})
    set(get(3,'Children'),'xTickLabel',{'filt-ord 1' 'fo 2' 'fo 3' 'fo 4'})
  end
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 6b) FFT + taper
  
  cfg=[];
  cfg.method='mtmconvol';
  cfg.output='fourier';
  cfg.taper='hanning';
  cfg.foi=freqtest(1):freqtest(end);
  cfg.keeptrials='yes';
  
  t_ftimwin{1}=timwin(1)*ones(size(cfg.foi)); % full length of data
  t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 2 periods 4 Hz
  t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  %   t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{5}=nan*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{6}=4./cfg.foi; %
  t_ftimwin{7}=3./cfg.foi; %
  t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
  t_ftimwin{9}=1./cfg.foi; % one period for each frequency
  t_ftimwin{10}=0.5./cfg.foi; %
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angfreq=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqerp=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpn=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpnerp=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  
  angfreqspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqerpspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpnspec=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpnerpspec=nan(numtrials,length(toi),length(t_ftimwin));
  
  angfreq10=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqerp10=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpn10=nan(numtrials,length(toi),length(t_ftimwin));
  angfreqpnerp10=nan(numtrials,length(toi),length(t_ftimwin));
  
  for frequse=1:7
    for kcampind=1:length(kcamp)
      for ss=1:numtrials
        cfg.trials=ss;
        for tt=1:length(toi)
          cfg.toi=toi{tt};
          for tf=setdiff(1:length(t_ftimwin),5)
            cfg.t_ftimwin=t_ftimwin{tf};
            freq  =ft_freqanalysis(cfg, raw{frequse,kcampind});
            freqerp=ft_freqanalysis(cfg, rawerp{frequse,kcampind});
            freqpn  =ft_freqanalysis(cfg, rawpn{frequse,kcampind});
            freqpnerp=ft_freqanalysis(cfg, rawpnerp{frequse,kcampind});
            if tt==1
              %           angfreq(:,:,tt,pp,tf)=angle(squeeze(freq{tt,pp,tf}.fourierspctrm))/(2*pi)*360;
              angfreq(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freq.fourierspctrm))));
              angfreqerp(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqerp.fourierspctrm))));
              angfreqpn(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freqpn.fourierspctrm))));
              angfreqpnerp(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqpnerp.fourierspctrm))));
            elseif tt==2
              angfreq(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freq.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqerp(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqerp.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpn(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freqpn.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpnerp(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqpnerp.fourierspctrm(1,1,:,:)))+pi  )));
            end
          end % tf
        end % tt
        if 0
          % given that we know exactly which freq was used for a given trial
          angfreqspec(ss,:,:)   =angfreq(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqerpspec(ss,:,:)=angfreqerp(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqpnspec(ss,:,:)   =angfreqpn(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          angfreqpnerpspec(ss,:,:)=angfreqpnerp(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:);
          % use the freq that was closest to the Kc (1.25Hz)
          angfreq10(ss,:,:)   =angfreq(ss,1,:,:);
          angfreqerp10(ss,:,:)=angfreqerp(ss,1,:,:);
          angfreqpn10(ss,:,:)   =angfreqpn(ss,1,:,:);
          angfreqpnerp10(ss,:,:)=angfreqpnerp(ss,1,:,:);
          
          angfreqspecdiff(ss,:,:,:)   =anglediff(angfreqspec(ss,:,:,:),   trueang_prenoise(ss,halftimeind),1);
          angfreqerpspecdiff(ss,:,:,:)=anglediff(angfreqerpspec(ss,:,:,:),trueang_prenoise(ss,halftimeind),1);
          angfreqpnspecdiff(ss,:,:,:)   =anglediff(angfreqpnspec(ss,:,:,:),   trueang_prenoise(ss,halftimeind),1);
          angfreqpnerpspecdiff(ss,:,:,:)=anglediff(angfreqpnerpspec(ss,:,:,:),trueang_prenoise(ss,halftimeind),1);
          
          angfreq10diff(ss,:,:,:)     =anglediff(angfreq10(ss,:,:,:),     trueang_prenoise(ss,halftimeind),1);
          angfreqerp10diff(ss,:,:,:)  =anglediff(angfreqerp10(ss,:,:,:),  trueang_prenoise(ss,halftimeind),1);
          angfreqpn10diff(ss,:,:,:)     =anglediff(angfreqpn10(ss,:,:,:),     trueang_prenoise(ss,halftimeind),1);
          angfreqpnerp10diff(ss,:,:,:)  =anglediff(angfreqpnerp10(ss,:,:,:),  trueang_prenoise(ss,halftimeind),1);
        end
        
        % compute differences
        angfreqdiff(ss,:,:,:)       =anglediff(angfreq(ss,:,:,:),       trueang_prenoise(halftimeind,ss,frequse,kcampind),1);
        angfreqerpdiff(ss,:,:,:)    =anglediff(angfreqerp(ss,:,:,:),    trueang_prenoise(halftimeind,ss,frequse,kcampind),1);
        angfreqpndiff(ss,:,:,:)       =anglediff(angfreqpn(ss,:,:,:),       trueang_prenoise(halftimeind,ss,frequse,kcampind),1);
        angfreqpnerpdiff(ss,:,:,:)    =anglediff(angfreqpnerp(ss,:,:,:),    trueang_prenoise(halftimeind,ss,frequse,kcampind),1);
        
      end
      if 0
        angfrms  =squeeze(rms(angfreqdiff,1));
        angferms =squeeze(rms(angfreqerpdiff,1));
        angfprms  =squeeze(rms(angfreqpndiff,1));
        angfperms =squeeze(rms(angfreqpnerpdiff,1));
        
        angfsrms =squeeze(rms(angfreqspecdiff,1));
        angfesrms=squeeze(rms(angfreqerpspecdiff,1));
        angfpsrms =squeeze(rms(angfreqpnspecdiff,1));
        angfpesrms=squeeze(rms(angfreqpnerpspecdiff,1));
        
        angftrms =squeeze(rms(angfreq10diff,1));
        angfetrms=squeeze(rms(angfreqerp10diff,1));
        angfptrms =squeeze(rms(angfreqpn10diff,1));
        angfpetrms=squeeze(rms(angfreqpnerp10diff,1));
      end
      
      angfmean(:,:,:,frequse,kcampind)=  abs(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials);
      angfemean(:,:,:,frequse,kcampind)=abs(sum(exp(i*deg2rad(angfreqerpdiff)),1)/numtrials);
      angfpmean(:,:,:,frequse,kcampind)=abs(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials);
      angfpemean(:,:,:,frequse,kcampind)=abs(sum(exp(i*deg2rad(angfreqpnerpdiff)),1)/numtrials);
      angfang(:,:,:,frequse,kcampind)=  rad2deg(angle(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials));
      angfeang(:,:,:,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqerpdiff)),1)/numtrials));
      angfpang(:,:,:,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials));
      angfpeang(:,:,:,frequse,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpnerpdiff)),1)/numtrials));
      
    end % erpampind
  end % frequse
  
  if 0
    save(['D:\phase_estimation\anglermsrun6.mat'],'ang*rms*','frequse')
  else
    save(['D:\phase_estimation\anglermsrun6b.mat'],'ang*mean','ang*ang','-append')
  end
  
  load(['D:\phase_estimation\anglermsrun6b.mat'])
  
  for kcampind=1:3
    for tt=1:2
      figure(100+tt+(kcampind-1)*2);
      for tf=[1:4 6:10],
        subplot(4,length(t_ftimwin),tf);                     bar(1-squeeze(angfmean(:,tt,tf,:,kcampind))')    ;axis([-inf inf 0 0.5])
        if tf<5
          title(['Time Window length ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['Time Window ' num2str(t_ftimwin{tf}(3)*3) ' periods'])
        end
        if tf==1,ylabel('No PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+length(t_ftimwin));      bar(1-squeeze(angfpmean(:,tt,tf,:,kcampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+2*length(t_ftimwin));    bar(1-squeeze(angfemean(:,tt,tf,:,kcampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('No PN, ERP');end
        subplot(4,length(t_ftimwin),tf+3*length(t_ftimwin));    bar(1-squeeze(angfpemean(:,tt,tf,:,kcampind))')    ;axis([-inf inf 0 0.5])
        if tf==1,ylabel('PN, ERP');end
      end
      legend({'cfg.foi 1 Hz' '2 Hz' '3 Hz' '4 Hz'})
      set(get(100+tt+(kcampind-1)*2,'Children'),'xTickLabel',{'Sim 1' '2 Hz' '3' '4' '5' '6' '7'})
    end
  end
  
  for kcampind=1:3
    for tt=1:2
      figure(110+tt+(kcampind-1)*2);
      for tf=[1:4 6:10],
        subplot(4,length(t_ftimwin),tf);                     bar(squeeze(angfang(:,tt,tf,:,kcampind))')    ;axis([-inf inf -20 20])
        if tf<5
          title(['Time Window length ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['Time Window ' num2str(t_ftimwin{tf}(3)*3) ' periods'])
        end
        if tf==1,ylabel('No PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+length(t_ftimwin));      bar(squeeze(angfpang(:,tt,tf,:,kcampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('PN, No ERP');end
        subplot(4,length(t_ftimwin),tf+2*length(t_ftimwin));    bar(squeeze(angfeang(:,tt,tf,:,kcampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('No PN, ERP');end
        subplot(4,length(t_ftimwin),tf+3*length(t_ftimwin));    bar(squeeze(angfpeang(:,tt,tf,:,kcampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('PN, ERP');end
      end
      legend({'cfg.foi 1 Hz' '2 Hz' '3 Hz' '4 Hz'})
      set(get(110+tt+(kcampind-1)*2,'Children'),'xTickLabel',{'Sim 1' '2 Hz' '3' '4' '5' '6' '7'})
    end
  end
  
  for kcampind=1:3
    for tt=1:2
      figure(120+tt+(kcampind-1)*2);
      for tf=1:length(t_ftimwin),
        subplot(4,5,tf);       bar(1-squeeze(angfpemean(:,tt,tf,1:4,kcampind))')    ;axis([-inf inf 0 0.5])
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['TW ' num2str(3*t_ftimwin{tf}(3)) ' per.'])
        end
        if tf==1,ylabel('Circ-Var');end
        subplot(4,5,tf+10);    bar(squeeze(angfpeang(:,tt,tf,1:4,kcampind))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('Mean Angle');end
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}(3)) ' s'])
        else
          title(['TW ' num2str(3*t_ftimwin{tf}(3)) ' per.'])
        end
      end
      set(get(120+tt+(kcampind-1)*2,'Children'),'XTickLabel',{'Sim 1' '2 Hz' '3' '4'})
      legend({'cfg.foi 1 Hz' '2 Hz' '3 Hz' '4 Hz'})
    end
  end
  
  
  
  
  if 0
    figure(100);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfrms(:,:,tf));axis([-inf inf 0 160]);end
    figure(101);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angferms(:,:,tf));axis([-inf inf 0 160]);end
    figure(102);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfprms(:,:,tf));axis([-inf inf 0 160]);end
    figure(103);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar(cfg.foi,angfperms(:,:,tf));axis([-inf inf 0 160]);end
    
    figure(104);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([angfsrms(:,tf),angfesrms(:,tf),angftrms(:,tf),angfetrms(:,tf), angfpsrms(:,tf),angfpesrms(:,tf),angfptrms(:,tf),angfpetrms(:,tf)]');axis([-inf inf 0 160]);end
    set(get(104,'Children'),'xTickLabel',{'spec', 'Kc spec', '1hz', 'Kc 1hz', 'PN spec', 'PN Kc spec', 'PN 1hz', 'PN Kc 1hz'})
    
    % no point in using '1hz' as it's already in 1-4.
    figure(110);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 5 ],[angfrms(:,:,tf); angfsrms(:,tf)'; ]);axis([-inf inf 0 160]);end
    set(get(110,'Children'),'xTickLabel',{'1', '2', '3', '4', 'custom'})
    figure(111);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 5 ],[angferms(:,:,tf); angfesrms(:,tf)'; ]);axis([-inf inf 0 160]);end
    set(get(111,'Children'),'xTickLabel',{'1', '2', '3', '4', 'custom'})
    figure(112);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 5],[angfprms(:,:,tf); angfpsrms(:,tf)';]);axis([-inf inf 0 160]);end
    set(get(112,'Children'),'xTickLabel',{'1', '2', '3', '4', 'custom'})
    figure(113);
    for tf=1:length(t_ftimwin),subplot(2,5,tf);bar([cfg.foi 5 ],[angfperms(:,:,tf); angfpesrms(:,tf)'; ]);axis([-inf inf 0 160]);end
    set(get(113,'Children'),'xTickLabel',{'1', '2', '3', '4', 'custom'})
  end
  
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 6c) Wavelet
  
  foi=freqtest(1):freqtest(end);
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  cfg.foi=foi;
  cfg.pad=time(end);
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angwave   =nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwaveerp=nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwavepn   =nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  angwavepnerp=nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),3);
  
  angwavediff    = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwaveerpdiff = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwavepndiff    = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  angwavepnerpdiff = nan(numtrials,length(cfg.foi),length(toi),length(wavwidth),length(timwin),3);
  
  angwavespecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwaveerpspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepnspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepnerpspecdiff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  
  angwave10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwaveerp10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepn10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  angwavepnerp10diff = nan(numtrials,length(toi),length(wavwidth),length(timwin),3);
  
  
  for tf=1:length(timwin) % more important than gwidth?
    scfg=[];
    scfg.latency=time([halftimeind-fsample*timwin(tf)/2 halftimeind+fsample*timwin(tf)/2]);
    rawuse=ft_selectdata(scfg,raw);
    rawerpuse=ft_selectdata(scfg,rawerp);
    rawpnuse=ft_selectdata(scfg,rawpn);
    rawpnerpuse=ft_selectdata(scfg,rawpnerp);
    halftimeinduse(tf)=dsearchn(rawuse.time{1}',time(halftimeind));
    
    for ss=1:numtrials
      cfg.trials=ss;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoi=ft_freqanalysis(cfg,rawuse);
          waveerpfoi=ft_freqanalysis(cfg,rawerpuse);
          wavepnfoi=ft_freqanalysis(cfg,rawpnuse);
          wavepnerpfoi=ft_freqanalysis(cfg,rawpnerpuse);
          indfoi1=unique(dsearchn(foi',wavefoi.freq'));
          if tt==1
            angwave(ss,indfoi1,tt,ww,1)   = rad2deg(wrapToPi(angle(squeeze(wavefoi.fourierspctrm))));
            angwaveerp(ss,indfoi1,tt,ww,1)= rad2deg(wrapToPi(angle(squeeze(waveerpfoi.fourierspctrm))));
            angwavepn(ss,indfoi1,tt,ww,1)   = rad2deg(wrapToPi(angle(squeeze(wavepnfoi.fourierspctrm))));
            angwavepnerp(ss,indfoi1,tt,ww,1)= rad2deg(wrapToPi(angle(squeeze(wavepnerpfoi.fourierspctrm))));
          elseif tt==2
            angwave(ss,indfoi1,tt,ww,1)   = rad2deg(diag(wrapToPi(   angle(squeeze(wavefoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwaveerp(ss,indfoi1,tt,ww,1)= rad2deg(diag(wrapToPi(   angle(squeeze(waveerpfoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwavepn(ss,indfoi1,tt,ww,1)   = rad2deg(diag(wrapToPi(   angle(squeeze(wavepnfoi.fourierspctrm(1,1,:,:)))+pi   )));
            angwavepnerp(ss,indfoi1,tt,ww,1)= rad2deg(diag(wrapToPi(   angle(squeeze(wavepnerpfoi.fourierspctrm(1,1,:,:)))+pi   )));
          end
        end % ww
      end % tt
      % given that we know exactly which freq was used for a given trial
      angwavespec(ss,:,:)   =angwave(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwaveerpspec(ss,:,:)=angwaveerp(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwavepnspec(ss,:,:)   =angwavepn(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      angwavepnerpspec(ss,:,:)=angwavepnerp(ss,dsearchn(cfg.foi',round(frequse(ss))),:,:,1);
      % use the freq that was closest to the Kc (1.25Hz)
      angwave10(ss,:,:)   =angwave(ss,1,:,:,1);
      angwaveerp10(ss,:,:)=angwaveerp(ss,1,:,:,1);
      angwavepn10(ss,:,:)   =angwavepn(ss,1,:,:,1);
      angwavepnerp10(ss,:,:)=angwavepnerp(ss,1,:,:,1);
      
      % compute differences
      angwavediff(ss,:,:,:,tf,1)       =anglediff(angwave(ss,:,:,:,1),       trueang_prenoise(ss,halftimeind),1);
      angwaveerpdiff(ss,:,:,:,tf,1)    =anglediff(angwaveerp(ss,:,:,:,1),    trueang_prenoise(ss,halftimeind),1);
      angwavepndiff(ss,:,:,:,tf,1)       =anglediff(angwavepn(ss,:,:,:,1),       trueang_prenoise(ss,halftimeind),1);
      angwavepnerpdiff(ss,:,:,:,tf,1)    =anglediff(angwavepnerp(ss,:,:,:,1),    trueang_prenoise(ss,halftimeind),1);
      
      angwavespecdiff(ss,:,:,tf,1)   =anglediff(angwavespec(ss,:,:),   trueang_prenoise(ss,halftimeind),1);
      angwaveerpspecdiff(ss,:,:,tf,1)=anglediff(angwaveerpspec(ss,:,:),trueang_prenoise(ss,halftimeind),1);
      angwavepnspecdiff(ss,:,:,tf,1)   =anglediff(angwavepnspec(ss,:,:),   trueang_prenoise(ss,halftimeind),1);
      angwavepnerpspecdiff(ss,:,:,tf,1)=anglediff(angwavepnerpspec(ss,:,:),trueang_prenoise(ss,halftimeind),1);
      
      angwave10diff(ss,:,:,tf,1)     =anglediff(angwave10(ss,:,:),     trueang_prenoise(ss,halftimeind),1);
      angwaveerp10diff(ss,:,:,tf,1)  =anglediff(angwaveerp10(ss,:,:),  trueang_prenoise(ss,halftimeind),1);
      angwavepn10diff(ss,:,:,tf,1)     =anglediff(angwavepn10(ss,:,:),     trueang_prenoise(ss,halftimeind),1);
      angwavepnerp10diff(ss,:,:,tf,1)  =anglediff(angwavepnerp10(ss,:,:),  trueang_prenoise(ss,halftimeind),1);
    end % ss
  end % tf
  
  
  cfg=[];
  cfg.method='wavelet';
  cfg.output='fourier';
  cfg.foilim=freqtest;
  cfg.pad=time(end);
  
  for tf=1:length(timwin) % more important than gwidth
    scfg=[];
    scfg.latency=time([halftimeind-fsample*timwin(tf)/2 halftimeind+fsample*timwin(tf)/2]);
    rawuse=ft_selectdata(scfg,raw);
    rawerpuse=ft_selectdata(scfg,rawerp);
    rawpnuse=ft_selectdata(scfg,rawpn);
    rawpnerpuse=ft_selectdata(scfg,rawpnerp);
    halftimeinduse(tf)=dsearchn(rawuse.time{1}',time(halftimeind));
    
    for ss=1:numtrials
      cfg.trials=ss;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for ww=1:length(wavwidth)
          cfg.width=wavwidth(ww);
          wavefoilim=ft_freqanalysis(cfg,rawuse);
          waveerpfoilim=ft_freqanalysis(cfg,rawerpuse);
          wavepnfoilim=ft_freqanalysis(cfg,rawpnuse);
          wavepnerpfoilim=ft_freqanalysis(cfg,rawpnerpuse);
          % this produces .freq of length varying with timwin (spacing is 1/timwin)
          
          % one option: take the middle freq within the range and use that angle.
          indfoi2=dsearchn(foi', round(wavefoilim.freq(ceil(length(wavefoilim.freq)/2))));
          if tt==1
            angwave(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi(angle(squeeze(wavefoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwaveerp(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(waveerpfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwavepn(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi(angle(squeeze(wavepnfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
            angwavepnerp(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi(angle(squeeze(wavepnerpfoilim.fourierspctrm(:,:,ceil(length(wavefoilim.freq)/2)  )))));
          elseif tt==2
            angwave(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi( angle(squeeze(wavefoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwaveerp(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi( angle(squeeze(waveerpfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwavepn(ss,indfoi2,tt,ww,2)   =rad2deg(wrapToPi( angle(squeeze(wavepnfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
            angwavepnerp(ss,indfoi2,tt,ww,2)=rad2deg(wrapToPi( angle(squeeze(wavepnerpfoilim.fourierspctrm(1,1,ceil(length(wavefoilim.freq)/2) ,indfoi2)))+pi ));
          end
          % another option: take average over
          if tt==1
            angwave(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(angle(squeeze(wavefoilim.fourierspctrm)))));
            angwaveerp(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(angle(squeeze(waveerpfoilim.fourierspctrm)))));
            angwavepn(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(angle(squeeze(wavepnfoilim.fourierspctrm)))));
            angwavepnerp(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(angle(squeeze(wavepnerpfoilim.fourierspctrm)))));
          elseif tt==2
            angwave(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavefoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwaveerp(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(  diag(angle(squeeze(waveerpfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwavepn(ss,indfoi2,tt,ww,3)   =mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavepnfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
            angwavepnerp(ss,indfoi2,tt,ww,3)=mean(rad2deg(wrapToPi(  diag(angle(squeeze(wavepnerpfoilim.fourierspctrm(1,1,dsearchn(wavefoilim.freq', [floor(cfg.foilim(1)):ceil(cfg.foilim(2))]'),:))))+pi   )));
          end
          
        end % ww
      end % tt
      
      % compute differences
      for dd=2:3
        angwavediff(ss,indfoi2,:,:,tf,dd)       =anglediff(angwave(ss,indfoi2,:,:,dd),       trueang_prenoise(ss,halftimeind),1);
        angwaveerpdiff(ss,indfoi2,:,:,tf,dd)    =anglediff(angwaveerp(ss,indfoi2,:,:,dd),    trueang_prenoise(ss,halftimeind),1);
        angwavepndiff(ss,indfoi2,:,:,tf,dd)       =anglediff(angwavepn(ss,indfoi2,:,:,dd),       trueang_prenoise(ss,halftimeind),1);
        angwavepnerpdiff(ss,indfoi2,:,:,tf,dd)    =anglediff(angwavepnerp(ss,indfoi2,:,:,dd),    trueang_prenoise(ss,halftimeind),1);
      end
    end  % ss
  end % tf
  
  
  
  angwrms  =squeeze(rms(angwavediff,1));
  angwerms =squeeze(rms(angwaveerpdiff,1));
  angwprms  =squeeze(rms(angwavepndiff,1));
  angwperms =squeeze(rms(angwavepnerpdiff,1));
  
  angwsrms =squeeze(rms(angwavespecdiff,1));
  angwesrms=squeeze(rms(angwaveerpspecdiff,1));
  angwpsrms =squeeze(rms(angwavepnspecdiff,1));
  angwpesrms=squeeze(rms(angwavepnerpspecdiff,1));
  
  angwtrms =squeeze(rms(angwave10diff,1));
  angwetrms=squeeze(rms(angwaveerp10diff,1));
  angwptrms =squeeze(rms(angwavepn10diff,1));
  angwpetrms=squeeze(rms(angwavepnerp10diff,1));
  
  save(['D:\phase_estimation\anglermsrun6.mat'],'ang*rms*','frequse')
  
  % 3rd dim index 4 means wavwidth=7
  for dd=1:3 % 3 diff ways of computing
    figure(200+dd*10);
    for tf=1:length(timwin),subplot(1,length(timwin),tf);bar(foi,angwrms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(201+dd*10);
    for tf=1:length(timwin),subplot(1,length(timwin),tf);bar(foi,angwerms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(202+dd*10);
    for tf=1:length(timwin),subplot(1,length(timwin),tf);bar(foi,angwprms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
    figure(203+dd*10);
    for tf=1:length(timwin),subplot(1,length(timwin),tf);bar(foi,angwperms(:,:,4,tf,dd));axis([-inf inf 0 160]);end
  end
  figure(204);
  for tf=1:length(timwin),subplot(1,length(timwin),tf);bar([angwsrms(:,tf),angwesrms(:,tf),angwtrms(:,tf),angwetrms(:,tf), angwpsrms(:,tf),angwpesrms(:,tf),angwptrms(:,tf),angwpetrms(:,tf)]');axis([-inf inf 0 160]);end
  set(get(204,'Children'),'xTickLabel',{'spec', 'Kc spec', '1hz', 'Kc 1hz', 'PN spec', 'PN Kc spec', 'PN 1hz', 'PN Kc 1hz'})
  
end % end 6

%% 8) Random frequency, then Hilbert
% Create new data, with ERP just after time point of interest, with additive model
% One raw structure with 100 trials: each at slight diff phase but same
% random frequency within alpha band; then additive ERP (8a) or phase-reset ERP (8b)

if run8
  
  
  % params for whole simulation
  clearvars -except run* time wdr fsample wav*width
  close all;
  numtrials=100;
  halftimeind=round(length(time)/2);
  plotflag=1;
  foilim=[8 12];
  frequse = diff(foilim)*rand(1)+foilim(1);
  
  cd('D:\fieldtrip_svn\utilities\private');
  state=randomseed(13);
  cd(wdr);
  
  
  % 8a)
  clear raw
  erpamp=[0.5 1 2];
  
  
  for erpampind=1:length(erpamp)
    raw{erpampind}.label{1}='test';
    raw{erpampind}.dimord='chan_time';
    raw{erpampind}.frequse=frequse;
    [raw{erpampind}.time{1:numtrials}]=deal(time);
    rawpn{erpampind}=raw{erpampind};
    rawerp{erpampind}=raw{erpampind};
    rawpnerp{erpampind}=raw{erpampind};
    for tr=1:numtrials % simulate 100 trials
      % Always use 'cos' to generate signals
      % Each trial has a base frequency, pink noise, and ERP added
      phaseshift(erpampind,tr) = wrapToPi(2*pi*rand(1));
      trueang(:,erpampind,tr) = frequse*2*pi*time+phaseshift(erpampind,tr)*ones(size(time)); % random freq and phase
      raw{erpampind}.trial{tr}=cos(trueang(:,erpampind,tr))';
      trueang_prenoise(:,erpampind,tr)=rad2deg(wrapToPi(trueang(:,erpampind,tr)));
      rawpn{erpampind}.trial{tr}=raw{erpampind}.trial{tr}+10*pinknoise(length(time));
      trueang_postnoise(:,erpampind,tr)=rad2deg(wrapToPi(angle(hilbert(rawpn{erpampind}.trial{tr})))); % after noise added.
      % add sinusoid (starting at amplitude zero) from time 0-100ms, at 10.1Hz frequency
      rawpnerp{erpampind}.trial{tr}=rawpn{erpampind}.trial{tr}+[zeros(1,halftimeind-1), erpamp(erpampind)*gausswin(301,2)'.*sin(10.1*2*pi*time(halftimeind:halftimeind+300)-10.1*2*pi*time(halftimeind)), zeros(1,halftimeind-301)];
      rawerp{erpampind}.trial{tr}  =raw{erpampind}.trial{tr}  +[zeros(1,halftimeind-1), erpamp(erpampind)*gausswin(301,2)'.*sin(10.1*2*pi*time(halftimeind:halftimeind+300)-10.1*2*pi*time(halftimeind)), zeros(1,halftimeind-301)];
    end
    if plotflag
      cfg=[];
      tlock=ft_timelockanalysis(cfg,raw{erpampind});
      tlockpn=ft_timelockanalysis(cfg,rawpn{erpampind});
      tlockerp=ft_timelockanalysis(cfg,rawerp{erpampind});
      tlockpnerp=ft_timelockanalysis(cfg,rawpnerp{erpampind});
      figure;plot(tlock.time,tlock.avg)
      hold on;plot(tlockpn.time,tlockpn.avg,'g');axis([-inf inf -2 2])
      hold on;plot(tlockerp.time,tlockerp.avg,'r');axis([-inf inf -2 2])
      hold on;plot(tlockerp.time,tlockpnerp.avg,'k');axis([-inf inf -2 2])
    end
  end
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 8a) Bandpass filter + Hilbert
  
  close all
  
  % Filter prestim and assess peak frequency with Hilbert
  
  
  
  
  for erpampind=1:length(erpamp)
    
    cfg=[];
    cfg.latency=[time(1) time(halftimeind)];
    rawpnerpuse=ft_selectdata(cfg,rawpnerp{erpampind});
    halftimeinduse=dsearchn(rawpnerpuse.time{1}',time(halftimeind));
    
    cfg=[];
    cfg.bpfilter='yes';
    cfg.bpfreq=[7 13];
    cfg.bpfilttype='fir';
    cfg.plotfiltresp='yes';
    cfg.fouse=round([3*fsample/cfg.bpfreq(1)]); % 3* is default
    cfg.figind=10;
    cfg.plotflag=0;
    cfg.hilbert='complex';
    rawpnerp_bpfir=filter4phase_estim8(cfg,rawpnerpuse);
    for tr=1:numtrials
      diffang=diff(angle(rawpnerp_bpfir{2}.trial{tr}));  % For this, we trust two-pass more than one-pass zerophase.
      diffang(diffang<-6)=nan;
      mediandiffang(tr)=nanmedian(diffang);
    end
    freqest(erpampind)=mean(mediandiffang)/(.001*2*pi);
    
    
    % The peak frequency found is biased toward centre of bandpass range.
    % Thus, run again now centred on peak found.  Typically it adjusts
    % slightly in the right direction (in ideal case of perfect sinusoid)
    cfg.bpfreq=[freqest(erpampind)-2 freqest(erpampind)+2];
    rawpnerp_bpfir=filter4phase_estim8(cfg,rawpnerpuse);
    for tr=1:numtrials
      diffang=diff(angle(rawpnerp_bpfir{2}.trial{tr}));  % For this, we trust two-pass more than one-pass zerophase.
      diffang(diffang<-6)=nan;
      mediandiffang(tr)=nanmedian(diffang);
    end
    freqest2(erpampind)=mean(mediandiffang)/(.001*2*pi);
    
  end
  
  % 8a) FFT + taper
  
  timwin=[4 2 1 .5 .25]; % duration in seconds; must be ordered from longest to shortest
  
  cfg=[];
  cfg.method='mtmconvol';
  cfg.output='fourier';
  cfg.taper='hanning';
  cfg.keeptrials='yes';
  cfg.foi=mean(freqest2);
  
  
  
  t_ftimwin{1}=timwin(1)*ones(size(cfg.foi)); % full length of data
  t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 2 periods 4 Hz
  t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{6}=4./cfg.foi; %
  t_ftimwin{7}=3./cfg.foi; %
  t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
  t_ftimwin{9}=1./cfg.foi; % one period for each frequency
  t_ftimwin{10}=0.5./cfg.foi; %
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angfreq=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqerp=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpn=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpnerp=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  
  for erpampind=1:length(erpamp)
    for ss=1:numtrials
      cfg.trials=ss;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for tf=1:length(t_ftimwin)
          cfg.t_ftimwin=t_ftimwin{tf};
          freq  =ft_freqanalysis(cfg, raw{erpampind});
          freqerp=ft_freqanalysis(cfg, rawerp{erpampind});
          freqpn  =ft_freqanalysis(cfg, rawpn{erpampind});
          freqpnerp=ft_freqanalysis(cfg, rawpnerp{erpampind});
          if tt==1
            %           angfreq(:,:,tt,pp,tf)=angle(squeeze(freq{tt,pp,tf}.fourierspctrm))/(2*pi)*360;
            angfreq(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freq.fourierspctrm))));
            angfreqerp(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqerp.fourierspctrm))));
            angfreqpn(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freqpn.fourierspctrm))));
            angfreqpnerp(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqpnerp.fourierspctrm))));
          elseif tt==2
            angfreq(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freq.fourierspctrm(1,1,:,:)))+pi  )));
            angfreqerp(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqerp.fourierspctrm(1,1,:,:)))+pi  )));
            angfreqpn(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freqpn.fourierspctrm(1,1,:,:)))+pi  )));
            angfreqpnerp(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqpnerp.fourierspctrm(1,1,:,:)))+pi  )));
          end
        end % tf
      end % tt
      
      % compute differences
      angfreqdiff(ss,:,:,:)       =anglediff(angfreq(ss,:,:,:),       trueang_prenoise(halftimeind,erpampind,ss),1);
      angfreqerpdiff(ss,:,:,:)    =anglediff(angfreqerp(ss,:,:,:),    trueang_prenoise(halftimeind,erpampind,ss),1);
      angfreqpndiff(ss,:,:,:)       =anglediff(angfreqpn(ss,:,:,:),       trueang_prenoise(halftimeind,erpampind,ss),1);
      angfreqpnerpdiff(ss,:,:,:)    =anglediff(angfreqpnerp(ss,:,:,:),    trueang_prenoise(halftimeind,erpampind,ss),1);
      
    end % ss
    
    angfmean(:,:,:,erpampind)=  abs(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials);
    angfemean(:,:,:,erpampind)=abs(sum(exp(i*deg2rad(angfreqerpdiff)),1)/numtrials);
    angfpmean(:,:,:,erpampind)=abs(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials);
    angfpemean(:,:,:,erpampind)=abs(sum(exp(i*deg2rad(angfreqpnerpdiff)),1)/numtrials);
    angfang(:,:,:,erpampind)=  rad2deg(angle(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials));
    angfeang(:,:,:,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqerpdiff)),1)/numtrials));
    angfpang(:,:,:,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials));
    angfpeang(:,:,:,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpnerpdiff)),1)/numtrials));
    
  end % erpampind
  
  try
    save(['D:\phase_estimation\anglermsrun8a.mat'],'ang*mean','ang*ang','freqest2','frequse','-append')
  catch
    save(['D:\phase_estimation\anglermsrun8a.mat'],'ang*mean','ang*ang','freqest2','frequse')
  end
  
  load(['D:\phase_estimation\anglermsrun8a.mat'])
  
  
  for tt=1:2
    figure(120+tt);
    for tf=1:length(t_ftimwin),
      subplot(4,5,tf);       bar(1-squeeze(angfpemean(:,tt,tf,:))')    ;axis([-inf inf 0 0.5])
      if tf<6
        title(['TW ' num2str(t_ftimwin{tf}) ' s'])
      else
        title(['TW ' num2str(t_ftimwin{tf}*cfg.foi) ' per.'])
      end
      if tf==1,ylabel('Circ-Var');end
      subplot(4,5,tf+10);    bar(squeeze(angfpeang(:,tt,tf,:))')    ;axis([-inf inf -20 20])
      if tf==1,ylabel('Mean Angle');end
      if tf<6
        title(['TW ' num2str(t_ftimwin{tf}) ' s'])
      else
        title(['TW ' num2str(t_ftimwin{tf}*cfg.foi) ' per.'])
      end
    end
    %       set(get(120+tt+(erpampind-1)*2,'Children'),'xTickLabel',{'Sim 8' '9' '10' '11' '12'})
    set(get(120+tt,'Children'),'XTickLabel',{'ERPamp1' 'ERPamp2' 'ERPamp3'})
    %legend({'cfg.foi 8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  end
  
  %% 8b
  % inspired by
  % phasereset (frames, epochs, srate, minfr, maxfr, position, tjitter)
  % from % Implemented by: Rafal Bogacz and Nick Yeung, September 2006
  
  clearvars -except run* time wdr fsample wav*width
  addpath('D:\phase_estimation\phasereset_modelling')
  halftimeind=round(length(time)/2);
  numtrials=100;
  foilim=[8 12];
  frequse = diff(foilim)*rand(1)+foilim(1);
  
  numones=[1 4 7];
  resetwin=[40 20 0];
  
  % Model 1:  frequency shift to meet phase at later time; return to new frequency
  % Model 12: Same effectively as 1, except modelled as linear phase procession (at new frequency not determinted explicitly); return to new frequency
  % Model 2:  frequency shift to meet phase at later time; return to original frequency
  % Model 22: Same effectively as 2, except modelled as linear phase procession (at new frequency not determinted explicitly); return to original frequency
  % Model 3:  short time to meet new phase via line;
  % Model 4:  short time to meet new phase via line; new phase determined by 'phase-response-curve'
  modeluse=12;
  
  for erpampind=1:length(resetwin)
    cd('D:\fieldtrip_svn\utilities\private');
    state=randomseed(13);
    cd(wdr);
    
    % control with no phase reset
    raw{erpampind}.label{1}='test';
    raw{erpampind}.dimord='chan_time';
    [raw{erpampind}.time{1:numtrials}]=deal(time);
    rawpn{erpampind}=raw{erpampind};
    rawpr{erpampind}=raw{erpampind};
    rawpnpr{erpampind}=raw{erpampind};
    for tr=1:numtrials
      switch modeluse
        case 1
          [raw{erpampind}.trial{tr}] = phasereset_jz (time, fsample, frequse, foilim(1), foilim(2), length(time)-1, length(time), 1);
        case 12
          [raw{erpampind}.trial{tr}] = phasereset_jz_linphase_postjitter (time, frequse, length(time)-1, length(time), 0 );
        case 3
          %       [raw{erpampind}.trial{tr}] = phasereset_jz_jump (time, frequse, length(time)-1, length(time), numones(erpampind) );
          [raw{erpampind}.trial{tr}] = phasereset_jz_jump (time, frequse, length(time)-1, length(time), 0 );
        case 2
          [raw{erpampind}.trial{tr}] = phasereset_jz_steadybasefreq (time, fsample, frequse, foilim(1), foilim(2), length(time)-1, length(time), 1);
        case 22
          [raw{erpampind}.trial{tr}] = phasereset_jz_linphase (time, frequse, length(time)-1, length(time), 0 );
        case 4 % phase response curve
          [raw{erpampind}.trial{tr}] = phasereset_jz_jump_prc (time, frequse, length(time)-1, length(time), 0 );
      end
      rawpn{erpampind}.trial{tr}=raw{erpampind}.trial{tr}+10*pinknoise(length(time));
    end
    
    % set back again to exact same state, thus same phase at time points until phase reset
    cd('D:\fieldtrip_svn\utilities\private');
    state=randomseed(13);
    cd(wdr);
    
    % do proper reset at halfway point
    for tr=1:numtrials
      switch modeluse
        case 1
          [rawpr{erpampind}.trial{tr},trueang_prenoise(tr,erpampind),newfreq(tr,erpampind)] = phasereset_jz (time, fsample, frequse, foilim(1), foilim(2), halftimeind, halftimeind+.1*fsample, numones(erpampind));
        case 12
          [rawpr{erpampind}.trial{tr},trueang_prenoise(tr,erpampind),resetang(tr,erpampind)] = phasereset_jz_linphase_postjitter (time, frequse, halftimeind, halftimeind+20, resetwin(erpampind) );
        case 3
          %       [rawpr{erpampind}.trial{tr},trueang_prenoise(tr,erpampind)] = phasereset_jz_jump (time, frequse, halftimeind, halftimeind+20, numones(erpampind) );
          [rawpr{erpampind}.trial{tr},trueang_prenoise(tr,erpampind),resetang(tr,erpampind)] = phasereset_jz_jump (time, frequse, halftimeind, halftimeind+20, resetwin(erpampind) );
        case 2
          [rawpr{erpampind}.trial{tr},trueang_prenoise(tr,erpampind)] = phasereset_jz_steadybasefreq (time, fsample, frequse, foilim(1), foilim(2), halftimeind, halftimeind+.1*fsample, numones(erpampind));
        case 22
          [rawpr{erpampind}.trial{tr},trueang_prenoise(tr,erpampind),resetang(tr,erpampind)] = phasereset_jz_linphase (time, frequse, halftimeind, halftimeind+20, resetwin(erpampind) );
        case 4
      end
      rawpnpr{erpampind}.trial{tr}=rawpr{erpampind}.trial{tr}+10*pinknoise(length(time));
    end
    
    prc(:,erpampind)=resetang(:,erpampind)-trueang_prenoise(:,erpampind);
    trueang_prenoise(:,erpampind)=rad2deg(wrapToPi(trueang_prenoise(:,erpampind)));
    
    cfg=[];
    tlock=ft_timelockanalysis(cfg,raw{erpampind});
    tlockpr=ft_timelockanalysis(cfg,rawpr{erpampind});
    tlockpn=ft_timelockanalysis(cfg,rawpn{erpampind});
    tlockpnpr=ft_timelockanalysis(cfg,rawpnpr{erpampind});
    
    if 1
      figure;plot(time,[tlock.avg; tlockpr.avg; tlockpn.avg; tlockpnpr.avg]);axis([-inf inf -1 1])
    end
    
  end % erpampind
  
  
  cfg=[];
  for erpampind=1:3
    tlockpnpr=ft_timelockanalysis(cfg,rawpnpr{erpampind});
    figure;plot(time,[tlockpnpr.avg]);axis([-inf inf -1 1])
  end
  
  
  startang=resetang-prc;
  
  
  % 8b) bandpass and Hilbert to find peak frequency
  for erpampind=1:length(resetwin)
    
    cfg=[];
    cfg.latency=[time(1) time(halftimeind)];
    rawpnerpuse=ft_selectdata(cfg,rawpnpr{erpampind});
    halftimeinduse=dsearchn(rawpnerpuse.time{1}',time(halftimeind));
    
    cfg=[];
    cfg.bpfilter='yes';
    cfg.bpfreq=[7 13];
    cfg.bpfilttype='fir';
    cfg.plotfiltresp='yes';
    cfg.fouse=round([3*fsample/cfg.bpfreq(1)]); % 3* is default
    cfg.figind=10;
    cfg.plotflag=0;
    cfg.hilbert='complex';
    rawpnerp_bpfir=filter4phase_estim8(cfg,rawpnerpuse);
    for tr=1:numtrials
      diffang=diff(angle(rawpnerp_bpfir{2}.trial{tr}));  % For this, we trust two-pass more than one-pass zerophase.
      diffang(diffang<-6)=nan;
      mediandiffang(tr)=nanmedian(diffang);
    end
    freqest(erpampind)=mean(mediandiffang)/(.001*2*pi);
    
    
    % The peak frequency found is biased toward centre of bandpass range.
    % Thus, run again now centred on peak found.  Typically it adjusts
    % slightly in the right direction (in ideal case of perfect sinusoid)
    cfg.bpfreq=[freqest(erpampind)-2 freqest(erpampind)+2];
    rawpnerp_bpfir=filter4phase_estim8(cfg,rawpnerpuse);
    for tr=1:numtrials
      diffang=diff(angle(rawpnerp_bpfir{2}.trial{tr}));  % For this, we trust two-pass more than one-pass zerophase.
      diffang(diffang<-6)=nan;
      mediandiffang(tr)=nanmedian(diffang);
    end
    freqest2(erpampind)=mean(mediandiffang)/(.001*2*pi);
    
  end
  
  % 8b) FFT at freqest
  
  timwin=[4 2 1 .5 .25]; % duration in seconds; must be ordered from longest to shortest
  cfg=[];
  cfg.method='mtmconvol';
  cfg.output='fourier';
  cfg.taper='hanning';
  cfg.keeptrials='yes';
  cfg.foi=mean(freqest2);
  
  
  t_ftimwin{1}=timwin(1)*ones(size(cfg.foi)); % full length of data
  t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 2 periods 4 Hz
  t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{6}=4./cfg.foi; %
  t_ftimwin{7}=3./cfg.foi; %
  t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
  t_ftimwin{9}=1./cfg.foi; % one period for each frequency
  t_ftimwin{10}=0.5./cfg.foi; %
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angfreq=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpr=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpn=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpnpr=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  
  for erpampind=1:length(resetwin)
    for ss=1:numtrials
      cfg.trials=ss;
      for tt=1:length(toi)
        cfg.toi=toi{tt};
        for tf=1:length(t_ftimwin)
          cfg.t_ftimwin=t_ftimwin{tf};
          freq  =ft_freqanalysis(cfg, raw{erpampind});
          freqpr=ft_freqanalysis(cfg, rawpr{erpampind});
          freqpn  =ft_freqanalysis(cfg, rawpn{erpampind});
          freqpnpr=ft_freqanalysis(cfg, rawpnpr{erpampind});
          if tt==1
            %           angfreq(:,:,tt,pp,tf)=angle(squeeze(freq{tt,pp,tf}.fourierspctrm))/(2*pi)*360;
            angfreq(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freq.fourierspctrm))));
            angfreqpr(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqpr.fourierspctrm))));
            angfreqpn(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freqpn.fourierspctrm))));
            angfreqpnpr(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqpnpr.fourierspctrm))));
          elseif tt==2
            angfreq(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freq.fourierspctrm(1,1,:,:)))+pi  )));
            angfreqpr(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqpr.fourierspctrm(1,1,:,:)))+pi  )));
            angfreqpn(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freqpn.fourierspctrm(1,1,:,:)))+pi  )));
            angfreqpnpr(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqpnpr.fourierspctrm(1,1,:,:)))+pi  )));
          end
        end % tf
      end % tt
      
      % compute differences
      angfreqdiff(ss,:,:,:)       =anglediff(angfreq(ss,:,:,:),       trueang_prenoise(ss,erpampind),1);
      angfreqprdiff(ss,:,:,:)    =anglediff(angfreqpr(ss,:,:,:),    trueang_prenoise(ss,erpampind),1);
      angfreqpndiff(ss,:,:,:)       =anglediff(angfreqpn(ss,:,:,:),       trueang_prenoise(ss,erpampind),1);
      angfreqpnprdiff(ss,:,:,:)    =anglediff(angfreqpnpr(ss,:,:,:),    trueang_prenoise(ss,erpampind),1);
      
      
    end
    
    angfmean(:,:,:,erpampind)=  abs(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials);
    angfemean(:,:,:,erpampind)=abs(sum(exp(i*deg2rad(angfreqprdiff)),1)/numtrials);
    angfpmean(:,:,:,erpampind)=abs(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials);
    angfpemean(:,:,:,erpampind)=abs(sum(exp(i*deg2rad(angfreqpnprdiff)),1)/numtrials);
    angfang(:,:,:,erpampind)=  rad2deg(angle(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials));
    angfeang(:,:,:,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqprdiff)),1)/numtrials));
    angfpang(:,:,:,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials));
    angfpeang(:,:,:,erpampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpnprdiff)),1)/numtrials));
    
  end % erpampind
  
  try
    save(['D:\phase_estimation\anglermsrun8b_model' num2str(modeluse) '.mat'],'ang*mean','ang*ang','freqest2','frequse','-append')
  catch
    save(['D:\phase_estimation\anglermsrun8b_model' num2str(modeluse) '.mat'],'ang*mean','ang*ang','freqest2','frequse')
  end
  
  load(['D:\phase_estimation\anglermsrun8b_model' num2str(modeluse) '.mat'])
  
  
  for tt=1:2
    figure(130+tt);
    for tf=1:length(t_ftimwin),
      subplot(4,5,tf);       bar(1-squeeze(angfpemean(1,tt,tf,:))')    ;axis([-inf inf 0 0.5])
      if tf<6
        title(['TW ' num2str(t_ftimwin{tf}) ' s'])
      else
        title(['TW ' num2str(t_ftimwin{tf}*cfg.foi) ' per.'])
      end
      if tf==1,ylabel('Circ-Var');end
      subplot(4,5,tf+10);    bar(squeeze(angfpeang(1,tt,tf,:))')    ;axis([-inf inf -20 20])
      if tf==1,ylabel('Mean Angle');end
      if tf<6
        title(['TW ' num2str(t_ftimwin{tf}) ' s'])
      else
        title(['TW ' num2str(t_ftimwin{tf}*cfg.foi) ' per.'])
      end
    end
    %       set(get(120+tt+(erpampind-1)*2,'Children'),'xTickLabel',{'Sim 8' '9' '10' '11' '12'})
    set(get(130+tt,'Children'),'XTickLabel',{'ERP PR1' 'ERP PR2' 'ERP PR3'})
    %     legend({'cfg.foi 8 Hz' '9 Hz', '10 Hz', '11 Hz' '12 Hz'})
  end
  
  
  
  
end % end 8



%% 9) Create new data, with Kc evoked at time point of interest, with additive model
% nearly same as (6)

if run9
  
  % params for whole simulation
  clearvars -except run* time wdr fsample wav*width
  close all;
  cd('D:\fieldtrip_svn\utilities\private');
  state=randomseed(13);
  cd(wdr);
  
  clear raw
  numtrials=100;
  halftimeind=round(length(time)/2);
  kcamp=[3 6 9];
  
  foilim=[1 4];
  frequse = diff(foilim)*rand(1)+foilim(1);
  
  for kcampind=1:length(kcamp)
    raw{kcampind}.label{1}='test';
    raw{kcampind}.dimord='chan_time';
    [raw{kcampind}.time{1:numtrials}]=deal(time);
    rawpn{kcampind}=raw{kcampind};
    rawerp{kcampind}=raw{kcampind};
    rawpnerp{kcampind}=raw{kcampind};
    
    for tr=1:numtrials % simulate 100 trials
      % Always use 'cos' to generate signals
      % Each trial has a base frequency, pink noise, and ERP added
      phaseshift(tr,kcampind) = wrapToPi(2*pi*rand(1));
      trueang(:,tr,kcampind) = frequse*2*pi*time+phaseshift(tr,kcampind)*ones(size(time)); % random phase
      raw{kcampind}.trial{tr}=cos(trueang(:,tr,kcampind)');
      trueang_prenoise(:,tr,kcampind)=rad2deg(wrapToPi(trueang(:,tr,kcampind)));
      rawpn{kcampind}.trial{tr}=raw{kcampind}.trial{tr}+10*pinknoise(length(time));
      trueang_postnoise(:,tr,kcampind)=rad2deg(wrapToPi(angle(hilbert(rawpn{kcampind}.trial{tr})))); % after noise added.
      % add sinusoid (starting at amplitude zero) from time 0-100ms, at 1.25Hz frequency
      rawpnerp{kcampind}.trial{tr}=rawpn{kcampind}.trial{tr}+[zeros(1,halftimeind-1), kcamp(kcampind)*gausswin(801,1)'.*sin(1.25*2*pi*time(halftimeind:halftimeind+800)-1.25*2*pi*time(halftimeind)), zeros(1,halftimeind-801)];
      rawerp{kcampind}.trial{tr}  =raw{kcampind}.trial{tr}  +[zeros(1,halftimeind-1), kcamp(kcampind)*gausswin(801,1)'.*sin(1.25*2*pi*time(halftimeind:halftimeind+800)-1.25*2*pi*time(halftimeind)), zeros(1,halftimeind-801)];
    end
    cfg=[];
    tlock=ft_timelockanalysis(cfg,raw{kcampind});
    tlockpn=ft_timelockanalysis(cfg,rawpn{kcampind});
    tlockerp=ft_timelockanalysis(cfg,rawerp{kcampind});
    tlockpnerp=ft_timelockanalysis(cfg,rawpnerp{kcampind});
    
    
    if 1
      figure;plot(tlock.time,raw{kcampind}.trial{1})
      hold on;plot(tlockpn.time,rawpn{kcampind}.trial{1},'g');axis([-inf inf -10 10])
      hold on;plot(tlockerp.time,rawerp{kcampind}.trial{1},'r');axis([-inf inf -10 10])
      hold on;plot(tlockerp.time,rawpnerp{kcampind}.trial{1},'k');axis([-inf inf -10 10])
    end
  end
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 9) Bandpass filter + Hilbert
  
  close all
  
  timwin=[4 2 1 .5]; % duration in seconds; must be ordered from longest to shortest
  
  for kcampind=1:length(kcamp)
    
    cfg=[];
    cfg.latency=[time(1) time(halftimeind)];
    rawpnerpuse=ft_selectdata(cfg,rawpnerp{kcampind});
    
    
    % Bandpass filter: FIR (Matlab 'fir1')
    rawpn_bpfir{3,2}=[];
    raw_bpfir{3,2}=[];
    rawerp_bpfir{3,2}=[];
    rawpnerp_bpfir{3,2}=[];
    cfg=[];
    cfg.bpfilter='yes';
    cfg.bpfreq=[0.5 5];
    cfg.bpfilttype='fir';
    cfg.plotfiltresp='yes';
    cfg.fouse=[3*fsample/cfg.bpfreq(1)]; % 3* is default
    cfg.figind=10;
    cfg.plotflag=0;
    cfg.hilbert='complex';
    rawpnerp_bpfir=filter4phase_estim8(cfg,rawpnerpuse);
    
    
    for tr=1:numtrials
      diffang=diff(angle(rawpnerp_bpfir{2}.trial{tr}));  % For this, we trust two-pass more than one-pass zerophase.
      diffang(diffang<-6)=nan;
      mediandiffang(tr)=nanmedian(diffang);
    end
    freqest(kcampind)=mean(mediandiffang)/(.001*2*pi);
    
    
    % The peak frequency found is biased toward centre of bandpass range.
    % Thus, run again now centred on peak found.  Typically it adjusts
    % slightly in the right direction (in ideal case of perfect sinusoid)
    cfg.bpfreq=[freqest(kcampind)-2 freqest(kcampind)+2];
    if cfg.bpfreq(1)<0.3, cfg.bpfreq(1)=0.3; end
    rawpnerp_bpfir=filter4phase_estim8(cfg,rawpnerpuse);
    for tr=1:numtrials
      diffang=diff(angle(rawpnerp_bpfir{2}.trial{tr}));  % For this, we trust two-pass more than one-pass zerophase.
      diffang(diffang<-6)=nan;
      mediandiffang(tr)=nanmedian(diffang);
    end
    freqest2(kcampind)=mean(mediandiffang)/(.001*2*pi);
    
    
  end % kcampind
  clear angdiff*
  
  
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 9) FFT + taper
  
  cfg=[];
  cfg.method='mtmconvol';
  cfg.output='fourier';
  cfg.taper='hanning';
  cfg.foi=mean(freqest2);
  cfg.keeptrials='yes';
  
  t_ftimwin{1}=timwin(1)*ones(size(cfg.foi)); % full length of data
  t_ftimwin{2}=timwin(2)*ones(size(cfg.foi)); % to match Hilbert calculations  % 2 periods 4 Hz
  t_ftimwin{3}=timwin(3)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{4}=timwin(4)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  %   t_ftimwin{5}=timwin(5)*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{5}=nan*ones(size(cfg.foi)); % to match Hilbert calculations % 1 period 4 Hz
  t_ftimwin{6}=4./cfg.foi; %
  t_ftimwin{7}=3./cfg.foi; %
  t_ftimwin{8}=2./cfg.foi; % two periods for each frequency
  t_ftimwin{9}=1./cfg.foi; % one period for each frequency
  t_ftimwin{10}=0.5./cfg.foi; %
  
  % first, centre on time of interest
  toi{1}=time(halftimeind);
  % second, centre half period before time of interest and add pi
  toi{2}=time(halftimeind)-0.5./cfg.foi;
  
  angfreq=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqerp=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpn=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  angfreqpnerp=nan(numtrials,length(cfg.foi),length(toi),length(t_ftimwin));
  
    for kcampind=1:length(kcamp)
      for ss=1:numtrials
        cfg.trials=ss;
        for tt=1:length(toi)
          cfg.toi=toi{tt};
          for tf=setdiff(1:length(t_ftimwin),5)
            cfg.t_ftimwin=t_ftimwin{tf};
            freq  =ft_freqanalysis(cfg, raw{kcampind});
            freqerp=ft_freqanalysis(cfg, rawerp{kcampind});
            freqpn  =ft_freqanalysis(cfg, rawpn{kcampind});
            freqpnerp=ft_freqanalysis(cfg, rawpnerp{kcampind});
            if tt==1
              %           angfreq(:,:,tt,pp,tf)=angle(squeeze(freq{tt,pp,tf}.fourierspctrm))/(2*pi)*360;
              angfreq(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freq.fourierspctrm))));
              angfreqerp(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqerp.fourierspctrm))));
              angfreqpn(ss,:,tt,tf)   =rad2deg(wrapToPi(angle(squeeze(freqpn.fourierspctrm))));
              angfreqpnerp(ss,:,tt,tf)=rad2deg(wrapToPi(angle(squeeze(freqpnerp.fourierspctrm))));
            elseif tt==2
              angfreq(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freq.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqerp(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqerp.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpn(ss,:,tt,tf)   =rad2deg(diag(wrapToPi(angle(squeeze(freqpn.fourierspctrm(1,1,:,:)))+pi  )));
              angfreqpnerp(ss,:,tt,tf)=rad2deg(diag(wrapToPi(angle(squeeze(freqpnerp.fourierspctrm(1,1,:,:)))+pi  )));
            end
          end % tf
        end % tt

        % compute differences
        angfreqdiff(ss,:,:,:)       =anglediff(angfreq(ss,:,:,:),       trueang_prenoise(halftimeind,ss,kcampind),1);
        angfreqerpdiff(ss,:,:,:)    =anglediff(angfreqerp(ss,:,:,:),    trueang_prenoise(halftimeind,ss,kcampind),1);
        angfreqpndiff(ss,:,:,:)       =anglediff(angfreqpn(ss,:,:,:),       trueang_prenoise(halftimeind,ss,kcampind),1);
        angfreqpnerpdiff(ss,:,:,:)    =anglediff(angfreqpnerp(ss,:,:,:),    trueang_prenoise(halftimeind,ss,kcampind),1);
        
      end
      
      angfmean(:,:,:,kcampind)=  abs(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials);
      angfemean(:,:,:,kcampind)=abs(sum(exp(i*deg2rad(angfreqerpdiff)),1)/numtrials);
      angfpmean(:,:,:,kcampind)=abs(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials);
      angfpemean(:,:,:,kcampind)=abs(sum(exp(i*deg2rad(angfreqpnerpdiff)),1)/numtrials);
      angfang(:,:,:,kcampind)=  rad2deg(angle(sum(exp(i*deg2rad(angfreqdiff)),1)/numtrials));
      angfeang(:,:,:,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqerpdiff)),1)/numtrials));
      angfpang(:,:,:,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpndiff)),1)/numtrials));
      angfpeang(:,:,:,kcampind)=rad2deg(angle(sum(exp(i*deg2rad(angfreqpnerpdiff)),1)/numtrials));
      
    end % kcampind
    
    try
      save(['D:\phase_estimation\anglermsrun9.mat'],'ang*mean','ang*ang','frequse','freqest2','-append')
    catch
      save(['D:\phase_estimation\anglermsrun9.mat'],'ang*mean','ang*ang','frequse','freqest2')
    end
    
    load(['D:\phase_estimation\anglermsrun9.mat'])
    
    for tt=1:2
      figure(140+tt);
      for tf=1:length(t_ftimwin),
        subplot(4,5,tf);       bar(1-squeeze(angfpemean(:,tt,tf,:))')    ;axis([-inf inf 0 0.5])
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}) ' s'])
        else
          title(['TW ' num2str(t_ftimwin{tf}*cfg.foi) ' per.'])
        end
        if tf==1,ylabel('Circ-Var');end
        subplot(4,5,tf+10);    bar(squeeze(angfpeang(:,tt,tf,:))')    ;axis([-inf inf -20 20])
        if tf==1,ylabel('Mean Angle');end
        if tf<6
          title(['TW ' num2str(t_ftimwin{tf}) ' s'])
        else
          title(['TW ' num2str(t_ftimwin{tf}*cfg.foi) ' per.'])
        end
      end
      set(get(140+tt,'Children'),'XTickLabel',{'Kc 1' 'Kc 2' 'Kc 3' })
      %       legend({'cfg.foi 1 Hz' '2 Hz' '3 Hz' '4 Hz'})
    end
  
end % run9