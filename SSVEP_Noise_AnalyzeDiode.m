% script to try out some data analysis techniques for broadband-noise-flicker and classical SSVEP signals
%
% broadband noise as described in
% Zhigalov, A. & Jensen, O. Alpha oscillations do not implement gain control in early visual cortex but rather gating in 
% parieto-occipital regions. Hum. Brain Mapp. 2020.04.03.021485 (2020) doi:10.1101/2020.04.03.021485

clearvars
%% initial parameters
p.path_bdf = 'O:\AllgPsy\experimental_data\2022_SSVEP_Noise_FShift\diode_data\DiodeTest.bdf';
p.path_diodeTiming = 'O:\AllgPsy\experimental_data\2022_SSVEP_Noise_FShift\diode_data\diode_timing.mat';

p.trig.trial_start = 88;
p.trig.trial_end = 89;
p.trig.cue = [[1 2 3 4], [1 2 3 4]+100,  [1 2 3 4]+200,  [1 2 3 4]+110,  [1 2 3 4]+120,  [1 2 3 4]+210,  [1 2 3 4]+220];

p.SSVEP_freq = 63;
p.BRBF_freq = 65;

p.fftres = 2^14;

p.filt_SSVEP = [60 65];
p.filt_Noise = [45 85];
p.plv_lagrange = [-200 200];
%% read in data
DIODE_Meas = pop_readbdf(p.path_bdf, [] ,2,[]);
% pop_eegplot(DIODE_Meas,1,1,1)

% read in stimulation logfile
t.file = open(p.path_diodeTiming);
DIODE_Stimlog = [t.file.resp.experiment{:}];

%% check timing and stuff
% extract timing (end MINUS start trigger) for each trial from measured data
timing.trialtime_meas = ([DIODE_Meas.event([DIODE_Meas.event.type] == p.trig.trial_end).latency] - ...
    [DIODE_Meas.event([DIODE_Meas.event.type] == p.trig.trial_start).latency]) ./ DIODE_Meas.srate*1000; % what's the measured trial time in ms

% extract timing (end MINUS start trigger) for each trial from log file ()
timing.trialtime_def = ([DIODE_Stimlog.pre_cue_times] + [DIODE_Stimlog.post_cue_frames]).*1000; % cave: bug in stimulation script

% plot results
figure; plot(timing.trialtime_meas-timing.trialtime_def)
hold on; plot([0 numel(timing.trialtime_meas)], [1 1]*1/DIODE_Meas.srate*1000,'r:')
plot([0 numel(timing.trialtime_meas)], [-1 -1]*1/DIODE_Meas.srate*1000,'r:')
xlabel('trial number')
ylabel('diff of triallength t_m_e_a_s - t_d_e_f in ms')

%% epoch data
% whole trial
DIODE_Meas_EP_wt = pop_epoch( DIODE_Meas, num2cell(p.trig.trial_start), [0 (max(timing.trialtime_meas)+100)/1000], 'epochinfo', 'yes');
% pop_eegplot(DIODE_Meas_EP_wt,1,1,1)

DIODE_Meas_EP_cue = pop_epoch( DIODE_Meas, num2cell(p.trig.cue), [-1 2.1], 'epochinfo', 'yes');
% pop_eegplot(DIODE_Meas_EP_cue,1,1,1)

%% some initial plotting
% index flicker type
idx.flickertype = [DIODE_Stimlog.flickertype]; idx.flickertype = idx.flickertype(1:2:end);

figure;
set(gcf,'Position',[100 100 1200 300],'PaperPositionMode','auto')
% plot single-trial stimulation signals
% first SSVEP trial
subplot(1,2,1)
pl.ydata = DIODE_Stimlog(find(strcmp(idx.flickertype,'SSVEP'),1,'first')).lummat;
pl.xdata = ((1:size(pl.ydata,2))-DIODE_Stimlog(find(strcmp(idx.flickertype,'SSVEP'),1,'first')).pre_cue_frames-1)./480*1000;
% plot(pl.xdata,pl.ydata+[0; 0.9])
h.pl = plot(pl.xdata,pl.ydata);
xlim([-100 400])
ylim([-0.1 1.1])
h.pl(1).Color = [1 0.4 0]; h.pl(2).Color = [0 0.4 1];
ylabel('luminance')
xlabel('time in ms relative to cue')
legend({'SSVEP 1';'SSVEP 2'},'Location','SouthOutside','Orientation','horizontal')

subplot(1,2,2)
pl.ydata = DIODE_Stimlog(find(strcmp(idx.flickertype,'BRBF'),1,'first')).lummat;
pl.xdata = ((1:size(pl.ydata,2))-DIODE_Stimlog(find(strcmp(idx.flickertype,'BRBF'),1,'first')).pre_cue_frames-1)./480*1000;
% plot(pl.xdata,pl.ydata+[0; 0.9])
h.pl = plot(pl.xdata,pl.ydata);
xlim([-100 400])
ylim([-0.1 1.1])
h.pl(1).Color = [1 0.4 0]; h.pl(2).Color = [0 0.4 1];
ylabel('luminance')
xlabel('time in ms relative to cue')
legend({'BRBF 1';'BRBF 2'},'Location','SouthOutside','Orientation','horizontal')

% SaveCurrentFigure([pwd '\figure' ], 'single_trial_luminance')


% plot recorded data
figure;
set(gcf,'Position',[100 100 1200 400],'PaperPositionMode','auto')
subplot(1,2,1)
t.idx = strcmp(idx.flickertype,'SSVEP');
pl.data = squeeze(DIODE_Meas_EP_cue.data(:,:,t.idx));
h.pl1 = plot(DIODE_Meas_EP_cue.times,pl.data','Color', [0.7 0.7 0.7]);
hold on;
h.pl2 = plot(DIODE_Meas_EP_cue.times,mean(pl.data,2),'Color', [0.1 0.1 0.1]);
xlabel('time in ms relative to cue')
ylabel('amplitude in a.u.')
xlim(DIODE_Meas_EP_cue.times([1 end]))
title(sprintf('measured SSVEP signals for %1.0f trials', sum(t.idx)))
legend([h.pl1(1) h.pl2], {'single trial';'average'},'Location','SouthOutside','Orientation','horizontal')

subplot(1,2,2)
t.idx = strcmp(idx.flickertype,'BRBF');
pl.data = squeeze(DIODE_Meas_EP_cue.data(:,:,t.idx));
h.pl1 = plot(DIODE_Meas_EP_cue.times,pl.data','Color', [0.7 0.7 0.7]);
hold on;
h.pl2 = plot(DIODE_Meas_EP_cue.times,mean(pl.data,2),'Color', [0.1 0.1 0.1]);
xlabel('time in ms relative to cue')
ylabel('amplitude in a.u.')
xlim(DIODE_Meas_EP_cue.times([1 end]))
title(sprintf('measured BRBF signals for %1.0f trials', sum(t.idx)))
legend([h.pl1(1) h.pl2], {'single trial';'average'},'Location','SouthOutside','Orientation','horizontal')

% SaveCurrentFigure([pwd '\figure' ], 'measured_luminance')

%% do FFT transforms + plotting
DIODE_Meas_EP_cue_stim = pop_select(DIODE_Meas_EP_cue, 'time', [-1 2]);

% index flicker type
idx.flickertype = [DIODE_Stimlog.flickertype]; idx.flickertype = idx.flickertype(1:2:end);

% induced
t.data = squeeze(DIODE_Meas_EP_cue_stim.data);
res.fft_ind =  squeeze(abs(fft(detrend(t.data),p.fftres,1))*2/size(t.data,1));
res.fft_xscale = ((0:size(res.fft_ind,1)-1)/size(res.fft_ind,1)) * DIODE_Meas_EP_cue_stim.srate;

% evoked
t.idx = strcmp(idx.flickertype,'SSVEP');
t.data = squeeze(mean(DIODE_Meas_EP_cue_stim.data(:,:,t.idx),3))';
res.fft_evo =  squeeze(abs(fft(detrend(t.data),p.fftres,1))*2/size(t.data,1));

t.idx = strcmp(idx.flickertype,'BRBF');
t.data = squeeze(mean(DIODE_Meas_EP_cue_stim.data(:,:,t.idx),3))';
res.fft_evo(:,2) =  squeeze(abs(fft(detrend(t.data),p.fftres,1))*2/size(t.data,1));

% plot data
figure;
set(gcf,'Position',[100 100 1200 400],'PaperPositionMode','auto')
subplot(1,2,1)
t.idx = strcmp(idx.flickertype,'SSVEP');
h.pl1 = plot(res.fft_xscale,res.fft_ind(:,t.idx),'Color', [0.7 0.7 0.7]);
hold on;
h.pl2 = plot(res.fft_xscale,mean(res.fft_ind(:,t.idx),2),'Color', [0.1 0.1 0.1]);
h.pl3 = plot(res.fft_xscale,res.fft_evo(:,1),'Color', [0.8 0.1 0.1]);
xlabel('frequency in Hz')
ylabel('amplitude in a.u.')
xlim([0 130])
title(sprintf('measured SSVEP spectra for %1.0f trials', sum(t.idx)))
legend([h.pl1(1) h.pl2 h.pl3],{'single trial';'average single trial';'evoked'})

% plot data
subplot(1,2,2)
t.idx = strcmp(idx.flickertype,'BRBF');
h.pl1 = plot(res.fft_xscale,res.fft_ind(:,t.idx),'Color', [0.7 0.7 0.7]);
hold on;
h.pl2 = plot(res.fft_xscale,mean(res.fft_ind(:,t.idx),2),'Color', [0.1 0.1 0.1]);
h.pl3 = plot(res.fft_xscale,res.fft_evo(:,2),'Color', [0.8 0.1 0.1]);
xlabel('frequency in Hz')
ylabel('amplitude in a.u.')
xlim([0 130])
title(sprintf('measured BRBF spectra for %1.0f trials', sum(t.idx)))
legend([h.pl1(1) h.pl2 h.pl3],{'single trial';'average single trial';'evoked'})

% SaveCurrentFigure([pwd '\figure' ], 'amplitude_spectra_measured_luminance')

%% calculate single trial coherence/PLV
% resample EEG signal to have same sampling rate as Projector?
DIODE_Meas_rs = pop_resample(DIODE_Meas,480);

% filter data for SSVEP
% DIODE_Meas_rsfS = pop_eegfiltnew(DIODE_Meas_rs, p.filt_SSVEP(1), p.filt_SSVEP(2), 16*DIODE_Meas_rs.srate, 0, [], 0);
DIODE_Meas_rsfS = DIODE_Meas_rs;
DIODE_Meas_rsfS.data = filterFGx(DIODE_Meas_rsfS.data,DIODE_Meas_rsfS.srate,p.SSVEP_freq,5,0);
% pop_eegplot(DIODE_Meas_rsfS,1,1,1)

% filter data for Noise
% DIODE_Meas_rsfN = pop_eegfiltnew(DIODE_Meas_rs, p.filt_Noise(1), p.filt_Noise(2), 16*DIODE_Meas_rs.srate, 0, [], 0);
DIODE_Meas_rsfN = DIODE_Meas_rs;
DIODE_Meas_rsfN.data = filterFGx(DIODE_Meas_rsfN.data,DIODE_Meas_rsfS.srate,p.BRBF_freq,30,0);

% epoch data
DIODE_Meas_rsfSep = pop_epoch(DIODE_Meas_rsfS, num2cell(p.trig.cue), [p.plv_lagrange(1)/1000 (2000+p.plv_lagrange(2))/1000], 'epochinfo', 'yes');
DIODE_Meas_rsfNep = pop_epoch(DIODE_Meas_rsfN, num2cell(p.trig.cue), [p.plv_lagrange(1)/1000 (2000+p.plv_lagrange(2))/1000], 'epochinfo', 'yes');
% pop_eegplot(DIODE_Meas_rsfSep,1,1,1)

% based on these calculate PLVs for different lags between measured signal and define diode signal
res.PLV.lag_range = [-200 200];
res.PLV.lag_rangeidx = dsearchn(DIODE_Meas_rsfSep.times',res.PLV.lag_range');
res.PLV.lag_idx = res.PLV.lag_rangeidx(1):res.PLV.lag_rangeidx(2);

% index flicker type
idx.flickertype = [DIODE_Stimlog.flickertype]; idx.flickertype = idx.flickertype(1:2:end);
idx.SSVEP = find(strcmp(idx.flickertype,'SSVEP'));
idx.BRBF = find(strcmp(idx.flickertype,'BRBF'));

% loop across lags
res.lagged_PLV_SSVEP_data = nan(DIODE_Stimlog(idx.SSVEP(1)).post_cue_times, numel(res.PLV.lag_idx),2);
res.lagged_PLV_SSVEP_time = nan(DIODE_Stimlog(idx.SSVEP(1)).post_cue_times, numel(res.PLV.lag_idx));
fprintf('\ncalculating phase locking value | SSVEP data ...')
for i_lag = 1:numel(res.PLV.lag_idx)
    % do SSVEP phase coherence
    t.PLV_tr = nan(DIODE_Stimlog(idx.SSVEP(1)).post_cue_times,numel(idx.SSVEP),2); % timepoint X trials X real/control
    for i_tr = 1:numel(idx.SSVEP)
        % diode signal from cue onwards
        t.ydata = DIODE_Stimlog(idx.SSVEP(i_tr)).lummat(1,DIODE_Stimlog(idx.SSVEP(i_tr)).pre_cue_frames+1:end);
        % filter signal
        t.ydataf =filterFGx(t.ydata,DIODE_Meas_rsfS.srate,p.SSVEP_freq,5,0);
        t.ydata_hilb = hilbert(t.ydataf);
        % figure; plot(angle(t.ydata_hilb))
        
        t.xdata = DIODE_Meas_rsfSep.data(1,(1:numel(t.ydata_hilb))+i_lag-1,idx.SSVEP(i_tr));
        t.xdata_hilb = hilbert(t.xdata);
        % figure; plot(angle(t.xdata_hilb))
        
        % calculate PLV for real signal
        t.PLV_tr(:,i_tr,1) = exp(1i*(angle(t.xdata_hilb) - angle(t.ydata_hilb)));
        
        % trouble shooting data
%         figure; plot(angle(t.ydata_hilb)); hold on; plot(angle(t.xdata_hilb))
%         figure; plot((t.ydataf)*0.6*10^6); hold on; plot(t.xdata)
%         figure; plot((angle(t.xdata_hilb) - angle(t.ydata_hilb)))
        
        % pseudo PLV witg BRBF (no locking expected!)
        % diode signal from cue onwards for random BRBF trial
        t.idx = randsample(idx.BRBF,1);
        t.ydata_r = DIODE_Stimlog(t.idx).lummat(1,DIODE_Stimlog(t.idx).pre_cue_frames+1:end);
        % filter signal
        t.ydataf_r =filterFGx(t.ydata_r,DIODE_Meas_rsfS.srate,p.BRBF_freq,30,0);
        t.ydata_hilb_r = hilbert(t.ydataf_r);
        % figure; plot(angle(t.ydata_hilb_r))
        
        % calculate PLV for pseudo signal
        t.PLV_tr(:,i_tr,2) = exp(1i*(angle(t.xdata_hilb) - angle(t.ydata_hilb_r)));
%         figure; plot((angle(t.xdata_hilb) - angle(t.ydata_hilb_r)))
    end
    % extract time value
    res.lagged_PLV_time(:,i_lag) = DIODE_Meas_rsfSep.times((1:numel(t.ydata_hilb))+i_lag-1);
    
    % extract data
    res.lagged_PLV_SSVEP_data(:,i_lag,:) = abs(sum(t.PLV_tr,2))/size(t.PLV_tr,2);
    % figure; plot(res.lagged_PLV_time(:,i_lag), squeeze(res.lagged_PLV_SSVEP_data(:,i_lag,:)))
end
fprintf('...done!\n')



% for BRBF signals
% loop across lags
res.lagged_PLV_BRBF_data = nan(DIODE_Stimlog(idx.BRBF(1)).post_cue_times, numel(res.PLV.lag_idx),3);
res.lagged_PLV_BRBF_time = nan(DIODE_Stimlog(idx.BRBF(1)).post_cue_times, numel(res.PLV.lag_idx));
fprintf('\ncalculating phase locking value | BRBF data ...')
for i_lag = 1:numel(res.PLV.lag_idx)
    % do SSVEP phase coherence
    t.PLV_tr = nan(DIODE_Stimlog(idx.BRBF(1)).post_cue_times,numel(idx.BRBF),3); % timepoint X trials X real/control1/control2
    for i_tr = 1:numel(idx.BRBF)
        % diode signal from cue onwards
        t.ydata = DIODE_Stimlog(idx.BRBF(i_tr)).lummat(1,DIODE_Stimlog(idx.BRBF(i_tr)).pre_cue_frames+1:end);
        % filter signal
        t.ydataf =filterFGx(t.ydata,DIODE_Meas_rsfN.srate,p.BRBF_freq,30,0);
        t.ydata_hilb = hilbert(t.ydataf);
        % figure; plot(angle(t.ydata_hilb))
        
        t.xdata = DIODE_Meas_rsfNep.data(1,(1:numel(t.ydata_hilb))+i_lag-1,idx.BRBF(i_tr));
        t.xdata_hilb = hilbert(t.xdata);
        % figure; plot(angle(t.xdata_hilb))
        
        % calculate PLV for real signal
        t.PLV_tr(:,i_tr,1) = exp(1i*(angle(t.xdata_hilb) - angle(t.ydata_hilb)));
        
        % trouble shooting data
%         figure; plot(angle(t.ydata_hilb)); hold on; plot(angle(t.xdata_hilb))
%         figure; plot((t.ydataf)*0.6*10^6); hold on; plot(t.xdata)
%         figure; plot((angle(t.xdata_hilb) - angle(t.ydata_hilb)))
        
        % control 1
        % pseudo PLV with other trial BRBF (no locking expected!)
        % diode signal from cue onwards for random BRBF trial
        t.idx = randsample(idx.BRBF,1);
        t.ydata_r = DIODE_Stimlog(t.idx).lummat(1,DIODE_Stimlog(t.idx).pre_cue_frames+1:end);
        % filter signal
        t.ydataf_r =filterFGx(t.ydata_r,DIODE_Meas_rsfN.srate,p.BRBF_freq,30,0);
        t.ydata_hilb_r = hilbert(t.ydataf_r);
        % figure; plot(angle(t.ydata_hilb_r))
        
        % calculate PLV for pseudo signal
        t.PLV_tr(:,i_tr,2) = exp(1i*(angle(t.xdata_hilb) - angle(t.ydata_hilb_r)));
%         figure; plot((angle(t.xdata_hilb) - angle(t.ydata_hilb_r)))

        % control 2
        % pseudo PLV with other trial SSVEP (no locking expected!)
        % diode signal from cue onwards for random SSVEP trial
        t.idx = randsample(idx.SSVEP,1);
        t.ydata_r2 = DIODE_Stimlog(t.idx).lummat(1,DIODE_Stimlog(t.idx).pre_cue_frames+1:end);
        % filter signal
        t.ydataf_r2 =filterFGx(t.ydata_r2,DIODE_Meas_rsfS.srate,p.SSVEP_freq,5,0);
        t.ydata_hilb_r2 = hilbert(t.ydataf_r2);
        % figure; plot(angle(t.ydata_hilb_r2))
        
        % calculate PLV for pseudo signal
        t.PLV_tr(:,i_tr,3) = exp(1i*(angle(t.xdata_hilb) - angle(t.ydata_hilb_r2)));
%         figure; plot((angle(t.xdata_hilb) - angle(t.ydata_hilb_r2)))
    end
    % extract time value
    res.lagged_PLV_time(:,i_lag) = DIODE_Meas_rsfSep.times((1:numel(t.ydata_hilb))+i_lag-1);
    
    % extract data
    res.lagged_PLV_BRBF_data(:,i_lag,:) = abs(sum(t.PLV_tr,2))/size(t.PLV_tr,2);
    % figure; plot(res.lagged_PLV_time(:,i_lag), squeeze(res.lagged_PLV_BRBF_data(:,i_lag,:)))
end
fprintf('...done!\n')

% plot results | SSVEP data
figure;
set(gcf,'Position',[100 100 700 700],'PaperPositionMode','auto')
subplot(2,1,1)
plot(res.lagged_PLV_time(1,:), squeeze(mean(res.lagged_PLV_SSVEP_data(:,:,:),1)))
xlabel('lag in ms')
ylabel('PLV')
ylim([0 1.1])
legend({'real effect';'random effect'}, 'Location','SouthOutside','Orientation','horizontal')
title(sprintf('SSVEP | avg. PLVs between presented and recorded signal for different lags | N = %1.0f trials', numel(idx.SSVEP)))

% plot results | BRBF data
subplot(2,1,2)
plot(res.lagged_PLV_time(1,:), squeeze(mean(res.lagged_PLV_BRBF_data(:,:,:),1)))
[t.m t.idx] = max(squeeze(mean(res.lagged_PLV_BRBF_data(:,:,1),1)));
text(res.lagged_PLV_time(1,t.idx), t.m,sprintf('| max at %1.3f ms', res.lagged_PLV_time(1,t.idx)))
xlabel('lag in ms')
ylabel('PLV')
ylim([0 1.1])
legend({'real effect';'random effect | other BRBF'; 'random effect | SSVEP'}, 'Location','SouthOutside','Orientation','horizontal')
title(sprintf('BRBF | avg. PLVs between presented and recorded signal for different lags | N = %1.0f trials', numel(idx.BRBF)))

% SaveCurrentFigure([pwd '\figure' ], 'PLV_measured_stimulated_luminance')