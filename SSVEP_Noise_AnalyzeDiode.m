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

p.fftres = 2^14;
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

t.idx = strcmp(idx.flickertype,'SSVEP');
pl.data = squeeze(DIODE_Meas_EP_cue.data(:,:,t.idx));
figure;
plot(DIODE_Meas_EP_cue.times,pl.data','Color', [0.7 0.7 0.7])
hold on;
plot(DIODE_Meas_EP_cue.times,mean(pl.data,2),'Color', [0.1 0.1 0.1])
xlabel('time in ms')
ylabel('amplitude in a.u.')
xlim(DIODE_Meas_EP_cue.times([1 end]))
title(sprintf('measured SSVEP signals for %1.0f trials', sum(t.idx)))


t.idx = strcmp(idx.flickertype,'BRBF');
pl.data = squeeze(DIODE_Meas_EP_cue.data(:,:,t.idx));
figure;
plot(DIODE_Meas_EP_cue.times,pl.data','Color', [0.7 0.7 0.7])
hold on;
plot(DIODE_Meas_EP_cue.times,mean(pl.data,2),'Color', [0.1 0.1 0.1])
xlabel('time in ms')
ylabel('amplitude in a.u.')
xlim(DIODE_Meas_EP_cue.times([1 end]))
title(sprintf('measured BRBF signals for %1.0f trials', sum(t.idx)))

%% do FFT transforms
DIODE_Meas_EP_cue_stim = pop_select(DIODE_Meas_EP_cue, 'time', [-1 2]);

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
t.idx = strcmp(idx.flickertype,'SSVEP');
figure;
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
t.idx = strcmp(idx.flickertype,'BRBF');
figure;
h.pl1 = plot(res.fft_xscale,res.fft_ind(:,t.idx),'Color', [0.7 0.7 0.7]);
hold on;
h.pl2 = plot(res.fft_xscale,mean(res.fft_ind(:,t.idx),2),'Color', [0.1 0.1 0.1]);
h.pl3 = plot(res.fft_xscale,res.fft_evo(:,2),'Color', [0.8 0.1 0.1]);
xlabel('frequency in Hz')
ylabel('amplitude in a.u.')
xlim([0 130])
title(sprintf('measured BRBF spectra for %1.0f trials', sum(t.idx)))
legend([h.pl1(1) h.pl2 h.pl3],{'single trial';'average single trial';'evoked'})
