% Script for analyzing EEG data from SSVEP_Noise_FShift
%
%
% C. Gundlach 2015, 2022, 2023



%% General Definitions 
clearvars
p.path=             'O:\AllgPsy\experimental_data\2022_SSVEP_Noise_FShift\';
p.bdf_path=         [p.path 'eeg\raw\'];
p.set_path=         [p.path 'eeg\set\'];
p.epoch_path=       [p.path 'EEG\epoch\'];
p.scads_path=       [p.path 'eeg\SCADS\'];
p.chanlocs_path=    ['C:\Users\psy05cvd\Dropbox\work\matlab\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_1020.epf'];
p.mean_path=        [p.path 'eeg\mean\'];
p.exp_name=         'SSVEP_Noise_FShift';
% p.subs=             arrayfun(@(x) sprintf('%02.0f',x),1:40,'UniformOutput',false)';
p.subs=             arrayfun(@(x) sprintf('%02.0f',x),101,'UniformOutput',false)'; % pilot
p.subs2use=         [1];%
p.part=             {'_2';'_3'};
p.events =          {[1 101 111 121 201 211 221];[2 102 112 122 202 212 222];...
                    [3 103 113 123 203 213 223];[4 104 114 124 204 214 224]};
% p.events=           {100; 200};
p.epoch=            [-2.5 2.5];
p.epoch2an=         [-1 2];
p.resample=         512;

p.AnaScriptName=    'SSVEP_Noise_FShift';

% flags
ImportFlag=         1; % set to 1 if files have ti be importet from raw .bdf files
EpochFlag=          1; % set to 1 if data has to be epoched

if ImportFlag==0 && EpochFlag==0 help(AnaScriptName),return, end

%% Main Script 
% loop for subjects
for i_sub=1:numel(p.subs2use)
    FileName=sprintf('VP%s',p.subs{p.subs2use(i_sub)});
        
    %% import
    if ImportFlag    %% Import Files and merge %%
        % load file
%         EEG=pop_biosig([p.bdf_path FileName '.bdf']);
        % load data
        temp.files = dir(sprintf('%sVP%s*.bdf',p.bdf_path,p.subs{p.subs2use(i_sub)}));
        
        for i_fi = 1:numel(temp.files)
            EEG(i_fi)=pop_readbdf(sprintf('%s%s',p.bdf_path,temp.files(i_fi).name),[],73,[]);
        end
        if numel(EEG)>1,EEG=pop_mergeset(EEG,1:numel(EEG),0); end
        %pop_eegplot(EEG,1,1,1)
        
        % replace boundary trigger by 999999
        if sum(strcmp({EEG.event.type},'boundary'))~=0
            [EEG.event(strcmp({EEG.event.type},'boundary')).type] = deal('9999999');
            for i_ev = 1:numel(EEG.event)
                EEG.event(i_ev).type = str2num(EEG.event(i_ev).type);
            end
        end
        
%         % troubleshooting
%         if any(unique([EEG.event.type])>16128)
%             t.triggernew = num2cell([EEG.event.type]-16128);
%             [EEG.event.type]=t.triggernew{:};
%             [EEG.urevent.type]=t.triggernew{:};
%         end
        
                
        % resample
        EEG=pop_resample(EEG,p.resample);
        
        %
%         figure; pwelch(EEG.data(29,:),EEG.srate*1,EEG.srate/2,256,EEG.srate)
%         figure; pop_eegplot(EEG,1,1,1)
%         figure; pwelch(EEG.data(45,1320*256:1324*256),EEG.srate*1,EEG.srate/2,256,EEG.srate)
        
        % check trigger
        trigg.all = unique(cell2mat({EEG.event.type}));
        trigg.freq = arrayfun(@(x) sum(cell2mat({EEG.event.type})==x), trigg.all);
        trigg.sum = [trigg.all; trigg.freq];
        trigg.events = cellfun(@(x) sum(ismember(cell2mat({EEG.event.type}), x)), p.events);
        
        
        % save in new format
        if ~exist(p.set_path); mkdir(p.set_path);end
        EEG = pop_saveset(EEG,'filename',[FileName '.set'], 'filepath', p.set_path);
        clear ('EEG');
    end
    
    %% epoch Bipolarize, Remove Baseline, Detrend, Blinks + EyeMovements Threshold, SCADS%%
    if EpochFlag
        % load
        EEG = pop_loadset('filename',[FileName '.set'], 'filepath', p.set_path);
        % pop_eegplot(EEG,1,1,1,1)
        
        % epoch data
        [EEG t.idx] = pop_epoch( EEG, num2cell(unique(cell2mat(p.events))), p.epoch, 'epochinfo', 'yes');
        
        % keep track of discarded trials
        PreProc.trial_nr = t.idx;
        PreProc.trial_con= nan(1,numel(EEG.epoch));
        PreProc.trial_blink=true(1,numel(EEG.epoch));
        PreProc.trial_eyemov=true(1,numel(EEG.epoch));
        PreProc.trial_SCADS=true(1,numel(EEG.epoch));
        for i_ep = 1:numel(EEG.epoch)
            if numel(EEG.epoch(i_ep).eventtype) == 1
                if isstr(cell2mat(EEG.epoch(i_ep).eventtype))
                    PreProc.trial_con(i_ep)=str2num(cell2mat(EEG.epoch(i_ep).eventtype));
                else
                    PreProc.trial_con(i_ep)=cell2mat(EEG.epoch(i_ep).eventtype);
                end
            else
                if isstr(EEG.epoch(i_ep).eventtype{1})
                    t.t=cellfun(@str2num,EEG.epoch(i_ep).eventtype);
                else
                    t.t=cell2mat(EEG.epoch(i_ep).eventtype);
                end
                if isstr(t.t(cell2mat(EEG.epoch(i_ep).eventlatency)==0))
                    PreProc.trial_con(i_ep)=str2num(t.t(cell2mat(EEG.epoch(i_ep).eventlatency)==0));
                else
                    PreProc.trial_con(i_ep)=t.t(cell2mat(EEG.epoch(i_ep).eventlatency)==0);
                end
            end
        end
        
        % create single horizontal and single vertical eye channel by
        % subtracting channels (e.g. VEOG1-VEOG2) from each other to
        % increase SNR of blinks and eye movements
        EEG = eegF_Bipolarize(EEG);
        
        % remove linear drift and offset of EEG signals
        EEG = eegF_Detrend(EEG,p.epoch2an); % 
        
        % index trials with blinks and eye movements for later rejection
        [EEG_temp,Trials2Remove1] = eegF_RemoveBlinks(EEG,65,p.epoch2an,1); % letzte Parameter '1': Artefakte nur markiert
        PreProc.trial_blink(Trials2Remove1)=false;
        [EEG_temp,Trials2Remove2] = eegF_RemoveEyeMovements(EEG,[65 66],p.epoch2an,25,1);
        PreProc.trial_eyemov(Trials2Remove2)=false;
        EEG = pop_select( EEG,'trial',find(PreProc.trial_eyemov& PreProc.trial_blink));
        
        % discard non-brain electrodes
        EEG = pop_select(EEG,'channel',1:64);
        % load channel info
        EEG.chanlocs = pop_chanedit(EEG.chanlocs,'load',{p.chanlocs_path,'filetype','besa'});
        % pop_eegplot(EEG,1,1,1)
        
        % filter before SCADS (is this a good idea?)
%         EEG = pop_eegfiltnew(EEG, 0.7, 0, 16*EEG.srate, 0, [], 0);
        
        % run SCADS to statistically detect artifacts
        [EEG Trials2Del SumData] = SCADS_Pass1(EEG,p.epoch2an,[6 6 15],1);%,[6 5 15]);  % ansonsten SCADS_Pass1()
        t.index = find(PreProc.trial_blink & PreProc.trial_eyemov);
        PreProc.trial_SCADS(t.index(Trials2Del))=false;        
        
        % save trials marked as artifacts for potential later inspection
        EEG_SCADS = pop_select(EEG, 'trial', Trials2Del);
        % discard trials indexed as containing artifacts
        EEG = pop_rejepoch(EEG, EEG.reject.rejmanual, 0);
        % rereference to average reference
        EEG=pop_reref(EEG,[],'refstate',0);   %
        
        if ~exist(p.scads_path); mkdir(p.scads_path);end
        if ~exist(p.epoch_path); mkdir(p.epoch_path);end
        save([p.scads_path FileName '_Preprocess_summary.mat'],'SumData','PreProc')
        pop_saveset(EEG_SCADS,[FileName '_SCADS.set'],p.scads_path);% Save Data
        pop_saveset(EEG, 'filename', [FileName '_e.set'], 'filepath', p.epoch_path);
        clear ('EEG')
    end
end

