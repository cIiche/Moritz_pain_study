clear all
clc
% https://www.fieldtriptoolbox.org/faq/how_can_i_import_my_own_dataformat/ 

% this script uses the fieldtrip toolbox to calculate coherence 
% for the Moritz Lab Pain Study 

% finding fieldtrip toolbox ; It is most convenient to have the addpath and ft_defaults in a script with the name startup.m, which is located in your own MATLAB directory. See this information from MathWorks.
restoredefaultpath
addpath 'C:\Users\Henry\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\FieldTrip'
ft_defaults

% Preprocessing options that you should only use when you are calling FT_PREPROCESSING with
% also the second input argument "data" are
%   cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')

% filePaths = {'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\8.18.21 5 min\', 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 1\8.24.21 r1 5 min', 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 2\8.24.21 r2 5 min'};
% fileNames = {'after stim 0', 'after stim 1', 'after stim 2', 'after stim 3', 'after stim 4', 'after stim 5', 'after stim 6', 'after stim 7', 'after stim 8', 'after stim 9'};
% cohmat = zeros(3,10);  
% 
% whichrat = input("What rat would you like to run? ('1' = 8/18 r1, '2' = 8/24 r2, '3' = 8/24 r3): ") ; 
% filePath = filePaths{whichrat} ;
% % for 8/18 
% set_channels=[1 2 3 4 5 6 7] ;
% RIL=set_channels(2);LBLA=set_channels(4);RBLA=set_channels(1);LIL=set_channels(3);stim=set_channels(5) ;
% fs = 20000 ;
% for i = 1:3
%     fileName = fileNames{i} ;
    
%% opening file 
% filePath = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\8.18.21 5 min\';
% filePath = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 1\8.24.21 r1 5 min\';
filePath = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 2\8.24.21 r2 5 min\'; 
fileName= 'last 5 min baseline' ;
load([filePath, fileName]);

set_channels=[1 2 3 4 5 6 7] ;
% % 8/18  
% RIL=set_channels(2);LBLA=set_channels(4);RBLA=set_channels(1);LIL=set_channels(3);stim=set_channels(5) ;
% % 8/24 r1 
% RIL=set_channels(4);LBLA=set_channels(1);RBLA=set_channels(2);LIL=set_channels(3);stim=set_channels(5) ;
% 8/24 r2
RIL=set_channels(2);LBLA=set_channels(1);RBLA=set_channels(4);LIL=set_channels(3);stim=set_channels(5) ;
fs = 20000 ;

timeax=1:dataend(1); %set time axis
time=timeax/fs/60;%frames to seconds to minutes (these are the time values for each data point)
timesec=timeax./fs;
tottime=length(timeax)./fs./60; % total experiment block length in minutes 
%% Organize data into structure array

alldata=[]; %initialize structure array (alldata is a struct)

alldata.RILdata=data(datastart(RIL):dataend(RIL)); % Call different fields as StructName.fieldName-> Struct is alldata and field is S1dataR
alldata.LBLAdata=data(datastart(LBLA):dataend(LBLA));
alldata.RBLAdata=data(datastart(RBLA):dataend(RBLA));
alldata.LILdata=data(datastart(LIL):dataend(LIL));
alldata.stimdata=data(datastart(stim):dataend(stim));

%% making data structure fieldtrip compatible 
% data.label      % cell-array containing strings, Nchan*1
% data.fsample    % sampling frequency in Hz, single number
% data.trial      % cell-array containing a data matrix for each
%                 % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
% data.time       % cell-array containing a time axis for each
%                 % trial (1*Ntrial), each time axis is a 1*Nsamples vector
% data.trialinfo  % this field is optional, but can be used to store
%                 % trial-specific information, such as condition numbers,
%                 % reaction times, correct responses etc. The dimensionality
%                 % is Ntrial*M, where M is an arbitrary number of columns.
% data.sampleinfo % optional array (Ntrial*2) containing the start and end
%                 % sample of each trial
eegdata=[];
eegdata.label = {'RBLA' ; 'RIL'; 'LIL'; 'LBLA'} ; 
eegdata.fsample = 20000;
% for 5 min chopped 
eegdatamat = [alldata.RBLAdata; alldata.RILdata; alldata.LILdata; alldata.LBLAdata];
eegdata.trial = {eegdatamat} ;
% for whole experiment
% fullsamp = length(RBLAdata) ;
% trialsamp = 
% data.trial
eegdata.time = {time}; % in minutes 
eegdata.sampleinfo = [1 length(alldata.RBLAdata)] ; 

%% downsampling 
% Use as
%   [data] = ft_resampledata(cfg, data)
%
% The data should be organised in a structure as obtained from the FT_PREPROCESSING
% function. The configuration should contain
%   cfg.resamplefs      = frequency at which the data will be resampled (default = 256 Hz)
%   cfg.detrend         = 'no' or 'yes', detrend the data prior to resampling (no default specified, see below)
%   cfg.demean          = 'no' or 'yes', whether to apply baseline correction (default = 'no')
%   cfg.baselinewindow  = [begin end] in seconds, the default is the complete trial (default = 'all')
%   cfg.feedback        = 'no', 'text', 'textbar', 'gui' (default = 'text')
%   cfg.trials          = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.sampleindex     = 'no' or 'yes', add a channel with the original sample indices (default = 'no')
% cfg.resamplefs = 108 ;
cfg.resamplefs = 30 ;
cfg.detrend = 'yes' ;
cfg.demean = 'no' ;
cfg.feedback = 'text' ;
cfg.trials = 'all';
cfg.sampleindex = 'no' ;
% cfg.sampleinfo = [1 length(alldata.RBLAdata)] ; 
data = ft_resampledata(cfg, eegdata) ;
data.sampleinfo = [1 length(data.time{1})] ; 

%% preproccesing  

% cfg = [];
% % bandpass filter 
% cfg.bpfilter = 'yes';
% % bandpass frequency range [lowFreq highFreq] in Hz
% cfg.bpfreq = [5 8]; 
% % cfg.datafile = 'after stim 0.mat' ;
% % cfg.headerfile = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\8.18.21 5 min\';
% % cfg.dataset = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\8.18.21 5 min\after stim 0.mat\' ;
% % cfg.dataset = alldata ;
% cfg.dataset = eegdata ;
% data_eeg = ft_preprocessing(cfg);

%% ft_freqanalysis for the computation of the cross- and power spectra

cfg            = [];
cfg.output     = 'powandcsd';
cfg.method     = 'mtmfft';
%   cfg.foilim     = [begin end], frequency band of interest
cfg.foilim     = [5 8];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'no';
% cfg.channel    = {'RBLA' 'RIL'};
cfg.channelcmb = {'RBLA' 'RIL'};
cfg.pad        = 'nextpow2';
freq           = ft_freqanalysis(cfg, data);

cfg            = [];
cfg.method     = 'coh';
cfg.channelcmb = {'RBLA' 'RIL'};
fd             = ft_connectivityanalysis(cfg, freq);

% cohmat(whichrat, i) = median(fd.cohspctrm) ;
% %% plotting coherence with ft_singleplotER 

% cfg                  = [];
% cfg.parameter        = 'cohspctrm';
% cfg.xlim             = 'maxmin';
% cfg.refchannel       = 'RIL';
% cfg.showlabels       = 'yes';
% % figure; ft_multiplotER(cfg, fd)
% % cfg.channel = 'MRC21';
% figure; 
% ft_singleplotER(cfg, fd);
% ylabel('Coherence')
% xlabel('Frequency (Hz)')
% title('Connectivity between RIL & RBLA with reference to RIL')
a = median(fd.cohspctrm);
a
% end 