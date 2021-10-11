%% Authors: Hannah Mach, Kat Floerchinger, Henry Tan

% This script plots CWTs (IN PROGRESS) 
% from pain study rat data 

%% 
clear 
close all
clc
% % %% load data

% filePath = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\8_11_21 rat\';
% fileName ='Trial 1' ;
% fileName ='4Hz_paincircuit10ms' ;
% fileName ='4Hz_paincircuit300mV' ;
% fileName ='4Hz_toBAsleep10ms' ;
% fileName ='4Hz_toBAsleep100ms' ;

filePath = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_16_2021 rat\';
% fileName ='Trial 1' ;
fileName= '4Hz_sleep_study_IL_before_BLA_200mV' ; 

% filePath = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\';
% fileName ='Trial 2' ;

% filePath = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 1\' ; 
% fileName ='Trial 1' ;

% filePath = 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 2\' ; 
% fileName ='Trial 1' ;

% file_list=dir([folder 'Trial 1.mat']);

load([filePath, fileName]);

%% set parameters
% for 8/11
% set_channels=[1 2 3 4 6] ;
% for 8/16, 8/18 
set_channels=[1 2 3 4 5 6 7] ;
% 8/24/21 rat 1 and rat 2
% set_channels=[1 2 3 4 5 6 7] ;

% set channel identities, stim should always be on index 5, no matter what channel 
decision = input('Do you want to manually input channels?(1 = yes, 0 = no):') ;
if decision == 1
    RIL = input('What channel is right IL?:') ;
    LBLA = input('What channel is left BLA?:') ;
    RBLA = input('What channel is right BLA?:') ;
    LIL = input('What channel is left IL?:'); 
    STIM = input('What channel is stimulus?:'); 
    RIL=set_channels(RIL);LBLA=set_channels(LBLA);RBLA=set_channels(RBLA);LIL=set_channels(LIL);stim=set_channels(STIM) ; 

end 
if decision == 0  
    % 8/11/21 
%     RIL=set_channels(1);LBLA=set_channels(4);RBLA=set_channels(3);LIL=set_channels(2);stim=set_channels(5) ;
    % 8/16/21
    RIL=set_channels(2);LBLA=set_channels(3);RBLA=set_channels(1);LIL=set_channels(4);stim=set_channels(5) ;
%     8/18/21
%     RIL=set_channels(2);LBLA=set_channels(4);RBLA=set_channels(1);LIL=set_channels(3);stim=set_channels(6) ;
%     8/24/21 rat 1
%     RIL=set_channels(4);LBLA=set_channels(1);RBLA=set_channels(2);LIL=set_channels(3);stim=set_channels(6) ;
%     8/24/21 rat 2
%     RIL=set_channels(2);LBLA=set_channels(1);RBLA=set_channels(4);LIL=set_channels(3);stim=set_channels(6) ;
end 

%% plot CWT or timeseries 

plotter = input("Plot CWT = '1', STA timeseries = '2': ") ; 
plotter2 = input("Coherence analysis? '1' = yes, '2' =no:  ") ; 
%% Set sampling rate
% fs = input('What is the tickrate/sampling rate?:') ;
fs = 20000 ;

%%
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


% % make bipolar channels
alldata.RILLILbipolardata=alldata.RILdata-alldata.LILdata;
alldata.RBLALBLAbipolardata=alldata.RBLAdata-alldata.LBLAdata;

% names={'RSdata','LSdata','RHdata','LHdata', 'stimdata', 'LHRSbipolardata', 'LHRHbipolardata'};
names={'RILdata','LBLAdata','RBLAdata','LILdata', 'stimdata', 'RILLILbipolardata', 'RBLALBLAbipolardata'};
%% plot raw data
figure
for i=1:length(names)
    subplot(length(names),1,i)
    plot(time,alldata.(char(names(i)))) % this is plotting time in minutes, but if you want seconds, use timesec instead of time
    title(names(i));
end
xlabel('time (minutes)')

% %%  Check frequency components of raw data 
% signal=(alldata.LHRSbipolardata);
% fftV1=fft(signal);
% figure
% plot(linspace(0,fs,length(fftV1)),abs(fftV1));
% xlim([0 100]);
% ylim([0 6]);
%% Filter the raw data with a lowpass butterworth filter and condition removing points > 4std from data 
%Note: it's better to filter the raw data before STAs because this limits
%the end effects that are inevitably created by the filter (you'd get end
%effects at the ends of each STA, vs end effects only at the beginning and
%end of the time series data

% lowEnd = 1; % Hz
% highEnd = 50; % Hz
lowEnd = 2; % Hz
highEnd = 100; % Hz
filterOrder = 3; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
[b, a] = butter(filterOrder, [lowEnd highEnd]/(fs/2)); % Generate filter coefficients
% [b, a] = butter(filterOrder, [lowEnd highEnd]/(fs/4)); % Generate filter coefficients

for ii=1:length(names)
filteredData.(char(names(ii))) = filtfilt(b, a,alldata.(char(names(ii)))); % Apply filter to data using zero-phase filtering
% is this correct here? before STAs? - in vis stim, for stats_analysis is
% filtered like this, after STA creation 
deviation=std(filteredData.(char(names(ii))));
trialmean = mean(filteredData.(char(names(ii))));
filteredData.(char(names(ii))) = filteredData.(char(names(ii)))(filteredData.(char(names(ii)))<trialmean+4*deviation);
end

%% detect stimuli
index_stim=[];%initialize reference array of stimulus onsets for STAs
index_allstim=[];%secondary array of stimuli for identifying first pulse of a train. 

X=alldata.stimdata;
X=X-min(X);
X=X/max(X);
Y=X>0.5;
%Y=X>0.04;used during debugging, works as well
Z=diff(Y);
% index_allstim=find(Z>0.5);
% index_allstim=index_allstim+1;
index_allstim=find(Z>0.04);
index_allstim=index_allstim+1;
% abs(Z>

%find first pulse of each train, if stimulation contains trains
index_trains=diff(index_allstim)>20000;
index_allstim(1)=[];
index_stim=index_allstim(index_trains);


%% Create STAs
for i=1:length(names) %initiate data array to hold STAs
stas.(char(names(i)))=[];
end

tb=1; %time before stim to start STA
ta=9; %time after stim to end STA

for i=1:length(names)
    for j=2:(length(index_stim)-1) %cycle through stimuli
        stas.(char(names(i)))=[stas.(char(names(i))); filteredData.(char(names(i)))((index_stim(j)-fs*tb):(index_stim(j)+fs*ta))];
    end
end

%% plot filtered STAS
figure
x2=1:length(stas.(char(names(1)))); % make an x axis
x2=x2/fs-1; % convert x axis from samples to time 
for i=1:(length(names)-1)
    subplot(length(names)-1,1,i)
    a=median(stas.(char(names(i))));
    a=a-median(stas.(char(names(i)))(1:fs));
    a=a/100*1000;
    asmooth=smoothdata(a); % smooth the data 
    plot(x2,asmooth,'linewidth',1.5)
    ylim([-0.01 0.01])
    ylabel(names(i))
end
subplot(4,1,4)
xlabel('time after stimulus onset (s)')

%% plot filtered CWTs of STAs or STA timeseries 
ticks=[0:.005:.1];
% ticks=[0:0.01:.1];
clear yticks
clear yticklabels

str = string(names) ;
    for i=1:length(names)    
%     for i=1:length(names)  
        figure
        caxis_track=[];
        %ylabels={'V1bipolar (Hz)';'S1A (Hz)';'S1V1(Hz)'; '40hzStim (Hz)'};
        xlabel('time after stimulus onset (s)');
        mediansig=median(stas.(char(names(i))));
        if plotter == 1
            [minfreq,maxfreq] = cwtfreqbounds(length(mediansig),fs); %determine the bandpass bounds for the signal
            cwt(mediansig,[],fs);
            ylabel('Frequency (Hz)')
            colormap(jet)
            title(names(i))
            ylim([.001, .1])
            yticks(ticks)
            
            yticklabels({  0    5.0000   10.0000   15.0000   20.0000   25.0000   30.0000   35.0000   40.0000   45.0000   50.0000  55.0000  60})
    
%             yticklabels({  0 10.0000 20.0000 30.0000 40.0000 50.0000 60})
            
            set(gca,'FontSize',15)
    %         caxis([.00008, .0002]);
            caxis([.00008, .001]);
%             caxis([.00008, .003]);
        else
            ts = timeseries(mediansig,x2) ;
            plot(ts)
%             n = 7500; 
            n = 5000; 
            ts2 = movmean(mediansig, n) ; 
           % median sta + moving average 
%             plot(x2, mediansig, 'y', x2, ts2, 'k')
%             legend(" " + n + " Pt. Average", "median STA");
           % just moving average 
%             plot(x2, ts2, 'DisplayName', title(names(i) + " mV & Moving Average"))
            plot(x2, ts2) 
            title(str(i))
            ylabel('mV') ;
            xlabel("Time (s)") ;
        end
        
         
%     pngFileName = sprintf('plot_%d.fig', i);
	%fullFileName = fullfile(folder, pngFileName);
		
	% Then save it
	%export_fig(fullFileName);
% 	   saveas(gcf, pngFileName)
    end
    
    %% movavg STA timeseries of IL and BLA 
% figure
% caxis_track=[];
% %ylabels={'V1bipolar (Hz)';'S1A (Hz)';'S1V1(Hz)'; '40hzStim (Hz)'};
% xlabel('time after stimulus onset (s)');
% mediansig1=median(stas.RILdata);
% mediansig2=median(stas.RBLAdata); 
% RILx=1:length(stas.RILdata);
% ts = timeseries(mediansig1,RILx) ;
% % moving average of 5000 points
% n = 5000; 
% ts1 = movmean(mediansig1, n) ; 
% ts2 = movmean(mediansig2, n) ; 
% % median sta + moving average 
% % plot(x2, mediansig, 'y', x2, ts2, 'k')
% % legend(" " + n + " Pt. Average", "median STA");
% % just moving average 
% plot(RILx, ts1, 'r')
% hold on 
% plot(RILx, ts2, 'b')
% legend("Right IL", "Right BLA");
% ylabel('mV') ;
% title(names(i) + " mV & Moving Average");
% xlabel("Time (s)") ;
