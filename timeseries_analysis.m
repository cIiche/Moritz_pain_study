%% Authors: Hannah Mach, Kat Floerchinger, Henry Tan

% This script plots CWTs (IN PROGRESS) 
% from pain study rat data 

clear all
close all
clc
% % %% load data
% 
% filePath = 'C:\Users\Henry Tan\MATLAB\Projects\Moritz Pain Study\Data\08_11_21_Rat\';
% fileName ='Trial 1' ;

filePath = 'C:\Users\Administrator\MATLAB\Projects\Pain_study\Data\08_16_2021 rat\';
fileName ='Trial 1' ;

% filePath = 'C:\Users\Administrator\MATLAB\Projects\Pain_study\Data\08_18_2021 rat\';
% fileName ='Trial 1' ;


load([filePath,fileName]);

%% set parameters

set_channels=[1 2 3 4 6] ;

% set channel identities, stim should always be on index 5, no matter what channel 
decision = input('Do you want to manually input channels?(1 = yes, 0 = no):') ;
if decision == 1
    RIL = input('What channel is right IL?:') ;
    LBLA = input('What channel is left BLA?:') ;
    RBLA = input('What channel is right BLA?:') ;
    LIL = input('What channel is left IL?:'); 
    RIL=set_channels(RIL);LBLA=set_channels(LBLA);RBLA=set_channels(RBLA);LIL=set_channels(LIL);stim=set_channels(5) ; 

end 
if decision == 0  
    % 8/11/21 Rat 
    RIL=set_channels(2);LBLA=set_channels(3);RBLA=set_channels(1);LIL=set_channels(4);stim=set_channels(5) ;
end 

%% Set sampling rate
% fs = input('What is the tickrate/sampling rate?:') ;
%Bobola Protocol sampling rate = 10k
% fs=10000 ;
%Eguchi Protocal sampling rate = 20k
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
    % only for 5/18 eguchi data lightstim data is too long
%     if i == 5
%         alldata.stimdata(3280001:6560000) = [] ;
%     end 
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
%% Filter the raw data with a lowpass butterworth filter 
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

%% arrange data for statistical analysis
names = {'RILdata', 'LBLAdata', 'RBLAdata', 'LILdata'} ;
for i=1:4 
    for_stats_analysis.(i) = alldata.(char(names(i))) ;
end 
% % after dividing data into normal distribution: 
% % >mean+4*stddev)(dataset) = [] 
% v1=std(for_stats_new.(conc));
% t1m = mean(for_stats_new.(conc));
% for_stats_analysis.(conc) = for_stats_new.(conc)(for_stats_new.(conc)<t1m+4*v1);

%% for kruskal-wallis 4 pairings 

first_second_third_forth_vector=[for_stats_analysis.RILdata for_stats_analysis.LBLAdata for_stats_analysis.RBLAdata for_stats_analysis.LILdata ];
str1=strings(1,length(for_stats_analysis.Trial_1));

% creating string of names for KW 

str1=strings(1,length(for_stats_analysis.RILdata));
for ii=1:length(for_stats_analysis.1)
    str1(ii)='Right IL';
end

str2=strings(1,length(for_stats_analysis.LBLAdata));
for ii=1:length(for_stats_analysis.2)
    str2(ii)='Left BLA';
end

str3=strings(1,length(for_stats_analysis.RBLAdata));
for ii=1:length(for_stats_analysis.3)
    str3(ii)='Right BLA';
end

str4=strings(1,length(for_stats_analysis.LILdata));
for ii=1:length(for_stats_analysis.4)
    str3(ii)='Left IL';
end

% creating a list of group labels corresponding to the data
first_vs_second_vs_third_vs_forth=[str1 str2 str3 str4];

%Kruskal-wallis tests between trials 1&2, 1&3, 2&3
run_stats_tests(first_second_third_forth_vector, first_vs_second_vs_third_vs_forth); 
% run_stats_tests(first_third_vector, first_vs_third);
% run_stats_tests(second_third_vector, second_vs_third);

% Mann-Whitney U test / Wilcoxon rank sum test significant if Kruskal-Wallis p < 0.05 
MWp1 = ranksum(for_stats_analysis.RILdata,for_stats_analysis.RBLAdata); % pairing 1 R IL to R BLA 
% MWp2 = ranksum(for_stats_analysis.Trial_1,for_stats_analysis.Trial_3); % pairing 2 1st. LO vs. 2nd LO 
% MWp3 = ranksum(for_stats_analysis.Trial_2,for_stats_analysis.Trial_3); % pairing 3 L+US vs. 2nd LO 


%% plot filtered CWTs of STAs
% ticks=[0:.005:.1];
ticks=[0:0.01:.1];
clear yticks
clear yticklabels


% filterbank initialization
% fb = cwtfilterbank('WaveletParameters',[2,5]);

    for i=1:length(names)    
%     for i=1:length(names)  
        figure
        caxis_track=[];
        %ylabels={'V1bipolar (Hz)';'S1A (Hz)';'S1V1(Hz)'; '40hzStim (Hz)'};
        xlabel('time after stimulus onset (s)');
        mediansig=median(stas.(char(names(i))));
        
        
        
        
        [minfreq,maxfreq] = cwtfreqbounds(length(mediansig),fs); %determine the bandpass bounds for the signal
        cwt(mediansig,[],fs);
%         cwt(mediansig,fb,fs)
        ylabel('Frequency (Hz)')
        colormap(jet)
        title(names(i))
        ylim([.001, .1])
        yticks(ticks)
%         yticklabels({  0    5.0000   10.0000   15.0000   20.0000   25.0000   30.0000   35.0000   40.0000   45.0000   50.0000  55.0000  60})
        yticklabels({  0 10.0000 20.0000 30.0000 40.0000 50.0000 60})
        set(gca,'FontSize',15)
%         caxis([.00008, .0002]);
        caxis([.00008, .001]);
        
        
%         pngFileName = sprintf('plot_%d.fig', i);
	%fullFileName = fullfile(folder, pngFileName);
		
	% Then save it
	%export_fig(fullFileName);
% 	   saveas(gcf, pngFileName)
	

    end
