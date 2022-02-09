%% Authors: Hannah Mach, Kat Floerchinger, Henry Tan

% This script plots CWTs (IN PROGRESS) 
% from pain study rat data 

clear all
close all
clc
% % %% load data

% filePath = 'C:\Users\Henry Tan\MATLAB\Projects\Pain_study\Data\08_11_21_Rat\';
% fileName ='Trial 1' ;
% 
% filePath = 'C:\Users\Henry Tan\MATLAB\Projects\Pain_study\Data\08_18_2021 rat\';
% fileName ='Trial 1' ; 

filePath = 'C:\Users\Henry Tan\MATLAB\Projects\Pain_study\Data\08_24_2021 rat 1\';
% fileName ='Trial 1' ;
fileName ='Refractory part 4' ;

% filePath = 'C:\Users\Henry Tan\MATLAB\Projects\Pain_study\Data\08_24_2021 rat 2\';
% % fileName ='Trial 1' ;
% fileName ='Refractory part 4' ;


load([filePath,fileName]);

%% set parameters
% for 8/11
% set_channels=[1 2 3 4 6] ;
% for 8/16, 8/18 
set_channels=[1 2 3 4 5 6 7] ;

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
    % 8/11/21 
%     RIL=set_channels(1);LBLA=set_channels(4);RBLA=set_channels(3);LIL=set_channels(2);stim=set_channels(5) ;
    % 8/16/21
%     RIL=set_channels(2);LBLA=set_channels(3);RBLA=set_channels(1);LIL=set_channels(4);stim=set_channels(5) 
    % 8/24/21 rat 1 
%     RIL=set_channels(4);LBLA=set_channels(1);RBLA=set_channels(2);LIL=set_channels(3);stim=set_channels(5) ;
    % 8/24/21 rat 2 
    RIL=set_channels(2);LBLA=set_channels(1);RBLA=set_channels(4);LIL=set_channels(3);stim=set_channels(5) ;
end 

%% running stim data or refractory 

runstas = input("Running stim data:'1', or refractory data:'2'?:") ; 

%% plot CWT or timeseries 

plotter = input("Plot CWT = '1', STA timeseries = '2': ") ;

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
if runstas == 1 
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
else
    index_stim=[];%initialize reference array of stimulus onsets for STAs
    index_allstim=[];%secondary array of stimuli for identifying first pulse of a train. 0
% implement artifical index_trains to plot avg of every second for 10 min refractories   
% 0s and then 1s for every 60 units is generally index_train == 1 light every 10 sec. 
% -- using every 100 units here
% train_length = int(length(alldata.(char(names(1))))/(3423.2)); 
%     train_length = length(alldata.(char(names(1))))/(2000) ;
% % 8/24/21 rat 1 lightstim = ~6,000,000 indexes long train becomes 
% % ~ 2000 indexes long. (lightstim data length)/train_length ~= 3423.2 
%     artificial_train=zeros(1,train_length) ; 
%     for k = 1:length(alldata.(char(names(1)))/100) 
%         artificial_train(60*k) = 1; 
%     end 
% %     index_allstim(1)=[];
%     index_stim=index_allstim(artificial_train);

% 8/24/21 rat 1 index_stim data 
    if filePath == 'C:\Users\Henry Tan\MATLAB\Projects\Pain_study\Data\08_24_2021 rat 1\';
    index_stim = [327677,527677,727677,927677,1127677,1327677,1527677,1727677,1927677,2127677,2327677,2527677,2727677,2927677,3127677,3327677,3527677,3727677,3927677,4127677,4327677,4527677,4727677,4927677,5127677,5327677,5527677,5727677,5927677,6127677,6327677,6527677,6727677];
    end
% 8/24/21 rat 1 index_stim data 
    if filePath == 'C:\Users\Henry Tan\MATLAB\Projects\Pain_study\Data\08_24_2021 rat 2\';
    index_stim = [312677,512677,712677,912677,1112677,1312677,1512677,1712677,1912677,2112677,2312677,2512677,2712677,2912677,3112677,3312677,3512677,3712677,3912677,4112677,4312677,4512677,4712677,4912677,5112677,5312677,5512677,5712677,5912677];
    end 
end

    %% Create STAs
    
for i=1:length(names) %initiate data array to hold STAs
    stas.(char(names(i)))=[];
end

tb=1; %time before stim to start STA
ta=9; %time after stim to end STA
% FOR refractory data plotting, this means the CWT output will be average
% over 10 seconds for stim fake stim event 

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
% ticks=[0:.005:.1];
ticks=[0:0.01:.1];
clear yticks
clear yticklabels

str = string(names) ;
    for i=1:length(names)-2    
%     for i=1:length(names)  
        figure
        caxis_track=[];
        %ylabels={'V1bipolar (Hz)';'S1A (Hz)';'S1V1(Hz)'; '40hzStim (Hz)'};
        xlabel('time after stimulus onset (s)');
        if plotter == 1
            if runstas == 1 
                % stimdata STA 
                mediansig=median(stas.(char(names(i))));
                [minfreq,maxfreq] = cwtfreqbounds(length(mediansig),fs); %determine the bandpass bounds for the signal
                cwt(mediansig,[],fs);
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
            else 
                % refractory cwt 
                mediansig=median(stas.(char(names(i))));
                [minfreq,maxfreq] = cwtfreqbounds(length(mediansig),fs); %determine the bandpass bounds for the signal
                cwt(mediansig,[],fs);
                ylabel('Frequency (Hz)')
                colormap(jet)
                title(names(i))
                ylim([.001, .1])
                yticks(ticks)
        %         yticklabels({  0    5.0000   10.0000   15.0000   20.0000   25.0000   30.0000   35.0000   40.0000   45.0000   50.0000  55.0000  60})
                yticklabels({  0 10.0000 20.0000 30.0000 40.0000 50.0000 60})
                set(gca,'FontSize',15)
%                 caxis([.00008, .001]);
                caxis([.00001, .0002]);
            end 
        else
            % Timeseries movavg 
            mediansig=median(stas.(char(names(i))));
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