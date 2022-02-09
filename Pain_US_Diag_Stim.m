%% Authors: Devon Griggs, John Kucewicz, Nels Schimek, Kat Floerchinger, Hannah Mach, Alissa Phutirat, Henry Tan
%this code is meant to set the index_stim equal to 60 so that we can
%standardize the length to calculate reliable p values 
fs=tickrate(1);
time1=1:dataend(1);
time=time1/fs/60;

%% change the bandpass for filtering pls
    % default done 
%     [bb,aa]=butter(3,[3,100]/(fs/2)); %trying to get the us noise out, 3 to 200
% [bb,aa]=butter(2,[3,100]/(fs/2)); %trying to get the us noise out, 3 to 200
%     low gamma 
%    [bb,aa]=butter(2,[30,59]/(fs/2)); 
    % beta
%   [bb,aa]=butter(2,[12,29]/(fs/2)); 
    % alpha 
%     [bb,aa]=butter(2,[8,11]/(fs/2));
    % theta
%   [bb,aa]=butter(2,[4,7]/(fs/2));


bandpasses = [3, 100; 30, 59; 12, 29; 8, 11; 4, 7;] ; 
    [bb,aa]=butter(2,bandpasses(brain_wave,:)/(fs/2));

%Organize data into structure array
alldata=[]; %initialize structure array

alldata.V1Ldata=filtfilt(bb,aa,data(datastart(V1L):dataend(V1L))')';
alldata.S1Ldata=filtfilt(bb,aa,data(datastart(S1L):dataend(S1L))')';
alldata.S1Rdata=filtfilt(bb,aa,data(datastart(S1R):dataend(S1R))')';
alldata.V1Rdata=filtfilt(bb,aa,data(datastart(V1R):dataend(V1R))')';
% alldata.lightstimdata=filtfilt(bb,aa,data(datastart(lightstim):dataend(lightstim))')';
alldata.lightstimdata=data(datastart(lightstim):dataend(lightstim));
%  alldata.lightstimdata=filtfilt(bb,aa,data(datastart(5):dataend(5))')'; % for 5/29 Hypothesis:
% this works if not all channels are imported; channel 7 = channel 5 vs.
% channel 1... 5, 6, 7 => channel 1-4, channel 7 
% what changes in experiments is whether channels 5-6 are included in
% export

%create names to access fields of 'alldata' for plotting loops
names={'V1Ldata','S1Ldata','S1Rdata','V1Rdata','lightstimdata'}; 
foranalysis={'V1Ldata','S1Ldata','S1Rdata','V1Rdata'}; 

%Filter out noise
% alldata.V1Ldata=alldata.V1Ldata(abs(alldata.V1Ldata)<0.02); %hardcoded filtering
%% detect stimuli

X = alldata.lightstimdata;
X=X-min(X);
X=X/max(X);
Y=X>0.5;
Z=diff(Y);
index_allstim=find(Z>0.5);index_allstim=index_allstim+1;

%find first pulse of each train, if stimulation contains trains
index_trains=diff(index_allstim)>2*fs; 
index_allstim(1)=[];
index_stim=index_allstim(index_trains);


%% Create STAs
for i=1:4 %initiate data array to hold STAs
stas.(char(names(i)))=[];
       % title(file_list(z).name,'interpreter','none');
end

if  time_series ==3 
    tb=1; %time before stim to start STA
    ta = 3; % first 3 seconds
else 
    tb=1;
    ta=10; %time after stim to end STA
end 

%% Data conditioning (cont.) 
% prevents errors based on discrepency between V1Ldata and
% lightstimdata length

for i=1:4  
%     %inner_loop_size = 0;
% %     if (fix(length(alldata.V1Ldata)/tickrate/10)) < length(index_stim)
% %         inner_loop_size = (fix(length(alldata.V1Ldata)/tickrate/10)-1);
% %     else
% %         inner_loop_size = length(index_stim)-1;
% %     end  
% %     for j=2:inner_loop_size %(length(index_stim)-1) %cycle through stimuli

%        for j=2:(length(index_stim)-1)
for j=2:(length(index_stim)-2) % to compensate for data chopping so data
% vectors are long enough (supposed to be 60 entries) for
% j=2:(length(index_stim)-3) % for 6/25/20 mouse experiment 1,  6/24/20
% experiment 1
%         for j=2:(length(index_stim)-4)
%      for j=2:(length(index_stim)-5) % 6/24/20 experiment 3 
        stas.(char(names(i)))=[stas.(char(names(i))); alldata.(char(names(i)))((index_stim(j)-fs*tb):(index_stim(j)+fs*ta))];

    end
end

% index_stim_conditioner(fs, tb, ta, names, index_stim, stas, alldata); 
% stas.(char(names(1))) = sta(1) ;
% stas.(char(names(2))) = sta(2) ;
% stas.(char(names(3))) = sta(3) ;
% stas.(char(names(4))) = sta(4) ;

% for i=1:4
% try 
%     for j =2:(length(index_stim)-1) 
%         stas.(char(names(i)))=[stas.(char(names(i))); alldata.(char(names(i)))((index_stim(j)-fs*tb):(index_stim(j)+fs*ta))];
%     end
% catch 
%     warning('Index exceeds the number of array elements. Trying j=2:(length(index_stim)-2)') 
%     for j =2:(length(index_stim)-2) 
%         stas.(char(names(i)))=[stas.(char(names(i))); alldata.(char(names(i)))((index_stim(j)-fs*tb):(index_stim(j)+fs*ta))];
%     end
% try
%    warning('Index exceeds the number of array elements. Trying j=2:(length(index_stim)-2)') 
%     for j =2:(length(index_stim)-2) 
%         stas.(char(names(i)))=[stas.(char(names(i))); alldata.(char(names(i)))((index_stim(j)-fs*tb):(index_stim(j)+fs*ta))];
%     end 
% catch 
%     for j =2:(length(index_stim)-3) 
%         stas.(char(names(i)))=[stas.(char(names(i))); alldata.(char(names(i)))((index_stim(j)-fs*tb):(index_stim(j)+fs*ta))];
%     end 
% end
% end
% end 
%% plot STAS

responseWindowEnd=0.4;

%figure
x2=1:length(stas.(char(names(1))));
x2=x2/fs-1;%time axis
%ylabels={'S1 (mV)';'A1 (mV)';'V1R (mV)'; 'V1L (mV)';'A1R (Hz)';'S1R (Hz)'};

minvalarray=zeros(length(foranalysis),1);
minidxarray=zeros(length(foranalysis),1);
maxvalarray=zeros(length(foranalysis),1);
maxidxarray=zeros(length(foranalysis),1);
RMSvalbarray=zeros(length(foranalysis),1);
RMSvalaarray=zeros(length(foranalysis),1);

%%
% for i=1:length(foranalysis)
%     figure(i+1)
%     ylabels={'V1L (Hz)';'S1L (Hz)';'S1R (Hz)'; 'V1R (Hz)'};
%     subplot(4,1,i);
%     a=mean(stas.(char(names(i))));
%     a=a-mean(stas.(char(names(i)))(1:fs));
%     a=a/100*1000;

%%
    

%     a=filtfilt(bb,aa,a); % used for the butter filter
%     %[minval,minidx]=min(a(fs*(tb):fs*(tb+responseWindowEnd))); %identifying min max within response window 
%     [maxval,maxidx]=max(a(fs*(tb):fs*(tb+responseWindowEnd)));
%     RMSvalb=rms(a(fs*(tb-0.25):fs*(tb)));%RMS before zero
%     RMSvala=rms(a(fs*(tb):fs*(tb+0.25)));%RMS after zero
%     
%     RMSvalbarray(i)=RMSvalb; %declares RMS arrays for print
%     RMSvalaarray(i)=RMSvala;
%     
%     %minvalarray(i)=minval; %declares minmax arrays for print 
%     %minidxarray(i)=minidx;
%     maxvalarray(i)=maxval;
%     maxidxarray(i)=maxidx; %maxtime array 
%     
% %     Plot STA's. IMPORTANT TO KEEP
%     figure;
%     plot(x2,a,'linewidth',1);hold on 
%    plot(x2(maxidx+fs*(tb)),a(maxidx+fs*(tb)),'o');hold on
%     xlim([-0.5 1.5]);% seconds that will be shown in the plot, stim is on at time 0    
%     ylim([-0.05 0.05]);%mV range on the plot, edit to get the entire signal to show
%    ylabel(ylabels(i));
%     
%     
% %     
% %     
% %     arrays of RMS values, one for each second of a 10 sec segment
% %     may want overlap in timeframes eventually
% %     This is gross and needs to be edited to reduce redundancy
%     all_points(i).RMSvalsb=rms(d(:,fs*(tb-0.25):fs*(tb))');
%     all_points(i).RMSvals_1=rms(d(:,fs*(tb):fs*(tb+1.00))');
%     all_points(i).RMSvals_2=rms(d(:,fs*(tb+1.00):fs*(tb+2.00))');
%     all_points(i).RMSvals_3=rms(d(:,fs*(tb+2.00):fs*(tb+3.00))');
%     all_points(i).RMSvals_4=rms(d(:,fs*(tb+3.00):fs*(tb+4.00))');
%     all_points(i).RMSvals_5=rms(d(:,fs*(tb+4.00):fs*(tb+5.00))');
%     all_points(i).RMSvals_6=rms(d(:,fs*(tb+5.00):fs*(tb+6.00))');
%     all_points(i).RMSvals_7=rms(d(:,fs*(tb+6.00):fs*(tb+7.00))');
%     all_points(i).RMSvals_8=rms(d(:,fs*(tb+7.00):fs*(tb+8.00))');
%     all_points(i).RMSvals_9=rms(d(:,fs*(tb+8.00):fs*(tb+9.00))');
%     all_points(i).RMSvals_10=rms(d(:,fs*(tb+9.00):fs*(ta))'); % I think this is redundant
% %     
% %    
% % 
%     disp('name:')
%     disp(names(i))
%     disp('Individual BEFORE rms values')
%     all_points(i).RMSvalsb'
%     disp('mean')
%     mean(all_points(i).RMSvalsb')
%     disp('stddev')
%     std(all_points(i).RMSvalsb')
%     disp('name:')
%     disp(names(i))
%     disp('Individual AFTER rms values')
%     all_points(i).RMSvals_1'
%     disp('mean')
%     mean(all_points(i).RMSvals_1')
%     disp('stddev')
%     std(all_points(i).RMSvals_1')
% 
% 
%    disp('Paired t-test');
%     [H,P,CI,STATS] = ttest(all_points(i).RMSvalsb',all_points(i).RMSvalsa');
%     disp('Individual Percentage Changes');
%    temp=(all_points(i).RMSvalsa' - all_points(i).RMSvalsb')./all_points(i).RMSvalsb'.*100;
%    disp('mean');
%    mean(temp);
%    disp('stddev');
%    std(temp);
%    disp('Individual Percentage Changes (absolute value)');
%    temp=abs(all_points(i).RMSvalsa' - all_points(i).RMSvalsb')./all_points(i).RMSvalsb'.*100;
%    disp('mean');
%    mean(temp);
%    disp('stddev');
%    std(temp);
%   
% end

%% collect all of the individual points of data for waterfalls

all_points(1).name=names(1);
d=stas.(char(names(1)));


%filter data
% d=filtfilt(bb,aa,d')';

% max_rms=0;

for k=1:10
       concat=['RMSvals_' num2str(k)];
%        all_points(1).(concat)=rms(alldata.(char(names))(:,fs*(tb+k-1):fs*(tb+k))');
       all_points(1).(concat)=rms(d(:,fs*(tb+k-1):fs*(tb+k))');
%        plot(d(:,fs*(tb+k-1):fs*(tb+k)))
%        maxrms=max(all_points(1).(concat));
%        if maxrms > max_rms
%            max_rms=maxrms;
%        end
end 

      
 % disp(max_rms);

% Create matrix of V1L RMS vals for plotting
% matrix=[all_points(1).RMSvals_9; all_points(1).RMSvals_8; all_points(1).RMSvals_7; all_points(1).RMSvals_6; all_points(1).RMSvals_5;
%     all_points(1).RMSvals_4; all_points(1).RMSvals_3; all_points(1).RMSvals_2; all_points(1).RMSvals_1];

% first 3 sec

    matrix=[all_points(1).RMSvals_1; all_points(1).RMSvals_2; all_points(1).RMSvals_3; all_points(1).RMSvals_4; all_points(1).RMSvals_5;
    all_points(1).RMSvals_6; all_points(1).RMSvals_7; all_points(1).RMSvals_8; all_points(1).RMSvals_9; all_points(1).RMSvals_10];


%% normalize the data using the baseline RMS
matrix=matrix/rms_baseline;

%% plotting waterfalls 

figure
%imagesc plot
% subplot(2,3,z);
imagesc(matrix')
ylim=[0 0.3];
colorbar

% naming waterfall plots based on 'z'
names = {'1st Light Only', 'This shouldnt be plotted', 'Light + US', '2nd Light Only'} ;
title(names(z)) % z = 1:4 trials in loopy

% setting waterfall axes 
ylim=[0 0.3];
ylabel('Stimulus event #'); 
ticks = 0:5:60 ; 
yticks(ticks) ; 
xlabel('Time after stimulus (s)') 

% for 3 second analysis 
if time_series == 3
    set(gca,'XTick',[1 2 3 4 5 6] ); %This is going to be the only values affected. 
    set(gca,'XTickLabel',[0.5 1 1.5 2 2.5 3] ); %This is what it's going to appear in those places.
else
    xticks(ticks)
end  
    
colorbar
%to set the magnitude for the color bar, change accordingly?????
caxis([0.00015 0.05])


%% calculate z-scores
[z_scores,mu,sigma]=find_zscores(matrix, baseline_rms);

%% arrange data for statistical analysis

% reshape rms matrix into a vector to be used to statistical testing
s=size(matrix);
S=s(1)*s(2);
for_stats=reshape(matrix,1,S);
conc=['Trial_' num2str(z)];
for_stats_analysis.(conc)=for_stats;

%% additional hardcoded filtering 
% removes outlier data points 4 standard deviations from mean 

deviation=std(for_stats_analysis.(conc));
trialmean = mean(for_stats_analysis.(conc));
for_stats_analysis.(conc) = for_stats_analysis.(conc)(for_stats_analysis.(conc)<trialmean+4 *deviation);

%% plotting the data

% subplot(2,3,z)
% histogram(matrix)
% title('RMS Distribution','interpreter','none');
% xlabel('RMS');
% ylim([0 200]);
% 
% subplot(2,3,z+3)
% histogram(z_scores)
% title('Z score Distribution','interpreter','none');
% xlabel('Z score');
% ylim([0 150]);


% plot raw EEG data
% figure
% plot(alldata.V1Ldata)
% title(trial_names(z),'interpreter','none');
% ylabel('V1L response (mV)');
% xlabel('Time (ms)');

%xlabel('time after stimulus onset (s)');
%prints minmax arrays declared above 
%disp('minvals') 
%minvalarray
%disp('minidx')
%minidxarray
% disp('maxvals')
% maxvalarray
%disp('maxtime');
%(maxidxarray/fs);
% disp('maxtime')
% x2(maxidxarray)' %where we think we need Devon for fixing time values
% a(maxidxarray)' %where we think we need Devon for fixing time values
%disp('mintime')
%x2(minidxarray)'
%disp('RMSval before zero');
%RMSvalbarray;;
%disp('RMSval after zero');
%RMSvalaarray;;

%disp('t-test on the before and after RMS values')
%[H,P,CI]=ttest(RMSvalbarray,RMSvalaarray)

% 
% %% plot CWTs of STAs
% 
% if plot_cwt==1 
% %     figure
% %     caxis_track=[];
% %     ylabels={'S1 (hz)';'A1 (hz)';'V1R (hz)'; 'V1L (hz)'};
% %      xlabel('time after stimulus onset (s)');
%     for i=1:length(foranalysis)
%         figure
%         caxis_track=[];
%         ylabels={'S1L (Hz)';'A1L (Hz)';'V1R (Hz)'; 'V1L (Hz)';'A1R (Hz)';'S1R (Hz)'};
%         xlabel('time after stimulus onset (s)');
%         a=mean(stas.(char(names(i))));
%         cwt(a,[],fs);
%         ylim([0.0005 0.032]);
%         ylabel(ylabels(i));
%         colormap(jet);
%         %caxis_track=[caxis_track;caxis]
% yticks([0.0005,0.001,0.002,0.008,0.032,0.1]);
% yticklabels({0.5,1,2,8,32,100});
% xticks([0,1,2,3,4,5,6,7,8,9]);
%  title(file_list(z).name,'interpreter','none');
% 
%     end
%      
% end