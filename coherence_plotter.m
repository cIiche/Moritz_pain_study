%%  Authors: Kat Floerchinger, Hannah Mach, Henry Tan
%This code plots dcoherence/dt for the Moritz Rat Pain Study 
clear all
close all
clc
%% reading rat experiment files 

% whichrat = input("Which experiment to run? '1' = 8/18/21 rat, '2' = 8/24/21 rat 1, '3' = 8/24/21 rat 2: "); 
% 
% if whichrat == 1 
%     MainDirectory = 'C:\Users\Administrator\MATLAB\Projects\Pain_study\Data\08_18_2021 rat\';
%     file1={'last 5 min baseline\', '5 min stim\', 'after stim 1\','after stim 2\','after stim 3\', 'after stim 4\', 'after stim 5\', 'after stim 6\', 'after stim 7\', 'after stim 8\', 'after stim 9\', 'after stim 10\','after stim 11\','after stim 12\'};
%     file1={'last 5 min baseline\', '5 min stim\', 'after stim 1\','after stim 2\','after stim 3\', 'after stim 4\', 'after stim 5\', 'after stim 6\', 'after stim 7\', 'after stim 8\', 'after stim 9\'};
% elseif whichrat == 2 
%     MainDirectory = 'C:\Users\Administrator\MATLAB\Projects\Pain_study\Data\08_24_2021 rat 1\';
%     file1={'last 5 min baseline\', '5 min stim\', 'after stim 1\','after stim 2\','after stim 3\', 'after stim 4\', 'after stim 5\', 'after stim 6\', 'after stim 7\', 'after stim 8\', 'after stim 9\'};
% elseif whichrat == 3 
%     MainDirectory = 'C:\Users\Administrator\MATLAB\Projects\Pain_study\Data\08_24_2021 rat 2\';
%     file1={'last 5 min baseline\', '5 min stim\', 'after stim 1\','after stim 2\','after stim 3\', 'after stim 4\', 'after stim 5\', 'after stim 6\', 'after stim 7\', 'after stim 8\', 'after stim 9\'};
% end 


%% creating coherence Matrix to store values 

%% Reading experiment dates 


file1= {'8.18.21 5 min renamed\', '8.24.21 r1 5 min\', '8.24.21 r2 5 min\'} ;
str=string(file1);
Directories = {'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\', 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 1\', 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 2\'} ; 
dirs = string(Directories) ; 
% MainDirectory = {'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\'} ; 

% for i = 1:length(dirs) 
for rat = 1
% for f=1:length(str) 
for f=1
MainDirectory = dirs{rat, rat} ;
folder = fullfile(MainDirectory,str{f});

%Change what is in the string depending on which file\files you want to run
file_list = dir([folder 'Trial*.mat']);
baseline = dir([folder 'Trial 1.mat']); 

set_channels=[1 2 3 4 5]; 
ch_names={'RIL','LBLA','RBLA', 'LIL', 'stimdata'}; 
trial_names={' FIRST LIGHT ONLY' 'LIGHT + US' ' SECOND LIGHT ONLY'};

%% channel configuration

% if folder == 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\8.18.21 5 min\' 
if rat == 1 
% if i == 1 
    RIL=set_channels(2);LBLA=set_channels(4);RBLA=set_channels(1);LIL=set_channels(3);stim=set_channels(5) ;
% end 
% if folder == 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 1\8.24.21 r1 5 min\'  
%     RIL=set_channels(4);LBLA=set_channels(1);RBLA=set_channels(2);LIL=set_channels(3);stim=set_channels(6) ;
% end 
% if folder == 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 2\8.24.21 r2 5 min\'  
%     RIL=set_channels(2);LBLA=set_channels(1);RBLA=set_channels(4);LIL=set_channels(3);stim=set_channels(6) ;
end 

%% calculate pre-stim RMS for normalization 

% baseline_rms=[];
% disp(baseline.name);
% load([folder baseline.name])
% calc_baseline;

% create matrix to hold data for statistical testing
for_stats_new = [];
for_stats_analysis=[];
%insert thresholding here

coh_matrix = zeros(4,11) ;

for z=1:11
%      if isequal(file_list(z).name,"Trial 2.mat"), continue, end %for 12-23
    disp(z)% display the number that the code is on in the terminal, do not put a ';' after it 
    disp(file_list(z).name);% displays the name of the file in the terminal
    load([folder file_list(z).name]);% bringing the file data into matlab so that the code can run
    calc_coherence;
    coh_matrix(rat, z) = cohmed ; 
end
end 
end 

min5 = 1:11 ;
plot(min5, coh_matrix)



