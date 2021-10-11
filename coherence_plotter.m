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
    file1={'last 5 min baseline\', '5 min stim\', 'after stim 1\','after stim 2\','after stim 3\', 'after stim 4\', 'after stim 5\', 'after stim 6\', 'after stim 7\', 'after stim 8\', 'after stim 9\'};
% elseif whichrat == 2 
%     MainDirectory = 'C:\Users\Administrator\MATLAB\Projects\Pain_study\Data\08_24_2021 rat 1\';
%     file1={'last 5 min baseline\', '5 min stim\', 'after stim 1\','after stim 2\','after stim 3\', 'after stim 4\', 'after stim 5\', 'after stim 6\', 'after stim 7\', 'after stim 8\', 'after stim 9\'};
% elseif whichrat == 3 
%     MainDirectory = 'C:\Users\Administrator\MATLAB\Projects\Pain_study\Data\08_24_2021 rat 2\';
%     file1={'last 5 min baseline\', '5 min stim\', 'after stim 1\','after stim 2\','after stim 3\', 'after stim 4\', 'after stim 5\', 'after stim 6\', 'after stim 7\', 'after stim 8\', 'after stim 9\'};
% end 


%% creating coherence Matrix to store values 

%% Reading experiment dates 


file1= {'8.18.21 5 min\', '8.24.21 r1 5 min\', '8.24.21 r2 5 min\'} ;
str=string(file1);
Directories = {'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\', 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 1\', 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_24_2021 rat 2\'} ; 
dirs = string(Directories) ; 
% MainDirectory = {'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\'} ; 

% for i = 1:length(dirs) 
for i = 1
for f=1:length(str) 
MainDirectory = dirs{i} ;
folder = fullfile(MainDirectory,str{f});

%Change what is in the string depending on which file\files you want to run
file_list = dir([folder 'after stim*.mat']);
baseline = dir([folder 'last 5 min of baseline.mat']); 

set_channels=[1 2 3 4 5]; 
ch_names={'RIL','LBLA','RBLA', 'LIL', 'stimdata'}; 
trial_names={' FIRST LIGHT ONLY' 'LIGHT + US' ' SECOND LIGHT ONLY'};
%% channel configuration

if folder == 'C:\Users\Henry\MATLAB\Projects\Pain_Study\Data\08_18_2021 rat\8.18.21 5 min\' 
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

baseline_rms=[];
 disp(baseline.name);
load([folder baseline.name])
calc_baseline;

% create matrix to hold data for statistical testing
for_stats_new = [];
for_stats_analysis=[];
%insert thresholding here

coh_matrix = zeros(length(dirs),9) ;

for z=1:9
%      if isequal(file_list(z).name,"Trial 2.mat"), continue, end %for 12-23
    disp(z)% display the number that the code is on in the terminal, do not put a ';' after it 
    disp(file_list(z).name);% displays the name of the file in the terminal
    load([folder file_list(z).name]);% bringing the file data into matlab so that the code can run
    calc_coherence;
%     coh_matrix(i, z) = cohmed ; 
end
end 
end 
%% filling matrix 
%     
%     PEN_MATRIX(f, :) = [PENY] ;
%     plot(1:4, PENY, 'o-', 'DisplayName','PEN DATA')
%     title('PEN DATA') 

%% for plotting each experiment data normalized its median of 1st LO 

% for ii = 1:3 
%     concat=['Trial_' num2str(ii)];
%     PEN_MATRIX.(f)(ii, : ) = for_stats_analysis.(concat)(1:550) ; 
% end 
% 
% counter = 0 ; 
% row = (counter*3);
% for ii = 1:3 
%     concat=['Trial_' num2str(ii)];
%     PEN_MATRIX.(ii) = for_stats_analysis.(concat) ; 
% end 
% counter = counter+f ;

% 
% for ii = 1:3 
% concat=[counter 'Trial_' num2str(ii)];
% concat2=['Trial_' num2str(ii)];
% PEN_MATRIX.(concat) = for_stats_analysis.(concat2) ;
% end
% counter = counter + 1 ;
% 
% PEN_MATRIX = for_stats_analysis ;

% dummy = [ for_stats_analysis.Trial_1(1:550);  for_stats_analysis.Trial_2(1:550);  for_stats_analysis.Trial_3(1:550)] ;

% for ii = 1+(3*counter):3*f 
%     concat=['Trial_' num2str(i)];
%     dummy(ii, :) = for_stats_analysis.(concat) ;
% end 
% counter = counter + 1 ;
% 
% if f == 1 
%     PEN_MATRIX1 = [for_stats_analysis.Trial_1(1:550); for_stats_analysis.Trial_2(1:550); for_stats_analysis.Trial_3(1:550)] ;
% elseif f == 2 
%     PEN_MATRIX2.a = [for_stats_analysis.Trial_1(1:329)] ;
%     PEN_MATRIX2.b = [for_stats_analysis.Trial_2(1:550)] ;
%     PEN_MATRIX2.c = [for_stats_analysis.Trial_3(1:550)] ; 
% elseif f == 3 
%     PEN_MATRIX3.a = [for_stats_analysis.Trial_1(1:416) ] ;
%     PEN_MATRIX3.b = [for_stats_analysis.Trial_2(1:550) ];
%     PEN_MATRIX3.c = [for_stats_analysis.Trial_3(1:550)] ; 
% elseif f == 4 
%     PEN_MATRIX4.a = [for_stats_analysis.Trial_1(1:529) ];
%     PEN_MATRIX4.b = [for_stats_analysis.Trial_2(1:529)] ;
%     PEN_MATRIX4.c = [for_stats_analysis.Trial_3(1:550)] ; 
% elseif f == 5 
%     PEN_MATRIX5 = [for_stats_analysis.Trial_1(1:550); for_stats_analysis.Trial_2(1:550); for_stats_analysis.Trial_3(1:550)] ; 
% elseif f == 6
%     PEN_MATRIX6 = [for_stats_analysis.Trial_1(1:550); for_stats_analysis.Trial_2(1:550); for_stats_analysis.Trial_3(1:550)] ; 
% elseif f == 7
%     PEN_MATRIX7 = [for_stats_analysis.Trial_1(1:550); for_stats_analysis.Trial_2(1:550); for_stats_analysis.Trial_3(1:550)] ; 
% end 

% for ii = 1+(3*counter):3*f 
%     for i= 1:3 
%         concat=['Trial_' num2str(i)];
%         PEN_MATRIX.(ii) = for_stats_analysis.(concat) ; 
%     end 
% end
% counter = counter + 1 ;

