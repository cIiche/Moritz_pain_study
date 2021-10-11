%% Authors: Kat Floerchinger, Hannah Mach, Henry Tan
%this code is meant to set the index_stim equal to 60 so that we can
%standardize the length to calculate reliable p values 

%Organize data into structure array
alldata=[]; 

alldata.RILdata=data(datastart(RIL):dataend(RIL)); % Call different fields as StructName.fieldName-> Struct is alldata and field is S1dataR
alldata.LBLAdata=data(datastart(LBLA):dataend(LBLA));
alldata.RBLAdata=data(datastart(RBLA):dataend(RBLA));
alldata.LILdata=data(datastart(LIL):dataend(LIL));
if z == 1 
alldata.stimdata=data(datastart(stim):dataend(stim));
end 

%create names to access fields of 'alldata' for plotting loops
names={'RILdata','LBLAdata','RBLAdata','LILdata', 'stimdata'};
foranalysis={'RILdata','LBLAdata','RBLAdata','LILdata'}; 

%% Filter the raw data with a lowpass butterworth filter and condition removing points > 4std from data 
%Note: it's better to filter the raw data before STAs because this limits
%the end effects that are inevitably created by the filter (you'd get end
%effects at the ends of each STA, vs end effects only at the beginning and
%end of the time series data

% lowEnd = 1; % Hz
% highEnd = 50; % Hz
lowEnd = 5; % Hz
highEnd = 8; % Hz
filterOrder = 3; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
[b, a] = butter(filterOrder, [lowEnd highEnd]/(fs/2)); % Generate filter coefficients

for ii=1:length(names)-1
filteredData.(char(names(ii))) = filtfilt(b, a,alldata.(char(names(ii)))); % Apply filter to data using zero-phase filtering
% is this correct here? before STAs? - in vis stim, for stats_analysis is
% filtered like this, after STA creation 

% for some reason, std and mean create NaN. This is commented out for now
deviation=std(filteredData.(char(names(ii))), 'omitnan');
trialmean = mean(filteredData.(char(names(ii))), 'omitnan');
filteredData.(char(names(ii))) = filteredData.(char(names(ii)))(filteredData.(char(names(ii)))<trialmean+4*deviation);
end

% name pages of wavs
% to match sizes 
if length(filteredData.RILdata) < length(filteredData.RBLAdata) 
    b = length(filteredData.RILdata);
else 
    b = length(filteredData.RBLAdata);
end

% is this too big to plot? 
% figure 
% plot(filteredData.RILdata)
% hold on 
% plot(filteredData.RBLAdata)

wav1 = filteredData.RILdata(1:1000);
wav2 = filteredData.RBLAdata(1:1000); 
cxy = mscohere(wav1, wav2) ;

%   wav(:,:,1) = filteredData.RILdata(1:b);
%   wav(:,:,2) = filteredData.RBLAdata(1:b);
%   
% % get the FFT of the waves
% % fftwav.RIL = fft(filteredData.RILdata) ; 
% % fftwav.RBLA = fft(filteredData.RBLAdata) ; 
% fftwav(:,:,1) = fft(filteredData.RILdata(1:6080000));
% fftwav(:,:,2) = fft(filteredData.RBLAdata(1:6080000));
% 
% % calculate the power-spectral densities (psd) and the cross-spectral
% % densities (csd) and sum them over repetitions
% numsmp = length(wav);
% psd = 2.*abs(fftwav).^2./(numsmp.^2);
% csd = 2.*(fftwav(:,:,1).*conj(fftwav(:,:,2)))./(numsmp.^2);
% sumpsd = squeeze(sum(psd,2));
% sumcsd = squeeze(sum(csd,2));
% 
% % calculate coherence
% coh = abs(sumcsd ./ sqrt(sumpsd(:,1) .* sumpsd(:,2)));
% 
% figure;
% plot(squeeze(wav(:,:,1)));
% figure;
% plot(squeeze(wav(:,:,2)));
% figure;
% plot(coh);
% 
% cohmed = median(coh) ;
