function [spect_wind1_sync, spect_wind2_sync, Fs_sync, fmax_filt_sync] = ...
    syncSpectrogram(spect_wind1, spect_wind2, Fs2, timeRC, fmax_filt)

Lf     = size(spect_wind1,1);

%% Resample for each cardiac cycle
timeRCech          = round(timeRC.*Fs2); % Sample of the beginning and the end
i_zero = find(timeRCech == 0);
timeRCech(i_zero) = 1;
nbPointAnalyseFreq = 1024;

spect_wind1_sync = zeros(length(timeRC),Lf,nbPointAnalyseFreq);
spect_wind2_sync = zeros(length(timeRC),Lf,nbPointAnalyseFreq);
Fs_sync          = zeros(length(timeRC),1);
fmax_filt_sync   = zeros(nbPointAnalyseFreq,length(timeRC));
parfor k=1:length(timeRC)
    t1           = linspace(timeRC(k,1),timeRC(k,2),timeRCech(k,2)-timeRCech(k,1)+1);
    t1_sync      = linspace(timeRC(k,1),timeRC(k,2),nbPointAnalyseFreq);
    Fs_sync(k,1) = 1/(t1_sync(2)-t1_sync(1));
    spect_wind1_sync(k,:,:) = interp1(t1,spect_wind1(:,timeRCech(k,1):timeRCech(k,2))',t1_sync,'spline')';
    spect_wind2_sync(k,:,:) = interp1(t1,spect_wind2(:,timeRCech(k,1):timeRCech(k,2))',t1_sync,'spline')';
    fmax_filt_sync(:,k)     = interp1(t1,fmax_filt(timeRCech(k,1):timeRCech(k,2))',t1_sync,'spline')';
end