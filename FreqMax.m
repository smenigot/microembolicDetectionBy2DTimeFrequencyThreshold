function [fmax, spect_windsat] = FreqMax(spect_wind1, iFMax, Fs, fcArt,bande)


%% Max frequency
Ncp    = round(Fs/2/fcArt);

Nbande = 3;
Lf     = size(spect_wind1,1);

spect_wind1_smooth = sgolayfilt(spect_wind1,1,21); 
amp_moy = zeros(1,Nbande);
parfor k=1:Nbande
    amp_moy(:,k) = mean(mean(spect_wind1_smooth(Ncp+(k-1)*round((iFMax-Ncp)/Nbande):Ncp+k*round((iFMax-Ncp)/Nbande),:)));
end


% Thresholding
spect_windsat = zeros(size(spect_wind1_smooth));
spect_windsat(spect_wind1_smooth>=amp_moy(bande)) = 1;

% Detection of maximal frequency
freq = (0:Lf-1)/Lf*Fs;
fmax = zeros(size(spect_windsat,2),1);
for k=1:size(spect_windsat,2)
    iF = find(spect_windsat(16:iFMax-8,k)==1);
    if ~isempty(iF)
        fmax(k,1) = freq(iF(end)+16);
    else
        if k>1
            fmax(k,1) = fmax(k-1,1);
        else
            fmax(k,1) = 0;
        end
    end
end
