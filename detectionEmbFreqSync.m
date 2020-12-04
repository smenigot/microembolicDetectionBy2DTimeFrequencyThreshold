% function [Embole, timeRCfinal, spect_wind1_moy] = detectionEmbFreqSync(X_complex_full, Fs, pEmbole, pArtefact)
function [Embole, timeRCfinal] = detectionEmbFreqSync(X_complex_full, Fs)

%% Settings
fcArt = 150;% [Hz] Low-Pass filter to remove artefacts

OverLap     = 80/100;
Lwind       = 64;
Lf          = 512;

pArtefact = erf(2/sqrt(2))*100;
pourcentile = erf(5/sqrt(2))*100;

windAnal = round(5*60/(1/Fs));
nbWindAnal = round(length(X_complex_full)/windAnal);
if nbWindAnal == 0
    nbWindAnal = 1;
end

for k=1:nbWindAnal
    k
    clear spect_wind1_sync spect_wind2_sync Fs_sync fmax_filt_sync timeRC
    fin = k*windAnal;
    if fin > length(X_complex_full)
        fin = length(X_complex_full);
    end
    X_complex = X_complex_full((k-1)*windAnal+1:fin);
    %% Area with usable signal
    [E,tE] = calculEnergie(X_complex, Fs, OverLap, Lwind, Lf, fcArt);
    iE = find(E>1e-5);
    X_complex = X_complex(iE);
    
    %% Heart rhythm
    [spect_wind1, spect_wind2, FsSpect] = spectrogramme2(X_complex, Fs, OverLap, Lwind, Lf, fcArt);
    [timeRC, fmax_filt]                 = detectionRC(spect_wind1, Fs, fcArt, FsSpect, spect_wind2);

    %% Embolus Detection
    [spect_wind1_sync, spect_wind2_sync, Fs_sync, fmax_filt_sync] = ...
        syncSpectrogram(spect_wind1, spect_wind2, FsSpect, timeRC, fmax_filt);
    
    %% Mean Spectrum
    if k==1
        spect_wind1_moy = squeeze(prctile(spect_wind1_sync,pourcentile,1));
    else
        spect_wind1_moy = spect_wind1_moy + squeeze(prctile(spect_wind1_sync,pourcentile,1));
    end
    
    save(['spect_wind' num2str(k)], 'spect_wind1_sync', 'spect_wind2_sync', 'Fs_sync', 'fmax_filt_sync', 'timeRC');
end   
clear spect_wind1_sync spect_wind2_sync Fs_sync fmax_filt_sync timeRC

spect_wind1_moy = spect_wind1_moy/nbWindAnal;% frequency threshold

seuil_artefact  = prctile(reshape(spect_wind2(1:128,:),1,numel(spect_wind2(1:128,:))),...
    pArtefact); % threshold for artefacts

%% detection
N=0;
timeRCfinal=[];
for kP=1:nbWindAnal
    load(['spect_wind' num2str(kP)], 'spect_wind1_sync', 'spect_wind2_sync', 'Fs_sync', 'fmax_filt_sync', 'timeRC');
    for k=19:21%:length(Fs_sync)
        EmboleCycle(k+N) = ...
            detectEmbole(squeeze(spect_wind1_sync(k,:,:)), squeeze(spect_wind2_sync(k,:,:)),...
            spect_wind1_moy, seuil_artefact, fmax_filt_sync(:,k), Fs, Fs_sync(k), k, 1);
        EmboleCycle(k+N).pos = EmboleCycle(k+N).pos+timeRC(k,1)+(kP-1)*5*60;
    end
    N=N+length(Fs_sync);
    
    timeRCfinal=[timeRCfinal;timeRC+(kP-1)*5*60];
    clear spect_wind1_sync spect_wind2_sync Fs_sync fmax_filt_sync timeRC
end
%%
kEmb=1;
Embole.RC     = [];
Embole.pos    = [];
Embole.length = [];
Embole.freq   = [];
Embole.bw     = [];
Embole.Amp    = [];
Embole.AmpMax = [];
for k=1:length(EmboleCycle)
    for k2=1:length(EmboleCycle(k).RC)
        Embole.RC(kEmb)     = EmboleCycle(k).RC(k2);
        Embole.pos(kEmb)    = EmboleCycle(k).pos(k2);
        Embole.length(kEmb) = EmboleCycle(k).length(k2);
        Embole.freq(kEmb)   = EmboleCycle(k).freq(k2);
        Embole.bw(kEmb)     = EmboleCycle(k).bw(k2);
        Embole.Amp(kEmb)    = EmboleCycle(k).Amp(k2);
        Embole.AmpMax(kEmb) = EmboleCycle(k).AmpMax(k2);
        kEmb=kEmb+1;
    end
end
%%
for k=1:nbWindAnal
    delete(['spect_wind' num2str(k) '.mat'])
end

