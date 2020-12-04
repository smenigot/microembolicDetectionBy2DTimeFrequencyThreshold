function [timeRC, fmax_smooth] = detectionRC(spect_wind1, Fs, fcArt, FsSpect, spect_wind2)

Lf     = size(spect_wind1,1);

%%
probFreq = prctile(spect_wind1',erf(3/sqrt(2))*100);
probFreq = probFreq-mean(probFreq(490-64:490));
probFreq(probFreq<0) = 0;
probFreq = probFreq/sum(probFreq);

iprobFreq = find(cumsum(probFreq) >=erf(3/sqrt(2)));

%% Max Frequency
[fmax, spect_windsat]   = FreqMax(spect_wind1, iprobFreq(1), Fs, fcArt,2);

% Smooth
fmax_smooth = sgolayfilt(fmax,1,21); 

%% Max detection
fmax_med  = sgolayfilt(smooth(fmax,51,'median'),1,21);
u_prctile = round(prctile(fmax_med,[1-erf(1/sqrt(2)) erf(1/sqrt(2))]*100));
maxtab    = peakdet(fmax_med, diff(u_prctile)/2);
indmax    = maxtab(:,1);

%%
Ncp    = round(Fs/2/fcArt);

Lf     = size(spect_wind1,1);

spect_wind2_smooth = sgolayfilt(spect_wind2,1,21); 
Nbande = 4;
fmax_bande = median(spect_wind2_smooth(512-(Ncp+Nbande*round((Lf-Ncp)/Nbande)-1):Lf-(Ncp+(Nbande-1)*round((Lf-Ncp)/Nbande)),:));
Fc = round(1/(64/FsSpect)*length(spect_wind2_smooth)/2/FsSpect);
Filtre          = [ones(1,Fc) zeros(1,length(fmax_bande)-2*Fc) ones(1,Fc)];
fmax_bande2 = real(ifft(fft(fmax_bande).*Filtre));

%%   
fmax_bande=fmax_bande2;
Fc = round(length(fmax_bande)/1024);
Filtre          = [ones(1,Fc) zeros(1,length(fmax_bande)-2*Fc) ones(1,Fc)];
fmax_bande_filt = real(ifft(fft(fmax_bande).*Filtre));
seuil_artefact  = fmax_bande_filt+diff(prctile(fmax_bande2,[50 erf(2/sqrt(2))*100]));

posArtefact = find(fmax_bande>=seuil_artefact);

indmaxOutside=[];
for k=1:length(indmax(:,1))
    ic  = find(indmax(k,1)>=posArtefact-20 & indmax(k,1)<=posArtefact+20);
    if (~isempty(ic) || fmax_med(indmax(k,1))<median(fmax_med))
        indmaxOutside=[indmaxOutside; indmax(k,1)];
    end
end


%% Time of cadiac cycle
timeE = (0:size(spect_windsat,2)-1)'./FsSpect;
c=1;

for k=1:length(indmax)-1
    if isempty(find(indmax(k)==indmaxOutside)) && isempty(find(indmax(k+1)==indmaxOutside)) % exclusion de cycle cardiaque
        timeRC(c,1)    = timeE(indmax(k));
        timeRC(c,2)    = timeE(indmax(k+1));
        c=c+1;
    end
end
