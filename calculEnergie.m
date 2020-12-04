function [E_interp,t_interp] = calculEnergie(X_complex, Fs, OverLap, Lwind, Lf, fcArt)

OverLap = 0;
Delay   = round(Lwind-(Lwind*OverLap));
Nwind   = round((length(X_complex)-Lwind)/Delay);
Fs2     = Fs/Delay;

Ncp   = round(Fs/2/fcArt); 

%% Compute energy and spectrogram

x_wind = zeros(Lwind,1);
spect = zeros(2*Lf,1);
E     = zeros(Nwind,1);

parfor k=1:Nwind
    x_wind = X_complex((k-1)*Delay+1:(k-1)*Delay+Lwind,:);
    spect  = fft(x_wind,2*Lf);
    Subband11=Ncp:length(spect)/2;
    E(k,1) = sum((abs(spect(Subband11)).*hamming(length(Subband11))).^2)/length(Subband11);
end
t=(0:length(E)-1)*1/Fs2;
t_interp=(0:length(X_complex)-1)*1/Fs;
E_interp=interp1(t,E,t_interp);

