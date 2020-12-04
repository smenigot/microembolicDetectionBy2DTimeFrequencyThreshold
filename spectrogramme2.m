function [spect_wind1, spect_wind2, Fs2] = spectrogramme2(X_complex, Fs, OverLap, Lwind, Lf, fcArt)

Delay   = round(Lwind-(Lwind*OverLap));
Nwind   = round((length(X_complex)-Lwind)/Delay);
Fs2     = Fs/Delay;

Ncp   = ceil(Fs/2/fcArt);  
Ncn   = ceil(Fs/2/fcArt);

%% Compute energy and spectrogram

x_wind      = zeros(Lwind,Nwind);
spect_wind1 = zeros(Lf,Nwind);
spect_wind2 = zeros(Lf,Nwind);
spect       = zeros(2*Lf,1);
ESub1       = zeros(Nwind,1);
ESub2       = zeros(Nwind,1);
f150        = round(150/Fs*Lf);

parfor k=1:Nwind
    x_wind(:,k)      = X_complex((k-1)*Delay+1:(k-1)*Delay+Lwind,:);
    spect            = fft(x_wind(:,k),2*Lf);
    
%     Subband11=Ncp:length(spect)/2;
%     ESub1(k,1)=sum((abs(spect(Subband11)).*hamming(length(Subband11))).^2)/length(Subband11);
    spect_wind1(:,k) = abs(spect(1:Lf) - flipud(spect(Lf+1:end)));
    
    
%     Subband12=length(spect)/2+1 : length(spect)-Ncn;
%     ESub2(k,1)=sum((abs(spect(Subband12)).*hamming(length(Subband12))).^2)/length(Subband12);
    spect_wind2(:,k) = abs(spect(Lf+1:end));
    
end

spect_wind1(1:f150,:) = 0;
spect_wind2(end-f150:end,:) = 0;


% ESub = [ESub1 ESub2];
% 
% %% substraction of the Energy of negative frequency
% ESubNorm(:,1)           = ESub(:,1)-ESub(:,2);
% ESubNorm(ESubNorm<0,1)  = 0;
% ESubNorm(:,2)           = ESub(:,2);
% 
% vSmooth = 25;
% ESubNormF(:,1)=smooth(ESubNorm(:,1)',vSmooth,'mean');
% ESubNormF(:,2)=smooth(ESubNorm(:,2)',vSmooth,'mean');
% 
% [n2,bin2] = hist(ESubNormF(:,2),length(ESubNorm));
% [~,b2]    = max(n2(10:end));
% Mpv       = bin2(b2+10);
% 
% %%
% Sat(1)  = mean(ESubNormF(:,1))+5*std(ESubNormF(:,1));
% Sat(2)  = 10*Mpv;
% 
% ESubNormLim(:,1)                    = ESubNorm(:,1);
% ESubNormLim(ESubNorm(:,1)>Sat(1),1) = Sat(1);
% ESubNormF(:,1) = smooth(ESubNormLim(:,1)',vSmooth,'mean');
% 
% ESubNormLim(:,2)                    = ESubNorm(:,2);
% ESubNormLim(ESubNorm(:,2)>Sat(2),2) = Sat(2);
% ESubNormF(:,2) = smooth(ESubNormLim(:,2)',vSmooth,'mean');
% 
% %%
% Fc = round((32/length(ESubNormF)*Fs)*4096/Fs);
% Filtre        = [ones(1,Fc) zeros(1,4096-2*Fc) ones(1,Fc)];
% 
% Spectre1(1,:)  = fft(ESubNormF(:,1),4096);
% temp           = real(ifft(Spectre1.*Filtre));
% ESubNormF(:,1) = temp(1:length(ESubNormF));
% 
% Spectre2(1,:)  = fft(ESubNormF(:,2),4096);
% temp           = real(ifft(Spectre2.*Filtre));
% ESubNormF(:,2) = temp(1:length(ESubNormF));
% 
% STD(1)            = std(ESubNormF(:,1));
% StdBruit(1)       = std(ESubNormLim(:,1)-ESubNormF(:,1));
% VarThreshold(:,1) = 2*(2*STD(1)+StdBruit(1));
% 
% STD(2)            = std(ESubNormF(:,2));
% StdBruit(2)       = std(ESubNormLim(:,2)-ESubNormF(:,2));
% VarThreshold(:,2) = 2*(2*STD(2)+StdBruit(2));
% %% nouveau
% D=ESubNorm-ESubNormF;
% 
% [a11,b11] = hist(D(D(:,1)<0,1),100);
% c11       = find( cumsum(a11/sum(a11))/2 < 0.1/100);
% if (isempty(c11) == 1)
%     std_flux(1) = abs(b11(1));
% else
%     std_flux(1) = abs(b11(c11(end)));
% end
% 
% [a22,b22] = hist(D(D(:,2)<0,2),100);
% c22       = find(cumsum(a22/sum(a22))/2<0.1/100);
% if (isempty(c22) == 1) 
%     std_flux(2) = abs(b22(1));
% else
%     std_flux(2) = abs(b22(c22(end))); 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     
%     