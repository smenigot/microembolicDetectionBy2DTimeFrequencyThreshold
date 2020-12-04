function EmboleCycle = ...
    detectEmbole(spect_wind1_sync, spect_wind2_sync, spect_wind1_moy, seuil_artefact, fmax_filt_sync, Fs, Fs_sync, kRC, display)

Lf                 = size(spect_wind1_sync,1);
nbPointAnalyseFreq = size(spect_wind1_sync,2);

%% init 
EmboleCycle.RC     = [];
EmboleCycle.pos    = [];
EmboleCycle.length = [];
EmboleCycle.freq   = [];
EmboleCycle.bw     = [];
EmboleCycle.Amp    = [];
EmboleCycle.AmpMax = [];

%% Event Detection
% timeOld = 0;

vFreq = (0:Lf-1)*Fs/Lf;
mFreq = repmat(vFreq,nbPointAnalyseFreq,1)';

timeE = (0:nbPointAnalyseFreq-1)'./Fs_sync;

mFmax_filt = repmat(fmax_filt_sync,1,Lf)';
temp0_1    = spect_wind1_sync;
temp0_1(mFreq >= mFmax_filt) = 0;
temp0_2    = spect_wind2_sync;
temp0_2(mFreq <= mFmax_filt) = 0;

spect_detect1 = imclose(temp0_1 >= spect_wind1_moy,strel('disk',64));
spect_detect2 = imclose(temp0_2 >= seuil_artefact,strel('rectangle',[512 4]));

temp_detect                  = spect_detect1 - 2*flipud(spect_detect2);
temp_detect(temp_detect==-2) = -1;
spect_detect_full            = temp_detect;

temp_detect(temp_detect~=0) = 1;
spect_detect_object = imdilate(medfilt2(temp_detect, [5 5]),strel('disk',40));
evenement           = bwconncomp(spect_detect_object, 8);

kEmb=1;
for kEv=1:length(evenement.PixelIdxList)
    grain = zeros(size(spect_detect_object));
    grain(evenement.PixelIdxList{kEv}) = 1;
    if (sum(sum((grain .* spect_detect_full)==-2)) == 0) && (sum(sum((grain .* spect_detect_full)==-1)) == 0)
        EmboleCycle.RC(kEmb) = kRC;
        
        grain2=spect_detect_full.*grain;
        grain2(grain~=0)=1;
        
        [ix,iy]=find(grain2==1);
        EmboleCycle.pos(kEmb)    = timeE(round(mean(iy)));
        if max(iy) >= length(timeE)
            if min(iy)-1 < 1
                EmboleCycle.length(kEmb) = timeE(max(iy))-timeE(min(iy));
            else
                EmboleCycle.length(kEmb) = timeE(max(iy))-timeE(min(iy)-1);
            end
        else
            EmboleCycle.length(kEmb) = timeE(max(iy)+1)-timeE(min(iy));
        end
        EmboleCycle.freq(kEmb) = round(mean(ix))*Fs/Lf;
        if max(iy) >= Lf
            EmboleCycle.bw(kEmb) = (max(ix)-min(ix)+1)*Fs/Lf;
        else
            EmboleCycle.bw(kEmb) = ((max(ix)+1)-min(ix))*Fs/Lf;
        end
        
        temp = grain2 .* spect_wind1_sync;
        EmboleCycle.Amp(kEmb)    = mean(mean(temp(temp~=0)));
        EmboleCycle.AmpMax(kEmb) = max(max(temp(temp~=0)));
        
        kEmb = kEmb + 1;
    end
end

