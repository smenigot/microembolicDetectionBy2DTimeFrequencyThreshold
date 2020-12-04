dir = './dataExample/';
filename = ['DopplerSignal'];
[y,Fs] = audioread([dir filename '.wav']);

% analytical signal
fin=length(y);
Xreal      =  y(1:fin,1);
Ximag      = -y(1:fin,2);
X_complex  = Xreal + 1i*Ximag;

% detection
[EmboleFreq, timeRC] = detectionEmbFreqSync(X_complex, Fs);


%% Figure
Lt      = 64;
Lf      = 256;
Yf1     = spectrogram(X_complex,Lt,round(Lt*0.8),Lf,Fs);
timeYf1 = linspace(0,length(X_complex)/Fs,length(X_complex)/round(Lt*0.8)-1);
freq    = (0:Lf-1)/Lf*Fs;

fig = figure(1);

ax(1) = subplot(5,3,1:2);
imagesc(timeYf1, freq, abs(Yf1(1:Lf/2,:))), axis xy
caxis([0 3.5545])
%
hold on

% Suppression des zones de non-detection
c=1;
for k=1:length(timeRC)-1
    if (timeRC(k,2) ~= timeRC(k+1,1))
        timeNonDetect(c,1) = timeRC(k,2);
        timeNonDetect(c,2) = timeRC(k+1,1);
        
        harea = area([timeNonDetect(c,1) timeNonDetect(c,2)], [Fs Fs]);
        set(harea,'FaceColor',[0 0 0],'linewidth',1)
        set(gca,'Layer','top')
        
        c=c+1;
    end
end

% affichage des detections
if ~isempty(EmboleFreq.pos);  plot(EmboleFreq.pos, EmboleFreq.freq,'rx');end
hold off
xlabel('Temps (s)'), ylabel('Frequence (Hz)')
xlim([0 timeYf1(end)])
ylim([0 Fs])
%
if ~isempty(EmboleFreq.pos)
    subplot(5,3,3);
    [nelements5,centers5] = hist(abs(EmboleFreq.freq),10);
    barh(centers5, nelements5);
    ylim([0 Fs])
    xlabel('Frequence (Hz)')
    
    % autres carateristiques
    ax(2) = subplot(5,3,4:5);
    hold on
    plot(EmboleFreq.pos, abs(EmboleFreq.length),'.');
    plot(repmat(EmboleFreq.pos(1):EmboleFreq.pos(end),11,1)',...
        repmat(prctile(EmboleFreq.length,[0:10:100]),length(EmboleFreq.pos(1):EmboleFreq.pos(end)),1));
    hold off
    xlabel('Temps (s)'), ylabel('Durée (s)')
    xlim([0 timeYf1(end)])
    ylim([0 max(abs(EmboleFreq.length))])
    hold on
    for k=1:length(timeNonDetect)
        harea = area([timeNonDetect(k,1) timeNonDetect(k,2)], [max(abs(EmboleFreq.length)) max(abs(EmboleFreq.length))]);
        set(harea,'FaceColor',[0 0 0],'linewidth',1)
        set(gca,'Layer','top')
    end
    hold off
    
    subplot(5,3,6);
    [nelements1,centers1] = hist(abs(EmboleFreq.length),10);
    barh(centers1, nelements1);
    ylim([0 max(abs(EmboleFreq.length))])
    xlabel('Durée (s)')
    
    %
    ax(3) = subplot(5,3,7:8);
    hold on
    plot(EmboleFreq.pos, abs(EmboleFreq.bw),'.');
    plot(repmat(EmboleFreq.pos(1):EmboleFreq.pos(end),11,1)',...
        repmat(prctile(EmboleFreq.bw,[0:10:100]),length(EmboleFreq.pos(1):EmboleFreq.pos(end)),1));
    hold off
    xlabel('Temps (s)'), ylabel('Bande passante (Hz)')
    xlim([0 timeYf1(end)])
    ylim([0 max(abs(EmboleFreq.bw))])
    hold on
    for k=1:length(timeNonDetect)
        harea = area([timeNonDetect(k,1) timeNonDetect(k,2)], [max(abs(EmboleFreq.bw)) max(abs(EmboleFreq.bw))]);
        set(harea,'FaceColor',[0 0 0],'linewidth',1)
        set(gca,'Layer','top')
    end
    hold off
    
    subplot(5,3,9);
    [nelements2,centers2] = hist(abs(EmboleFreq.bw),10);
    barh(centers2, nelements2);
    ylim([0 max(abs(EmboleFreq.bw))])
    xlabel('Bande passante (Hz)')
    
    %
    ax(4) = subplot(5,3,10:11);
    hold on
    plot(EmboleFreq.pos, EmboleFreq.Amp,'.');
    plot(repmat(EmboleFreq.pos(1):EmboleFreq.pos(end),11,1)',...
        repmat(prctile(EmboleFreq.Amp,[0:10:100]),length(EmboleFreq.pos(1):EmboleFreq.pos(end)),1));
    hold off
    xlabel('Temps (s)'), ylabel('Amplitude moyenne (u.a.')
    xlim([0 timeYf1(end)])
    ylim([0 max(abs(EmboleFreq.Amp))])
    hold on
    for k=1:length(timeNonDetect)
        harea = area([timeNonDetect(k,1) timeNonDetect(k,2)], [max(abs(EmboleFreq.Amp)) max(abs(EmboleFreq.Amp))]);
        set(harea,'FaceColor',[0 0 0],'linewidth',1)
        set(gca,'Layer','top')
    end
    hold off
    
    subplot(5,3,12);
    [nelements3,centers3] = hist(abs(EmboleFreq.Amp),10);
    barh(centers3, nelements3);
    ylim([0 max(abs(EmboleFreq.Amp))])
    xlabel('Amplitude moyenne (u.a.)')
    
    %
    ax(5) = subplot(5,3,13:14);
    hold on
    plot(EmboleFreq.pos, EmboleFreq.AmpMax,'.');
    plot(repmat(EmboleFreq.pos(1):EmboleFreq.pos(end),11,1)',...
        repmat(prctile(EmboleFreq.AmpMax,[0:10:100]),length(EmboleFreq.pos(1):EmboleFreq.pos(end)),1));
    hold off
    xlabel('Temps (s)'), ylabel('Amplitude maximale (u.a.')
    xlim([0 timeYf1(end)])
    ylim([0 max(abs(EmboleFreq.AmpMax))])
    hold on
    for k=1:length(timeNonDetect)
        harea = area([timeNonDetect(k,1) timeNonDetect(k,2)], [max(abs(EmboleFreq.AmpMax)) max(abs(EmboleFreq.AmpMax))]);
        set(harea,'FaceColor',[0 0 0],'linewidth',1)
        set(gca,'Layer','top')
    end
    hold off
    
    subplot(5,3,15);
    [nelements4,centers4] = hist(abs(EmboleFreq.AmpMax),10);
    barh(centers4, nelements4);
    ylim([0 max(abs(EmboleFreq.AmpMax))])
    xlabel('Amplitude maximale (u.a.')
    
    linkaxes(ax,'x');
end
