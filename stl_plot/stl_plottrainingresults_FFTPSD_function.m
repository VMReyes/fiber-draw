function [FY, FY_filt, FYpredict, FYpredict_filt, frq_discretes,...
          PY, PY_filt, PYpredict, PYpredict_filt, WW ] = ... 
          stl_plottrainingresults_FFTPSD_function(Y, Y_filt, Ypredict, Ypredict_filt, ...
                nWhichBatch, nWhichSub, FIGNUMBEROFFSET, subplot_mnpIN, strTitlePart )
bUseSubplot = 0;
subplot_mnp = [];

if(nargin < 7)
    FIGNUMBEROFFSET = 0;
end

if(nargin < 8)
    bUseSubplot = 0;
else
   bUseSubplot  = 1; 
   subplot_mnp = subplot_mnpIN;
end
if(isempty(subplot_mnp))
    bUseSubplot = 0;
end

if(nargin < 9)
    strTitlePart = 'Training';
end

%Edits for using real sampling rate
dTsampled = .5;

%% FFT and PSD
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

    %take the FFT
    FY = fft(Y);
    FY_s = fftshift(FY);
    FY_filt = fft(Y_filt);
    %take the FFT
    FYpredict = fft(Ypredict);
    FYpredict_s = fftshift(FYpredict);
    FYpredict_filt = fft(Ypredict_filt);
    %the radial frequencies
    frq_discretes = 2/length(Y).*([0:(length(Y)-1)]-length(Y)/2);
    
    %
    %find the periodgram
    [PY,WW] = pwelch(Y);
    PY_filt = pwelch(Y_filt);
    %take the FFT
    PYpredict = pwelch(Ypredict);
    PYpredict_filt = pwelch(Ypredict_filt);
    %the radial frequencies
    %frq_discretes = 2*pi/length(Y).*([0:(length(Y)-1)]-length(Y)/2);
    
    %
    % Plot resutls of FFT
    strlistlegend = [];
    figure(3000+FIGNUMBEROFFSET)
    if(bUseSubplot)
        subplot(subplot_mnp(1), subplot_mnp(2), subplot_mnp(3))
    end
    plot(frq_discretes,log10(abs(FY_s).^2),'k');            strlistlegend{length(strlistlegend)+1} = 'Data';
    hold on
  %  stem(frq_discretes,log10(abs(FY_filt_s).^2),'g');       strlistlegend{length(strlistlegend)+1} = 'FY_filt_s';
    plot(frq_discretes,log10(abs(FYpredict_s).^2),'r');     strlistlegend{length(strlistlegend)+1} = 'Prediction';
  %  stem(frq_discretes,log10(abs(FYpredict_filt_s).^2),'b');        strlistlegend{length(strlistlegend)+1} = 'FYpredict_filt_s';
    hold off
    legend(strlistlegend)
    ylabel('Log Magnitude')
    xlabel('Frequency')
  	%axis([0 1 -.1 max(log10(abs(FY_s).^2))])
  	axis([0 1 -.1  max([   max(log10(abs(FY_s).^2))   max(log10(abs(FYpredict_s).^2))   ])   ])
    % title('Power Spectrum:  FFT^2  vs Freq (-pi to pi)')
    title('Power Spectrum');
   
    %
    % Plot resutls of Welch
    strlistlegend = [];
    figure(3001+FIGNUMBEROFFSET)
    if(bUseSubplot)
        subplot(subplot_mnp(1), subplot_mnp(2), subplot_mnp(3))
    end
    plot(WW/pi,fmag1(PY),'k');            strlistlegend{length(strlistlegend)+1} = 'Data';
    hold on
  %  plot(WW,fmag1(PY_filt),'g');       strlistlegend{length(strlistlegend)+1} = 'PY_filt';
    plot(WW/pi,fmag1(PYpredict),'r');     strlistlegend{length(strlistlegend)+1} = 'Prediction';
  %  plot(WW,fmag1(PYpredict_filt),'b');        strlistlegend{length(strlistlegend)+1} = 'PYpredict_filt';
    hold off
    legend(strlistlegend)
    ylabel('Normalized to Max Magnitude')
    xlabel('Frequency')
  	axis([0 1 0 max(fmag1(PY))])
  	%axis([-.25 .25 -.1 max(log10(abs(FY_s).^2))])
    % title('Power Spectrum:  FFT^2  vs Freq (-pi to pi)')
    title('Welch Power Spectrum');