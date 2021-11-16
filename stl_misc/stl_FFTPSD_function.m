function [FY, FY_filt, FYpredict, FYpredict_filt, frq_discretes,...
          PY, PY_filt, PYpredict, PYpredict_filt, WW ] = ... 
          stl_FFTPSD_function(Y, Y_filt, Ypredict, Ypredict_filt)


%% FFT and PSD


    %take the FFT
    FY = fft(Y);
    FY_s = fftshift(FY);
    FY_filt = fft(Y_filt);
    FY_filt_s = fftshift(FY_filt);
    %take the FFT
    FYpredict = fft(Ypredict);
    FYpredict_s = fftshift(FYpredict);
    FYpredict_filt = fft(Ypredict_filt);
    FYpredict_filt_s = fftshift(FYpredict_filt);
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
    %see fix also in stl_plottrainingresults_FFTPSD_function
    
  
