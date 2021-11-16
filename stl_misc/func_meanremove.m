function [Yout, meanY] = func_meanremove(Y)

% mean remove on Y and Y_filt
meanY   = mean(Y);
Yout    = Y- meanY;
