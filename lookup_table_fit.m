%% demo
clc; clear; close all;
slope = [0 0.01 0.2 0.3 0.5 0.8 1 1.5 2 3 4];
speed100 = [4.2 4 2.5 1.6 0 -3 -5 -8 -10 -10 -10];
speed130 = [2.5 2 1 0.5 -1 -2 -3 -4 -8 -10 -10];

slopefit100 = slope(1:end-2);
speedfit100 = speed100(1:end-2);
slopefit130 = slope(1:end-1);
speedfit130 = speed130(1:end-1);

% [P100,S100,MU100] = polyfit(slopefit100, speedfit100, 1);
% [Y100,DELTA100] = polyval(P100,slopefit100,S100,MU100);

% [P130,S130,MU130] = polyfit(slopefit130, speedfit130, 1);
% [Y130,DELTA130] = polyval(P100,slopefit130,S130,MU130);
ft = fittype('a*exp(-b*x)+c','independent','x');
[F100, G100] =  fit(slopefit100', speedfit100', ft);
[F130, G130] =  fit(slopefit130', speedfit130', ft);

y100 = max(F100(slope), -10);
y130 = max(F130(slope), -10);

figure; hold on; plot(slope,speed100, 'ro'); plot(slope, speed130, 'k^');
plot(slope, y100, 'r'); plot(slope, y130,'k');
grid minor; title('Exponential Regression of Feed Speed Lookup Table');
xlabel('Slope of Capstan Speed'); ylabel('Feed Speed @ 2700 mm/min');
legend({'100mm Preform', '130mm Preform'}); ylim([-12 5])
customPlot