function visualize_model(x_test, y_test, net)

reset_net = net.resetState();
y_test_net_pred = reset_net.predict(x_test, "MiniBatchSize",1);

% tiledlayout(3,2);
% for i = 1:3
tiledlayout(2,1)
batch = 5
for i = batch:batch

%     nexttile;
%     plot(y_test{i}); hold on;
%     plot(y_test_net_pred{i}); hold off;
%     title('Actual Vs. Predicted BFD on Testing Data');
%     ylim([-0.2 0.2]);
    
    nexttile;
    plot(smoothdata(y_test{i}, 'movmean', 20),'b'); hold on;
    plot(y_test_net_pred{i},'r'); hold off;
    title('BFD as Predicted by Side-by-Side LSTM Network');
    legend('Actual (smoothed)', 'Prediction')
    ylim([-0.075 0.075]);
end

% stl_plottrainingresults_FFTPSD_function(y_test{1}, y_test{1}, ...
%     y_test_net_pred{1}, y_test_net_pred{1})
nexttile;
plot_fft(y_test{1}, y_test{1}, ...
    y_test_net_pred{1}, y_test_net_pred{1})
end

function plot_fft(Y, Y_filt, Ypredict, Ypredict_filt)
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

    plot(frq_discretes,log10(abs(FY_s).^2),'b');
    hold on
    plot(frq_discretes,log10(abs(FYpredict_s).^2),'r'); 
    hold off
    
    ylabel('Log Magnitude')
    xlabel('Frequency')
    legend('Data', 'Prediction')
    axis([0 1 -.1  max([   max(log10(abs(FY_s).^2))   max(log10(abs(FYpredict_s).^2))   ])   ])
    title('Power Spectrum');
    ylim([-2 3.5])
end