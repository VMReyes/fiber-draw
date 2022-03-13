function visualize_model(x_test, y_test, net)

reset_net = net.resetState();
y_test_net_pred = reset_net.predict(x_test);

tiledlayout(3,2);
for i = 1:3

    nexttile;
    plot(y_test{i}); hold on;
    plot(y_test_net_pred{i}); hold off;
    title('Actual Vs. Predicted BFD on Testing Data');
    ylim([-0.2 0.2]);
    
    nexttile;
    plot(smoothdata(y_test{i}, 'movmean', 20)); hold on;
    plot(y_test_net_pred{i}); hold off;
    title('Actual (smoothed) Vs. Predicted BFD on Testing Data');
    ylim([-0.075 0.075]);
end

stl_plottrainingresults_FFTPSD_function(y_test{1}, y_test{1}, ...
    y_test_net_pred{1}, y_test_net_pred{1})