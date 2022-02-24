function visualize_model(x_test, y_test, net)

reset_net = net.resetState();
y_test_net_pred = reset_net.predict(x_test);

tiledlayout(2,1);

nexttile;
plot(y_test{1}); hold on;
plot(y_test_net_pred{1}); hold off;
title(net_name + ' - Actual Vs. Predicted BFD on Testing Data');
ylim([-0.2 0.2]);

nexttile;
plot(abs(y_test{1} - y_test_net_pred{1}));
title(net_name + ' - Actual minus Predicted BFD on Testing Data');
ylim([0 0.2]);