load("run_results\architecture_experiment.mat")
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat")

net = deep_lstm.resetState();

combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});
y_pred = net.predict(combined_Xdata, "MiniBatchSize", 1);

errors = {};
for b = 1:length(y_pred)
    errors{end+1} = sum((y_pred{b} - combined_Ydata{b}).^2);
end
plot(errors);