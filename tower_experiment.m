load("all_file_data.mat", 'all_file_data');
load("architecture_experiment.mat", "deep_lstm");
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat")

net = deep_lstm.resetState();

combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});

y_pred = net.predict(combined_Xdata, "MiniBatchSize", 1);

errors = zeros(1,length(y_pred));
for b = 1:length(y_pred)
    errors(b) = sqrt(sum((y_pred{b} - combined_Ydata{b}).^2/length(y_pred{b})));
end

figure; 
plot(errors); title('Model RMSE over Time'); 
ylabel('RMSE'); xlabel('Subbatch')
xlim([0 length(errors)])
latexify_plot;


