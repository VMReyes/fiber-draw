load("run_results\architecture_experiment.mat")
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat")

net = deep_lstm.resetState();

% Xdata = [(num_subbatches, 4, 8000), (num_subbatches, 4, 800), .. for
% every file)]
combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});
% combined_Xdata.shape = (num_subbatches, 4 (inputs), 8000 (seq. length))
% bode_example = combined_Xdata{1}
% bode_example.shape = (4, 8000)
% bode_example{2} = sine_wave(amp, freq, len=8000)
% y_pred = net.predict(bode_example, "MiniBatchSize", 1)
y_pred = net.predict(combined_Xdata, "MiniBatchSize", 1);

errors = zeros(1,length(y_pred));
for b = 1:length(y_pred)
    errors(b) = sqrt(sum((y_pred{b} - combined_Ydata{b}).^2/length(y_pred{b})));
%     errors(b) = sum((y_pred{b} - combined_Ydata{b}).^2);
end

figure; 
plot(errors); title('Model RMSE over Time'); 
ylabel('RMSE'); xlabel('Subbatch')
xlim([0 length(errors)])
latexify_plot;