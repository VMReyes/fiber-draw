load("alldatatrain\all_data_processed_4in_1out_yremove125_DT48_on100ms.mat", "Xdata", "Ydata", ...
     "x_train", "y_train", "x_test", "y_test", "filenames")   

% % combine data
% combined_Xdata = cat(2, Xdata{:});
% combined_Ydata = cat(2, Ydata{:});
% 
% [~, num_batches] = size(combined_Xdata);
% [train_ind, test_ind] = dividerand(num_batches, 0.9, 0.1, 0.0);
% 
% x_train = combined_Xdata(train_ind);
% y_train = combined_Ydata(train_ind);
% x_test = combined_Xdata(test_ind);
% y_test = combined_Ydata(test_ind);

% for every architecture, train;

% simple_lstm = create_simple_lstm(x_train, y_train);
% [simple_lstm, simple_lstm_info] = train_lstm(simple_lstm, x_train, y_train, x_test, y_test, 8);

downsampled_x_train = cell(size(x_train));
downsampled_y_train = cell(size(y_train));
downsampled_x_test = cell(size(x_test));
downsampled_y_test = cell(size(y_test));

for i = 1:length(downsampled_x_train)
    downsampled_x_train{i} = downsample(x_train{i}',5)';
    downsampled_y_train{i} = downsample(y_train{i}',5)';
end
for i = 1:length(downsampled_x_test)
    downsampled_x_test{i} = downsample(x_test{i}',5)';
    downsampled_y_test{i} = downsample(y_test{i}',5)';
end

% downsampled_x_train = cellfun(downsample(x_train, 5), x_train);
% downsampled_y_train = downsample(y_train, 5);
% downsampled_x_test = downsample(x_test, 5);
% downsampled_y_test = downsample(y_test, 5);
downsampled_deep_lstm = create_deep_lstm(downsampled_x_train, downsampled_y_train);
[deep_lstm_downsampled, deep_lstm_downsampled_info] = train_lstm(downsampled_deep_lstm, downsampled_x_train, downsampled_y_train, downsampled_x_test, downsampled_y_test, 16);


deep_lstm = create_deep_lstm(x_train, y_train);
[deep_lstm, deep_lstm_info] = train_lstm(deep_lstm, x_train, y_train, x_test, y_test, 2);

deep_lstm_batchnorm = create_deep_lstm_batchnorm(x_train, y_train);
[deep_lstm_batchnorm, deep_lstm_batchnorm_info] = train_lstm(deep_lstm_batchnorm, x_train, y_train, x_test, y_test, 2);

% simple_bilstm = create_simple_bilstm(x_train, y_train);
% [simple_bilstm, simple_bilstm_info] = train_lstm(simple_bilstm, x_train, y_train, x_test, y_test, 2);
% 
% side_by_side_3_lstm = create_3_side_by_side_lstm(x_train, y_train);
% [side_by_side_3_lstm, side_by_side_3_lstm_info] = train_lstm(side_by_side_3_lstm, x_train, y_train, x_test, y_test, 4);

save("results/architecture_experiment_4in_1out_DT48_on100ms.mat", "deep_lstm", "deep_lstm_info", "deep_lstm_batchnorm", "deep_lstm_batchnorm_info", "deep_lstm_downsampled", "deep_lstm_downsampled_info");
% save("results/architecture_experiment_4in_1out_DT48_on100ms.mat", "simple_lstm", "simple_lstm_info", ...
%      "deep_lstm", "deep_lstm_info", "simple_bilstm", "simple_bilstm_info", "side_by_side_3_lstm", ...
%      "side_by_side_3_lstm_info");