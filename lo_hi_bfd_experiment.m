load("experiments/low_hi_bfd/low_hi_bfd_all_x_y_data.mat", "all_Xdata", "all_Ydata")

nets = {};
infos = {};
target_x_test = 0;
target_y_test = 0;
for i = 1:10
    % get the corresponding data
    Xdata = all_Xdata{i};
    Ydata = all_Ydata{i};
    % concat
    combined_Xdata = cat(2, Xdata{:});
    combined_Ydata = cat(2, Ydata{:});

    % remove 125
    combined_Ydata = cellfun(@(x) x - 125, combined_Ydata, 'un', 0);

    % separate into train / test
    [~, num_batches] = size(combined_Xdata);
    [train_ind, test_ind] = dividerand(num_batches, 0.9, 0.1, 0.0);
    
    x_train = combined_Xdata(train_ind);
    y_train = combined_Ydata(train_ind);
    x_test = combined_Xdata(test_ind);
    y_test = combined_Ydata(test_ind);

    % train model
    deep_lstm = create_deep_lstm(x_train, y_train);
    [deep_lstm, deep_lstm_info] = train_lstm(deep_lstm, x_train, y_train, x_test, y_test, 16);

    nets{end+1} = deep_lstm;
    infos{end+1} = deep_lstm_info;
    if i == 1
        target_x_test = x_test;
        target_y_test = y_test;
    end
end

x_test = target_x_test;
y_test = target_y_test;
error_matrix = zeros(1,10);
for i = 1:10
    net = nets{i};
    net = net.resetState();
    y_pred = net.predict(x_test, 'MiniBatchSize', 1);
    mse = 0;
    for b = 1:length(y_pred)
        squared_loss = sum((y_pred{b} - y_test{b}).^2);
        mean_squared_loss = squared_loss / length(y_pred{b});
        mse = mse + mean_squared_loss;
    end
    mse = mse / length(y_pred);  
    rmse = sqrt(mse);
    error_matrix(1,i) = rmse;
end


save("run_results/lohi_bfd_experiment.mat")