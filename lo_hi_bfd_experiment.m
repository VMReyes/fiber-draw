load("experiments/low_hi_bfd/low_hi_bfd_all_x_y_data.mat", "all_Xdata", "all_Ydata")

nets = {};
infos = {};
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
    simple_bilstm = create_simple_bilstm(x_train, y_train);
    [simple_bilstm, simple_bilstm_info] = train_lstm(simple_bilstm, x_train, y_train, x_test, y_test, 16);

    nets{end+1} = simple_bilstm;
    infos{end+1} = simple_bilstm_info;
end

