load("alldatatrain\all_data_processed_4in_1out_yremove125.mat", "Xdata", "Ydata", ...
     "filenames")   

% combine data
combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});

[~, num_batches] = size(combined_Xdata);
[train_ind, test_ind] = dividerand(num_batches, 0.9, 0.1, 0.0);

x_train = combined_Xdata(train_ind);
y_train = combined_Ydata(train_ind);
x_test = combined_Xdata(test_ind);
y_test = combined_Ydata(test_ind);

% for every architecture, train);

simple_lstm = create_simple_lstm(x_train, y_train);
[simple_lstm, simple_lstm_info] = train_lstm(simple_lstm, x_train, y_train, x_test, y_test, 64);

simple_bilstm = create_simple_bilstm(x_train, y_train);
[simple_bilstm, simple_bilstm_info] = train_lstm(simple_bilstm, x_train, y_train, x_test, y_test, 16);

side_by_side_3_lstm = create_3_side_by_side_lstm(x_train, y_train);
[side_by_side_3_lstm, side_by_side_3_lstm_info] = train_lstm(side_by_side_3_lstm, x_train, y_train, x_test, y_test, 32);


save("alldatatrain\architecture_models.mat", "simple_lstm", "simple_lstm_info", ...
     "simple_bilstm", "simple_bilstm_info", "side_by_side_3_lstm", ...
     "side_by_side_3_lstm_info");

visualize_model(x_test, y_test, simple_bilstm);