% parameters
strDataPath         = 'C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\';
strOutputPath       = 'C:\Users\Victor\Desktop\fiber-draw\alldatatrain\';

% BatchInfo Parameters
bXLSLoad = 1;
bPlotAll = 0;
bPlot_each_preform_on_subplot = 1;
bPlot_each_preform_on_subplot_with_inrangesubbatches_ = 1;
loBFD = 124;
hiBFD = 126;
nTower = 48; % The tower number 
subbatchMinLen 	= 2000; % a batch is the same as a preform, multiple 
subbatchMaxLen  = 8000; % batches (or preforms) are run, one after the 
                        % other in the tower. a subbatch is defined as a
                        % contiguous region of production
x_columns = ["cpspdactval", "frnpwrmv", "hetubetemp", "pfspdactval"]
y_columns = ["barefibrediadisplay"]

% Subbatch Parameters
fltLEN = 21; 
bPlot = 0; % Plot batch
PrefltLEN = 1;
limit_subbatches = 0;
yRemove125 = 1;

dataFiles = dir(fullfile(strDataPath, "*.csv"));
y_data_with_varying_filter = {};

filenames = {};

% Create sets of Y data
for PrefltLEN = [1 3 5 7]
    Xdata = {};
    Ydata = {};
    for i = 1:14 % Only tower 48 data
        strDataFilename = dataFiles(i).name
        % get batch info
        [BatchInfo, STRDEF] = stl_load_batchinfo(bXLSLoad, strDataPath, ...
            strDataFilename, nTower, bPlotAll, ...
            bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
            bPlot_each_preform_on_subplot, loBFD, hiBFD, subbatchMinLen, subbatchMaxLen, ...
            x_columns, y_columns);
        
        % turn it into a train / test array
        [XTrainTranspose, YTrainTranspose] = stl_prep_training_data(BatchInfo, ...
            STRDEF, x_columns, y_columns, fltLEN, PrefltLEN, bPlot, limit_subbatches, yRemove125);
        Xdata{end+1} = XTrainTranspose;
        Ydata{end+1} = YTrainTranspose;
        filenames{end+1} = strDataFilename;
    end
    y_data_with_varying_filter{end+1} = Ydata;
end

[train_ind, test_ind] = dividerand(424, 0.9, 0.1, 0.0);

% Use data to train an equal number of models
nets = {};
infos = {};
for i = 1:4
    Ydata = y_data_with_varying_filter{i};
    combined_Xdata = cat(2, Xdata{:});
    combined_Ydata = cat(2, Ydata{:});
    
    
    x_train = combined_Xdata(train_ind);
    y_train = combined_Ydata(train_ind);
    x_test = combined_Xdata(test_ind);
    y_test = combined_Ydata(test_ind);

    deep_lstm = create_deep_lstm(x_train, y_train);
    [deep_lstm, deep_lstm_info] = train_lstm(deep_lstm, x_train, y_train, x_test, y_test, 16);
    nets{end+1} = deep_lstm;
    infos{end+1} = deep_lstm_info;
end

unfiltered_y = y_data_with_varying_filter{1};
unfiltered_y = cat(2, unfiltered_y{:});
unfiltered_y_test = unfiltered_y(test_ind);

visualize_model(x_test, unfiltered_y_test, nets{1});
visualize_model(x_test, unfiltered_y_test, nets{2});
visualize_model(x_test, unfiltered_y_test, nets{3});
visualize_model(x_test, unfiltered_y_test, nets{4});

error_matrix = zeros(4,4);
for i = 1% 1:4 % model filter
    net = nets{i}; 
    for j = 1%1:4 % data filter
        input_data = x_test;
        output_data = y_data_with_varying_filter{j};
        output_data = cat(2, output_data{:});
        output_data = output_data(test_ind);
        mse = 0;
        for b = 1:length(input_data)
            net = net.resetState();
            model_prediction = net.predict(input_data{b});
            squared_loss = sum((output_data{b} - model_prediction).^2);
            mse = mse + squared_loss;
        end
        squared_loss_total = mse;
        error_matrix(i,j) = squared_loss_total;
    end
end

save("run_results/filter_len_experiment.mat");

%% scratch
for i = 2% 1:4 % model filter
    net = nets{i}; 
    for j = 4%1:4 % data filter
        input_data = x_test;
        output_data = y_data_with_varying_filter{j};
        output_data = cat(2, output_data{:});
        output_data = output_data(test_ind);
        figure(1); plot(output_data{b}); hold on;
        plot(model_prediction); hold off; 
        latexify_plot
    end
end