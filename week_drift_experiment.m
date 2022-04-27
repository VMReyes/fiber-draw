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
x_columns = ["cpspdactval", "frnpwrmv", "hetubetemp", "pfspdactval"];
y_columns = ["barefibrediadisplay"]; % keep bfd as first output for batch
                                     % creation code to work

% Subbatch Parameters
fltLEN = 21; 
bPlot = 0; % Plot batch
PrefltLEN = 1; 

dataFiles = dir(fullfile(strDataPath, "*.csv"));

Xdata = {};
Ydata = {};
for i = 1:15 % Only tower 48 data, 2 months
    strDataFilename = dataFiles(i).name;

    % get batch info
    [BatchInfo, STRDEF] = stl_load_batchinfo(bXLSLoad, strDataPath, ...
        strDataFilename, nTower, bPlotAll, ...
        bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
        bPlot_each_preform_on_subplot, loBFD, hiBFD, subbatchMinLen, subbatchMaxLen, ...
        x_columns, y_columns);
    
    % turn it into a train / test array
    limitSubbatches = 1;
    [XTrainTranspose, YTrainTranspose] = stl_prep_training_data(BatchInfo, ...
        STRDEF, x_columns, y_columns, fltLEN, PrefltLEN, bPlot, limitSubbatches);
    
    Xdata{end+1} = XTrainTranspose;
    Ydata{end+1} = YTrainTranspose;
end

nets = {};
for i = 1:length(Xdata)
    nets{end+1} = train_simple_lstm(Xdata{i}, Ydata{i});
end

error_matrix = zeros(length(nets));

for i = 1:length(nets)
    curr_net = nets{i};
    for j = 1:length(nets)
        input_data = Xdata{j};
        output_data = Ydata{j};
        mse = 0;
        for b = 1:length(input_data)
            curr_net = curr_net.resetState();
            model_prediction = curr_net.predict(input_data{b});
            squared_loss = sum((output_data{b} - model_prediction).^2);
            mse = mse + squared_loss;
        end
        mse = mse / length(input_data);
        error_matrix(i,j) = mse;
    end
end

save("alldatatrain\week_drift.mat", "nets", "Xdata", "Ydata", "error_matrix");


