% parameters
strDataPath         = 'C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\';
strOutputPath       = 'C:\Users\Victor\Desktop\fiber-draw\alldatatrain\';
% BatchInfo Parameters
bXLSLoad = 1;
bPlotAll = 0;
bPlot_each_preform_on_subplot = 1;
bPlot_each_preform_on_subplot_with_inrangesubbatches_ = 1;
loBFD = 115;
hiBFD = 135;
nTower = 48; % The tower number 
subbatchMinLen 	= 2000 * 5; % a batch is the same as a preform, multiple 
subbatchMaxLen  = 8000 * 5; % batches (or preforms) are run, one after the 
                        % other in the tower. a subbatch is defined as a
                        % contiguous region of production
x_columns = ["cpspdactval", "frnpwrmv", "hetubetemp", "pfspdactval"]
y_columns = ["barefibrediadisplay", "tenncmv"]

% Subbatch Parameters
fltLEN = 21; 
bPlot = 0; % Plot batch
PrefltLEN = 1;
limit_subbatches = 0;
yRemove125 = 1;

dataFiles = dir(fullfile(strDataPath, "*.csv"));
Xdata = {};
Ydata = {};
filenames = {};
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

combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});

[~, num_batches] = size(combined_Xdata);
[train_ind, test_ind] = dividerand(num_batches, 0.9, 0.1, 0.0);

x_train = combined_Xdata(train_ind);
y_train = combined_Ydata(train_ind);
x_test = combined_Xdata(test_ind);
y_test = combined_Ydata(test_ind);


save("alldatatrain\all_data_processed_4in_2out_yremove125_DT48_100ms.mat", "Xdata", "Ydata", ...
     "x_train", "y_train", "x_test", "y_test", "filenames")