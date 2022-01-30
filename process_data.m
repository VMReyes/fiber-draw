% parameters
strDataPath         = 'C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\';
strOutputPath       = 'C:\Users\Victor\Desktop\fiber-draw\alldatatrain\';
strDataFilename = 'DrawData_Tower48_2020-12-01_to2020-12-08.csv';

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
        STRDEF, x_columns, y_columns, fltLEN, PrefltLEN, bPlot, 0);
    Xdata{end+1} = XTrainTranspose;
    Ydata{end+1} = YTrainTranspose;
    filenames{end+1} = strDataFilename;
end
save("alldatatrain\all_data_processed_4in_1out.mat", "Xdata", "Ydata", "filenames")