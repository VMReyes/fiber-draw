# fiber-draw

This repository contains various experiments that explore how to use an LSTM network to model the fiber draw process. This codebase contains legacy code written by Brian Anthony. 

## First Steps

Clone or download this repository locally. Navigate to this folder in MATLAB. Before running any code, go to the file list panel on the left of the MATLAB window, select all folders and right click > "Add to Path" > "Selected folders and subfolders".

## Machine Learning Experiment Workflow

An experiment workflow typically consists of two parts: data preparation and model training.
First, we will examine "`process_data.m`".
Second, we will examine "`architecture_experiment.m`"

First, the data must be prepared and the script "`process_data.m`" can accomplish this.
The script contains many options and parameters to set.
For a new user, the main options is `strDataPath`.
This path should point to the root directory of Sterlite data.

Next, we have some parameters that affect the data we extract.
Currently, we define "good data" as "steady state" data and we mathematically define such data with the `loBFD` and `hiBFD` parameters.
Consecutive datapoints that lie between `loBFD` and `hiBFD`  are treated as a sequence of steady state.
At the moment, if the sequence is longer than 2000 timesteps, we save it as a subbatch.
If the sequence is longer than 8000, we make only one subbatch of length 8000.
`x_columns` and `y_columns` select what columns we want to use as input and output values, respectively.
`PrefltLEN` controls how much we smooth the output data.
`yRemove125` removes 125 from the output column.

The script then iterates through the datafiles in the folder specified by `strDataPath`.
It goes through every file, extracts the `BatchInfo` (a data structure containing relevant subbatch metadata.
It then uses this to create a `XTrainTranspose` and `YTrainTranspose` matrix and stores it into `Xdata`.
By the end of the main loop, `Xdata` and `Ydata` contain an element for every file processed.
Next, we combine the data and split it into train and test sets to produce `x_train`, `y_train`, `x_test`, and `y_test`.
We save all of these variables to be consumed by experiments.
This allows us to ensure all experiments use a standard set of data.

Now, we will go over "`architecture_experiment.m`" this script analyzes the performance of 4 different LSTM architectures on the fiber draw data.
First, the script loads "`alldatatrain\all_data_processed_4in_1out_yremove125.mat`" which contains the prepared data for the experiment.
In particular, this data file has been prepared so it only has 4 inputs of importance and a single output (bare fiber diameter).
Additionally, it has had the value 125 removed from the output data, to get the data to a range of around -1 to 1.

It creates the model architectures using the helper functions "`create_simple_lstm`", "`create_deep_lstm`", etc.
Then, it invokes "`train_lstm`" which takes in the split data to train a model.

`train_lstm` sets up the training environment and has many parameters.
At the moment, we only expose `miniBatchSize` but there are many to tune here such as learning rate, epochs, and gradient threshold.
It is currently configured to output the model with the lowest error on the validation data seen during training.
It also outputs the info data structure which contains training performance history and the error of the model it output.

To visualize models, we are currently building up `visualize_model.m`.
It takes in data and a model and performs inference and plots against the true values.

## System Identification Pipeline

In `sys_id_pipeline.m`, the main section contains code for the grid search process and merging models. Helper functions are defined below, implementing the algorithm mentioned in my thesis. 

`load_data()` iterates through all files in the folder, creates subbatches and stores it in a cell array called `all_file_data`. For the first time, this process takes hours to run, so storing it as a `.mat` file would be time-saving. 

`grid_search_oe()` and `grid_search_armax()` functions perform the grid search process in nested for loops. The arguments are controller name (either`'kd'` or `'kt'`), upper bound of the order (set to be 10 in this study), and a boolean indicating if the time- and frequency- domain plots are displayed (useful for debugging but slows the process). These functions return the number of "good" subbatches and the average accuracy for each fixed order. 

The `merge_models()` function similarly takes in the controller name (either`'kd`' or `'kt`'), model structure name (either `'oe'` or `'armax'`), the selected model order after reviewing all models, and the boolean flag for plotting (see above). Of particular importance is the `all_iddata` it produces, which is a data structure containing all subbatches as experiments which `oe()` and `armax()` can fit to. A `.mat` file containing the merged model is saved to the current folder. 

If you'd prefer to see the above in a script-like code instead of wrapped functions, see `sys_id_pipeline_old.m`. It's less clean as portions of the code are copy-pasted many times, but it's easier to debug. 