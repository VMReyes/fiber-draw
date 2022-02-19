function [XTrain,YTrain] = stl_prep_training_data(BatchInfo, STRDEF, ...
    x_columns, y_columns, fltLEN, PrefltLEN, bPlot, limit_subbatches, yRemove125)
%STL_PREP_TRAINING_DATA Summary of this function goes here
%   Detailed explanation goes here
% generate nImportantColumns
num_batches = length(BatchInfo);

if limit_subbatches
    num_batches = min(num_batches, 12);
end

for nWhichBatch = 1:num_batches
    %select and order columns of importance.
    nNumColumns = length(x_columns) + length(y_columns);
    nImportantColumns        = 1:nNumColumns; 
    bMeanRemove             = 0;
    XTrainTRANSPOSE = {};
    YTrainTRANSPOSE = {};

    numsubbatch = length(BatchInfo(nWhichBatch).subbatchIndsSTORE);
    if(numsubbatch>0)
        bSubBatchCounter = 0;
        for iinnWhichSub = 1:numsubbatch
            [Y_sub,~,x_sub, ~,  ~, ~]  = stl_prep_trainingdata_allt_function_Wbfdref_2lstm_Wprefilt( BatchInfo, ...
                                                                STRDEF,...
                                                                nWhichBatch, ...
                                                                iinnWhichSub,...
                                                                nImportantColumns,...
                                                                fltLEN,...
                                                                bPlot, ...
                                                                0, ... %bPlotAllSelectedColumns, ...
                                                                bMeanRemove, ...
                                                                PrefltLEN, ...
                                                                x_columns, ...
                                                                y_columns, ...
                                                                yRemove125);
            bSubBatchCounter = bSubBatchCounter + 1;
            XTrainTRANSPOSE{bSubBatchCounter}=x_sub';
            YTrainTRANSPOSE{bSubBatchCounter}=Y_sub';
        end
    end
    XTrainTRANSPOSE_ARRAY{nWhichBatch} = XTrainTRANSPOSE;
    YTrainTRANSPOSE_ARRAY{nWhichBatch} = YTrainTRANSPOSE;
end

%all data from file
XTrainTRANSPOSE_fromONEfile = [XTrainTRANSPOSE_ARRAY{:}];
YTrainTRANSPOSE_fromONEfile = [YTrainTRANSPOSE_ARRAY{:}];
XTrain = XTrainTRANSPOSE_fromONEfile;
YTrain = YTrainTRANSPOSE_fromONEfile;

