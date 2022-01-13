function [XTrain,YTrain] = stl_prep_training_data(BatchInfo, ...
    x_columns, y_columns, fltLEN, bPlot)
%STL_PREP_TRAINING_DATA Summary of this function goes here
%   Detailed explanation goes here
% generate nImportantColumns

for nWhichBatch = 1:length(BatchInfo) 
    %select and order columns of importance.
    nImportantColumns        = [1 2 3 4 5 6 7]; 
    bMeanRemove             = 0;
    XTrainTRANSPOSE = {};
    YTrainTRANSPOSE = {};

    numsubbatch = length(BatchInfo(nWhichBatch).subbatchIndsSTORE);
    if(numsubbatch>0)
        bSubBatchCounter = 0;
        for iinnWhichSub = 1:numsubbatch
            [Y_sub,Y_filt_sub,x_sub, sample_indexfilt_sub,  meanY_sub, meanX_sub]  = stl_prep_trainingdata_allt_function_Wbfdref_2lstm_Wprefilt( BatchInfo, ...
                                                                STRDEF,...
                                                                nWhichBatch, ...
                                                                iinnWhichSub,...
                                                                nImportantColumns,...
                                                                fltLEN,...
                                                                bPlot, ...
                                                                0, ... %bPlotAllSelectedColumns, ...
                                                                bMeanRemove, ...
                                                                PrefltLEN);
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

