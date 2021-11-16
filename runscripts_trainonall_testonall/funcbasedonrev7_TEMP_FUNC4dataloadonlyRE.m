function [XTrainTRANSPOSE_ARRAY, YTrainTRANSPOSE_ARRAY] = funcbasedonrev7_TEMP_FUNC4dataloadonlyRE(filenamebase,PrefltLEN)


%% stl_load_parse_allt
bXLSLoad                                                = 1;
bPlotAll                                                = 0;
bPlot_each_preform_on_subplot                           = 0;
bPlot_each_preform_on_subplot_with_inrangesubbatches_   = 0;

loBFD           = 124;
hiBFD           = 126;
% a batch is the same as a preform, multiple batches (or preforms) are run,
% one after the other in in the tower.
% a subbatch is defined as a contiguous region of production
subbatchMinLen 	= 2000; 
subbatchMaxLen  = 8000;

strDataPath     = 'E:\Dropbox (SquareCircleMITtoo)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
% %xls_file        = 'DrawData_Tower51_2021-01-05_to2021-01-12.csv';
% filenamebase =  'DrawData_Tower51_2020-12-22_to2020-12-29';
xls_file        = [filenamebase '.csv'];
 nTower          =  0;



%data00, dataOTHER, uniPreformID to file
[BatchInfo, STRDEF,~, data00, dataOTHER, uniPreformID ]     =  ...
    stl_load_parse_allt_function_rev4(...
        bXLSLoad, strDataPath, xls_file, ...
        nTower, ...
        bPlotAll, bPlot_each_preform_on_subplot, bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
        loBFD, hiBFD , subbatchMinLen, subbatchMaxLen);

%%
% bRNNnarx = 1;
% 
% PrefltLEN = 1;


bSubBatchCounter = 0;

for nWhichBatch = 1:length(BatchInfo)
%for nWhichBatch = [3 5 19] % for 48 1:length(BatchInfo)
%for nWhichBatch = 1:4
%for nWhichBatch = 1:10
%for nWhichBatch = 2
    nWhichSub   = -1; 
    
    %select and order columns of importance.
    %nImportantColumns        = [6 16 15 10   22    26 29]; %removing 11 it is an output
    nImportantColumns        = [1 2 3 4 5 6 7]; 
    fltLEN                  = 21;
    bPlot                   = 0;
    bPlotAllSelectedColumns = 0;
    bMeanRemove             = 1;
    
    XTrainTRANSPOSE = {};
    YTrainTRANSPOSE = {};
    
    dsiddata = [];
    dsiddata_filt = [];
     fitsys1 = [];
     fitsys1_filt =[];     
     fitnet = [];
     fitnet_filt =[];
    
    numsubbatch = length(BatchInfo(nWhichBatch).subbatchIndsSTORE);
    if(numsubbatch>0)
        [Y,Y_filt,x, sample_indexfilt,  meanY, meanX]  = stl_prep_trainingdata_allt_function_Wbfdref_2lstm_Wprefilt( BatchInfo, ...
                                                            STRDEF,...
                                                            nWhichBatch, ...
                                                            nWhichSub,...
                                                            nImportantColumns,...
                                                            fltLEN,...
                                                            bPlot, ...
                                                            bPlotAllSelectedColumns, ...
                                                            bMeanRemove,...
                                                            PrefltLEN);
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

        %[rrr,ccc] = size(XTrainTRANSPOSE{1});
        [rrr,ccc] = size(x_sub');
        %% 
        % Similarly, create a shorter validation signal to use during network training.

        xval = [];%idinput(valsignalLength,signalType);
        yval = [];%lsim(fourthOrderMdl,xval,trgs(1:valsignalLength));



    end

    XTrainTRANSPOSE_ARRAY{nWhichBatch} = XTrainTRANSPOSE;
    YTrainTRANSPOSE_ARRAY{nWhichBatch} = YTrainTRANSPOSE;
    

end

