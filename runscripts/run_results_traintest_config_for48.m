
%TestSetName = 'ttresults_a100_000_TRAINON_51AllofEachWeek_TESTON_AllOfOneweek_Tower48_2021-01-05_to2021-01-12__filter5';
%TestDataFile = 'rev7_newtest51_DrawData_Tower48_2021-01-05_to2021-01-12__filter5.mat'; 

% TestSetName = 'ttresults_a100_000_TRAINON_51AllofEachWeek_TESTON_AllOfOneweek_Tower48_2020-12-22_to2020-12-29__filter5';
% TestDataFile = 'rev7_newtest51_DrawData_Tower48_2020-12-22_to2020-12-29__filter5.mat'; 

% FileSet{1} = 'Tower48_2021-01-05_to2021-01-12';
% FileSet{2} = 'Tower48_2020-12-22_to2020-12-29';
% FileSet{3} = 'Tower48_2020-12-15_to2020-12-22';
% FileSet{4} = 'Tower48_2020-12-08_to2020-12-15';
% FileSet{5} = 'Tower48_2020-12-01_to2020-12-08';

FileSet{1} = 'Tower48_2020-12-01_to2020-12-08';  
FileSet{2} = 'Tower48_2020-12-08_to2020-12-15';  
FileSet{3} = 'Tower48_2020-12-15_to2020-12-22'; 
FileSet{4} = 'Tower48_2020-12-22_to2020-12-29'; 
FileSet{5} = 'Tower48_2020-12-29_to2021-01-05'; 
FileSet{6} = 'Tower48_2021-01-05_to2021-01-12';  
FileSet{7} = 'Tower48_2021-01-12_to2021-01-19';  
FileSet{8} = 'Tower48_2021-02-02_to2021-02-09';  
FileSet{9} = 'Tower48_2021-02-23_to2021-03-02';  
FileSet{10} = 'Tower48_2021-03-02_to2021-03-09';  
FileSet{11} = 'Tower48_2021-03-09_to2021-03-16';  
FileSet{12} = 'Tower48_2021-03-16_to2021-03-23';  


FILTERlen = 5; %and '_filter5'
%[XTrainTRANSPOSE_ARRAY, YTrainTRANSPOSE_ARRAY] = funcbasedonrev7_TEMP_FUNC4dataloadonlyRE('DrawData_Tower48_2020-12-01_to2020-12-08',1);
%[XTrainTRANSPOSE_ARRAY, YTrainTRANSPOSE_ARRAY] = funcbasedonrev7_TEMP_FUNC4dataloadonlyRE('DrawData_Tower48_2021-01-05_to2021-01-12',1);

%maxE250_mx8000_drop200_lstm350lstm0    is     b100
%maxE150_mx8000_drop100_lstm100lstm50   is    a100


for iifs = 1:length(FileSet)
    
    TestSetName = ['ttresults_b100_000_TRAINON_48AllofEachWeek_TESTON_AllOfOneweek_' FileSet{iifs} '__filter5']
    %TestDataFile = ['rev7_newtest51_DrawData_' FileSet{iifs} '__filter5.mat'];

    if(1)
        %TestSetName = 'tt_a100_000_TRAINON_51AllofEachWeek_TESTON_AllOfOneweek_Tower48_2021-01-05_to2021-01-12__filter5';

        % For each single ”week” train a model on all data
        % Train on each week (here only one) separately, 
        % Eval: Compare the avg / std for prediction from each week-model on other week

        % pdf ex in rev7_newtest51_DrawData_Tower48_2021-01-05_to2021-01-12__filter5_ERR_TrALLvsALL.pdf


        %
        % test
        % XTrainTRANSPOSE_ARRAY = [];
        % load rev7_newtest51_DrawData_Tower48_2021-01-05_to2021-01-12__filter5.mat
        % % rewrite of the immediate above
        XTrainTRANSPOSE_ARRAY = [];
        YTrainTRANSPOSE_ARRAY = [];
%         %ssLoad = load('rev7_newtest51_DrawData_Tower48_2021-01-05_to2021-01-12__filter5.mat','XTrainTRANSPOSE_ARRAY', 'YTrainTRANSPOSE_ARRAY');
%         ssLoad = load(TestDataFile,'XTrainTRANSPOSE_ARRAY', 'YTrainTRANSPOSE_ARRAY');
%         XTrainTRANSPOSE_ARRAY = ssLoad.XTrainTRANSPOSE_ARRAY;
%         YTrainTRANSPOSE_ARRAY = ssLoad.YTrainTRANSPOSE_ARRAY;
%         %
        [XTrainTRANSPOSE_ARRAY, YTrainTRANSPOSE_ARRAY] = funcbasedonrev7_TEMP_FUNC4dataloadonlyRE(['DrawData_' FileSet{iifs}],FILTERlen);

        %
        % trained model from each week of 51
        trainedNetworkModel_ARRAY = [];
        matfilename ='rev7_DrawData_Tower48_2020-12-01_to2020-12-08__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{1} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2020-12-08_to2020-12-15__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{2} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2020-12-15_to2020-12-22__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{3} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2020-12-22_to2020-12-29__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{4} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2021-01-05_to2021-01-12__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{5} = tmpP{1};
        
        matfilename ='rev7_DrawData_Tower48_2020-12-29_to2021-01-05__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{6} = tmpP{1};
        
        matfilename ='rev7_DrawData_Tower48_2021-01-12_to2021-01-19__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{7} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2021-02-02_to2021-02-09__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{8} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2021-02-23_to2021-03-02__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{9} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2021-03-02_to2021-03-09__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{10} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2021-03-09_to2021-03-16__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{11} = tmpP{1};
        matfilename ='rev7_DrawData_Tower48_2021-03-16_to2021-03-23__trALL_maxE250_mx8000_drop200_lstm350lstm0_filter5.mat';
        tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
        trainedNetworkModel_ARRAY{12} = tmpP{1};
        %
        %
        bPlotAllFreqOverride = 1;

        run_results_traintest;
        save([TestSetName '.mat'])
        
        close all
    end


end