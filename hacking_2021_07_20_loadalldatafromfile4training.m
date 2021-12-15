function [model, XTrain, YTrain] = hacking_2021_07_20_loadalldatafromfile4training(dataFilePath, outPath, prefltLEN)
    close all; 
    [path, name, ext] = fileparts(dataFilePath);
    strDataPath         = path;
    strOutputPath       = outPath;
  
    
    filenamebasearray{1} = name+ext;  
    %filenamebasearray{2} = 'DrawData_Tower48_2020-12-08_to2020-12-15';  
    %filenamebasearray{1} = 'DrawData_Tower48_2020-12-15_to2020-12-22'; 
    %filenamebasearray{2} = 'DrawData_Tower48_2020-12-22_to2020-12-29'; 
    % filenamebasearray{1} = 'DrawData_Tower48_2020-12-29_to2021-01-05'; 
    % filenamebasearray{2} = 'DrawData_Tower48_2021-01-05_to2021-01-12';  
    % filenamebasearray{1} = 'DrawData_Tower48_2021-01-12_to2021-01-19';  
    % filenamebasearray{2} = 'DrawData_Tower48_2021-02-02_to2021-02-09';  
    % filenamebasearray{3} = 'DrawData_Tower48_2021-02-23_to2021-03-02';  
    %filenamebasearray{2} = 'DrawData_Tower48_2021-03-02_to2021-03-09';  
    %filenamebasearray{1} = 'DrawData_Tower48_2021-03-09_to2021-03-16';  
    %filenamebasearray{1} = 'DrawData_Tower48_2021-03-16_to2021-03-23';  
    
    % PrefltLEN_array = [1 5 7]; 
    prefltLEN
    PrefltLEN = prefltLEN;
    %% stl_load_parse_allt
    bXLSLoad                                                = 1;
    bPlotAll                                                = 0;
    bPlot_each_preform_on_subplot                           = 1;
    bPlot_each_preform_on_subplot_with_inrangesubbatches_   = 1;

    loBFD           = 124;
    hiBFD           = 126;
    % a batch is the same as a preform, multiple batches (or preforms) are run,
    % one after the other in in the tower.
    % a subbatch is defined as a contiguous region of production
    subbatchMinLen 	= 2000; 
    subbatchMaxLen  = 8000;

    % strDataPath     = 'E:\Dropbox (SquareCircleMITtoo)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
    %xls_file        = 'DrawData_Tower51_2021-01-05_to2021-01-12.csv';
    %filenamebase =  'DrawData_Tower51_2021-01-05_to2021-01-12';
    %filenamebase =  'DrawData_Tower51_2020-12-22_to2020-12-29';
    xls_file        = name+ext;
    nTower          =  48;

    nnin_LearnRateDropPeriod    = 200;
    nnin_maxEpochs              = 250 ; %150; %500; %1000
    nnin_miniBatchSize          = 200;%200;

    clear BatchInfo;
    [BatchInfo, STRDEF,~, data00, dataOTHER, uniPreformID ]     =  ...
        stl_load_parse_allt_function_rev4(...
            bXLSLoad, strDataPath, xls_file, ...
            nTower, ...
            bPlotAll, bPlot_each_preform_on_subplot, bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
            loBFD, hiBFD , subbatchMinLen, subbatchMaxLen);


    tempstr = ['rev7_' name '__trALL' '_maxE' num2str(nnin_maxEpochs) '_mx' num2str(subbatchMaxLen) '_drop' num2str(nnin_LearnRateDropPeriod) '_lstm350lstm0' '_filter' num2str(PrefltLEN)];
    filenameLock = [tempstr '.txt']; 
    fullfileout = [strOutputPath filenameLock];

%     bLocked = 0;
%     if(exist(fullfileout,'file') == 2)
%         sprintf('lock exists');
%         bLocked = 1;
%     else
%         save(fullfileout,'-ascii','nTower');
%     end
%         if(~bLocked)
    if (1)

        bSubBatchCounter = 0;

        clear XTrainTRANSPOSE_fromONEfile;
        clear YTrainTRANSPOSE_fromONEfile;
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
            bPlotAllSelectedColumns = 1;
            bMeanRemove             = 0;

            XTrainTRANSPOSE = {};
            YTrainTRANSPOSE = {};

            XTrainTRANSPOSEWithMean = {};
            YTrainTRANSPOSEWithMean = {};


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


        %% Train on all data

        [rrr,ccc] = size(XTrainTRANSPOSE_fromONEfile{1});

        %% 
        % Similarly, create a shorter validation signal to use during network training.

        xval = [];%idinput(valsignalLength,signalType);
        yval = [];%lsim(fourthOrderMdl,xval,trgs(1:valsignalLength));
        %% Create and Train Network
        % The following network architecture was determined by using a Bayesian optimization 
        % routine where the Bayesian optimization cost function uses independent validation 
        % data (see the accompanying |bayesianOptimizationForLSTM.mlx| for the details). 
        % Although multiple architectures may work, this optimization provides the most 
        % computationally efficient one. The optimization process also showed that as 
        % the complexity of the transfer function increases when applying LSTM to other 
        % linear transfer functions, the architecture of the network does not change significantly. 
        % Rather, the number of epochs needed to train the network increases. The number 
        % of hidden units required for modeling a system is related to how long the dynamics 
        % take to damp out. In this case there are two distinct parts to the response: 
        % a high frequency response and a low frequency response. A higher number of hidden 
        % units are required to capture the low frequency response. If a lower number 
        % of units are selected the high frequency response is still modeled. However, 
        % the estimation of the low frequency response deteriorates. 
        % 
        % Create the network architecture.

        numResponses = 4; % TODO(gcfchen): this was originally 1, modified to be 4
        featureDimension = rrr; %1;
        %
        %numHiddenUnits = 100;
        maxEpochs       = nnin_maxEpochs; %150 ; %150; %500; %1000
        miniBatchSize   = nnin_miniBatchSize; %200;%200;
                     
%             lgraph = layerGraph
%             tempLayers = sequenceInputLayer(featureDimension, "Name", "sequence");
%             lgraph = addLayers(lgraph, tempLayers);
%             
%             tempLayers = lstmLayer(128, "Name", "lstm_2");
%             lgraph = addLayers(lgraph, tempLayers);
%             
%             tempLayers = lstmLayer(128, "Name", "lstm_3");
%             lgraph = addLayers(lgraph, tempLayers);
%             
%             tempLayers = lstmLayer(128, "Name", "lstm_1");
%             lgraph = addLayers(lgraph, tempLayers);
%             
%             tempLayers = [
%                 concatenationLayer(1,3, "Name", "concat")
%                 fullyConnectedLayer(numResponses, "Name", "fc")
%                 regressionLayer("Name", "regressionoutput")];
%             lgraph = addLayers(lgraph,tempLayers);
%             clear tempLayers;
%             lgraph = connectLayers(lgraph,"sequence","lstm_2");
%             lgraph = connectLayers(lgraph,"sequence","lstm_3");
%             lgraph = connectLayers(lgraph,"sequence","lstm_1");
%             lgraph = connectLayers(lgraph,"lstm_2","concat/in2");
%             lgraph = connectLayers(lgraph,"lstm_3","concat/in3");
%             lgraph = connectLayers(lgraph,"lstm_1","concat/in1");

       Networklayers = [sequenceInputLayer(featureDimension) ...
           lstmLayer(100) ...
           fullyConnectedLayer(numResponses) ...
           regressionLayer];

        options = trainingOptions('adam', ...
            'MaxEpochs',maxEpochs, ...
            'MiniBatchSize',miniBatchSize, ...
            'GradientThreshold',10, ...
            'SequenceLength','longest', ...             % padding
            'SequencePaddingDirection', 'right', ...    %
            'SequencePaddingValue',0, ...               %
            'Shuffle','once', ...
            'Plots','training-progress',...
            'LearnRateSchedule','piecewise',...
            'LearnRateDropPeriod',nnin_LearnRateDropPeriod,... %was 100
            'Verbose',1);%,...
            %'ValidationData',[{xval'} {yval'}]);

        % options = trainingOptions('adam', ...
        %     'MaxEpochs',maxEpochs, ...
        %     'MiniBatchSize',miniBatchSize, ...
        %     'GradientThreshold',10, ...
        %     'Shuffle','once', ...
        %     'Plots','training-progress',...
        %     'ExecutionEnvironment','gpu',...
        %     'LearnRateSchedule','piecewise',...
        %     'LearnRateDropPeriod',100,...
        %     'Verbose',0,...
        %     'ValidationData',[{xval'} {yval'}]);


        %poolobj = parpool;
        trainedNetworkModel = trainNetwork(XTrainTRANSPOSE_fromONEfile,YTrainTRANSPOSE_fromONEfile,Networklayers,options);
        model = trainedNetworkModel;
        XTrain = XTrainTRANSPOSE_fromONEfile
        YTrain = YTrainTRANSPOSE_fromONEfile

        %trainedNetworkModel_ARRAY{1}    = trainedNetworkModel;
        %XTrainTRANSPOSE_ARRAY{1}        = XTrainTRANSPOSE_fromONEfile;
        %YTrainTRANSPOSE_ARRAY{1}        = YTrainTRANSPOSE_fromONEfile;

%         filenameRslt= [tempstr '.mat']; 
%         save([strOutputPath filenameRslt]);
%         save("alldatatrainfromonefile\model.mat", "trainedNetworkModel")
    end
%     %% plot
%     close all;
%     set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
%     set(groot, 'defaultLegendInterpreter','latex');
%     set(groot, 'defaultTextInterpreter','latex');
%     
%     y_pred = trainedNetworkModel.predict(x_sub');
%     y_pred = y_pred';
%     figure; plot(Y_sub(:,1)); hold on; plot(y_pred(:,1),'r'); 
%     legend({'Data','Prediction'});
%     
%     stl_plottrainingresults_FFTPSD_function(Y_sub(:,1), Y_sub(:,1), y_pred(:,1), y_pred(:,1));
%     
%     figure(3000); xlim([0 1]); ylim([-5 4]); 
%     figure(3001); grid minor;  

