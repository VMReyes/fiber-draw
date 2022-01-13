strDataPath         = 'C:\Users\Victor\Desktop\fiber-draw\MIT_DrawData_48and51\';

strOutputPath       = 'C:\Users\Victor\Desktop\fiber-draw\alldatatrain\';

% filenamebasearray{1} = 'DrawData_Tower51_2020-12-01_to2020-12-08';  
% filenamebasearray{2} = 'DrawData_Tower51_2020-12-08_to2020-12-15';  
% filenamebasearray{3} = 'DrawData_Tower51_2020-12-15_to2020-12-22';  
% filenamebasearray{4} = 'DrawData_Tower51_2020-12-22_to2020-12-29';  
% filenamebasearray{5} = 'DrawData_Tower51_2021-01-05_to2021-01-12';  

filenamebasearray{1} = 'DrawData_Tower48_2020-12-01_to2020-12-08';  
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

PrefltLEN_array = [1]; 

for fn = 1:length(filenamebasearray)
    filenamebase = filenamebasearray{fn};
    for PrefltLEN = PrefltLEN_array
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

        xls_file        = [filenamebase '.csv'];
        nTower          =  48;



        nnin_LearnRateDropPeriod    = 200;
        nnin_maxEpochs              = 250 ; %150; %500; %1000
        nnin_miniBatchSize          = 200;%200;

        clear BatchInfo;
        %data00, dataOTHER, uniPreformID to file
        [BatchInfo, STRDEF,~, data00, dataOTHER, uniPreformID ]     =  ...
            stl_load_parse_allt_function_rev4(...
                bXLSLoad, strDataPath, xls_file, ...
                nTower, ...
                bPlotAll, bPlot_each_preform_on_subplot, bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
                loBFD, hiBFD , subbatchMinLen, subbatchMaxLen);


        tempstr = ['rev7_' filenamebase '__trALL' '_maxE' num2str(nnin_maxEpochs) '_mx' num2str(subbatchMaxLen) '_drop' num2str(nnin_LearnRateDropPeriod) '_lstm350lstm0' '_filter' num2str(PrefltLEN)];
        filenameLock = [tempstr '.txt']; 
        fullfileout = [strOutputPath filenameLock];

        bLocked = 0;
        if(exist(fullfileout,'file') == 2)
            sprintf('lock exists');
            bLocked = 1;
        else
            save(fullfileout,'-ascii','nTower');
        end

        if(~bLocked)

            bSubBatchCounter = 0;

            clear XTrainTRANSPOSE_fromONEfile;
            clear YTrainTRANSPOSE_fromONEfile;
            for nWhichBatch = 1:length(BatchInfo)
                nWhichSub   = -1; 

                %select and order columns of importance.
                nImportantColumns        = [1 2 3 4 5 6 7]; 
                fltLEN                  = 21;
                bPlot                   = 0;
                bPlotAllSelectedColumns = 1;
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

            numResponses = 1;
            featureDimension = rrr; %1;
            %
            %numHiddenUnits = 100;
            maxEpochs       = nnin_maxEpochs; %150 ; %150; %500; %1000
            miniBatchSize   = nnin_miniBatchSize; %200;%200;

            Networklayers = [sequenceInputLayer(featureDimension, "Normalization","zscore") ...
                lstmLayer(100) ...
                fullyConnectedLayer(numResponses) ...
                regressionLayer];

            %% 
            % The initial learning rate impacts the success of the network. Using an initial 
            % learning rate that is too high results in high gradients, which lead to longer 
            % training times. Longer training times can lead to saturation of the fully connected 
            % layer of the network. When the network saturates, the outputs diverge and the 
            % network outputs a |NaN| value. Hence, use the default value of 0.01, which is 
            % a relatively low initial learning rate. This results in a monotonically decreasing 
            % residual and loss curves. Use a piecewise rate schedule to keep the optimization 
            % algorithm from getting trapped in local minima at the start of the optimization 
            % routine.

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

            trainedNetworkModel = trainNetwork(XTrainTRANSPOSE_fromONEfile,YTrainTRANSPOSE_fromONEfile,Networklayers,options);


            clear XTrainTRANSPOSE_ARRAY;
            clear YTrainTRANSPOSE_ARRAY;
            clear trainedNetworkModel_ARRAY;
            trainedNetworkModel_ARRAY{1}    = trainedNetworkModel;
            XTrainTRANSPOSE_ARRAY{1}        = XTrainTRANSPOSE_fromONEfile;
            YTrainTRANSPOSE_ARRAY{1}        = YTrainTRANSPOSE_fromONEfile;

            filenameRslt= [tempstr '.mat']; 

            save([strOutputPath filenameRslt]);

        end
    end

end
    
%% plot
close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
trainedNetworkModel = trainedNetworkModel.resetState();
y_pred = trainedNetworkModel.predict(x_sub');
Y_sub = Y_sub;
y_pred = y_pred';
figure; plot(Y_sub(:,1)); hold on; plot(y_pred(:,1),'r'); 
legend({'Data','Prediction'});

stl_plottrainingresults_FFTPSD_function(Y_sub(:,1), Y_sub(:,1), y_pred(:,1), y_pred(:,1));


figure(3000); xlim([0 1]); ylim([-5 4]); 
figure(3001); grid minor;      

