%%
% Sterlite overview script
% Various sections/ block to demonstrate functionality.
% September 16, 2021

%% CURRENT EXPERIMENTS - approximate current set of active expirmentation for training and test is at end of this file


%%  First, here, after directory changes, should run to demo some basic functionality

%% Data parse and load and plot
%
%
% BA notes for context:   
% run_trainon_allpreforms_rev7_2lstm_withnewdataload_TEMP_auto.m
%
%
%

% New style load

strDataPath         = 'E:\Dropbox (SquareCircleMITtoo)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
filenamebase = 'DrawData_Tower48_2021-03-16_to2021-03-23';

        bXLSLoad                                                = 1;
        bPlotAll                                                = 1;
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

        %PrefltLEN = 1;
                
        nnin_LearnRateDropPeriod    = 200;
        nnin_maxEpochs              = 250 ; %150; %500; %1000
        nnin_miniBatchSize          = 200;%200;

%         nnin_LearnRateDropPeriod    = 100;
%         nnin_maxEpochs              = 150 ; %150; %500; %1000
%         nnin_miniBatchSize          = 200;%200;
        
        %data00, dataOTHER, uniPreformID to file
        [BatchInfo, STRDEF,~, data00, dataOTHER, uniPreformID ]     =  ...
            stl_load_parse_allt_function_rev4(...
                bXLSLoad, strDataPath, xls_file, ...
                nTower, ...
                bPlotAll, bPlot_each_preform_on_subplot, bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
                loBFD, hiBFD , subbatchMinLen, subbatchMaxLen);

%% for new data coming from from stl_load_parse_allt_function_rev4 
% these are the important columns nImportantColumns        = [1 2 3 4 5 6 7]; 
% in old styleloading the important columns were not presorted

%%
% Ba notes for context:
% run_scripts_allt_evaleachsb.m
% run_scripts_allt.m
% run_scripts.m

%% some optional figures / plotting
%
%    stl_plot_parsed_subbatches_on_subplots
%    stl_plot_parsed_PSDs
nImportantColumns        = [1 2 3 4 5 6 7];
%
% below is from   stl_plot_parsed_subbatches_on_subplots    -- also see stl_plot_parsed_filtered_subbatches (old)
%  

    %was dev_plot_filtered_subbatches.m
    %
    %run after stl_load_parse (or other older files)

    numBatches = length(BatchInfo);
    BatchRange = [1:numBatches];
    fltLEN = 21;

    m6or8 = floor(sqrt(numBatches));
    n5or4 = ceil(numBatches/m6or8);

    figure;

    % use dataSTORE and BatchInfo to plot each Batch with subBatch regions
    % on seperate figures (or subplot) for each batch
    for i = 1:1
        strlistclr{1} = 'r'; strlistclr{2} = 'g'; strlistclr{3} = 'b'; strlistclr{4} = 'm';
        strlistlegend = [];
        for ii = BatchRange
            %BFD = dataSTORE{ii}(:,6);

            dataaB = BatchInfo(ii).data;
            %select and re-order for ONLY the columns of importance. 
            % (Batch Inportant Reordered)
            dataaBIR = dataaB(:,nImportantColumns);
            BFD = dataaBIR(:,1);

            %filter data 
            [BFD,sample_indexfilt] = func_simplefilter(BFD,fltLEN,'same');
            %figure
            % - or
            subplot(m6or8,n5or4,ii)
            %
            plot(BFD);
            hold on

            strlistlegend   = [];
            tit             = ['Batch ' num2str(ii)]; strlistlegend{1} = tit;

            numsubbatch = length(BatchInfo(ii).subbatchIndsSTORE);
            for jj = 1:numsubbatch
                tempInds = BatchInfo(ii).subbatchIndsSTORE{jj};
                clrstr = strlistclr{func_clrstr(strlistclr,jj)};
                plot(tempInds,BFD(tempInds),clrstr)

                tit  = ['Sub: ' num2str(jj) ' Len: ' num2str(length(tempInds))]; strlistlegend{1+jj} = tit;

            end
            hold off
            title(['Batch ' num2str(ii)])
            legend(strlistlegend, 'Location', 'Best');    
        end

        sgtitle([xls_file '   ' num2str(subbatchMinLen) '   ' num2str(subbatchMaxLen)], 'Interpreter', 'none')

    end

%% Example of putting data (new style) into structs for training
%
% Ba notes for context:
% hacking_2021_07_20_loadalldatafromfile4training_2021_09_15.m

    PrefltLEN  = 1;

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
        bMeanRemove             = 1;

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
                                                                    bPlotAllSelectedColumns, ...
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

    
%% stl_plottrainingresults_FFTPSD_function
% 
% %from ___
% [FY, FY_filt, FYpredict, FYpredict_filt, frq_discretes, ...
%  PY, PY_filt, PYpredict, PYpredict_filt, WW, ...
% ]           = stl_plottrainingresults_FFTPSD_function(Y_sub, Y_sub, modelYOutResponse, modelYOutResponse, -1, -1, 0, [], 'str');
% %]           = stl_plottrainingresults_FFTPSD_function(Y, Y, modelYOutResponse, modelYOutResponse, nWhichBatch, nWhichSub, nWhichSub*100000, [], 'str');

[FY, FY_filt, frq_discretes,...
          PY, PY_filt, WW ] = ... 
          stl_plot_FFTPSD_function(Y_sub,Y_filt_sub, ...
                -1, -1, 0, [], 'str' );
            
            
            
            
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDummyFlag4allowingdepednincies = 0;
if(bDummyFlag4allowingdepednincies)
    %% __trALL
    % hacking_2021_07_20_loadalldatafromfile4training.m
    % hacking_2021_07_20_loadalldatafromfile4training_2021_09_15.m
    % runs training - (on entire week worth of data)
    
    hacking_2021_07_20_loadalldatafromfile4training
    hacking_2021_07_20_loadalldatafromfile4training_2021_09_15
    
    %% plot cross trained/tested results
    % hacking_2021_06_02and09and23and30_and07_14and21.m
    % edited to become run_results_traintest
    % and see hack_subplot_reorg  for loading
    %       and func_tightsubplotreorg
    hacking_2021_06_02and09and23and30_and07_14and21
    
    % run_results_traintest_config will be  into run_results_traintest
    %
    %
    %
    % run_results_traintest_config (fcurrently for 51)
    % run_results_traintest_config_for48.m
    %   calls: 
    %   funcbasedonrev7_TEMP_FUNC4dataloadonlyRE to load the test data for
    %   entire week.
    %   - hacked fast / derived from
    %           run_trainon_allpreforms_rev7_2lstm_withnewdataload_TEMP_auto
    %           and
    %           run_trainon_allpreforms_rev7_2lstm_withnewdataload.m
    %   calls: 
    %   run_results_traintest to generate figures and saves mat files with figure handles
    %     -->> generates: ttresults_b100_000_TRAINON_... or 
    %                      ttresults_a100_000_TRAINON_4 ..
    %
    run_results_traintest_config
    run_results_traintest_config_for48
    
    % run_results_traintest_makepdffigs.m 
    %    loads the matfiles with saved figures - cleans up and saves pdfs to disk
    %    calls:
    %       hack_subplot_reorg
    %
    %
    run_results_traintest_makepdffigs
    
    
    %% Used in above
    % Stl_load \ stlload_trainedNetworkModel_ARRAY_fromfile.m

    %%
    % hacking_2021_07_14_errorfunctioncompareoftimeseries.m

    %% __trEACH
    % \runscripts_trainonall_testonall\
    % run_trainon_allpreforms_rev7_2lstm_withnewdataload_TEMP_auto.m
    run_trainon_allpreforms_rev7_2lstm_withnewdataload_TEMP_auto
end

%% Above - the approximate current set of active expirmentation for training and test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            