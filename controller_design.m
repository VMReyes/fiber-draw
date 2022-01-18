clc; clear; close all;

strDataPath         = 'C:\Users\georg\Dropbox (MIT)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
strOutputPath       = 'C:\Users\georg\Desktop\MS Research\fiber-draw\alldatatrainfromonefile\';

filenamebasearray{1} = 'DrawData_Tower48_2020-12-01_to2020-12-08';  

%% File to struct using Brian's code

% PrefltLEN_array = [1 5 7]; 
PrefltLEN_array = [7];

for fn = 1:length(filenamebasearray)
    filenamebase = filenamebasearray{fn};
    for PrefltLEN = PrefltLEN_array
        PrefltLEN

        bXLSLoad                                                = 1;
        bPlotAll                                                = 0;
        bPlot_each_preform_on_subplot                           = 0;
        bPlot_each_preform_on_subplot_with_inrangesubbatches_   = 0;

        loBFD           = 124;
        hiBFD           = 126;

        nnin_LearnRateDropPeriod    = 200;
        nnin_maxEpochs              = 250 ; %150; %500; %1000
        nnin_miniBatchSize          = 200;%200;


        % a batch is the same as a preform, multiple batches (or preforms) are run,
        % one after the other in in the tower.
        % a subbatch is defined as a contiguous region of production
        subbatchMinLen 	= 2000; 
        subbatchMaxLen  = 8000;
        xls_file        = [filenamebase '.csv'];
        nTower          =  48;

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

%         if(~bLocked)
        if (1)
            bSubBatchCounter = 0;

            for nWhichBatch = 1:length(BatchInfo)
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


                numsubbatch = length(BatchInfo(nWhichBatch).subbatchIndsSTORE);
                if(numsubbatch>0)
                    bSubBatchCounter = 0;
                    for iinnWhichSub = 1:numsubbatch
                        x_cols = [4,1]; y_cols = [5,3]; % according to STRDEF
                        
                        % following adapted from
                        % stl_prep_trainingdata_allt_function_Wbfdref_2lstm_Wprefilt.m
                        [Y_sub,Y_filt_sub,x_sub, sample_indexfilt_sub]  = controller_design_parse(BatchInfo, ...
                                                                            nWhichBatch, ...
                                                                            iinnWhichSub,...
                                                                            nImportantColumns,...
                                                                            bPlot, ...
                                                                            bMeanRemove, ...
                                                                            PrefltLEN, ...
                                                                            x_cols, ...
                                                                            y_cols);
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
 
        end
    end
end

%% From X Y ONEfile, Kt Controller
for i = 1:1 %length(YTrainTRANSPOSE_fromONEfile)
    xs = cell2mat(XTrainTRANSPOSE_fromONEfile(i));
    ys = cell2mat(YTrainTRANSPOSE_fromONEfile(i));
    tension = xs(1,:); bfd = xs(2,:);
    furnace_power = ys(1,:); capstan_speed = ys(2,:);

    % experiment with [1/2 2 1/2 1]
    
    sys_kt = armax(iddata(furnace_power', tension', 1), [1 2 2 1]);
    sys_kd = armax(iddata(capstan_speed', bfd',     1), [1 2 2 1]);

    % back out PID gains

    
end
disp('Done!')