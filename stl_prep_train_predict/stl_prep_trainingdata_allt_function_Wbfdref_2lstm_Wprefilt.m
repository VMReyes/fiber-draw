function [Y,Y_filt,x,sample_indexfilt, meanY, meanX]  = stl_prep_trainingdata_allt_function_Wbfdref_2lstm_Wprefilt(BatchInfo, ...
                                                    STRDEF,...
                                                    nWhichBatch, ...
                                                    nWhichSub,...
                                                    nImportantColumns,...
                                                    fltLEN,...
                                                    bPlot, ...
                                                    bPlotAllSelectedColumns, ...
                                                    bMeanRemove, ... 
                                                    PrefltLEN, ...
                                                    x_columns, ...
                                                    y_columns, ...
                                                    yRemove125)


% time windor for all extended regions within UL bounds
[timewindow] = func_stl_alltimewindow(BatchInfo,nWhichBatch);

%for i = 1:1
    %total number of numsubbatch for the selected batch
    numsubbatch = length(BatchInfo(nWhichBatch).subbatchIndsSTORE);
    
    %get data array for selected batch from correct dataSTORE
    % (Batch)
        %dataaB    = dataSTORE{nWhichBatch};
    dataaB = BatchInfo(nWhichBatch).data;
    %select and re-order for ONLY the columns of importance. 
    % (Batch Inportant Reordered)
    dataaBIR = dataaB(:,nImportantColumns);
    
    %PrefltLEN = %1 does nothing, should be odd
    dataaBIRpadd = [dataaBIR(1,:)'*ones(1,floor(PrefltLEN/2))   dataaBIR'      dataaBIR(end,:)' * ones(1,floor(PrefltLEN/2))]';
    [dataaBIR,foooo] = func_simplefilter(dataaBIRpadd,PrefltLEN,'valid');
    
    %create data matrix for training
    dataaFIN = [];
    % select subbatch rows within continuous bounds
    if(nWhichSub >= 0)
        tempInds = BatchInfo(nWhichBatch).subbatchIndsSTORE{nWhichSub};
        dataaFIN = dataaBIR(tempInds,:);
    end
    %if nWhichSub < 0 then we want all of the data with gaps filled in with
    %the timewindow
    if(nWhichSub <  0)
        for cc = 1:length(nImportantColumns)
            %dataaFIN(:,cc) = dataaBIR(:,cc);
            dataaFIN(:,cc) = dataaBIR(:,cc).*timewindow;
        end
    end
    
    if(nWhichSub <  0)
        %here construct dataFIN from subbatches in order to have more control
        %over how deal with gaps between subbatches
        dataaFIN2 = [];
        for sb = 1:(numsubbatch-1)
            tempInds                = BatchInfo(nWhichBatch).subbatchIndsSTORE{sb};
            dataa4oneSUB            = dataaBIR(tempInds,:);
            dataaFIN2(tempInds,:)   = dataa4oneSUB;

            %fill in the gaps between....
            tempIndsPLUS    = BatchInfo(nWhichBatch).subbatchIndsSTORE{sb+1};
            tempIndsGAP     = (tempInds(end)+1):(tempIndsPLUS(1)-1);

            %here hold mean
            %INTERPvals = mean(dataa4oneSUB);
            %INTERPvals = mean(dataa4oneSUB(  floor(length(tempInds)/10*8):end ,:));
            INTERPvals = mean(dataa4oneSUB(  floor(length(tempInds)/10*8):floor(length(tempInds)/10*9) ,:));
            %INTERPvals = mean(dataa4oneSUB(  (end-100):end ,:));
            for cc = 1:length(nImportantColumns)
                %dataaFIN(:,cc) = dataaBIR(:,cc);
                dataaFIN2(tempIndsGAP,cc) = ones(length(tempIndsGAP),1)*INTERPvals(cc);
            end
        end
        %lat subbatch
        sb = numsubbatch
            tempInds                = BatchInfo(nWhichBatch).subbatchIndsSTORE{sb};
            dataa4oneSUB            = dataaBIR(tempInds,:);
            dataaFIN2(tempInds,:)   = dataa4oneSUB;

        %here hold mean - to end    
        tempIndsGAP     = ( tempInds(end)+1):length(dataaB) ;
        INTERPvals = mean(dataa4oneSUB(  floor(length(tempInds)/10*8):floor(length(tempInds)/10*9) ,:));
        for cc = 1:length(nImportantColumns)
            %dataaFIN(:,cc) = dataaBIR(:,cc);
            dataaFIN2(tempIndsGAP,cc) = ones(length(tempIndsGAP),1)*INTERPvals(cc);
        end
        dataaFIN = dataaFIN2;
    end
    
    for ii = 1:0
        %here construct dataFIN from subbatchesd in order to have more control
        %over how deal with gaps between 
        dataaFIN2 = [];
        for sb = 1:numsubbatch
            tempInds = BatchInfo(nWhichBatch).subbatchIndsSTORE{sb};
            dataa4oneSUB = dataaBIR(tempInds,:);
            dataaFIN2(tempInds,:) = dataa4oneSUB;
        end

        %assuming that smoothing is turned off on the time window, 
        %this should be zero:
        if(0)
            plot(dataaFIN-dataaFIN2)
        end
    end
    

    %size of data matrix
    [mRR,nRR] = size(dataaFIN);

    %plot all selected columns for selected batch and subbatch
    if(bPlot || bPlotAllSelectedColumns)
        figure
        for ii = 1:nRR
            subplot(nRR,1,ii)
            plot(dataaFIN(:,ii))
            tit = STRDEF{nImportantColumns(ii)};
            title(tit)
        end
    end

    %% OUTPUT - first column is the output (i.e. BFD)
    Y = dataaFIN(:,1:length(y_columns));
    %low pass filter on the noisey data
    %Ypadd   = [Y(1)*ones(1,floor(fltLEN/2)) Y' Y(end)*ones(1,floor(fltLEN/2))]';
    Ypadd = 0;
    %[Y_filt,sample_indexfilt] = func_simplefilter(Ypadd,fltLEN,'valid');
    Y_filt = 0;
    sample_indexfilt = 0;

    if(bPlot)
        figure
        plot(Y)
        hold on 
        plot(sample_indexfilt,Y_filt,'g');
        hold off
        title('BFD and LPF of BFD')
    end
    
    % mean remove on Y and Y_filt
    meanY       = mean(Y);
    meanY_filt  = mean(Y_filt);
    %
    if(bMeanRemove)
        Y           = Y         - meanY;
        Y_filt      = Y_filt    - meanY;
    end

    if(yRemove125)
        Y = Y - 125;
    end
    %% INPUTS
    x = dataaFIN(:,[length(y_columns)+1:end]);
    % mean remove on all inputs
    meanX = mean(x);
    if(bMeanRemove)
        x = x - meanX;
    end
    
    %% remove begin and end time - to deal with transients
    n2trim      = fltLEN;
%     x       = x(n2trim:(end-n2trim),:);
%     Y       = Y(n2trim:(end-n2trim),:);
%     Y_filt  = Y_filt(n2trim:(end-n2trim),:);
    sample_indexfilt = sample_indexfilt(n2trim:(end-n2trim));

    %% Create Data Structure - for training  Y and X
%     dsiddata = iddata(Y,x,0.005);
%     dsiddata_filt = iddata(Y_filt,x,0.005);
%end
