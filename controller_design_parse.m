function [Y,Y_filt,x,sample_indexfilt]  = controller_design_parse(BatchInfo, ...
    nWhichBatch, ...
    nWhichSub,...
    nImportantColumns,...
    bPlot, ...
    bMeanRemove, ...
    PrefltLEN, ...
    x_cols, ...
    y_cols)

% adapted from stl_prep_trainingdata_allt_function_Wbfdref_2lstm_Wprefilt.m

% first attempt at auto regress fit in a good range of data -- lots to do here

% time window for all extended regions within UL bounds
[timewindow] = func_stl_alltimewindow(BatchInfo,nWhichBatch);


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

        INTERPvals = mean(dataa4oneSUB(  floor(length(tempInds)/10*8):floor(length(tempInds)/10*9) ,:));
        for cc = 1:length(nImportantColumns)
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

    %
    %
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

end

%size of data matrix
[mRR,nRR] = size(dataaFIN);

%% OUTPUT - first column is the output (i.e. BFD)
% Y = dataaFIN(2:length(dataaFIN), y_cols);
Y = dataaFIN(1:length(dataaFIN), y_cols);


Y_filt = Y;

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
%% INPUTS
% x = dataaFIN(1:length(dataaFIN)-1,x_cols);
x = dataaFIN(1:length(dataaFIN),x_cols);

% mean remove on all inputs
meanX = mean(x);
if(bMeanRemove)
    x = x - meanX;
end
sample_indexfilt = 0;
%% Create Data Structure - for training  Y and X
%     dsiddata = iddata(Y,x,0.005);
%     dsiddata_filt = iddata(Y_filt,x,0.005);
%end
