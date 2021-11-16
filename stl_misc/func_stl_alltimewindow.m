function [timewindow] = func_stl_alltimewindow(BatchInfo,nWhichBatch)

%BatchInfo(1).subbatchIndsSTORE{1}

numBatches = length(BatchInfo);

ii          = nWhichBatch; %31;  %which Batch (ie preform)_

% Get the length of batch (of the indicated
lengthbatch = length(BatchInfo(ii).origiinds);
% Get the number of subbatches (number of long reqions within UL limits
numsubbatch = length(BatchInfo(ii).subbatchIndsSTORE);

%create a timevector same length as preform/batch
timewindow = zeros(lengthbatch,1);
% for the location where subbatches are within UL limits 
% let time window be 1
for jj = 1:numsubbatch
    iindssubbatch = BatchInfo(ii).subbatchIndsSTORE{jj};
    timewindow(iindssubbatch) = 1;
end

%smooth the edges of the time window
flt = ones(5,1);
flt = flt / length(flt);

%timewindow = conv(timewindow,flt,'same');

if(0)
    XX = BatchInfo(ii).data(:,6);
    plot(XX.*timewindow,'.')
end   
%    origiinds(BatchInfo(1).origiinds(BatchInfo(1).subbatchIndsSTORE{2});