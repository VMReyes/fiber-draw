function [BatchInfo, STRDEF, dataa, data00, dataOTHER, uniPreformID  ]  = ...
    stl_load_parse_allt_function_rev4(...
        bXLSLoad, strDataPath, xls_file, ...
        nTower, ...
        bPlotAll, bPlot_each_preform_on_subplot, bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
        loBFD, hiBFD , subbatchMinLen, subbatchMaxLen)

if(0)
    bXLSLoad = 0;

    bPlotAll                                                = 0;
    bPlot_each_preform_on_subplot                           = 1;
    bPlot_each_preform_on_subplot_with_inrangesubbatches_   = 1;

    loBFD = 124.8;
    hiBFD = 125.2; %this HERE
    % loBFD = 124.7;
    % hiBFD = 125.3; %this HERE
    % loBFD = 124;
    % hiBFD = 126;
    % 
    % loBFD = 124.8;
    % hiBFD = 125.2;

    % loBFD = 115;
    % hiBFD = 135;
    %loBFD = 50;
    %hiBFD = 500;

    %subbatchMinLen = 500; %100 %HERE

    subbatchMinLen = 200; %100 %HERE

%     [BatchInfo, STRDEF]  = stl_load_parse_allt_function(bXLSLoad, bPlotAll, bPlot_each_preform_on_subplot, ...
%         bPlot_each_preform_on_subplot_with_inrangesubbatches_, loBFD, hiBFD , subbatchMinLen);

end


%main output is BatchInfo

%%
STRDEF{1} = 'barefibrediadisplay';
STRDEF{2} = 'rampupslopevalflt';
STRDEF{3} = 'cpspdactval';
STRDEF{4} = 'tenncmv';
STRDEF{5} = 'frnpwrmv';
STRDEF{6} = 'pfspdactval';
STRDEF{7} = 'hetubetemp';

%%
%check to see of mat files exist
% bMatFileExists = 0;
% iindperiod = strfind(xls_file,'.');
% if(~isempty(iindperiod))
%     matfilename = [xls_file(1:(iindperiod-1)) '.mat'];
%     nexist = exist(matfilename);
%     if(nexist == 2)
%         bMatFileExists = 1
%     end
% end

%%

if(1)
    %%

    % strDataPath     = 'E:\Dropbox (SquareCircleMITtoo)\minigroup_mit_sterlite\from Sterlite\data\';
    % xls_file        = 'Draw_Tower_48.xlsx';

%     strDataPath     = 'E:\Dropbox (SquareCircleMITtoo)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
%     %xls_file        = 'DrawData_Tower48_2020-12-01_to2020-12-08.csv';
%     xls_file        = 'DrawData_Tower48_2021-03-16_to2021-03-23.csv';


    strfullpth      = strDataPath + filesep + xls_file;
    headers         = readcell(strfullpth,'Range','1:1');

    %which columns
    nImportantColumns = [ find(strcmpi(headers,'barefibrediadisplay')==1)
                            find(strcmpi(headers,'rampupslopevalflt')==1)
                            find(strcmpi(headers,'cpspdactval')==1)
                            find(strcmpi(headers,'tenncmv')==1)
                            find(strcmpi(headers,'frnpwrmv')==1)
                            find(strcmpi(headers,'pfspdactval')==1)
                            find(strcmpi(headers,'hetubetemp')==1)]';

%     nImportantColumns = [ find(strcmpi(headers,'barefibrediadisplay')==1)
%                             find(strcmpi(headers,'tenncmv')==1) % TODO(gcfchen): why was this commented
%                             find(strcmpi(headers,'rampupslopevalflt')==1)
%                             find(strcmpi(headers,'cpspdactval')==1)
%                             find(strcmpi(headers,'frnpwrmv')==1) % TODO(gcfchen): why was this commented
%                             find(strcmpi(headers,'pfspdactval')==1)
%                             find(strcmpi(headers,'hetubetemp')==1)]';

    disp('Using following important columns.')
    xlsColNum2Str(nImportantColumns)

    figure(100)
    data00 =[];
    for ccc = 1:length(nImportantColumns)

        ColLetterCell = xlsColNum2Str(nImportantColumns(ccc));
        ColLetter       = ColLetterCell{1};

        datcol0    = readmatrix(strfullpth,'Range',[ColLetter ':' ColLetter]);

        %trim off the first row which has the headers
        datcol0 = datcol0(2:end);

        subplot(length(nImportantColumns),1,ccc)
        plot(datcol0)

        data00 = [data00 datcol0];
    end

    %good fiber start
        nGoodFiberCol = [ find(strcmpi(headers,'GoodFibreStartState')==1) ];
        xlsColNum2Str(nGoodFiberCol)
        ColLetterCell = xlsColNum2Str(nGoodFiberCol);
        ColLetter = ColLetterCell{1};
        datcolGOOD    = readmatrix(strfullpth,'Range',[ColLetter ':' ColLetter]);
        datcolGOOD    = datcolGOOD(2:end);

    % Preform Process Position
        nCol = [ find(strcmpi(headers,'pfprocesspsn')==1) ];
        xlsColNum2Str(nCol)
        ColLetterCell = xlsColNum2Str(nCol);
        ColLetter = ColLetterCell{1};
        datcolPreformPos    = readmatrix(strfullpth,'Range',[ColLetter ':' ColLetter]);
        datcolPreformPos    = datcolPreformPos(2:end);

    % Fiber Draw length / position
        nCol = [ find(strcmpi(headers,'CpLenTareVal')==1) ];
        xlsColNum2Str(nCol)
        ColLetterCell = xlsColNum2Str(nCol);
        ColLetter = ColLetterCell{1};
        datcolDrawLength    = readmatrix(strfullpth,'Range',[ColLetter ':' ColLetter]);
        datcolDrawLength    = datcolDrawLength(2:end);

    % %this is slow
    % datcolDATE    = readcell(strfullpth,'Range','AE:AE');
    % datcolTIME    = readcell(strfullpth,'Range','AF:AF');
    %this is much faster
    % TT = readtimetable(strfullpth,'RowTimes','date');
    % T = readtable(strfullpth);

    %date
        nCol = [ find(strcmpi(headers,'date')==1) ];
        xlsColNum2Str(nCol)
        ColLetterCell = xlsColNum2Str(nCol);
        ColLetter = ColLetterCell{1};
        Tdate = readtable(strfullpth, 'Range',[ColLetter ':' ColLetter]);
        %Tdate = readtable(strfullpth, 'Range','AE:AE');

    %time
        nCol = [ find(strcmpi(headers,'time')==1) ];
        xlsColNum2Str(nCol)
        ColLetterCell = xlsColNum2Str(nCol);
        ColLetter = ColLetterCell{1};
        Ttime = readtable(strfullpth, 'Range',[ColLetter ':' ColLetter]);
        %Ttime = readtable(strfullpth, 'Range','AF:AF');
        % Tdate.Variables
        % Ttime.Variables


    datcolDATEasnumber = datenum(table2array(Tdate(:,1)));
    datcolTIMEasnumber = datenum(table2array(Ttime(:,1)));

    %ms time
        nCol = [ find(strcmpi(headers,'Millitm')==1) ];
        xlsColNum2Str(nCol)
        ColLetterCell = xlsColNum2Str(nCol);
        ColLetter = ColLetterCell{1};
        datcolMSTIMEasnumber = readmatrix(strfullpth,'Range',[ColLetter ':' ColLetter]);
        %datcolMSTIMEasnumber = readmatrix(strfullpth,'Range','AG:AG');
        datcolMSTIMEasnumber = datcolMSTIMEasnumber(2:end);


    %Spool ID
        nCol = [ find(strcmpi(headers,'Spool_ID')==1) ];
        xlsColNum2Str(nCol)
        ColLetterCell = xlsColNum2Str(nCol);
        ColLetter = ColLetterCell{1};
        Tspoolid = readtable(strfullpth, 'Range',[ColLetter ':' ColLetter]);

    temp = table2array(Tspoolid(:,1));
    strLenID = strlength(temp);
    indswithValidID = find(strLenID == 11);
    temp_valid = temp(indswithValidID);

    cellarrayPreform_valid    =  extractBetween(temp_valid,1,8);
    cellarraySpoolNum_valid   =  extractAfter(temp_valid,8);
    datacolSpoolNum_valid = str2double(cellarraySpoolNum_valid);

    cellarrayPreform = cell(size(temp));
    cellarrayPreform(:) = {'foo'};
    cellarrayPreform(indswithValidID) = cellarrayPreform_valid;

    cellarraySpoolNum = cell(size(temp));
    cellarraySpoolNum(indswithValidID) = cellarraySpoolNum_valid;

    %numericalID for spool
    datacolSpoolNum = -ones(size(temp));
    datacolSpoolNum(indswithValidID) = datacolSpoolNum_valid;

    %average time of each preform
    uniPreformID = unique(cellarrayPreform_valid)
    for uu = 1:length(uniPreformID)
        idflag         = strcmp(cellarrayPreform,uniPreformID{uu});
        timepreform(uu) = mean(datcolDATEasnumber(find(idflag)) + datcolTIMEasnumber(find(idflag)))
    end
    %sort uni ids by time
    [timepreform,iinds] = sort(timepreform);
    uniPreformID = uniPreformID(iinds)
    %
    %
    %(index) numericalID for preform
    datacolPreFormNum = -ones(size(temp));
    %uniPreformID = unique(cellarrayPreform_valid);
    for uu = 1:length(uniPreformID)
        idflag         = strcmp(cellarrayPreform,uniPreformID{uu});
        datacolPreFormNum(find(idflag)) = uu;
    end

    %% SORT

    dataOTHER = [datcolDATEasnumber datcolTIMEasnumber datcolMSTIMEasnumber   datacolPreFormNum datacolSpoolNum datcolGOOD datcolPreformPos datcolDrawLength ];

    %this isn't needed if excel file is sorted
    %sort by time
    [dataOTHER, iinds] = sortrows(dataOTHER,3,'ascend');
    data00    = data00(iinds,:);
    [dataOTHER, iinds] = sortrows(dataOTHER,2,'ascend');
    data00    = data00(iinds,:);
    [dataOTHER, iinds] = sortrows(dataOTHER,1,'ascend');
    data00    = data00(iinds,:);

else
   matfilename
   load(matfilename);
   
end
%%

% figure(200)
% [m,n] = size(data00)
% for ccc = 1:n
%     
%     subplot(n,1,ccc)
%     plot(data00(:,ccc))
% 
% end
% 
% figure(300)
% [m,n] = size(dataOTHER)
% for ccc = 1:n
%     
%     subplot(n,1,ccc)
%     plot(dataOTHER(:,ccc))
% 
% end

%%
%stick data together
dataa = [data00 dataOTHER];
%getindex range for each batch
datacolPreFormNumSORTED = dataOTHER(:,4);

% --- see for format stl_load_parse_allt_function_rev3 or stl_load_parse_allt_function_rev2

numBatches = length(uniPreformID);
m6or8 = floor(sqrt(numBatches));
n5or4 = ceil(numBatches/m6or8);


clear BatchInfo;
for uu = 1:numBatches
    inds         = find(datacolPreFormNumSORTED == uu);
    [min(inds) max(inds)] 
    %
    origiinds = min(inds):max(inds);
    data4pf = dataa(origiinds,:);
    %
    BatchInfo(uu).data = data4pf;
    BatchInfo(uu).origiinds = origiinds;
end


%% plot each batches (preforms) on a subplot
if(bPlotAll || bPlot_each_preform_on_subplot)
    figure(10)
    for ii = 1:numBatches
        subplot(m6or8,n5or4,ii)
        %plot(dataSTORE{ii}(:,6))
        plot(BatchInfo(ii).data(:,1))
        title(['Batch ' num2str(ii)]);
    end
end

%%  Plot each batches (preforms) on a subplot
%   show within reange of loBFD and hiBFD
%loBFD = ...;
%hiBFD = ... 
% plot batches
if(bPlotAll || bPlot_each_preform_on_subplot_with_inrangesubbatches_)
    figure(11)
    for ii = 1:numBatches
        %BFD = dataSTORE{ii}(:,6);
        BFD = BatchInfo(ii).data(:,1);
        %
        %find all indices with BFD in range lo hi for the entire batch
        iinds_inspecA = find((BFD >= loBFD));
        iinds_inspecB = find((BFD <= hiBFD));
        iinds_inspec = intersect(iinds_inspecA, iinds_inspecB);

        subplot(m6or8,n5or4,ii)
        %plot(dataSTORE{ii}(:,6))
        plot(BatchInfo(ii).data(:,1))
        hold on
        %plot(iinds_inspec, dataSTORE{ii}(iinds_inspec,6),'r.')
        plot(iinds_inspec, BatchInfo(ii).data(iinds_inspec,1),'r.')
        hold off
        title(['Batch ' num2str(ii) ' w/ in-bound']);
    end
end

%% parse for subbatch regions of MinLen that are within low BFD and high BFD bounds
%  find long time to find 'good' SubBatches.
%  store indices within a Batch for all SubBatches in __ BatchInfo __.
%
% subbatchMinLen = ...
for ii = 1:numBatches
%for ii = 13
    ['batch ' num2str(ii)]

    %BFD = dataSTORE{ii}(:,6);
    BFD = BatchInfo(ii).data(:,1);
    %
    %find all indices with BFD in range lo hi for the entire batch
    iinds_inspecA = find((BFD >= loBFD));
    iinds_inspecB = find((BFD <= hiBFD));
    iinds_inspec = intersect(iinds_inspecA, iinds_inspecB);
        
    % take diff on indices
    % non continuoous region will have indice diff greater than 1
    inddiff = diff(iinds_inspec);
    
    %[length(iinds_inspec) length(BFD)]
    
    indsOut = find(inddiff ~= 1);
    
    %catch case when 1 subbatch for entire or 1 long a beginning
    %(DO) this can be improved a little as it cuts off a valid pt at the end
    if(isempty(indsOut))
        indsOut = [indsOut iinds_inspec(end)-1];
    end
    
    subbatchIndsSTORE = [];
    subbatchcount=0;
    if(length(indsOut) >= 1)
        tempInds = iinds_inspec(1):iinds_inspec(indsOut(1));
        %  HACK >> fixed with subbatchMaxLen
        if(length(tempInds)>subbatchMaxLen)
            tempInds = tempInds(1:subbatchMaxLen);
        end
        %  << HACK
        if(length(tempInds)>subbatchMinLen)
            subbatchcount = subbatchcount+1;
            subbatchIndsSTORE{subbatchcount} = tempInds;
        end
        %
        for jj = 1:(length(indsOut)-1)
            tempInds = iinds_inspec(indsOut(jj)+1):iinds_inspec(indsOut(jj+1));
            %  HACK >>
            if(length(tempInds)>subbatchMaxLen)
                tempInds = tempInds(1:subbatchMaxLen);
            end
            %  << HACK
            if(length(tempInds)>subbatchMinLen)
                subbatchcount = subbatchcount+1;
                subbatchIndsSTORE{subbatchcount} = tempInds;
            end
        end
        %
        jj = length(indsOut);
        tempInds = iinds_inspec(indsOut(jj)+1):iinds_inspec(end);
        %  HACK >>
        if(length(tempInds)>subbatchMaxLen)
            tempInds = tempInds(1:subbatchMaxLen);
        end
        %  << HACK
        if(length(tempInds)>subbatchMinLen)
            subbatchcount = subbatchcount+1;
            subbatchIndsSTORE{subbatchcount} = tempInds;
        end
        
        %subbatchIndsSTORE

    end
    BatchInfo(ii).subbatchIndsSTORE = subbatchIndsSTORE;
    
    %subbatchIndsSTORE is a cell array with the indices, relative to 
    %the batch, of each subbatch.
    %Stored in BatchInfo, a structure indexed by batch number.
    
end


%% use dataSTORE and BatchInfo to plot each Batch with subBatch regions
% on seperate figures (or subplot) for each batch
if(bPlotAll || bPlot_each_preform_on_subplot_with_inrangesubbatches_)
    figure(12)
    BatchRange = [1:numBatches];
    %BatchRange = 25;

    for i = 1:1
        strlistclr{1} = 'r'; strlistclr{2} = 'g'; strlistclr{3} = 'b'; strlistclr{4} = 'm';
        strlistlegend = [];
        for ii = BatchRange
            %BFD = dataSTORE{ii}(:,6);
            BFD = BatchInfo(ii).data(:,1);
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

                tit             = ['Sub: ' num2str(jj) ' Len: ' num2str(length(tempInds))]; strlistlegend{1+jj} = tit;

            end
            hold off
            title(['Batch ' num2str(ii)])
            legend(strlistlegend, 'Location', 'Best');    
        end
    end
    sgtitle([xls_file '   ' num2str(subbatchMinLen) '   ' num2str(subbatchMaxLen)], 'Interpreter', 'none') 
end
