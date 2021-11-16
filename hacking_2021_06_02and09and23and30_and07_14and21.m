


% first 
% run run_trainon_allpreforms_re2qv5_2lstm.m


% layers 100 50
% run run_trainon_allpreforms_rev5_2lstm.m
% loBFD             = 124;
% hiBFD             = 126; 
% nTower            = 51
% nWhichBatch       =  1;
% subbatchMinLen 	= 2000;
% save rev5_temp10_all51_maxEpochs150
% save rev5_temp10_all51_maxEpochs80
%
% now with stl_prep_trainingdata_allt_function_Wbfdref_2lstm_Wprefilt
% filt with PrefltLEN= 5
% save rev5_temp10_all51_maxEpochs150_filter5
% save rev5_temp10_all48_maxEpochs150_filter5
%
% see dev_figs2_2018.m
%
%

%%% for 6/23/2021
%see run_trainon_allpreforms_rev6_2lstm_withnewdataload.m
% new data
% filt with PrefltLEN= 1
% subbatchMinLen 	= 2000;
% hack in load with top at 3000
%save rev6_temp10_newtest48_maxEpochs150_mx3000_filter1
%save rev6_temp10_newtest48_maxEpochs150_mx3000_filter5
%save rev6_temp10_newtest48_maxEpochs150_mx8000_filter5
%save rev6_temp10_newtest48_maxEpochs150_mx8000_filter1

%for 6/28
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_filter1.mat
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm100lstm0_filter1
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm150lstm0_filter1
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm200lstm0_filter1.mat
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm250lstm0_filter1
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm100lstm50lstm25_filter1

%for 6/29
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm200lstm100lstm50_filter1
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm350lstm0_filter1.mat
%
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm150lstm0_filter5.mat
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm200lstm0_filter5.mat
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm250lstm0_filter5.mat
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm350lstm0_filter5.mat
%
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm100lstm50lstm25_filter5
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm200lstm100lstm50_filter5
%
%load rev6_temp10_newtest48_maxEpochs250_mx3000_drop200_lstm500lstm0_filter5.mat

%for 6/30
%load rev6_temp10_newtest48_maxEpochs250_mx8000_drop200_lstm350lstm0_filter5
%
%load rev6_temp10_newtest48_maxEpochs250_mx8000_drop200_lstm100lstm50lstm25_filter7
%load rev6_temp10_newtest48_maxEpochs250_mx8000_drop200_lstm100lstm50lstm25_filter5

%
%

%load rev6_temp10_newtest48_maxEpochs250_mx8000_drop200_lstm350lstm0_filter1.mat
%load rev6_temp10_newtest48_maxEpochs250_mx8000_drop200_lstm350lstm0_filter5.mat
%load rev6_temp10_newtest48_maxEpochs250_mx8000_drop200_lstm350lstm0_filter7.mat
% 
%load rev7_newtest51_DrawData_Tower51_2020-12-22_to2020-12-29__filter1.mat
%load rev7_newtest51_DrawData_Tower51_2020-12-22_to2020-12-29__filter5.mat
% 
%load rev7_newtest51_DrawData_Tower51_2021-01-05_to2021-01-12__filter1.mat
%load  rev7_newtest51_DrawData_Tower51_2021-01-05_to2021-01-12__filter5.mat



% load rev7_newtest51_DrawData_Tower51_2020-12-22_to2020-12-29__filter1.mat
% trainedNetworkModel_ARRAY = [];
% matfilename ='rev7_DrawData_Tower51_2020-12-22_to2020-12-29__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter1.mat'
% [trainedNetworkModel_ARRAY] = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
%
% load rev7_newtest51_DrawData_Tower51_2020-12-22_to2020-12-29__filter5.mat
% trainedNetworkModel_ARRAY = [];
% matfilename ='rev7_DrawData_Tower51_2020-12-22_to2020-12-29__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter5.mat'
% [trainedNetworkModel_ARRAY] = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
%
% load rev7_newtest51_DrawData_Tower51_2021-01-05_to2021-01-12__filter1.mat
% trainedNetworkModel_ARRAY = [];
% matfilename ='rev7_DrawData_Tower51_2021-01-05_to2021-01-12__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter1.mat'
% [trainedNetworkModel_ARRAY] = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
%
% load rev7_newtest51_DrawData_Tower51_2021-01-05_to2021-01-12__filter5.mat
% trainedNetworkModel_ARRAY = [];
% matfilename = 'rev7_DrawData_Tower51_2021-01-05_to2021-01-12__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter5.mat'
% [trainedNetworkModel_ARRAY] = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);


% %load rev7_newtest51_DrawData_Tower51_2020-12-22_to2020-12-29__filter1.mat
% load rev7_newtest51_DrawData_Tower51_2021-01-05_to2021-01-12__filter1.mat
% trainedNetworkModel_ARRAY = [];
% matfilename ='rev7_DrawData_Tower51_2020-12-01_to2020-12-08__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter1.mat';
% tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% trainedNetworkModel_ARRAY{1} = tmpP{1};
% matfilename ='rev7_DrawData_Tower51_2020-12-08_to2020-12-15__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter1.mat';
% tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% trainedNetworkModel_ARRAY{2} = tmpP{1};
% matfilename ='rev7_DrawData_Tower51_2020-12-15_to2020-12-22__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter1.mat';
% tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% trainedNetworkModel_ARRAY{3} = tmpP{1};
% matfilename ='rev7_DrawData_Tower51_2020-12-22_to2020-12-29__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter1.mat';
% tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% trainedNetworkModel_ARRAY{4} = tmpP{1};
% matfilename ='rev7_DrawData_Tower51_2021-01-05_to2021-01-12__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter1.mat';
% tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% trainedNetworkModel_ARRAY{5} = tmpP{1};
%
%bPlotAllFreqOverride = 1;



% load rev7_DrawData_Tower48_2020-12-08_to2020-12-15__trEACH_maxE250_mx8000_drop200_lstm350lstm0_filter1_7.mat
% bPlotAllFreqOverride = 0;

% load rev7_DrawData_Tower48_2020-12-08_to2020-12-15__trEACH_maxE250_mx8000_drop200_lstm350lstm0_filter1_7.mat
% trainedNetworkModel_ARRAY = [];
% matfilename ='rev7_DrawData_Tower48_2020-12-08_to2020-12-15__trALL_maxE250_mx8000_drop200_lstm100lstm50_filter1.mat'
% [trainedNetworkModel_ARRAY] = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% bPlotAllFreqOverride = 1;

% load rev7_DrawData_Tower48_2020-12-08_to2020-12-15__trEACH_maxE250_mx8000_drop200_lstm350lstm0_filter1_7.mat
% trainedNetworkModel_ARRAY = [];
% matfilename ='rev7_DrawData_Tower48_2021-03-16_to2021-03-23__trALL_maxE250_mx8000_drop200_lstm100lstm50_filter1.mat'
% tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% trainedNetworkModel_ARRAY{1} = tmpP{1};
% matfilename ='rev7_DrawData_Tower48_2020-12-08_to2020-12-15__trALL_maxE250_mx8000_drop200_lstm100lstm50_filter1.mat'
% tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% trainedNetworkModel_ARRAY{2} = tmpP{1};
% bPlotAllFreqOverride = 1;


% load rev6_temp10_newtest48_maxEpochs250_mx8000_drop200_lstm350lstm0_filter1.mat
% trainedNetworkModel_ARRAY = [];
% matfilename = 'rev7_DrawData_Tower48_2021-03-16_to2021-03-23__trALL_maxE250_mx8000_drop200_lstm100lstm50_filter1.mat'
% [trainedNetworkModel_ARRAY] = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
% bPlotAllFreqOverride = 1;

load rev6_temp10_newtest48_maxEpochs250_mx8000_drop200_lstm350lstm0_filter1.mat
trainedNetworkModel_ARRAY = [];
matfilename ='rev7_DrawData_Tower48_2021-03-16_to2021-03-23__trALL_maxE250_mx8000_drop200_lstm100lstm50_filter1.mat'
tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
trainedNetworkModel_ARRAY{1} = tmpP{1};
matfilename ='rev7_DrawData_Tower48_2020-12-08_to2020-12-15__trALL_maxE250_mx8000_drop200_lstm100lstm50_filter1.mat'
tmpP = stlload_trainedNetworkModel_ARRAY_fromfile(matfilename);
trainedNetworkModel_ARRAY{2} = tmpP{1};
bPlotAllFreqOverride = 1;

%for7/21
%load rev7_DrawData_Tower51_2020-12-15_to2020-12-22__trALL_maxE150_mx8000_drop100_lstm100lstm50_filter1.mat


% on 9/18/2021
%moved to run_results_traintest

errorMtime = [];
errorMfreq = [];
        
for nWhichBatch2Test = 1:length(XTrainTRANSPOSE_ARRAY)
%for nWhichBatch2Test = length(trainedNetworkModel_ARRAY)

    Ypredict2Matrix = [];
    XTRANSPOSE_2Test =     XTrainTRANSPOSE_ARRAY{nWhichBatch2Test};
    YTRANSPOSE_2Test =     YTrainTRANSPOSE_ARRAY{nWhichBatch2Test};
        
    %for nWhichBatch2use4model = nWhichBatch2Test %just the diagonal
    for nWhichBatch2use4model = 1:length(trainedNetworkModel_ARRAY)  %
        
        trainedNetworkModel = trainedNetworkModel_ARRAY{nWhichBatch2use4model} ; %i.e. which batch to use for the training data that lead to model



        Y_4allsubbatchesinbatch = [];
        Ymodel_4allsubbatchesinbatch = [];
        MSEtimeArray_subs   = [];
        MSEfreqArray_subs   = [];
        MSElenThelp_subs    = [];
        MSElenFhelp_subs    = [];
        %
        nsubs = length(XTRANSPOSE_2Test);
        for nWhichSub = 1:nsubs
            x_sub = XTRANSPOSE_2Test{nWhichSub}';
            Y_sub = YTRANSPOSE_2Test{nWhichSub}';

            modelYOutResponse = predict(trainedNetworkModel,x_sub')';
            modelYOutResponse0 = modelYOutResponse;

            %from ___
            [FY, FY_filt, FYpredict, FYpredict_filt, frq_discretes, ...
             PY, PY_filt, PYpredict, PYpredict_filt, WW, ...
            ]           = stl_plottrainingresults_FFTPSD_function(Y_sub, Y_sub, modelYOutResponse, modelYOutResponse, -1, -1, 0, [], 'str');
            %]           = stl_plottrainingresults_FFTPSD_function(Y, Y, modelYOutResponse, modelYOutResponse, nWhichBatch, nWhichSub, nWhichSub*100000, [], 'str');


            Y_4allsubbatchesinbatch = [Y_4allsubbatchesinbatch' Y_sub']';
            Ymodel_4allsubbatchesinbatch = [Ymodel_4allsubbatchesinbatch' modelYOutResponse']';
            

            bMeanRemove4freqanalysis = 1;
            [MSEtime_OUT, MSEfreq_OUT, MSElenT_OUT, MSElenF_OUT] = func_errors_rev3_errsbetween2timeseries(Y_sub,modelYOutResponse, bMeanRemove4freqanalysis );
            MSEtimeArray_subs(nWhichSub)    = MSEtime_OUT;
            MSEfreqArray_subs(nWhichSub)    = MSEfreq_OUT;
            MSElenThelp_subs(nWhichSub)     = MSElenT_OUT;
            MSElenFhelp_subs(nWhichSub)     = MSElenF_OUT;
            

            
            if(1) %plotting

                for i2 = 1:1

                figure(nWhichBatch2Test*10000)
                subplot(4,nsubs,nWhichSub)
                title('response estimation')
                if(nWhichBatch2use4model == 1) %only first time though
                plot(Y_sub,'k')
                end
                hold on
                plot(modelYOutResponse,'r')
                grid on
                legend('System','Estimated')
                %title('Evaluation')
                title(['Eval. Batch ' num2str(nWhichBatch2Test) '.  Sub ' num2str(nWhichSub) '. ']);
                ylim([-.5 .5])
                xlim([1 length(Y_sub)-1])

                tinds = floor(length(Y_sub)/2):( floor(length(Y_sub)/2) + 200);
                subplot(4,nsubs,1*nsubs+nWhichSub)
                %title('response estimation')
                plot(tinds, Y_sub(tinds),'k')
                hold on
                plot(tinds, modelYOutResponse(tinds),'r')
                grid on
                legend('System','Estimated')
                title('Eval - time zoom')
                %ylim([-.5 .5])
                xlim([tinds(1) tinds(end-1)])


                if((nWhichBatch2use4model == nWhichBatch2Test) || bPlotAllFreqOverride)
                strTitlePart = 'str';
                %
                % Plot resutls of FFT
                FY_s = fftshift(FY);
                FYpredict_s = fftshift(FYpredict);
                strlistlegend = [];
                subplot(4,nsubs,2*nsubs+nWhichSub)
                plot(frq_discretes,log10(abs(FY_s).^2),'k');            strlistlegend{length(strlistlegend)+1} = 'FY_s';
                hold on
              %  stem(frq_discretes,log10(abs(FY_filt_s).^2),'g');       strlistlegend{length(strlistlegend)+1} = 'FY_filt_s';
                plot(frq_discretes,log10(abs(FYpredict_s).^2),'r');     strlistlegend{length(strlistlegend)+1} = 'FYpredict_s';
              %  stem(frq_discretes,log10(abs(FYpredict_filt_s).^2),'b');        strlistlegend{length(strlistlegend)+1} = 'FYpredict_filt_s';
                hold on %off
                legend(strlistlegend)
                ylabel('Log Magnitude')
                xlabel('Frequency')
                %axis([0 1 -.1 max(log10(abs(FY_s).^2))])
                axis([0 1 -.1  max([   max(log10(abs(FY_s).^2))   max(log10(abs(FYpredict_s).^2))   ])   ])
                % title('Power Spectrum:  FFT^2  vs Freq (-pi to pi)')
                title(['Power Spectrum. ' strTitlePart ': Batch ' num2str(nWhichBatch2Test) '.  Sub ' num2str(nWhichSub) '. ']);
                xlim([0 .6])

                %
                % Plot resutls of Welch
                strlistlegend = [];
                subplot(4,nsubs,3*nsubs+nWhichSub)
                plot(WW/pi,fmag1(PY),'k');            strlistlegend{length(strlistlegend)+1} = 'PY';
                hold on
              %  plot(WW,fmag1(PY_filt),'g');       strlistlegend{length(strlistlegend)+1} = 'PY_filt';
                plot(WW/pi,fmag1(PYpredict),'r');     strlistlegend{length(strlistlegend)+1} = 'PYpredict';
              %  plot(WW,fmag1(PYpredict_filt),'b');        strlistlegend{length(strlistlegend)+1} = 'PYpredict_filt';
                hold on %off
                legend(strlistlegend)
                ylabel('Normalized to Max Magnitude')
                xlabel('Frequency')
                axis([0 1 0 max(fmag1(PY))])
                %axis([-.25 .25 -.1 max(log10(abs(FY_s).^2))])
                % title('Power Spectrum:  FFT^2  vs Freq (-pi to pi)')
                title(['Welch Power Spectrum. ' strTitlePart ': Batch ' num2str(nWhichBatch2Test) '.  Sub ' num2str(nWhichSub) '. ']);
                xlim([0 .6])
                end
                
                end %end of for used just for collapsing code
            end

        end %end of subbatch loop
        
        Ypredict2Matrix =  [Ypredict2Matrix Ymodel_4allsubbatchesinbatch];

        MSEtimetotal = sum(MSEtimeArray_subs.*MSElenThelp_subs)/sum(MSElenThelp_subs);
        MSEfreqtotal = sum(MSEfreqArray_subs.*MSElenFhelp_subs)/sum(MSElenFhelp_subs);
        
        errorMtime(nWhichBatch2use4model,nWhichBatch2Test ) = MSEtimetotal;
        errorMfreq(nWhichBatch2use4model,nWhichBatch2Test ) = MSEfreqtotal;

    end %end of nWhichBatch2use4model
rsltSTORE{nWhichBatch2Test}.Ypredict2Matrix = Ypredict2Matrix;
rsltSTORE{nWhichBatch2Test}.Y_4allsubbatchesinbatch = Y_4allsubbatchesinbatch;

end %end of nWhichBatch2Test

%%
nWhichBatch2Test_array = 1:length(XTrainTRANSPOSE_ARRAY);
%nWhichBatch2Test_array = length(XTrainTRANSPOSE_ARRAY);

for nWhichBatch2Test        = nWhichBatch2Test_array
    Ypredict2Matrix         = rsltSTORE{nWhichBatch2Test}.Ypredict2Matrix ;
    Y_4allsubbatchesinbatch = rsltSTORE{nWhichBatch2Test}.Y_4allsubbatchesinbatch ;
    
    Y2 = Y_4allsubbatchesinbatch;

    nWhichSub = -1;

    figure(nWhichBatch2Test*1000000)
    subplot 211
    Y2mean = mean(Ypredict2Matrix,2);
    Y2std  = std(Ypredict2Matrix,0,2);
    %fill([1:length(Y2mean) fliplr(1:length(Y2mean))], [[Y2mean-3*Y2std]' fliplr([Y2mean+3*Y2std]')], 'g')
    hold on
    plot(Y2,'k');
    hold on 
    h = plot(Y2mean,'r');
    set(h, 'LineWidth', 2)
    plot(Y2mean+1*Y2std,'b-.');
    plot(Y2mean-1*Y2std,'b-.');
    % %axis([1 length(Y2) loBFD hiBFD]);
    % axis([1 length(Y2) 121 129]);
    % if(bMeanRemove)
    %     axis([1 length(Y2) loBFD-meanY2 hiBFD-meanY2]);
    % end         
    legend(['Test: Batch ' num2str(nWhichBatch2Test) ],'Prediction avg from training on all other batches')
    hold off
    title(['Concatinated. Test: Batch ' num2str(nWhichBatch2Test) '.  Sub ' num2str(nWhichSub) '. with ' num2str(inf) ]);
        
    subplot 212
    Y2mean = mean(Ypredict2Matrix,2);
    Y2std  = std(Ypredict2Matrix,0,2);
    patch([1:length(Y2mean) fliplr(1:length(Y2mean))], [[Y2mean-1*Y2std]' fliplr([Y2mean+1*Y2std]')], 'g')
    hold on
    %plot(Y2,'k');
    hold on 
    h = plot(Y2mean,'r');
    set(h, 'LineWidth', 2)
    %plot(Y2mean+1*Y2std,'b-.');
    %plot(Y2mean-1*Y2std,'b-.');
    % %axis([1 length(Y2) loBFD hiBFD]);
    % axis([1 length(Y2) 121 129]);
    % if(bMeanRemove)
    %     axis([1 length(Y2) loBFD-meanY2 hiBFD-meanY2]);
    % end         
    legend(['Test: Batch ' num2str(nWhichBatch2Test) ',  +/- 1 stdev'],'Prediction avg from training on all other batches')
    hold off
    title(['Concatinated. Test: Batch ' num2str(nWhichBatch2Test) '.  Sub ' num2str(nWhichSub) '. with ' num2str(inf) ]);
    
    figure(nWhichBatch2Test*1000001)
    h = plot(Y2,'k');
    set(h, 'LineWidth', 2)
    hold on
    plot(Ypredict2Matrix);
    % %axis([1 length(Y2) loBFD hiBFD]);
    % axis([1 length(Y2) 121 129]);
    % if(bMeanRemove)
    %     axis([1 length(Y2) loBFD-meanY2 hiBFD-meanY2]);
    % end
    legend(num2str([0 nWhichBatch2Test_array]'))
    hold off
    title(['Concatinated. Test: Batch ' num2str(nWhichBatch2Test) '.  Sub ' num2str(nWhichSub) '. with ' num2str(inf) ]);



end



% 
% 
% %https://www.mathworks.com/matlabcentral/answers/99917-how-can-i-label-each-pixel-of-the-output-of-imagesc-with-its-value
% 
% % Generate Random Data
% N = 5;
% M = rand(N);
% x = repmat(1:N,N,1); % generate x-coordinates
% y = x'; % generate y-coordinates
% % Generate Labels
% t = num2cell(M); % extact values into cells
% t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% % Draw Image and Label Pixels
% imagesc(M)
% text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
% 



%%


figure(777777700)
%see https://stackoverflow.com/questions/3942892/how-do-i-visualize-a-matrix-with-colors-and-values-displayed
mat = errorMtime           % A 5-by-5 matrix of random values from 0 to 1
imagesc(mat);            % Create a colored plot of the matrix values
% colormap(flipud(gray));  % Change the colormap to gray (so higher values are
%                          %   black and lower values are white)

textStrings = num2str(mat(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
%[x, y] = meshgrid(1:16);  % Create x and y coordinates for the strings
[szx,szy] = size(mat');
[x, y] = meshgrid(1:szx,1:szy);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(mat(:) > midValue, 1, 3);  % Choose white or black for the
                                               %   text color of the strings so
                                               %   they can be easily seen over
                                               %   the background color
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors

% set(gca, 'XTick', 1:5, ...                             % Change the axes tick marks
%          'XTickLabel', {'A', 'B', 'C', 'D', 'E'}, ...  %   and tick labels
%          'YTick', 1:5, ...
%          'YTickLabel', {'A', 'B', 'C', 'D', 'E'}, ...
%          'TickLength', [0 0]);
% 
title('MSE error matrix time')

figure(777777701)
%see https://stackoverflow.com/questions/3942892/how-do-i-visualize-a-matrix-with-colors-and-values-displayed
mat = errorMfreq           % A 5-by-5 matrix of random values from 0 to 1
imagesc(mat);            % Create a colored plot of the matrix values
% colormap(flipud(gray));  % Change the colormap to gray (so higher values are
%                          %   black and lower values are white)

textStrings = num2str(mat(:), '%0.2f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
%[x, y] = meshgrid(1:16);  % Create x and y coordinates for the strings
[szx,szy] = size(mat');
[x, y] = meshgrid(1:szx,1:szy);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(mat(:) > midValue, 1, 3);  % Choose white or black for the
                                               %   text color of the strings so
                                               %   they can be easily seen over
                                               %   the background color
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors

% set(gca, 'XTick', 1:5, ...                             % Change the axes tick marks
%          'XTickLabel', {'A', 'B', 'C', 'D', 'E'}, ...  %   and tick labels
%          'YTick', 1:5, ...
%          'YTickLabel', {'A', 'B', 'C', 'D', 'E'}, ...
%          'TickLength', [0 0]);
% 
title('MSE error matrix freq')

