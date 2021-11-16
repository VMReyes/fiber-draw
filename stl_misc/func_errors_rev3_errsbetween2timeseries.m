function [MSEtime_OUT, MSEfreq_OUT, MSElenT_OUT, MSElenF_OUT] = func_errors_rev3_errsbetween2timeseries(Y2,Ypredict2, ...
    bMeanRemove4freqanalysis )
% 

    %MeanRemove for the frequency analysis
    if(bMeanRemove4freqanalysis)
        [Y002, meanYa2]                = func_meanremove(Y2);
        [Ypredict002, meanYc2]         = func_meanremove(Ypredict2);
    else
        Y002                = Y2;
        Ypredict002         = Ypredict2;
    end


    [FY2, FY_filt2, FYpredict2, FYpredict_filt2, frq_discretes2, ...
     PY2, PY_filt2, PYpredict2, PYpredict_filt2, WW2, ...
    ]           = stl_FFTPSD_function(Y002, Y002, Ypredict002, Ypredict002);
                                                                              
    %single error measures
    MSEtime2    = (sum((Y2 - Ypredict2).^2))/length(Y2);
    %calculate error for low frequencies
    iindLOWFre  = find(WW2 <= 1);
    MSEfreq2    = (sum((PY2(iindLOWFre) - PYpredict2(iindLOWFre)).^2))/length(iindLOWFre);

    MSEtime_OUT = MSEtime2;
    MSEfreq_OUT = MSEfreq2;
    MSElenT_OUT = length(Y2);
    MSElenF_OUT = length(iindLOWFre);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [MSEtimeArray2, MSEfreqArray2, MSElenThelp2, MSElenFhelp2] = func_errors_rev3(BatchInfo,nWhichBatch2Test,nWhichSub2Test, Y2,Y_filt2, Ypredict2, Ypredict_filt2, bMeanRemove4freqanalysis )
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% here get a particular (all) subbatch to plot on the Test
% 
% numsubbatch2Test = length(BatchInfo(nWhichBatch2Test).subbatchIndsSTORE);
%     %if nWhichSub2Test is a pos number 
%     %then we are evaluating on just that one subbatch
%     if(nWhichSub2Test >= 0)
%         numsubbatch2Test = 1;
%     end
% 
% MSEtimeArray2 =[];
% MSEfreqArray2 =[];
% MSElenThelp2 = [];
% MSElenFhelp2 = [];
%     
% cnt = 0;
% for nWhichSub_2PLOT = 1:numsubbatch2Test
%     %if nWhichSub2Test is a pos number 
%     %then we are evaluating on just that one subbatch
%     if(nWhichSub2Test >= 0)
%         nWhichSub_2PLOT = nWhichSub2Test;
%     end
%     
%     cnt = cnt + 1;
%     
%     Y_2PLOT2             = Y2;
%     Y_filt_2PLOT2        = Y_filt2;
%     Ypredict_2PLOT2      = Ypredict2;
%     Ypredict_filt_2PLOT2 = Ypredict_filt2;
%     if(1)
%         %nWhichSub_2PLOT 	= 5;
%         tempInds_2PLOT      = BatchInfo(nWhichBatch2Test).subbatchIndsSTORE{nWhichSub_2PLOT};
%         %if nWhichSub2Test is a pos number 
%         %then we are evaluating on just that one subbatch
%         if(nWhichSub2Test >= 0)
%             tempInds_2PLOT = 1:length(Y2);
%         end
%         
%         %
%         Y_2PLOT2             = Y2(tempInds_2PLOT);
%         Y_filt_2PLOT2        = Y_filt2(tempInds_2PLOT);
%         Ypredict_2PLOT2      = Ypredict2(tempInds_2PLOT);
%         Ypredict_filt_2PLOT2 = Ypredict_filt2(tempInds_2PLOT);    
% 
%     end
% 
%     %MeanRemove for the frequency analysis
%     if(bMeanRemove4freqanalysis)
%         Y002                = Y2;
%         Y_filt002           = Y_filt2;
%         Ypredict002         = Ypredict2;
%         Ypredict_filt002    = Ypredict_filt2;
%     else
%         [Y002, meanYa2]                = func_meanremove(Y2);
%         [Y_filt002, meanYb2]           = func_meanremove(Y_filt2);
%         [Ypredict002, meanYc2]         = func_meanremove(Ypredict2);
%         [Ypredict_filt002, meanYd2]    = func_meanremove(Ypredict_filt2);
%     end
%     % and subbatch for 2PLOT ( if analyzing all and get a particular subbatch
%     % to plot)
%         Y00_2FREQ2             = Y002(tempInds_2PLOT);
%         Y_filt00_2FREQ2        = Y_filt002(tempInds_2PLOT);
%         Ypredict00_2FREQ2      = Ypredict002(tempInds_2PLOT);
%         Ypredict_filt00_2FREQ2 = Ypredict_filt002(tempInds_2PLOT);   
% 
%     [FY2, FY_filt2, FYpredict2, FYpredict_filt2, frq_discretes2, ...
%      PY2, PY_filt2, PYpredict2, PYpredict_filt2, WW2, ...
%     ]           = stl_FFTPSD_function(Y00_2FREQ2, Y_filt00_2FREQ2, Ypredict00_2FREQ2, Ypredict_filt00_2FREQ2);
%                                                                               
%     %single error measures
%     MSEtime2    = (sum((Y_2PLOT2 - Ypredict_2PLOT2).^2))/length(Y_2PLOT2);
%     %calculate error for low frequencies
%     iindLOWFre  = find(WW2 <= 1);
%     MSEfreq2    = (sum((PY2(iindLOWFre) - PYpredict2(iindLOWFre)).^2))/length(iindLOWFre);
% 
%     MSEtimeArray2(cnt) = MSEtime2;
%     MSEfreqArray2(cnt) = MSEfreq2;
%     MSElenThelp2(cnt) = length(Y_2PLOT2);
%     MSElenFhelp2(cnt) = length(iindLOWFre);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%