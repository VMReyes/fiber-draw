



% FileSet{1} = 'Tower51_2021-01-05_to2021-01-12';
% FileSet{2} = 'Tower51_2020-12-22_to2020-12-29';
% FileSet{3} = 'Tower51_2020-12-15_to2020-12-22';
% FileSet{4} = 'Tower51_2020-12-08_to2020-12-15';
% FileSet{5} = 'Tower51_2020-12-01_to2020-12-08';
% 
FileSet{1} = 'Tower48_2020-12-01_to2020-12-08';  
FileSet{2} = 'Tower48_2020-12-08_to2020-12-15';  
FileSet{3} = 'Tower48_2020-12-15_to2020-12-22'; 
FileSet{4} = 'Tower48_2020-12-22_to2020-12-29'; 
FileSet{5} = 'Tower48_2020-12-29_to2021-01-05'; 
FileSet{6} = 'Tower48_2021-01-05_to2021-01-12';  
% FileSet{7} = 'Tower48_2021-01-12_to2021-01-19';  
% FileSet{8} = 'Tower48_2021-02-02_to2021-02-09';  
% FileSet{9} = 'Tower48_2021-02-23_to2021-03-02';  
% FileSet{10} = 'Tower48_2021-03-02_to2021-03-09';  
% FileSet{11} = 'Tower48_2021-03-09_to2021-03-16';  
% FileSet{12} = 'Tower48_2021-03-16_to2021-03-23';  

for iifs = 1:6
%for iifs = 1:12

%strbase = ['ttresults_a100_000_TRAINON_51AllofEachWeek_TESTON_AllOfOneweek_' FileSet{iifs} '__filter7']; %HERE for different tests
%strbase = ['ttresults_b100_000_TRAINON_51AllofEachWeek_TESTON_AllOfOneweek_' FileSet{iifs} '__filter5']; %HERE for different tests
strbase = ['ttresults_b100_000_TRAINON_48AllofEachWeek_TESTON_AllOfOneweek_' FileSet{iifs} '__filter5']; %HERE for different tests

%strbase = 'ttresults_a100_000_TRAINON_51AllofEachWeek_TESTON_AllOfOneweek_Tower51_2020-12-01_to2020-12-08__filter1';

strpath = 'E:\Dropbox (SquareCircleMITtoo)\_y_Code\zzResultsFolder\';
strfull = [strpath strbase]

mkdir([strpath strbase])
cd([strpath strbase])

clear myfigHANDLEarray;
clear myfigNUMBERarray;

load([strbase '.mat'])
hack_subplot_reorg
dev_figs2_2018

close all

end

