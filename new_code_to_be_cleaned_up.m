%% Init
% install Statistics and Machine Learning Toolbox
clear; clc; close all;

% parameters
strDataPath = 'C:\Users\georg\Dropbox (MIT)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
all_files   = dir(strDataPath);
curr_path   = 'C:\Users\georg\Desktop\fiber-draw';
% curr_path   = 'D:\GEORGE\fiber-draw';

% BatchInfo Parameters
bXLSLoad = 1;
bPlotAll = 0;
bPlot_each_preform_on_subplot = 1;
bPlot_each_preform_on_subplot_with_inrangesubbatches_ = 1;
loBFD = 124;
hiBFD = 126;
nTower = 48; % The tower number
subbatchMinLen 	= 2000; % a batch is the same as a preform, multiple
subbatchMaxLen  = 16000; 
x_columns = ["cpspdactval", "frnpwrmv", "pfspdactval"];
y_columns = ["barefibrediadisplay", "tenncmv"];

% Subbatch Parameters
fltLEN = 21;
bPlot = 0; % Plot batch
PrefltLEN = 1;

dt = 0.5;
orders_kd_oe = [1 4 3];

if exist('all_file_data.mat', 'file') ~= 2
    all_file_data = cell(length(all_files),2);
    
    for i = 1:length(all_files)
        curr_file = all_files(i);
        if ~curr_file.isdir
            % get batch info
            [BatchInfo, STRDEF] = stl_load_batchinfo(bXLSLoad, strDataPath, ...
                    curr_file.name, nTower, bPlotAll, ...
                    bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
                    bPlot_each_preform_on_subplot, loBFD, hiBFD, subbatchMinLen, subbatchMaxLen, ...
                    x_columns, y_columns);
            
            % turn it into a train / test array
            [XTrainTranspose, YTrainTranspose] = stl_prep_training_data(BatchInfo, ...
                STRDEF, x_columns, y_columns, fltLEN, PrefltLEN, bPlot, 0, 0);
    
            all_file_data{i,1} = XTrainTranspose;
            all_file_data{i,2} = YTrainTranspose;
        end
        fprintf('%d/%d files loaded...\n', i, length(all_files))
    end

    save('all_file_data','all_file_data')
else
    load("all_file_data.mat")
end

disp('Done Loading!')

%% PLOT ANALYSIS

folder_name = "_analysis";
cd(curr_path)


if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir

        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd'-125,     dt);
            sys_kd = oe(iddata_bfd_to_capstan_speed, orders_kd_oe);

            [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
            fprintf('file %d/%d\t subbatch %d/%d %2.4f\n', ...
                    file_ind,length(all_files),i,length(XTrainTranspose), fit)

            if (80 < abs(fit) && abs(fit) < 100)
                y = capstan_speed';
                y_hat = y_hat.OutputData;
                res = y - y_hat;
                [xcres,lags] = xcorr(res,res,'normalized');

                fig1 = figure(1); set(gcf, 'Position', [221 65 762 888])
                subplot(5,2,[1, 2]); compareplot(iddata_bfd_to_capstan_speed, sys_kd);
                ax = gca; ax.Legend.Location = 'southeast';

                subplot(5,2,[3, 4]); plot(res);
                subplot(5,2,[5, 6]); plot_fft_comparison(y, y_hat);

                subplot(5,2,7); plot(y,y_hat, '.'); hold on; 
%                 plot(linspace(min(min(y,y_hat)),max(max(y,y_hat)),length(y_hat)),linspace(min(min(y,y_hat)),max(max(y,y_hat)),length(y_hat)),'--');
                plot(linspace(min(y), max(y), length(y)),linspace(min(y), max(y), length(y)), 'LineWidth',2)
                hold off; xlabel('Actual'); ylabel('Predicted'); 

                subplot(5,2,8); histfit(res,50, 'logistic'); title('Residual Histogram')
                subplot(5,2,9); normplot(res);
                subplot(5,2,10); plot(lags, xcres); title('Residual Autocorrelation')
                
                [m,p] = bode(sys_kd, {10^-8, 10});
                m = m(:); p = p(:); m = 20 * log10(m); p = p - p(1);
                if (97 > max(m) && max(m) > 50 && max(p) < 100)
                    fig2 = figure(2); set(gcf, 'Position', [1048 222 754 696])
                    hold on; 
                    bodeopt = bodeoptions("cstprefs");
                    bodeopt.PhaseMatching = 'on';
                    bodeopt.PhaseMatchingFreq = 10^-8;
                    bodeopt.PhaseMatchingValue = 0;
                    bode(sys_kd, {10^-8, 10}, bodeopt);
                end
                
                saveas(fig1, sprintf('%s\\%s\\%d,%d-1', curr_path, folder_name, file_ind, i),'png');
            end

        end
    end
end
saveas(fig2, sprintf('%s\\%s\\%d,%d-2', curr_path, folder_name, file_ind, i),'png');  
disp('Done!')
%% helper functions
function plot_fft_comparison(Y, Ypredict)
    %take the FFT
    FY = fft(Y);
    FY_s = fftshift(FY);
    %take the FFT
    FYpredict = fft(Ypredict);
    FYpredict_s = fftshift(FYpredict);

    %the radial frequencies
    frq_discretes = 2/length(Y).*([0:(length(Y)-1)]-length(Y)/2);

    plot(frq_discretes,log10(abs(FY_s).^2),'k');
    hold on
    plot(frq_discretes,log10(abs(FYpredict_s).^2),'r'); 
    hold off
    
    ylabel('Log Magnitude')
    xlabel('Frequency')
    legend('Data', 'Prediction')
    axis([0 1 -.1  max([   max(log10(abs(FY_s).^2))   max(log10(abs(FYpredict_s).^2))   ])   ])
    title('Power Spectrum');
end