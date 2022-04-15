%% Init
% install Statistics and Machine Learning Toolbox
clear; clc; close all;

% parameters
strDataPath = 'C:\Users\georg\Dropbox (MIT)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
all_files   = dir(strDataPath);
curr_path = 'C:\Users\georg\Desktop\MS Research\fiber-draw';
% curr_path   = 'C:\Users\georg\Desktop\fiber-draw';
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
orders_kd_armax = [4 3 1 3];

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

folder_name = "_analysis_armax";
cd(curr_path)

plot_residuals = false;
plot_bode = true;

if plot_bode 
    bodeopt = bodeoptions("cstprefs");
    bodeopt.PhaseMatching = 'on';
    bodeopt.PhaseMatchingFreq = 10^-8;
    bodeopt.PhaseMatchingValue = 0;
    fig2 = figure(2); set(gcf, 'Position', [4589.8 -225.4 809.6 1008])
    hold on; grid on; grid minor; latexify_plot
end

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
%             sys_kd = oe(iddata_bfd_to_capstan_speed, orders_kd_oe);
            sys_kd = armax(iddata_bfd_to_capstan_speed, orders_kd_armax);

            [y_hat, fit, ~] = compare(iddata_bfd_to_capstan_speed, sys_kd);
            fprintf('file %d/%d\t subbatch %d/%d %2.4f\n', ...
                    file_ind,length(all_files),i,length(XTrainTranspose), fit)

            if (80 < fit && fit < 100)
                if plot_residuals
                    y = capstan_speed';
                    y_hat = y_hat.OutputData;
                    res = y - y_hat;
                    [xcres,lags] = xcorr(res,res,'normalized');
    
                    fig1 = figure(1); set(gcf, 'Position', [4589.8 -225.4 809.6 1008]); latexify_plot;
                    subplot(5,2,[1, 2]); compareplot(iddata_bfd_to_capstan_speed, sys_kd);
                    ax = gca; ax.Legend.Location = 'southeast'; 
                    
                    t = linspace(0, length(res)/2, length(res));
                    subplot(5,2,[3, 4]); plot(t,res); hold on; plot(t, ones(1,length(res))*mean(res), '--'); hold off;
                    ylabel('Residual'); xlabel('Time (seconds)'); legend({'Residual', 'MSE'}); 
                    xlim([0 length(res)/2])
                    subplot(5,2,[5, 6]); plot_fft_comparison(y, y_hat);
                    
                    subplot(5,2,7); plot(y, y_hat, '.'); hold on; 
                    plot(0:max(max(y),max(y_hat)), 0:max(max(y),max(y_hat)),'--','LineWidth',2);
                    xlim([min(min(y),min(y_hat)), max(max(y),max(y_hat))]);
                    ylim([min(min(y),min(y_hat)), max(max(y),max(y_hat))]);
                    hold off; xlabel('Actual'); ylabel('Predicted'); 
                    
                    subplot(5,2,8); histogram(res, 50); title('Residual Histogram')
                    subplot(5,2,9); normplot(res);
                    subplot(5,2,10); plot(lags, xcres); title('Residual Autocorrelation')
    
                    saveas(fig1, sprintf('%s\\%s\\%d,%d', curr_path, folder_name, file_ind, i),'png');
                end

                if plot_bode
                    fig2 = figure(2); [m,p] = bode(sys_kd, {10^-8, 10});
                    m = m(:); p = p(:); m = 20 * log10(m); p = p - p(1);
                    bode(sys_kd, {10^-8, 10}, bodeopt);
                end
            end

        end
    end
end

if plot_bode
    fig2 = figure(2); hold off;
    saveas(fig2, sprintf('%s\\%s\\all_bode', curr_path, folder_name),'png');  
    saveas(fig2, sprintf('%s\\%s\\all_bode', curr_path, folder_name),'fig');  
end
disp('Done!')

%% Plotting for thesis
ind = [8 18]; %[7 39], [12 4], [8 18]
file_num = ind(1); subbatch_num = ind(2);
curr_file = all_files(file_num);
XTrainTranspose = all_file_data{file_num,1};
YTrainTranspose = all_file_data{file_num,2};

xs = cell2mat(XTrainTranspose(subbatch_num));
ys = cell2mat(YTrainTranspose(subbatch_num));
capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
bfd = ys(1,:); tension = ys(2,:);
orders_kd_armax = [4 3 1 3];

iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd'-125,     dt);
sys_kd = armax(iddata_bfd_to_capstan_speed, orders_kd_armax);
[y_hat, fit, ~] = compare(iddata_bfd_to_capstan_speed, sys_kd);
y = capstan_speed';
y_hat = y_hat.OutputData;
res = y - y_hat;
[xcres,lags] = xcorr(res,res,'normalized');

fig1 = figure(1); set(gcf, 'Position', [4589.8 -225.4 809.6 1008]); latexify_plot;
subplot(5,2,[1, 2]); compareplot(iddata_bfd_to_capstan_speed, sys_kd);
ax = gca; ax.Legend.Location = 'southeast'; 


subplot(5,2,[3, 4]); plot(t,res); hold on; plot(ones(1,length(res))*mean(res), '--'); hold off;
ylabel('Residual'); xlabel('Time (seconds)'); legend({'Residual', 'MSE'}); 
subplot(5,2,[5, 6]); plot_fft_comparison(y, y_hat);

subplot(5,2,7); plot(y, y_hat, '.'); hold on; 
plot(0:max(max(y),max(y_hat)), 0:max(max(y),max(y_hat)),'--','LineWidth',2);
xlim([min(min(y),min(y_hat)), max(max(y),max(y_hat))]);
ylim([min(min(y),min(y_hat)), max(max(y),max(y_hat))]);
hold off; xlabel('Actual'); ylabel('Predicted'); 

subplot(5,2,8); histogram(res, 50); title('Residual Histogram')
subplot(5,2,9); normplot(res);
subplot(5,2,10); plot(lags, xcres); title('Residual Autocorrelation')

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