clc; close all;
model48 = load("run_results\architecture_experiment.mat", "deep_lstm"); % load model
model48 = model48.deep_lstm;
model51 = load("run_results\architecture_experiment_4in_1out_tower51.mat", "deep_lstm"); % load model
model51 = model51.deep_lstm;
%% LOAD DATA
% parameters
strDataPath = 'C:\Users\georg\Dropbox (MIT)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
all_files   = dir(strDataPath);
curr_path   = 'C:\Users\georg\Desktop\fiber-draw';

% BatchInfo Parameters
bXLSLoad = 1;
bPlotAll = 0;
bPlot_each_preform_on_subplot = 1;
bPlot_each_preform_on_subplot_with_inrangesubbatches_ = 1;
loBFD = 124;
hiBFD = 126;
nTower = 48; % The tower number
subbatchMinLen 	= 2000; 
subbatchMaxLen  = 16000;
x_columns = ["cpspdactval", "frnpwrmv", "hetubetemp", "pfspdactval"];
y_columns = ["barefibrediadisplay"];

% Subbatch Parameters
fltLEN = 21;
bPlot = 0; 
PrefltLEN = 1;
dt = 0.5;
bfd_setpoint = 125;

if exist('all_file_data_4_in_1_out.mat', 'file') ~= 2
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

    save('all_file_data_4_in_1_out','all_file_data')
else
    load("all_file_data_4_in_1_out.mat", 'all_file_data');
end

disp('Done Loading!')

%% INFERENCE & RMSE ON MODEL48

errors_48 = [];
net = model48.resetState();

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_file_data{file_ind,1};
    if ~isempty(curr_file)

        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i)) - bfd_setpoint;
            net = model48.resetState();
            y_pred = net.predict(xs, "MiniBatchSize", 1);
            rmse = sqrt(sum((y_pred - ys).^2/length(y_pred)));
            errors_48(end+1) = rmse;
            fprintf('tower: 48\t file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_file_data),i,length(XTrainTranspose), rmse)
        end
    end
end


errors_51 = [];
net = model48.resetState();

for file_ind = 18:length(all_files) % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_file_data{file_ind,1};
    if ~isempty(curr_file)

        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i)) - bfd_setpoint;
            net = model48.resetState();
            y_pred = net.predict(xs, "MiniBatchSize", 1);
            rmse = sqrt(sum((y_pred - ys).^2/length(y_pred)));
            errors_51(end+1) = rmse;
            fprintf('tower: 51\t file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_file_data),i,length(XTrainTranspose), rmse)
        end
    end
end

save("model48_errors", "errors_48", "errors_51");

%% INFERENCE & RMSE ON MODEL51

errors_48 = [];
net = model51.resetState();

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_file_data{file_ind,1};
    if ~isempty(curr_file)

        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i)) - bfd_setpoint;
            net = model51.resetState();
            y_pred = net.predict(xs, "MiniBatchSize", 1);
            rmse = sqrt(sum((y_pred - ys).^2/length(y_pred)));
            errors_48(end+1) = rmse;
            fprintf('tower: 48\t file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_file_data),i,length(XTrainTranspose), rmse)
        end
    end
end


errors_51 = [];
net = model51.resetState();

for file_ind = 18:length(all_files) % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_file_data{file_ind,1};
    if ~isempty(curr_file)

        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i)) - bfd_setpoint;
            net = model51.resetState();
            y_pred = net.predict(xs, "MiniBatchSize", 1);
            rmse = sqrt(sum((y_pred - ys).^2/length(y_pred)));
            errors_51(end+1) = rmse;
            fprintf('tower: 51\t file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_file_data),i,length(XTrainTranspose), rmse)
        end
    end
end

save("model51_errors", "errors_48", "errors_51");

%% plots
rmse48 = load("model48_errors", "errors_48", "errors_51");
rmse51 = load("model51_errors", "errors_48", "errors_51");

figure(1); 
histogram(rmse48.errors_48(1:length(errors_51)), 50);
hold on;
% xline(mean(errors_48),'r')
histogram(rmse48.errors_51, 50);
% xline(mean(errors_51),'b')
hold off;
title('RMSE Distribution using Tower 48 Model'); 
ylabel('Density'); xlabel('RMSE')
legend({'Tower 48', 'Tower 51'})
latexify_plot;

figure(2); 
histogram(rmse51.errors_48(1:length(errors_51)), 50);
hold on;
% xline(mean(errors_48),'r')
histogram(rmse51.errors_51, 50);
% xline(mean(errors_51),'b')
hold off;

title('RMSE Distribution using Tower 51 Model'); 
ylabel('Density'); xlabel('RMSE')
legend({'Tower 48', 'Tower 51'})
latexify_plot;

%% Freq content


folder_name = 'tower_freq';
if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 18:length(all_files) % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_file_data{file_ind,1};
    if ~isempty(curr_file)

        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        for i = 1:length(XTrainTranspose)
            f = figure(2);
            xs = cell2mat(XTrainTranspose(i)); capstan_speed = xs(1,:);
            ys = cell2mat(YTrainTranspose(i)) - bfd_setpoint;
%             net = model48.resetState();
%             y_pred = net.predict(xs, "MiniBatchSize", 1);
            t = 0:0.5:floor(length(ys)/2);
            if (length(t) ~= length(ys)) t = t(1:length(ys)); end
            subplot(2,1,1); plot(t, ys); % hold on; plot(t, y_pred); hold off;
%             legend({'Data', 'Prediction'}); 
            title('Tower 51 BFD')
            xlabel('Time (s)'); ylabel('BFD Error ($\mu m$)')
            xlim([0 length(ys)/2])
            subplot(2,1,2);
            plot_fft_comparison(ys, ys);
            latexify_plot;
            saveas(f, sprintf('%s\\%i,%i', folder_name, file_ind, i),'png')
        end
    end
end

%%
function plot_fft_comparison(Y, Ypredict)
    %take the FFT
    FY = fft(Y);
    FY_s = fftshift(FY);
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
%     legend('Data', 'Prediction')
    axis([0 1 -1  max([   max(log10(abs(FY_s).^2))   max(log10(abs(FYpredict_s).^2))   ])   ])
    title('Power Spectrum');
end
