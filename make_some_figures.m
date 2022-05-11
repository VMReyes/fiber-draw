%% Init
% parameters
strDataPath         = 'C:\Users\georg\Dropbox (MIT)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
all_files = dir(strDataPath);
curr_path = 'C:\Users\georg\Desktop\fiber-draw';

% BatchInfo Parameters
bXLSLoad = 1;
bPlotAll = 0;
bPlot_each_preform_on_subplot = 1;
bPlot_each_preform_on_subplot_with_inrangesubbatches_ = 1;
loBFD = 124;
hiBFD = 126;
nTower = 48; % The tower number
subbatchMinLen 	= 2000; % a batch is the same as a preform, multiple
subbatchMaxLen  = 16000; % batches (or preforms) are run, one after the
% other in the tower. a subbatch is defined as a
% contiguous region of production
x_columns = ["cpspdactval", "frnpwrmv", "pfspdactval"];
y_columns = ["barefibrediadisplay", "tenncmv"];

% Subbatch Parameters
fltLEN = 21;
bPlot = 0; % Plot batch
PrefltLEN = 1;

if exist('all_file_data.mat', 'file') ~= 2
    all_file_data = cell(length(all_files),2);
    
    for i = 1:length(all_files)
        curr_file = all_files(i);
        if ~curr_file.isdir
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
    end

    save('all_file_data','all_file_data')
else
    load("all_file_data.mat")
end

disp('Done Loading!')

%% Plot Input / Output (specific subbatch)
file_num = 12;
subbatch_num = 17;

curr_file = all_files(file_num);
if ~curr_file.isdir

    XTrainTranspose = all_file_data{file_num,1}; 
    YTrainTranspose = all_file_data{file_num,2};

    fig = figure(1);
    set(gcf, 'Position', [229 222 754 696])

    xs = cell2mat(XTrainTranspose(subbatch_num));
    ys = cell2mat(YTrainTranspose(subbatch_num));
    capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
    bfd = ys(1,:); tension = ys(2,:);

    subplot(2,1,1); plot(bfd-125); title('Input: Bare Fiber Diameter Error'); ylabel('BFD Error ($\mu m$)')
    subplot(2,1,2); plot(capstan_speed); title('Output: Capstan Speed'); ylabel('Capstan Speed (m/min)');
    xlabel('Samples'); latexify_plot
end

%% Plot Input / Output (specific subbatch)

file_num = 4;
subbatch_num = 7;

% for file_num = 1:16
    curr_file = all_files(file_num);
    if ~curr_file.isdir

        XTrainTranspose = all_file_data{file_num,1};
        YTrainTranspose = all_file_data{file_num,2};

%         for subbatch_num = 1:length(XTrainTranspose)
            fig = figure(1);
            set(gcf, 'Position', [229 222 754 696])

            xs = cell2mat(XTrainTranspose(subbatch_num));
            ys = cell2mat(YTrainTranspose(subbatch_num));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            dt = 0.5;
            iddata_tension_to_power     = iddata(furnace_power', tension', dt);
            %             iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);
            fprintf('%d\t %d\n', file_num, subbatch_num)
            subplot(2,1,1); plot(tension - median(tension)); title('Input: Tension Error'); ylabel('Tension (g)')
            xlim([0 length(tension)])
            subplot(2,1,2); plot(furnace_power); title('Output: Furnace Power'); ylabel('Furnace Power (\%)');
            xlim([0 length(tension)])
            xlabel('Samples'); latexify_plot; pause(0.5)
            
%         end
    end
% end
%% data wrangling
disp('Started...')
for curr_file_ind = 3:3%16 %3 / 20;
% example_bfd = readmatrix([strDataPath all_files(curr_file_ind).name], 'Range', 'C:C');
% example_bfd(example_bfd > 126) = 125;
% example_bfd(example_bfd < 124) = 125;

controller_on = readmatrix([strDataPath all_files(curr_file_ind).name], 'Range', 'Q:Q');

% example_bfd(controller_on ~=1) = 125;

tension = readmatrix([strDataPath all_files(curr_file_ind).name], 'Range', 'G:G');
tension(controller_on ~= 1) = 120;

figure; set(gcf, 'Position', [96 558 1429 420]);
plot(example_bfd); ylim([124 126]); xlim([0 1e5])
title('Example Batch of BFD Data'); xlabel('Samples'); ylabel('Bare Fiber Diameter ($\mu m$)')
end

disp('Done!')
%% generate simulated bode (typical runtime: 2 hrs)
clc; clear; close all;
load("run_results\architecture_experiment.mat");
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat");
curr_path = 'C:\Users\georg\Desktop\fiber-draw';
folder_name = 'sim_bode_results';
if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

combined_Xdata = cat(2, Xdata{:});

% Xdata rows: capstan speed, furnace power, He temp, preform velocity

mean_matrix = zeros(length(combined_Xdata),4);
for i = 1:length(combined_Xdata)
    subbatch = combined_Xdata{i};
    mean_matrix(i,:) = [mean(subbatch(1,:)) mean(subbatch(2,:)) mean(subbatch(3,:)) mean(subbatch(4,:))];
end
median1 = median(mean_matrix(:,1));
median2 = median(mean_matrix(:,2));
median3 = median(mean_matrix(:,3));
median4 = median(mean_matrix(:,4));
% figure; plot(stats_matrix)

% visualize the distribution of ranges (max-min) of signal
% for row = 1:4
%     ranges = zeros(1, length(combined_Xdata));
%     for i = 1:length(combined_Xdata)
%         ranges(i) = (max(combined_Xdata{i}(row,:)) - min(combined_Xdata{i}(row,:)));
%     end
%     figure; 
%     
%     switch row
%         case 1
%             histogram(ranges,200,"BinLimits",[0 2700]); xline(700, 'LineWidth',2); % 1st row
%         case 2
%             histogram(ranges,100,"BinLimits",[0 15]); xline(10, 'LineWidth',2); % 2nd row
%         case 3
%             histogram(ranges,100,"BinLimits",[0 5]); xline(3, 'LineWidth',2); % 3rd row
%         case 4
%             histogram(ranges,100,"BinLimits",[0 35]); xline(15, 'LineWidth',2) % 4th row
%     end
%     latexify_plot;
% end

combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});

sampling_freq = 2; % Hz
nyq_freq = sampling_freq * 2*3 / 2;
avg_capstan_speed = 2500; % 2500 or 2685
subbatch_length = 8000;
t = 1:subbatch_length; 

all_w = logspace(-3,log10(nyq_freq),100);
all_max_A = [700 10 3 15]; % 1st, 2nd, 3rd, 4th row
fit_opt = fitoptions('Method','LinearLeastSquares' , 'Robust', 'Bisquare');

bPlot = 0;

for row = 1:4
    if row == 1
        all_A = 5:5:all_max_A(row);
    else
        all_A = 0:0.1:all_max_A(row);
    end

    bode_p2p = zeros(length(all_A), length(all_w));
    bode_fit_obj = cell(length(all_A), length(all_w));
    
    for i = 1:length(all_A)
        for j = 1:length(all_w)
            
            A = all_A(i); w = all_w(j);
            sine_wave = A * sin(w*t);
            
            % add sine wave to 1st, 2nd, 3rd, 4th row
            new_Xdata = [ones(1,subbatch_length)*avg_capstan_speed; 
                         ones(1,subbatch_length)*median2;
                         ones(1,subbatch_length)*median3; 
                         ones(1,subbatch_length)*median4];

            new_Xdata(row,:) = new_Xdata(row,:) + sine_wave;
    
            % see, it's similar to input!
%             figure(6); plot(combined_Xdata{1}'); hold on; plot(new_Xdata', 'k'); hold off;

            net = deep_lstm.resetState();
            y_pred = net.predict(new_Xdata, "MiniBatchSize", 1);
    
            eqn_str = sprintf('a*sin(%f*x+phi)+c',w);
            ft = fittype(eqn_str, 'independent', 'x', 'dependent', 'y');
            [fit_obj, goodness_info] = fit(t', y_pred', ft);
    
            % identify peaks
            peaks_ind = islocalmax(y_pred); troughs_ind = islocalmin(y_pred);
            all_peaks = y_pred(peaks_ind);
            all_peaks = all_peaks(all_peaks > (fit_obj.c + abs(fit_obj.a)/2));
            all_troughs = y_pred(troughs_ind);
            all_troughs = all_troughs(all_troughs < (fit_obj.c - abs(fit_obj.a)/2)); 
            A_p2p = (median(all_peaks) - median(all_troughs))/2;
            bode_fit_obj{i,j} = fit_obj; bode_p2p(i,j) = A_p2p;

%             if bPlot
%                 peaks = zeros(1,length(y_pred)); troughs = zeros(1,length(y_pred));
%                 for k = 1:length(peaks_ind)
%                     if peaks_ind(k) && y_pred(k) > (fit_obj.c + abs(fit_obj.a)/2)
%                         peaks(k) = y_pred(k);
%                     else
%                         peaks(k) = nan;
%                     end
%                 end
%                 for k = 1:length(troughs_ind)
%                     if troughs_ind(k) && y_pred(k) < (fit_obj.c + abs(fit_obj.a)/2)
%                         troughs(k) = y_pred(k);
%                     else
%                         troughs(k) = nan;
%                     end
%                 end
%     
%                 figure(7); set(gcf, 'Position', [31 500 1815 478]);
%                 plot(fit_obj, t, y_pred, '-'); hold on;
%                 plot(ones(1, length(subbatch))*median(all_peaks),'LineWidth',3); 
%                 plot(ones(1, length(subbatch))*median(all_troughs),'LineWidth',3);
%                 
%                 plot(t, peaks, '.', 'MarkerSize', 20);
%                 plot(t, troughs, '.', 'MarkerSize', 20);
%                 hold off; 
%                 if (j > 50) xlim([0 1e2]); end
% %                 saveas(gcf, sprintf('%s\\sim_bode_plots\\%d,%d,%d', curr_path, row, i, j),'fig');
% %                 saveas(gcf, sprintf('%s\\sim_bode_plots\\%d,%d,%d', curr_path, row, i, j),'png');
%             end

            fprintf('row: %d\t i: %d/%d\t j:%d/%d\n', row, i, length(all_A), j, length(all_w))
        end
    end
    if (i == length(all_A) && j == length(all_w))
        save(sprintf('%s\\%s\\simulated_bode_%d.mat', curr_path, folder_name, row),"bode_fit_obj",'bode_p2p');
    end
end
%% plot simulated bode
clc; close all;

row = 1; % 1, 2, 3, 4

load(sprintf('sim_bode_results\\simulated_bode_%d.mat',row));

sampling_freq = 2; % Hz
nyq_freq = sampling_freq * 2*3 / 2;
avg_capstan_speed = 2500; % 2500 or 2685
subbatch_length = 8000;
t = 1:subbatch_length; 
all_w = logspace(-3,log10(nyq_freq),100);
all_max_A = [700 10 3 20]; % 1st, 2nd, 3rd, 4th row
if row == 1
    all_A = 5:5:all_max_A(row);
else
    all_A = 0:0.1:all_max_A(row);
end

figure; axes('XScale', 'log', 'YScale', 'log')
hold on;
for A = 2:size(bode_fit_obj,1)
    if A~=12
        row = zeros(1,size(bode_fit_obj,2));
        for i = 1:length(row)
            row(i) = abs(bode_fit_obj{A,i}.a)/all_A(A);
        end
        plot(all_w, row) % mag2db(row)?
    end
%     pause(0.5)
end
xline(nyq_freq);
xlabel('Frequency (rad/s)'); ylabel('Gain');
title('Simulated Bode Plot of Fiber Drawing Plant')
grid on; grid minor;
xlim([1e-3 10])
latexify_plot;

figure; axes('XScale', 'log', 'YScale', 'log');
hold on;
for A = 2:size(bode_p2p,1)
    plot(all_w, bode_p2p(A,:)./all_A(A)) % mag2db()?
end
xline(nyq_freq);
xlabel('Frequency (rad/s)'); ylabel('Gain');
title('Simulated Bode Plot of Fiber Drawing Plant')
grid on; grid minor;
xlim([1e-3 10]) 
latexify_plot;