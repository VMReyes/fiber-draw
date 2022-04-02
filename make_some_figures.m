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

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

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

    dt = 1;
        iddata_tension_to_power     = iddata(furnace_power', tension', dt);
%             iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);

    subplot(2,1,1); plot(bfd); title('Input: Bare Fiber Diameter'); ylabel('BFD ($\mu m$)')
    subplot(2,1,2); plot(capstan_speed); title('Output: Capstan Speed'); ylabel('Capstan Speed (m/min)');
    xlabel('Samples')
end

%% data wrangling
disp('Started...')
for curr_file_ind = 3:3%16
% curr_file_ind = 3; %3 / 20;
example_bfd = readmatrix([strDataPath all_files(curr_file_ind).name], 'Range', 'C:C');
example_bfd(example_bfd > 126) = 125;
example_bfd(example_bfd < 124) = 125;

controller_on = readmatrix([strDataPath all_files(curr_file_ind).name], 'Range', 'Q:Q');

example_bfd(controller_on ~=1) = 125;

figure; set(gcf, 'Position', [96 558 1429 420]);
plot(example_bfd); ylim([124 126]); xlim([0 1e5])
title('Example Batch of BFD Data'); xlabel('Samples'); ylabel('Bare Fiber Diameter ($\mu m$)')
end

disp('Done!')
%% generate simulated / empirical bode 
clc; clear;
load("run_results\architecture_experiment.mat");
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat");
combined_Xdata = cat(2, Xdata{:});

% Xdata rows: capstan speed, furnace power, He temp, preform velocity

stats_matrix = zeros(length(combined_Xdata),4);
for i = 1:length(combined_Xdata)
    subbatch = combined_Xdata{i};
    stats_matrix(i,:) = [mean(subbatch(1,:)) mean(subbatch(2,:)) mean(subbatch(3,:)) mean(subbatch(4,:))];
end
% figure; plot(stats_matrix)

% alist = zeros(1, length(combined_Xdata));
% for i = 1:length(combined_Xdata)
%     alist(i) = (max(combined_Xdata{i}(1,:)) - min(combined_Xdata{i}(1,:)));
% end
% figure; histogram(alist,200,"BinLimits",[0 2700])

combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});

sampling_freq = 2; % Hz
nyq_freq = sampling_freq * 2*3 / 2;
avg_capstan_speed = 2500; % 2500 or 2685
subbatch_length = 8000;
t = 1:subbatch_length; % subbatch length

all_w = logspace(-3,log10(nyq_freq),100);
% all_A = 5:5:700; % for capstan speed
all_max_A = [700 10 3 20]; % 1st, 2nd, 3rd, 4th row

for row = 2:length(all_max_A) 
    all_A = 0:0.1:all_max_A(row);
    bode_p2p = zeros(length(all_A), length(all_w));
    bode_fit_obj = cell(length(all_A), length(all_w));
    
    for i = 1:length(all_A)
        for j = 1:length(all_w)
            
            A = all_A(i); w = all_w(j);
    
    %         sine_wave = (avg_capstan_speed-A/2) * sin(w*t);
        
            % modify 1st, 2nd, 3rd, 4th line to be sine wave
    %         new_Xdata = [ones(1,subbatch_length)*avg_capstan_speed + sine_wave; 
    %                      ones(1,subbatch_length)*median(stats_matrix(:,2)); 
    %                      ones(1,subbatch_length)*median(stats_matrix(:,3)); 
    %                      ones(1,subbatch_length)*median(stats_matrix(:,4))];
            
            sine_wave = A * sin(w*t);
            new_Xdata = [ones(1,subbatch_length)*avg_capstan_speed; 
                         ones(1,subbatch_length)*median(stats_matrix(:,2));
                         ones(1,subbatch_length)*median(stats_matrix(:,3)); 
                         ones(1,subbatch_length)*median(stats_matrix(:,4))];

            new_Xdata(row,:) = new_Xdata(row,:) + sine_wave;
    
            % see, it's similar to input!
            % figure(1); plot(combined_Xdata{1}'); hold on; plot(new_Xdata{1}', 'k'); hold off;
            net = deep_lstm.resetState();
            y_pred = net.predict(new_Xdata, "MiniBatchSize", 1);
    
            eqn_str = sprintf('a*sin(%f*x+phi)+c',w);
            ft = fittype(eqn_str, 'independent', 'x', 'dependent', 'y');
            fit_opt = fitoptions('Method','LinearLeastSquares' , 'Robust', 'Bisquare');
            [fit_obj, goodness_info] = fit(t', y_pred', ft);
    
            % identify peaks
            peaks_ind = islocalmax(y_pred); troughs_ind = islocalmin(y_pred);
            all_peaks = y_pred(peaks_ind);
            all_peaks = all_peaks(all_peaks > (fit_obj.c + abs(fit_obj.a)));
            all_troughs = y_pred(troughs_ind);
            all_troughs = all_troughs(all_troughs < (fit_obj.c)); 
            A_p2p = (mean(all_peaks(2:end)) - mean(all_troughs(2:end)))/2;
%             peaks = zeros(1,length(y_pred)); troughs = zeros(1,length(y_pred));
%             for k = 1:length(peaks_ind)
%                 if peaks_ind(k) && y_pred(k) > (fit_obj.c + abs(fit_obj.a))
%                     peaks(k) = y_pred(k);
%                 else
%                     peaks(k) = nan;
%                 end
%             end
%             for k = 1:length(troughs_ind)
%                 if troughs_ind(k) && y_pred(k) < (fit_obj.c)
%                     troughs(k) = y_pred(k);
%                 else
%                     troughs(k) = nan;
%                 end
%             end

                    
    
            figure(2); 
            plot(fit_obj, t, y_pred, '-'); hold on;
            plot(ones(1, length(subbatch))*max(all_peaks(2:end))); 
            plot(ones(1, length(subbatch))*min(all_troughs(2:end)));
            
%             plot(t, peaks, '.', 'MarkerSize', 20);
%             plot(t, troughs, '.', 'MarkerSize', 20);
            hold off; pause(2)

%             xlim([0 1e2]); 
            fprintf('row: %d\t i: %d/%d\t j:%d/%d\n', row, i, length(all_A), j, length(all_w))
        end
    end
    save(sprintf('simulated_bode_%d.mat',row),"bode_fit_obj",'bode_p2p');
end
%% plot simulated bode
clc; close all;
load simulated_bode.mat

figure; axes('XScale', 'log', 'YScale', 'log')
hold on;
for A = 1:size(bode_fit_obj,1)
    row = zeros(1,size(bode_fit_obj,2));
    for i = 1:length(row)
        row(i) = abs(bode_fit_obj{A,i}.a)/all_A(A);
    end
    plot(all_w, row)
end

figure; axes('XScale', 'log', 'YScale', 'log');
hold on;
for A = 1:10 %1:size(bode_p2p,1)
    plot(all_w, bode_p2p(A,:)./all_A(A))
end
xline(nyq_freq);
xlabel('Frequency (rad/s)'); ylabel('Gain');
title('Simulated Bode Plot of Fiber Drawing Plant')
xlim([1e-3 10]) %band-aid fix, need fix p2p calc & resample
