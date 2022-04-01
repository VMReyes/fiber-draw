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
%% simulated / empirical bode plot
load("run_results\architecture_experiment.mat");
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat");
net = deep_lstm.resetState();

% Xdata = [(num_subbatches, 4, 8000), (num_subbatches, 4, 800), .. for
% every file)]
% combined_Xdata.shape = (num_subbatches, 4 (inputs), 8000 (seq. length))
% bode_example = combined_Xdata{1}
% bode_example.shape = (4, 8000)
% bode_example{2} = sine_wave(amp, freq, len=8000)
% y_pred = net.predict(bode_example, "MiniBatchSize", 1)
% x_columns = ["cpspdactval", "frnpwrmv", "hetubetemp", "pfspdactval"]
% y_columns = ["barefibrediadisplay"]

% Xdata rows: capstan speed, furnace power, He temp, preform velocity

% stats_matrix = zeros(length(combined_Xdata),3);
% for i = 1:length(combined_Xdata)
%     subbatch = combined_Xdata{i};
%     stats_matrix(i,:) = [mean(subbatch(2,:)) mean(subbatch(3,:)) mean(subbatch(4,:))];
% end

% alist = zeros(1, length(combined_Xdata));
% for i = 1: length(combined_Xdata)
%     alist(i) = (max(combined_Xdata{i}(1,:)) - min(combined_Xdata{i}(1,:)));
% end
% figure; plot(alist,'.'); ylim([0, 2500])

combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});
new_Xdata = cell(length(combined_Xdata),1);

for i = 1:1 % length(combined_Xdata)

    subbatch = combined_Xdata{i};
    A = 250;
    f = 0.03;
    t = 1:length(subbatch);
    sine_wave = A*sin(2*pi*f*t);
    % figure(1); plot(t,combined_Ydata{1}); xlim([0 1e2]); hold on;
    % plot(t, sine_wave); hold off;

    new_Xdata{i} = [ones(1,length(subbatch))*mean(subbatch(1,:)) + sine_wave; 
                    ones(1,length(subbatch))*mean(subbatch(2,:)); 
                    ones(1,length(subbatch))*mean(subbatch(3,:)); 
                    ones(1,length(subbatch))*mean(subbatch(4,:))];
end

figure(1); plot(combined_Xdata{1}'); hold on; plot(new_Xdata{1}', 'k'); hold off;

y_pred = net.predict(new_Xdata{1}, "MiniBatchSize", 1);

eqn_str = sprintf('a*sin(2*pi*%s*x+c)+d',f);
ft = fittype(eqn_str, 'independent', 'x', 'dependent', 'y');
[fit_obj, goodness_info] = fit(t', y_pred', ft);
(mean(y_pred(islocalmax(y_pred))) - mean(y_pred(islocalmin(y_pred))))/2

figure(2); plot(y_pred); hold on; plot(fit_obj); xlim([0 1e3]); hold off;
