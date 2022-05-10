%% Init
clc; clear; close all;
curr_path = 'C:\Users\georg\Desktop\fiber-draw';
% curr_path = 'D:\GEORGE\fiber-draw';
cd(curr_path);
load("run_results\architecture_experiment.mat");
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat");
folder_name = 'sim_bode_results';
if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

combined_Xdata = cat(2, Xdata{:});

labels = {'Capstan Speed', 'Furnace Power', 'Helium Temperature', 'Preform Velocity'};
mean_matrix = zeros(length(combined_Xdata),4);
for i = 1:length(combined_Xdata)
    subbatch = combined_Xdata{i};
    mean_matrix(i,:) = [mean(subbatch(1,:)) mean(subbatch(2,:)) mean(subbatch(3,:)) mean(subbatch(4,:))];
end
median1 = median(mean_matrix(:,1));
median2 = median(mean_matrix(:,2));
median3 = median(mean_matrix(:,3));
median4 = median(mean_matrix(:,4));

combined_Xdata = cat(2, Xdata{:});
combined_Ydata = cat(2, Ydata{:});

sampling_freq = 2; % Hz
nyq_freq = sampling_freq * 2*3 / 2; % rad/s
avg_capstan_speed = 2500; % 2500 or 2685
subbatch_length = 8000;

T = 1/sampling_freq;   % Sampling period
L = subbatch_length;   % Length of signal
t = (0:subbatch_length-1)*T;        % Time vector
f = sampling_freq*(0:(L/2))/L; % frequency vector

all_w = logspace(-2,log10(nyq_freq),100);
all_ffts = cell(4, 1);
num_peaks_counted = 10;
significant_peak_ratio = 0.05;
fit_opt = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'Bisquare');

plot_progress = true;

if (plot_progress) figure(1); latexify_plot; end

%% generate simulated bode (typical runtime: 2 hrs)
typical_A = [500 8 2.5 10]; % 1st, 2nd, 3rd, 4th row

for row = 1:4
    A = typical_A(row);
    all_ffts{row} = zeros(length(f), length(all_w));

    for j = 1:length(all_w)
        w = all_w(j);
        sine_wave = A * sin(w*t);

        % add sine wave to 1st, 2nd, 3rd, 4th row
        new_Xdata = [ones(1,subbatch_length)*avg_capstan_speed;
            ones(1,subbatch_length)*median2;
            ones(1,subbatch_length)*median3;
            ones(1,subbatch_length)*median4];

        new_Xdata(row,:) = new_Xdata(row,:) + sine_wave;

        % see, it's similar to input!
        % figure(6); plot(combined_Xdata{1}'); hold on; plot(new_Xdata', 'k'); hold off;

        net = deep_lstm.resetState();
        y_pred = net.predict(new_Xdata, "MiniBatchSize", 1);

        y_pred_fft = fft(y_pred);
        P2 = abs(y_pred_fft/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);

        fft_peaks = P1(islocalmax(P1));
        max_peak = max(fft_peaks);
        highest_peaks_inds = find(fft_peaks > significant_peak_ratio*max_peak, num_peaks_counted, 'first');
        highest_peaks = fft_peaks(highest_peaks_inds);

        peaks_plot = zeros(size(P1));
        for ind = 1:length(highest_peaks_inds)
            peaks_plot(P1 == highest_peaks(ind)) = highest_peaks(ind);
        end

        all_ffts{row}(:,j) = peaks_plot';

        %         freq_peaks = zeros(1, length(highest_peaks));
        %         for k = 1:length(highest_peaks)
        %             freq_peaks(k) = f(P1 == highest_peaks(k));
        %         end

        if plot_progress
            peaks_plot_progress = ones(size(P1))*nan;
            for ind = 1:length(highest_peaks_inds)
                peaks_plot_progress(P1 == highest_peaks(ind)) = highest_peaks(ind);
            end
            subplot(2,1,1); plot(t,y_pred);
            title(sprintf('Plant Output with Sinusoidal %s', labels{row}));
            xlabel('Time (s)'); ylabel('BFD Error');
            subplot(2,1,2); plot(f,P1); hold on;
            plot(f, peaks_plot_progress, '.', 'MarkerSize', 10); hold off;
            xlimit_ind = find(P1 > 1e-4, 1, 'last');
            xlim([0, f(xlimit_ind)]);
            xlabel('Frequency (Hz)'); ylabel('FFT Amplitude');
            title('FFT of Plant Output');
            pause(0.01)
        end

        fprintf('row: %d\t freq: %d/%d\n', row, j, length(all_w))
    end
end

if (row == 4 &&j == length(all_w))
    save(sprintf('%s\\%s\\bode_fft.mat', curr_path, folder_name),"all_ffts");
end

%% plot 3d
cd(curr_path);
load([curr_path '\' folder_name '\' 'bode_fft.mat']);
for row = 1:4
    figure; mesh(all_w/(2*pi), f, all_ffts{row});
    latexify_plot;
    set(gca, 'ZLim', get(gca, 'ZLim') .* [0 1] + [7e-6 0]);
    title(sprintf('Frequency Response w.r.t. %s', labels{row}))
    xlabel('Input frequency (Hz)');
    ylabel('Output frequency (Hz)');
    zlabel('FFT Amplitude');
    axis square; grid on; grid minor; box on;
end

%% gather bode profile of many amplitudes (estimated runtime: 4 hrs)
all_max_A = [700 10 3 15]; % 1st, 2nd, 3rd, 4th row
all_bode_profile = cell(4,1);
for row = 1:4
    if row == 1
        all_A = 5:5:all_max_A(row);
    else
        all_A = 0:0.1:all_max_A(row);
    end

    all_bode_profile{row} = zeros(length(all_A), length(all_w));

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
            % figure(6); plot(combined_Xdata{1}'); hold on; plot(new_Xdata', 'k'); hold off;

            net = deep_lstm.resetState();
            y_pred = net.predict(new_Xdata, "MiniBatchSize", 1);

            y_pred_fft = fft(y_pred);
            P2 = abs(y_pred_fft/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);

            all_bode_profile{row}(i,j) = max(P1);

            if plot_progress
                peaks_plot_progress = ones(size(P1))*nan;
                peaks_plot_progress(P1==max(P1)) = max(P1);
                disp('here:')
                disp([f(find(P1==max(P1),1,'first')), w/(2*pi)])
                subplot(2,1,1); plot(t,y_pred)
                subplot(2,1,2); plot(f,P1); hold on;
                plot(f, peaks_plot_progress, '.', 'MarkerSize', 20);
                hold off;
                xlimit_ind = find(P1 > 1e-4, 1, 'last');
                xlim([0, f(xlimit_ind)])
                pause(0.1)
            end

            fprintf('row: %d\t A: %d/%d\t freq: %d/%d\n', row, i, length(all_A), j, length(all_w))
        end
    end
end

if (row == 4 && i == length(all_A) && j == length(all_w))
    save(sprintf('%s\\%s\\all_bode_profile.mat', curr_path, folder_name),"all_bode_profile");
    disp('Saved and done!')
else
    disp('Done.')
end

%% plot bode profile (fft amplitude is different from actual)
clc; close all;
load('sim_bode_results\all_bode_profile.mat');
colors = {'black', 'red', 'blue', 'magenta'};
for i = 1:length(all_bode_profile)
    bode_profile = all_bode_profile{i};
    figure; axes('XScale', 'log', 'YScale', 'log')
    latexify_plot;
    switch i
        case 2
            all_jp = [14 52];
        case 4
            all_jp = [14 80]; % [14 25 50 80 97 150];
        case 1
            all_jp = [75,];
        case 3
            all_jp = [25,];
    end
    
    for j = 1:size(bode_profile, 1)
        p = loglog(all_w, bode_profile(j,:));
        p.Color(4) = 0.2; % set transparency
        hold on;
    end

    for jp = 1:length(all_jp)
        ind = all_jp(jp);
        plot(all_w, bode_profile(ind,:), 'LineWidth', 2, 'Color', colors{jp});
%         disp(jp);
%         pause(0.5);
%         hold off;
    end
    hold off;
    title(sprintf('Bode Plot of %s Input to BFD', labels{i}));
    xlabel('Input Frequency (rad/s)'); ylabel('Amplitude Ratio');

end

%% plot amplitude slice (not used)
all_w = logspace(-2,log10(nyq_freq),100);
w = all_w(50);
all_max_A = [700 10 3 15]; % 1st, 2nd, 3rd, 4th row
figure;
for row = 1:4
    if row == 1
        all_A = 5:5:all_max_A(row);
    else
        all_A = 0:0.1:all_max_A(row);
    end
    for j = 1:length(all_A)
        
        A = all_A(j);
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
        
        subplot(2,1,1); plot(t,y_pred)
        y_pred_fft = fft(y_pred);
        P2 = abs(y_pred_fft/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = sampling_freq*(0:(L/2))/L;
        subplot(2,1,2); plot(f,P1)
    
        xlimit_ind = find(P1 > 1e-4, 1, 'last')+50;
        xlim([0, 0.5]); %f(xlimit_ind)])
        pause(0.1)
        
    
        fprintf('row: %d\t freq: %d/%d\n', row, j, length(all_A))
    end
end
