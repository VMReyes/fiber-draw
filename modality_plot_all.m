clc; clear; close all; 
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

%% Choose model / order (BY SINGLE SUBBATCH) (ARMAX)
diary on; dt = 0.5;
fit_matrix_bool = cell(length(all_files), 1);
fit_matrix = cell(length(all_files), 1);
folder_name = "_armax_all_2";

if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:length(all_files) 
    curr_file = all_files(file_ind);
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

        fit_matrix_bool{file_ind} = zeros(length(XTrainTranspose), 8*8*8*4);
        fit_matrix{file_ind} = zeros(length(XTrainTranspose), 8*8*8*4);

%         cd 'C:\Users\georg\Desktop\MS Research\fiber-draw'
%         cd 'C:\Users\georg\Dropbox (MIT)\fiber-draw'
        cd(curr_path)

        fig = figure(1);
        set(gcf, 'Position', [229 222 754 696])

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            %     sys_kt = armax(iddata_tension_to_power,  [2 2 2 1]);
            %     sys_kt = bj(iddata_tension_to_power,  [2 2 2 2 1]);
            %     sys_kt = oe(iddata_tension_to_power,  [2 2 1]);

            %     iddata_tension_to_power     = iddata(furnace_power', tension', dt);
            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);

            for a = 1:8 for b = 1:8 for c = 1:8 for k = 1:4
                    %     sys_kt = armax(iddata_tension_to_power,  [a b c 1]);
                    %     sys_kt = bj(iddata_tension_to_power,  [a b c d 1]);
                if ~(a < 3 && b < 3 && c < 3)
                    sys_kd = armax(iddata_bfd_to_capstan_speed,  [a b c k]);
                    [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
                    %     [a b c d fit]
                    fprintf('file %d/%d\t subbatch %d/%d \t %d \t %d \t %d \t %d \t %2.4f\n', ...
                            file_ind,length(all_files),i,length(XTrainTranspose), a, b, c, k, fit)

                    if (30 < abs(fit) && abs(fit) < 150)
                        fit_matrix_bool{file_ind}(i, (a-1)*(8*8*4)+(b-1)*(8*4)+(c-1)*4+k) = 1;
                        fit_matrix{file_ind}(i, (a-1)*(8*8*4)+(b-1)*(8*4)+(c-1)*4+k) = abs(fit);
                        subplot(2,1,1); compare(iddata_bfd_to_capstan_speed, sys_kd);
                        title(sprintf('%d, %d, (%d,%d,%d,%d)', file_ind,i,a,b,c,k));
                        ax = gca; ax.Legend.Location = 'southeast';

                        subplot(2,1,2);
                        plot_fft(capstan_speed, y_hat.OutputData)

                        saveas(fig, sprintf('%s\\%s\\%d,%d,(%d,%d,%d,%d)', curr_path, folder_name, file_ind,i,a,b,c,k),'png');
                    end
                end
            end; end; end; end
        end
    end
end
disp('Done!')
diary off;
save('fit_results'+folder_name, 'fit_matrix', 'fit_matrix_bool');

%% Choose model / order (BY SINGLE SUBBATCH) (OE)
max_order = 10;

diary on;
fit_matrix_bool = cell(length(all_files), 1);
fit_matrix = cell(length(all_files), 1);
folder_name = "_kd_oe_all";

if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:length(all_files) 
    curr_file = all_files(file_ind);
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

        fit_matrix_bool{file_ind} = zeros(length(XTrainTranspose), max_order*max_order*4);
        fit_matrix{file_ind} = zeros(length(XTrainTranspose), max_order*max_order*4);

%         cd 'C:\Users\georg\Desktop\MS Research\fiber-draw'
%         cd 'C:\Users\georg\Dropbox (MIT)\fiber-draw'
        cd(curr_path)

        fig = figure(1);
        set(gcf, 'Position', [229 222 754 696])

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            dt = 1;
                iddata_tension_to_power     = iddata(furnace_power', tension', dt);
%             iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);

            for b = 1:max_order for a = 1:b for k = 1:4
                    if ~(a<3 && b<3)
                        %     sys_kt = armax(iddata_tension_to_power,  [a b c 1]);
                        %     sys_kt = bj(iddata_tension_to_power,  [a b c d 1]);

%                         sys_kd = oe(iddata_bfd_to_capstan_speed,  [a b k]);
                        sys_kt = oe(iddata_tension_to_power,  [a b k]);
%                         [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
                        [y_hat, fit, x0] = compare(iddata_tension_to_power, sys_kt);
                        %     [a b c d fit]
                        fprintf('file %d/%d\t subbatch %d/%d \t %d \t %d \t %d \t %2.4f\n', ...
                                file_ind,length(all_files),i,length(XTrainTranspose), a, b, k, fit)
    
                        if (30 < abs(fit) && abs(fit) < 150)
                            fit_matrix_bool{file_ind}(i, (a-1)*max_order*4+(b-1)*4+k) = 1;
                            fit_matrix{file_ind}(i, (a-1)*max_order*4+(b-1)*4+k) = abs(fit);
                            subplot(2,1,1); 
%                             compare(iddata_bfd_to_capstan_speed, sys_kd);
                            compare(iddata_tension_to_power, sys_kt);
                            title(sprintf('%d, %d, (%d,%d,%d)', file_ind,i,a,b,k));
                            ax = gca; ax.Legend.Location = 'southeast';
    
                            subplot(2,1,2);
                            plot_fft(capstan_speed, y_hat.OutputData)
    
                            cd(folder_name);
                            saveas(fig, sprintf('%d,%d,(%d,%d,%d)',file_ind,i,a,b,k),'png');
                            cd ..;
                        end
                    end
            end; end; end

        end
    end
end
disp('Done!')
diary off;
save('fit_results'+folder_name, 'fit_matrix', 'fit_matrix_bool');
%% obtain coeffs for modality
diary on;
all_as = cell(length(all_files),1);
all_bs = cell(length(all_files),1);
all_cs = cell(length(all_files),1);
all_ds = cell(length(all_files),1);
all_fs = cell(length(all_files),1);

for file_ind = 1:length(all_files)
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        all_as{file_ind} = zeros(1,length(XTrainTranspose));
        all_bs{file_ind} = zeros(7,length(XTrainTranspose));
        all_cs{file_ind} = zeros(1,length(XTrainTranspose));
        all_ds{file_ind} = zeros(1,length(XTrainTranspose));
        all_fs{file_ind} = zeros(8,length(XTrainTranspose));

        % get batch info
        [BatchInfo, STRDEF] = stl_load_batchinfo(bXLSLoad, strDataPath, ...
            curr_file.name, nTower, bPlotAll, ...
            bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
            bPlot_each_preform_on_subplot, loBFD, hiBFD, subbatchMinLen, subbatchMaxLen, ...
            x_columns, y_columns);

        % turn it into a train / test array
        [XTrainTranspose, YTrainTranspose] = stl_prep_training_data(BatchInfo, ...
            STRDEF, x_columns, y_columns, fltLEN, PrefltLEN, bPlot, 0);

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);
            dt = 1;
%             iddata_tension_to_power     = iddata(furnace_power', tension', dt);
            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);

            % replace for different graph
            sys_kd = oe(iddata_bfd_to_capstan_speed,  [6 7 1]);
            all_as{file_ind}(i) = sys_kd.A;
            all_bs{file_ind}(:,i) = sys_kd.B;
            all_cs{file_ind}(i) = sys_kd.C;
            all_ds{file_ind}(i) = sys_kd.D;
            all_fs{file_ind}(:,i) = sys_kd.F;
            fprintf("file %d/%d \t subbatch %d/%d\n", file_ind, length(all_files), i, length(XTrainTranspose))
        end
    end
end
diary off;
%% plotting modality
bs = zeros(7,1);
fs = zeros(8,1);
for i = 3:24
    for j = 1:length(all_bs{i})
        if bs == zeros(7,1)
            bs = all_bs{i}(:,j);
        else
            bs = horzcat(bs, all_bs{i}(:,j));
        end
    end
    for j = 1:length(all_fs{i})
        if fs == zeros(8,1)
            fs = all_fs{i}(:,j);
        else
            fs = horzcat(fs, all_fs{i}(:,j));
        end
    end
end
% remove cols of zeros
bs=bs(:,any(bs));
fs=fs(:,any(fs));

figure; 
set(gcf, 'Position', [229 222 754 696])
for i = 2:7
    subplot(2,3,i-1); 
    h = histogram(bs(i,:),200,"BinLimits",[-50 50]);
    title(sprintf('B(%d) Distribution',i));
    ylabel('Frequency');
end

figure; 
set(gcf, 'Position', [229 222 754 696])
for i = 2:8
    subplot(2,4,i-1); 
    histogram(fs(i,:),100);
    histfit(fs(i,:),100);
    d = fitdist(fs(i,:)', 'Normal');
    fprintf('%1.6f\t%1.6f\n',d.mu, d.sigma)
    title(sprintf('F(%d) Distribution',i));
    ylabel('Frequency');
end

%% analyze - convert cell to matrix

fit_matrix_all = [];
fit_matrix_bool_all = [];

for i = 1:length(fit_matrix)
    fit_matrix_all = vertcat(fit_matrix_all,fit_matrix{i});
    fit_matrix_bool_all = vertcat(fit_matrix_bool_all,fit_matrix_bool{i});
end

%% PLOT ANALYSIS

folder_name = "_analysis";

% cd 'C:\Users\georg\Dropbox (MIT)\fiber-draw'
cd(curr_path)

if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:length(all_files) 
    curr_file = all_files(file_ind);
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

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);
            
            dt = 1; 
%             a = 5; b = 5; k = 3;
            a = 4; b = 2; c = 2; k = 2;
            %     iddata_tension_to_power     = iddata(furnace_power', tension', dt);
            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);
%             sys_kd = oe(iddata_bfd_to_capstan_speed,  [a b k]);
            sys_kd = armax(iddata_bfd_to_capstan_speed,  [a b c k]);

            [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
            fprintf('file %d/%d\t subbatch %d/%d %2.4f\n', ...
                    file_ind,length(all_files),i,length(XTrainTranspose), fit)

            if (30 < abs(fit) && abs(fit) < 150)
                y = capstan_speed';
                y_hat = y_hat.OutputData;
                res = y - y_hat;
                [xcres,lags] = xcorr(res,res,'normalized');

                fig1 = figure(1); set(gcf, 'Position', [221 65 762 888])
                subplot(5,2,[1, 2]); compareplot(iddata_bfd_to_capstan_speed, sys_kd);
                ax = gca; ax.Legend.Location = 'southeast';

                subplot(5,2,[3, 4]); plot(res);
                subplot(5,2,[5, 6]); plot_fft(y, y_hat);

                subplot(5,2,7); plot(y,y_hat, '.'); hold on; 
%                 plot(linspace(min(min(y,y_hat)),max(max(y,y_hat)),length(y_hat)),linspace(min(min(y,y_hat)),max(max(y,y_hat)),length(y_hat)),'--');
                plot(linspace(min(y), max(y), length(y)),linspace(min(y), max(y), length(y)), 'LineWidth',2)
                hold off; xlabel('Actual'); ylabel('Predicted'); 

                subplot(5,2,8); histfit(res,50, 'logistic'); title('Residual Histogram')
                subplot(5,2,9); normplot(res);
                subplot(5,2,10); plot(lags, xcres); title('Residual Autocorrelation')
                
                fig2 = figure(2); set(gcf, 'Position', [1048 222 754 696])
                bode(sys_kd, {10^-8, 10})
                if (80 < abs(fit) && abs(fit) < 100) 
                    saveas(fig1, sprintf('%s\\%s\\%d,%d-1', curr_path, folder_name, file_ind, i),'png');
                    saveas(fig2, sprintf('%s\\%s\\%d,%d-2', curr_path, folder_name, file_ind, i),'png');  
                end
            end

        end
    end
end
disp('Done!')

%% PLOT OVERLAP BODE

folder_name = "_bode";

% cd 'C:\Users\georg\Dropbox (MIT)\fiber-draw'
cd(curr_path)

if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:length(all_files) 
    curr_file = all_files(file_ind);
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

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);
            
            dt = 1; 
%             a = 5; b = 5; k = 3;
            a = 4; b = 2; c = 2; k = 2;
            %     iddata_tension_to_power     = iddata(furnace_power', tension', dt);
            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);
%             sys_kd = oe(iddata_bfd_to_capstan_speed,  [a b k]);
            sys_kd = armax(iddata_bfd_to_capstan_speed,  [a b c k]);

            [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
            fprintf('file %d/%d\t subbatch %d/%d %2.4f\n', ...
                    file_ind,length(all_files),i,length(XTrainTranspose), fit)

            if (80 < abs(fit) && abs(fit) < 100)
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
%                 if (80 < abs(fit) && abs(fit) < 100) 
%                     saveas(fig2, sprintf('%s\\%s\\%d,%d-2', curr_path, folder_name, file_ind, i),'png');  
%                 end
            end

        end
    end
end
disp('Done!')

%% plotting poles (4th order, not really good fit)
figure(2); hold on;
for file_ind = 1:length(all_files)
    curr_file = all_files(file_ind);
    if ~curr_file.isdir

        % get batch info
        [BatchInfo, STRDEF] = stl_load_batchinfo(bXLSLoad, strDataPath, ...
            curr_file.name, nTower, bPlotAll, ...
            bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
            bPlot_each_preform_on_subplot, loBFD, hiBFD, subbatchMinLen, subbatchMaxLen, ...
            x_columns, y_columns);

        % turn it into a train / test array
        [XTrainTranspose, YTrainTranspose] = stl_prep_training_data(BatchInfo, ...
            STRDEF, x_columns, y_columns, fltLEN, PrefltLEN, bPlot, 0);

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);
            dt = 1;
            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);

            sys = ssest(iddata_bfd_to_capstan_speed, 4);
            figure(1); compare(iddata_bfd_to_capstan_speed, sys)
            figure(2); pzmap(sys); hold on;
            figure(3); bode(sys);
            sys.Report.Fit.FitPercent
        end

    end
end


%% helper function
function [f] = plot_fft(Y, Ypredict)
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