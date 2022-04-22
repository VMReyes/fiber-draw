%% Main
clear; clc; close all;

[num_good_subbatches, average_accuracy] = grid_search_oe('kd', 10, false); 
writematrix([num_good_subbatches; average_accuracy], 'kd_oe_results.xlsx');

[num_good_subbatches, average_accuracy] = grid_search_oe('kt', 10, false);
writematrix([num_good_subbatches; average_accuracy], 'kt_oe_results.xlsx');

[num_good_subbatches, average_accuracy] = grid_search_armax('kd', 8, false); 
writematrix([num_good_subbatches; average_accuracy], 'kd_armax_results.xlsx');

[num_good_subbatches, average_accuracy] = grid_search_armax('kt', 8, false);
writematrix([num_good_subbatches; average_accuracy], 'kt_armax_results.xlsx');

[all_sys, all_iddata] = merge_models('kd', 'oe', [1 4 3], false);
[all_sys, all_iddata] = merge_models('kd', 'armax', [4 3 1 3], false);
[all_sys, all_iddata] = merge_models('kt', 'oe', [4 5 4], false);
[all_sys, all_iddata] = merge_models('kt', 'armax', [3 2 2 3], false);

%% Function Definitions
function [all_file_data, dt] = load_data()
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
    x_columns = ["cpspdactval", "frnpwrmv", "pfspdactval"];
    y_columns = ["barefibrediadisplay", "tenncmv"];
    
    % Subbatch Parameters
    fltLEN = 21;
    bPlot = 0; 
    PrefltLEN = 1;
    dt = 0.5;
    
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
        load("all_file_data.mat", 'all_file_data');
    end
    
    disp('Done Loading!')
end

function [num_good_subbatches, average_accuracy] = grid_search_oe(controller, max_order, plot_sysid_process)

    [all_file_data, dt] = load_data();
    
    switch controller
        case 'kd'
            all_nk = [3];
        case 'kt'
            all_nk = 1:4;
        otherwise
            error("invalid controller argument!")
    end
    
    fit_matrix = cell(length(all_file_data), 1);
    folder_name = sprintf("oe_%s_plots", controller);
    
    if plot_sysid_process && exist(folder_name, 'dir') ~= 7
        mkdir(folder_name);
    end
    
    for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
        curr_file = all_file_data{file_ind,1};
        if ~isempty(curr_file)
    
            % load from loaded data
            XTrainTranspose = all_file_data{file_ind,1};
            YTrainTranspose = all_file_data{file_ind,2};
    
            fit_matrix{file_ind} = zeros(length(XTrainTranspose), max_order^2*length(all_nk));
    
            for i = 1:length(XTrainTranspose)
                xs = cell2mat(XTrainTranspose(i));
                ys = cell2mat(YTrainTranspose(i));
                capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
                bfd = ys(1,:); tension = ys(2,:);
    
                switch controller
                    case 'kd'
                        outputs = capstan_speed';
                        inputs = bfd'-125;
    
                    case 'kt'
                        outputs = furnace_power';
                        inputs = tension'-median(tension);
                end
    
                iddata_subbatch = iddata(outputs, inputs, dt);
    
                for a = 1:max_order for b = 1:max_order for k = all_nk
                    if ~(a<4 && b<4)
                        sys_subbatch = oe(iddata_subbatch,  [a b k]);
    
                        [y_hat, fit, ~] = compare(iddata_subbatch, sys_subbatch);
    
                        fprintf('file %d/%d\t subbatch %d/%d \t %d \t %d \t %d \t %2.4f\n', ...
                            file_ind,length(all_file_data),i,length(XTrainTranspose), a, b, k, fit)
    
                        if (30 < abs(fit) && abs(fit) < 150)
                            switch controller
                                case 'kd'
                                    fit_matrix{file_ind}(i, (a-1)*max_order+b) = abs(fit);
                                case 'kt'
                                    fit_matrix{file_ind}(i, (a-1)*max_order*4+(b-1)*4+k) = abs(fit);
                            end
    
                            if plot_sysid_process
                                fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                                subplot(2,1,1);
                                compare(iddata_subbatch, sys_subbatch);
                                title(sprintf('%d, %d, (%d,%d,%d)', file_ind,i,a,b,k));
                                ax = gca; ax.Legend.Location = 'southeast';
    
                                subplot(2,1,2);
                                plot_fft_comparison(outputs', y_hat.OutputData)
                                latexify_plot;
                                saveas(fig, sprintf('%s\\%d,%d,(%d,%d,%d)',folder_name,file_ind,i,a,b,k),'png');
                            end
                        end
                    end
                end; end; end
            end
        end
    end
    
    fit_matrix_all = [];
    
    for i = 1:length(fit_matrix)
        fit_matrix_all = vertcat(fit_matrix_all,fit_matrix{i});
    end
    
    num_good_subbatches = zeros(1,size(fit_matrix_all,2));
    average_accuracy = zeros(1,size(fit_matrix_all,2));
    
    for i = 1:size(fit_matrix_all,2)
        col = fit_matrix_all(:,i);
        num_good_subbatches(i) = sum(col > 30 & col <= 100);
        average_accuracy(i) = mean(col(col > 30 & col <= 100));
    end
end

function [num_good_subbatches, average_accuracy] = grid_search_armax(controller, max_order, plot_sysid_process)

    [all_file_data, dt] = load_data();
    
    switch controller
        case 'kd'
            all_nk = [3];
        case 'kt'
            all_nk = 1:4;
        otherwise
            error("invalid controller argument!")
    end
    
    fit_matrix = cell(length(all_file_data), 1);
    folder_name = sprintf("armax_%s_plots", controller);
    
    if plot_sysid_process && exist(folder_name, 'dir') ~= 7
        mkdir(folder_name);
    end
    
    for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
        curr_file = all_file_data{file_ind,1};
        if ~isempty(curr_file)
    
            % load from loaded data
            XTrainTranspose = all_file_data{file_ind,1};
            YTrainTranspose = all_file_data{file_ind,2};
    
            fit_matrix{file_ind} = zeros(length(XTrainTranspose), max_order^3*length(all_nk));
    
            for i = 1:length(XTrainTranspose)
                xs = cell2mat(XTrainTranspose(i));
                ys = cell2mat(YTrainTranspose(i));
                capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
                bfd = ys(1,:); tension = ys(2,:);
    
                switch controller
                    case 'kd'
                        outputs = capstan_speed';
                        inputs = bfd'-125;
    
                    case 'kt'
                        outputs = furnace_power';
                        inputs = tension'-median(tension);
                end
    
                iddata_subbatch = iddata(outputs, inputs, dt);
    
                for a = 1:max_order for b = 1:max_order for c = 1:max_order for k = all_nk
                    if ~(a<4 && b<4 && c<4)
                        sys_subbatch = armax(iddata_subbatch,  [a b c k]);
    
                        [y_hat, fit, ~] = compare(iddata_subbatch, sys_subbatch);
                        fprintf('file %d/%d\t subbatch %d/%d \t %d \t %d \t %d \t %d \t %2.4f\n', ...
                            file_ind,length(all_file_data),i,length(XTrainTranspose), a, b, c, k, fit)
    
                        if (30 < abs(fit) && abs(fit) < 150)
                            switch controller
                                case 'kd'
                                    fit_matrix{file_ind}(i, (a-1)*max_order^2+(b-1)*max_order+c) = abs(fit);
                                case 'kt'
                                    fit_matrix{file_ind}(i, (a-1)*max_order^3+(b-1)*max_order^2+(c-1)*max_order+k) = abs(fit);
                            end
    
    
                            if plot_sysid_process
                                fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                                subplot(2,1,1);
                                compare(iddata_subbatch, sys_subbatch);
                                title(sprintf('%d, %d, (%d,%d,%d)', file_ind,i,a,b,c));
                                ax = gca; ax.Legend.Location = 'southeast';
    
                                subplot(2,1,2);
                                plot_fft_comparison(outputs', y_hat.OutputData)
                                latexify_plot;
                                saveas(fig, sprintf('%s\\%d,%d,(%d,%d,%d)',folder_name,file_ind,i,a,b,k),'png');
                            end
                        end
                    end
                end; end; end; end
            end
        end
    end
    
    fit_matrix_all = [];
    
    for i = 1:length(fit_matrix)
        fit_matrix_all = vertcat(fit_matrix_all,fit_matrix{i});
    end
    
    num_good_subbatches = zeros(1,size(fit_matrix_all,2));
    average_accuracy = zeros(1,size(fit_matrix_all,2));
    
    for i = 1:size(fit_matrix_all,2)
        col = fit_matrix_all(:,i);
        num_good_subbatches(i) = sum(col > 30 & col <= 100);
        average_accuracy(i) = mean(col(col > 30 & col <= 100));
    end

end


function [all_sys, all_iddata] = merge_models(controller, model_struc, selected_orders, plot_sysid_process)

    [all_file_data, dt] = load_data();
    
    if ~any(strcmp({'kd','kt'}, controller))
        error('invalid controller argument!')
    end
    
    if ~any(strcmp({'armax','oe'}, model_struc))
        error('invalid model structure argument!')
    end
    
    if ~((strcmpi(model_struc, 'armax') && length(selected_orders) == 4) || (strcmpi(model_struc, 'oe') && length(selected_orders) == 3))
        error('invalid specified orders!')
    end
    
    all_sys = {};
    all_iddata = {};
    folder_name = sprintf("%s_%s_merge", controller, model_struc);
    if plot_sysid_process && exist(folder_name, 'dir') ~= 7
        mkdir(folder_name);
    end
    
    for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
        curr_file = all_file_data(file_ind);
        if ~isempty(curr_file)
    
            % load from loaded data
            XTrainTranspose = all_file_data{file_ind,1};
            YTrainTranspose = all_file_data{file_ind,2};
    
            for i = 1:length(XTrainTranspose)
                xs = cell2mat(XTrainTranspose(i));
                ys = cell2mat(YTrainTranspose(i));
                capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
                bfd = ys(1,:); tension = ys(2,:);
    
                switch controller
                    case 'kd'
                        outputs = capstan_speed';
                        inputs = bfd'-125;
    
                    case 'kt'
                        outputs = furnace_power';
                        inputs = tension'-median(tension);
                end
    
                iddata_subbatch = iddata(outputs, inputs, dt);
    
                switch model_struc
                    case 'oe'
                        sys_subbatch = oe(iddata_subbatch,  selected_orders);
                    case 'armax'
                        sys_subbatch = armax(iddata_subbatch,  selected_orders);
                end
    
                [y_hat, fit, ~] = compare(iddata_subbatch, sys_subbatch);
    
                fprintf('file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_file_data),i,length(XTrainTranspose), fit)
    
                if (30 < abs(fit) && abs(fit) < 150)
                    if isempty(all_iddata)
                        all_iddata = {iddata_subbatch};
                    else
                        all_iddata = horzcat(all_iddata, {iddata_subbatch});
                    end
    
                    if isempty(all_sys)
                        all_sys = {sys_subbatch};
                    else
                        all_sys = horzcat(all_sys, {sys_subbatch});
                    end
    
                    if plot_sysid_process
                        fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                        subplot(2,1,1);
                        compare(iddata_subbatch, sys_subbatch);
                        title(sprintf('File %d, Subbatch %d', file_ind,i));
                        ax = gca; ax.Legend.Location = 'southeast';
    
                        subplot(2,1,2);
                        plot_fft_comparison(outputs', y_hat.OutputData)
                        latexify_plot;
                        saveas(fig, sprintf('%s\\%d,%d',folder_name,file_ind,i),'png');
                    end
                end
            end
        end
    end
    disp('Done Iterating!');
    
    [m, tv] = merge(all_sys{:});
    disp('Done with merge(sys)!')
    
    switch model_struc
        case 'oe'
            sys_total = oe(merge(all_iddata{:}), selected_orders);
            disp('Done with oe(merge(), [])!')
        case 'armax'
            sys_total = armax(merge(all_iddata{:}), selected_orders);
            disp('Done with armax(merge(), [])!')
    end
    
    save(sprintf('sys_%s_total_%s.mat', controller, model_struc), "all_sys", "all_iddata", "sys_total", "m", "tv");
end

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
    legend('Data', 'Prediction')
    axis([0 1 -.1  max([   max(log10(abs(FY_s).^2))   max(log10(abs(FYpredict_s).^2))   ])   ])
    title('Power Spectrum');
end