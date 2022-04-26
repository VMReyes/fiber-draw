%% Init
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

%% Grid Search Order (OE, kd)
max_order = 10;
plot_sysid_process = false;

diary on;
fit_matrix_bool = cell(length(all_files), 1);
fit_matrix = cell(length(all_files), 1);
folder_name = "_kd_oe_final";

if plot_sysid_process && exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        
        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        fit_matrix_bool{file_ind} = zeros(length(XTrainTranspose), max_order*max_order);
        fit_matrix{file_ind} = zeros(length(XTrainTranspose), max_order*max_order);

        cd(curr_path)

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            dt = 0.5;
            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd'-125,     dt);

            for b = 1:max_order for a = 1:max_order % for k = 1:4
                if ~(a<4 && b<4)
                    sys_kd = oe(iddata_bfd_to_capstan_speed,  [a b 3]);
%                         sys_kt = oe(iddata_tension_to_power,  [a b k]);

                    [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
%                         [y_hat, fit, x0] = compare(iddata_tension_to_power, sys_kt);
                    fprintf('file %d/%d\t subbatch %d/%d \t %d \t %d \t %2.4f\n', ...
                            file_ind,length(all_files),i,length(XTrainTranspose), a, b, fit)

                    if (30 < abs(fit) && abs(fit) < 150)
                        fit_matrix_bool{file_ind}(i, (a-1)*max_order+b) = 1;
                        fit_matrix{file_ind}(i, (a-1)*max_order+b) = abs(fit);

                        if plot_sysid_process
                            fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                            subplot(2,1,1); 
                            compare(iddata_bfd_to_capstan_speed, sys_kd);
%                                 compare(iddata_tension_to_power, sys_kt);
                            title(sprintf('%d, %d, (%d,%d)', file_ind,i,a,b));
                            ax = gca; ax.Legend.Location = 'southeast';
    
                            subplot(2,1,2);
                            plot_fft_comparison(capstan_speed, y_hat.OutputData)
    
                            cd(folder_name);
                            saveas(fig, sprintf('%d,%d,(%d,%d)',file_ind,i,a,b),'png');
                            cd ..;
                        end
                    end
                end
            end; end
        end
    end
end
disp('Done!')
diary off;
save('fit_results'+folder_name, 'fit_matrix', 'fit_matrix_bool');

%% Grid Search Order (ARMAX, kd)
max_order = 8;
plot_sysid_process = false;

diary on;
fit_matrix_bool = cell(length(all_files), 1);
fit_matrix = cell(length(all_files), 1);
folder_name = "_kd_armax_final";

if plot_sysid_process && exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        
        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        fit_matrix_bool{file_ind} = zeros(length(XTrainTranspose), max_order^3);
        fit_matrix{file_ind} = zeros(length(XTrainTranspose), max_order^3);

        cd(curr_path)

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            dt = 0.5;
%             iddata_tension_to_power     = iddata(furnace_power', tension', dt);
            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd'-125,     dt);

            for a = 1:max_order for b = 1:max_order for c = 1:max_order 
                if ~(a<4 && b<4 && c<4)
                    sys_kd = armax(iddata_bfd_to_capstan_speed,  [a b c 3]);
%                         sys_kt = oe(iddata_tension_to_power,  [a b k]);

                    [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
%                         [y_hat, fit, x0] = compare(iddata_tension_to_power, sys_kt);
                    fprintf('file %d/%d\t subbatch %d/%d \t %d \t %d \t %d \t %2.4f\n', ...
                            file_ind,length(all_files),i,length(XTrainTranspose), a, b, c, fit)

                    if (30 < abs(fit) && abs(fit) < 150)
                        fit_matrix_bool{file_ind}(i, (a-1)*max_order^2+(b-1)*max_order+c) = 1;
                        fit_matrix{file_ind}(i, (a-1)*max_order^2+(b-1)*max_order+c) = abs(fit);

                        if plot_sysid_process
                            fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                            subplot(2,1,1); 
                            compare(iddata_bfd_to_capstan_speed, sys_kd);
%                                 compare(iddata_tension_to_power, sys_kt);
                            title(sprintf('%d, %d, (%d,%d,%d)', file_ind,i,a,b,c));
                            ax = gca; ax.Legend.Location = 'southeast';
    
                            subplot(2,1,2);
                            plot_fft_comparison(capstan_speed, y_hat.OutputData)
    
                            cd(folder_name);
                            saveas(fig, sprintf('%d,%d,(%d,%d,%d)',file_ind,i,a,b,c),'png');
                            cd ..;
                        end
                    end
                end
            end; end; end
        end
    end
end
disp('Done!')
diary off;
save('fit_results'+folder_name, 'fit_matrix', 'fit_matrix_bool');

%% Grid Search Order (OE, kt)
max_order = 10;
plot_sysid_process = false;

diary on;
fit_matrix_bool = cell(length(all_files), 1);
fit_matrix = cell(length(all_files), 1);
folder_name = "_kt_oe_final";

if plot_sysid_process && exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        
        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        fit_matrix_bool{file_ind} = zeros(length(XTrainTranspose), max_order*max_order);
        fit_matrix{file_ind} = zeros(length(XTrainTranspose), max_order*max_order);

        cd(curr_path)
        
        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            dt = 0.5;
            iddata_tension_to_power     = iddata(furnace_power', tension'-median(tension), dt);
            for b = 1:max_order for a = 1:max_order for k = 1:4
                if ~(a<4 && b<4)
                    sys_kt = oe(iddata_tension_to_power,  [a b k]);

                    [y_hat, fit, x0] = compare(iddata_tension_to_power, sys_kt);
                    fprintf('file %d/%d\t subbatch %d/%d \t %d \t %d \t %d\t %2.4f\n', ...
                            file_ind,length(all_files),i,length(XTrainTranspose), a, b, k, fit)

                    if (30 < abs(fit) && abs(fit) < 150)
                        fit_matrix_bool{file_ind}(i, (a-1)*max_order*4+(b-1)*4+k) = 1;
                        fit_matrix{file_ind}(i, (a-1)*max_order*4+(b-1)*4+k) = abs(fit);

                        if plot_sysid_process
                            fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                            subplot(2,1,1); 
                            compare(iddata_tension_to_power, sys_kt);
                            title(sprintf('%d, %d, (%d,%d,%d)', file_ind,i,a,b,k));
                            ax = gca; ax.Legend.Location = 'southeast';
    
                            subplot(2,1,2);
                            plot_fft_comparison(furnace_power, y_hat.OutputData)
    
                            cd(folder_name);
                            saveas(fig, sprintf('%d,%d,(%d,%d,%d)',file_ind,i,a,b,k),'png');
                            cd ..;
                        end
                    end
                end
            end; end; end
            end
    end
end
disp('Done!')
diary off;
save('fit_results'+folder_name, 'fit_matrix', 'fit_matrix_bool');

%% Grid Search Order (ARMAX, kt)
max_order = 8;
plot_sysid_process = false;
diary on;
fit_matrix_bool = cell(length(all_files), 1);
fit_matrix = cell(length(all_files), 1);
folder_name = "_kt_armax_final";

if plot_sysid_process && exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        
        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        fit_matrix_bool{file_ind} = zeros(length(XTrainTranspose), max_order^3);
        fit_matrix{file_ind} = zeros(length(XTrainTranspose), max_order^3);

        cd(curr_path)

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            dt = 0.5;
            iddata_tension_to_power     = iddata(furnace_power', tension'-median(tension), dt);
            for a = 1:max_order for b = 1:max_order for c = 1:max_order for k = 1:4
                if ~(a<3 && b<3 && c<3)
                    sys_kt = armax(iddata_tension_to_power,  [a b c k]);

                    [y_hat, fit, x0] = compare(iddata_tension_to_power, sys_kt);
                    fprintf('file %d/%d\t subbatch %d/%d \t %d \t %d \t %d \t %d \t %2.4f\n', ...
                            file_ind,length(all_files),i,length(XTrainTranspose), a, b, c, k, fit)

                    if (30 < abs(fit) && abs(fit) < 150)
                        fit_matrix_bool{file_ind}(i, (a-1)*max_order^2+(b-1)*max_order+c) = 1;
                        fit_matrix{file_ind}(i, (a-1)*max_order^2+(b-1)*max_order+c) = abs(fit);

                        if plot_sysid_process
                            fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                            subplot(2,1,1); 
                            compare(iddata_tension_to_power, sys_kt);
                            title(sprintf('%d, %d, (%d,%d,%d)', file_ind,i,a,b,c));
                            ax = gca; ax.Legend.Location = 'southeast';
    
                            subplot(2,1,2);
                            plot_fft_comparison(furnace_power, y_hat.OutputData)
    
                            cd(folder_name);
                            saveas(fig, sprintf('%d,%d,(%d,%d,%d,%d)',file_ind,i,a,b,c,k),'png');
                            cd ..;
                        end
                    end
                end
            end; end; end; end
        end
    end
end
disp('Done!')
diary off;
save('fit_results'+folder_name, 'fit_matrix', 'fit_matrix_bool');

%% Merge Models (OE, kd)

selected_order = [1 4 3];

all_sys = {};
all_iddata = {};
dt = 0.5;
plot_sysid_process = false;
folder_name = "_kd_oe_final_select";
if plot_sysid_process && exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        
        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        cd(curr_path)

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd'-125,     dt);

            sys_kd = oe(iddata_bfd_to_capstan_speed,  selected_order);

            [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);

            fprintf('file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_files),i,length(XTrainTranspose), fit)

            if (30 < abs(fit) && abs(fit) < 150)
                if isempty(all_iddata)
                    all_iddata = {iddata_bfd_to_capstan_speed};
                else
                    all_iddata = horzcat(all_iddata, {iddata_bfd_to_capstan_speed});
                end

                if isempty(all_sys)
                    all_sys = {sys_kd};
                else
                    all_sys = horzcat(all_sys, {sys_kd});
                end

                if plot_sysid_process
                    fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                    subplot(2,1,1); 
                    compare(iddata_bfd_to_capstan_speed, sys_kd);
                    title(sprintf('File %d, Subbatch %d', file_ind,i));
                    ax = gca; ax.Legend.Location = 'southeast';

                    subplot(2,1,2);
                    plot_fft_comparison(capstan_speed, y_hat.OutputData)

                    cd(folder_name);
                    saveas(fig, sprintf('%d,%d',file_ind,i),'png');
                    cd ..;
                end
            end
        end
    end
end
disp('Done Iterating!');
[m, tv] = merge(all_sys{:});
disp('Done with merge(sys)!')
sys_kd_total = oe(merge(all_iddata{:}), selected_order);
disp('Done with oe(merge(), [])!')
save("sys_kd_total_oe.mat", "all_sys", "all_iddata", "sys_kd_total", "m", "tv");


%% Merge Models (ARMAX, kd)

selected_order = [4 3 1 3];

all_sys = {};
all_iddata = {};
dt = 0.5;
plot_sysid_process = false;
folder_name = "_kd_armax_final_select";
if plot_sysid_process && exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        
        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        cd(curr_path)

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd'-125,     dt);

            sys_kd = armax(iddata_bfd_to_capstan_speed,  selected_order);

            [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);

            fprintf('file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_files),i,length(XTrainTranspose), fit)

            if (30 < abs(fit) && abs(fit) < 150)
                if isempty(all_iddata)
                    all_iddata = {iddata_bfd_to_capstan_speed};
                else
                    all_iddata = horzcat(all_iddata, {iddata_bfd_to_capstan_speed});
                end

                if isempty(all_sys)
                    all_sys = {sys_kd};
                else
                    all_sys = horzcat(all_sys, {sys_kd});
                end

                if plot_sysid_process
                    fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                    subplot(2,1,1); 
                    compare(iddata_bfd_to_capstan_speed, sys_kd);
                    title(sprintf('File %d, Subbatch %d', file_ind,i));
                    ax = gca; ax.Legend.Location = 'southeast';

                    subplot(2,1,2);
                    plot_fft_comparison(capstan_speed, y_hat.OutputData)

                    cd(folder_name);
                    saveas(fig, sprintf('%d,%d',file_ind,i),'png');
                    cd ..;
                end
            end
        end
    end
end
disp('Done Iterating!');
[m, tv] = merge(all_sys{:});
disp('Done with merge(sys)!')
sys_kd_total = armax(merge(all_iddata{:}), selected_order);
disp('Done with oe(merge(), [])!')
save("sys_kd_total_armax.mat", "all_sys", "all_iddata", "sys_kd_total", "m", "tv");

%% Merge Models (OE, kt)

selected_order = [4 5 4];

all_sys = {};
all_iddata = {};
dt = 0.5;
plot_sysid_process = false;
folder_name = "_kt_oe_final_select";
if plot_sysid_process && exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        
        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        cd(curr_path)

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            iddata_tension_to_power     = iddata(furnace_power', tension'-median(tension), dt);
            
            sys_kd = oe(iddata_tension_to_power,  selected_order);

            [y_hat, fit, x0] = compare(iddata_tension_to_power, sys_kd);

            fprintf('file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_files),i,length(XTrainTranspose), fit)

            if (30 < abs(fit) && abs(fit) < 150)
                if isempty(all_iddata)
                    all_iddata = {iddata_tension_to_power};
                else
                    all_iddata = horzcat(all_iddata, {iddata_tension_to_power});
                end

                if isempty(all_sys)
                    all_sys = {sys_kd};
                else
                    all_sys = horzcat(all_sys, {sys_kd});
                end

                if plot_sysid_process
                    fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                    subplot(2,1,1); 
                    compare(iddata_tension_to_power, sys_kd);
                    title(sprintf('File %d, Subbatch %d', file_ind,i));
                    ax = gca; ax.Legend.Location = 'southeast';

                    subplot(2,1,2);
                    plot_fft_comparison(capstan_speed, y_hat.OutputData)

                    cd(folder_name);
                    saveas(fig, sprintf('%d,%d',file_ind,i),'png');
                    cd ..;
                end
            end
        end
    end
end
disp('Done Iterating!');
[m, tv] = merge(all_sys{:});
disp('Done with merge(sys)!')
sys_kt_total = oe(merge(all_iddata{:}), selected_order);
disp('Done with oe(merge(), [])!')
save("sys_kt_total_oe.mat", "all_sys", "all_iddata", "sys_kt_total", "m", "tv");

%% Merge Models (ARMAX, kt)

selected_order = [3 2 2 3];

all_sys = {};
all_iddata = {};
dt = 0.5;
plot_sysid_process = false;
folder_name = "_kt_oe_final_select";
if plot_sysid_process && exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for file_ind = 1:16 % 1:16 for tower 48, 18:length(all_files) for tower 51
    curr_file = all_files(file_ind);
    if ~curr_file.isdir
        
        % load from loaded data
        XTrainTranspose = all_file_data{file_ind,1};
        YTrainTranspose = all_file_data{file_ind,2};

        cd(curr_path)

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            iddata_tension_to_power     = iddata(furnace_power', tension'-median(tension), dt);
            sys_kt = armax(iddata_tension_to_power,  selected_order);
            [y_hat, fit, x0] = compare(iddata_tension_to_power, sys_kt);

            fprintf('file %d/%d\t subbatch %d/%d \t %2.4f\n', ...
                    file_ind,length(all_files),i,length(XTrainTranspose), fit)

            if (30 < abs(fit) && abs(fit) < 150)
                if isempty(all_iddata)
                    all_iddata = {iddata_tension_to_power};
                else
                    all_iddata = horzcat(all_iddata, {iddata_tension_to_power});
                end

                if isempty(all_sys)
                    all_sys = {sys_kt};
                else
                    all_sys = horzcat(all_sys, {sys_kt});
                end

                if plot_sysid_process
                    fig = figure(1); set(gcf, 'Position', [229 222 754 696])
                    subplot(2,1,1); 
                    compare(iddata_tension_to_power, sys_kt);
                    title(sprintf('File %d, Subbatch %d', file_ind,i));
                    ax = gca; ax.Legend.Location = 'southeast';

                    subplot(2,1,2);
                    plot_fft_comparison(capstan_speed, y_hat.OutputData)

                    cd(folder_name);
                    saveas(fig, sprintf('%d,%d',file_ind,i),'png');
                    cd ..;
                end
            end
        end
    end
end
disp('Done Iterating!');
[m, tv] = merge(all_sys{:});
disp('Done with merge(sys)!')
sys_kt_total = armax(merge(all_iddata{:}), selected_order);
disp('Done with oe(merge(), [])!')
save("sys_kt_total_armax_5.mat", "all_sys", "all_iddata", "sys_kt_total", "m", "tv");

%% plot noise variances
% load('sys_kd_total_oe.mat');
% load('sys_kd_total_armax.mat');
load('sys_kt_total_oe.mat')

all_vars = zeros(1, length(all_sys));
for i = 1:length(all_sys)
    all_vars(i) = all_sys{i}.NoiseVariance;
end
figure(2); 
% kd_oe
% histogram(all_vars, 100, 'BinLimits', [0 5e3]) 
% title('Noise Variances of OE Models for $K_d$ Controller')

% kd_armax
% histogram(all_vars, 100);
% title('Noise Variances of ARMAX Models for $K_d$ Controller')

% kt_oe
histogram(all_vars, 100, 'BinLimits', [0 2]);
title('Noise Variances of OE Models for $K_t$ Controller')
latexify_plot

%% merged models misc
clc;

% load sys_kd_total_oe.mat
% % tv = 0.0121


load sys_kd_total_armax.mat
% % tv = 499.5298

% load sys_kt_total_oe.mat
% % tv = 0.0122

% m
% [A,B,C,D,F,dA,dB,dC,dD,dF] = polydata(m)
% [PVEC, DPVEC]  = getpvec(m)
% covariance = getcov(m)
% 
% sys_kd_total
% [A,B,C,D,F,dA,dB,dC,dD,dF] = polydata(sys_kd_total)
% [PVEC, DPVEC]  = getpvec(sys_kd_total)
% covariance = getcov(sys_kd_total)

% sys_kd_total.Report
% sys_kd_total.Report.Fit
mean(abs(sys_kd_total.Report.Fit.FitPercent))
mean(sys_kd_total.Report.Fit.MSE)

%% analyze - convert cell to matrix

fit_matrix_all = [];
fit_matrix_bool_all = [];

for i = 1:length(fit_matrix)
    fit_matrix_all = vertcat(fit_matrix_all,fit_matrix{i});
    fit_matrix_bool_all = vertcat(fit_matrix_bool_all,fit_matrix_bool{i});
end

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

