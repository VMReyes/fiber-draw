clear; clc; close all;
%% Init
% parameters
strDataPath = 'C:\Users\georg\Dropbox (MIT)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
all_files   = dir(strDataPath);
% curr_path   = 'C:\Users\georg\Desktop\fiber-draw';
curr_path   = 'D:\GEORGE\fiber-draw';

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
end

disp('Done Loading!')

%% Grid Search Order (OE)
max_order = 10;
plot_sysid_process = false;

diary on;
fit_matrix_bool = cell(length(all_files), 1);
fit_matrix = cell(length(all_files), 1);
folder_name = "_kd_oe_final";

if exist(folder_name, 'dir') ~= 7
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

        fig = figure(1);
        set(gcf, 'Position', [229 222 754 696])

        for i = 1:length(XTrainTranspose)
            xs = cell2mat(XTrainTranspose(i));
            ys = cell2mat(YTrainTranspose(i));
            capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
            bfd = ys(1,:); tension = ys(2,:);

            dt = 0.5;
%             iddata_tension_to_power     = iddata(furnace_power', tension', dt);
            iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);

            for b = 1:max_order for a = 1:b % for k = 1:4
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
                            subplot(2,1,1); 
                            compare(iddata_bfd_to_capstan_speed, sys_kd);
%                                 compare(iddata_tension_to_power, sys_kt);
                            title(sprintf('%d, %d, (%d,%d)', file_ind,i,a,b));
                            ax = gca; ax.Legend.Location = 'southeast';
    
                            subplot(2,1,2);
                            plot_fft(capstan_speed, capstan_speed, y_hat.OutputData, y_hat.OutputData)
    
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

%% analyze - convert cell to matrix

fit_matrix_all = [];
fit_matrix_bool_all = [];

for i = 1:length(fit_matrix)
    fit_matrix_all = vertcat(fit_matrix_all,fit_matrix{i});
    fit_matrix_bool_all = vertcat(fit_matrix_bool_all,fit_matrix_bool{i});
end