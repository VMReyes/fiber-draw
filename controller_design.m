clc; clear; close all;
%% Init
% parameters
strDataPath         = 'C:\Users\georg\Dropbox (MIT)\minigroup_mit_sterlite\from Sterlite\data\MIT_DrawData_48and51\';
strDataFilename     = 'DrawData_Tower48_2020-12-01_to2020-12-08.csv';

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

% get batch info
[BatchInfo, STRDEF] = stl_load_batchinfo(bXLSLoad, strDataPath, ...
    strDataFilename, nTower, bPlotAll, ...
    bPlot_each_preform_on_subplot_with_inrangesubbatches_, ...
    bPlot_each_preform_on_subplot, loBFD, hiBFD, subbatchMinLen, subbatchMaxLen, ...
    x_columns, y_columns);

% turn it into a train / test array
[XTrainTranspose, YTrainTranspose] = stl_prep_training_data(BatchInfo, ...
    STRDEF, x_columns, y_columns, fltLEN, PrefltLEN, bPlot, 0);


%% visualize data

cd 'C:\Users\georg\Desktop\MS Research\fiber-draw'
fig = figure; 

folder_name = '_viz';
if exist(folder_name) == 7
    rmdir("C:\Users\georg\Desktop\MS Research\fiber-draw\"+folder_name, 's');
end
mkdir(folder_name);

for i = 1:length(XTrainTranspose)
    xs = cell2mat(XTrainTranspose(i));
    ys = cell2mat(YTrainTranspose(i));
    capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
    bfd = ys(1,:); tension = ys(2,:);
    
    subplot(2,2,1); plot(tension); title 'Tension'; ylim([80 160]);
    subplot(2,2,2); plot(furnace_power); title 'Furnace Power'; ylim([48 70]);
    subplot(2,2,3); plot(bfd); title 'BFD'; ylim([124 126]);
    subplot(2,2,4); plot(capstan_speed); title 'Capstan Speed'; ylim([1600 2800])
    if i == 1
        pause(10) % drag the window
    end
    cd(folder_name);
    saveas(fig, num2str(i), 'png');
    cd ..;
end

%% Choose model / order (BY SINGLE SUBBATCH)

cd 'C:\Users\georg\Desktop\MS Research\fiber-draw'

fig = figure;
set(gcf, 'Position',[133.8000 95.4000 1.2952e+03 658.4000])

if exist("_oe") == 7
    rmdir("C:\Users\georg\Desktop\MS Research\fiber-draw\_oe", 's');
end
mkdir('_oe');

max_order = 8;
num_params = 2;
fit_matrix_bool = zeros(length(XTrainTranspose), max_order^num_params);
fit_matrix = zeros(length(XTrainTranspose), max_order^num_params);

for i = 1:length(XTrainTranspose)
    xs = cell2mat(XTrainTranspose(i));
    ys = cell2mat(YTrainTranspose(i));
    capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
    bfd = ys(1,:); tension = ys(2,:);

%     sys_kt = armax(iddata_tension_to_power,  [2 2 2 1]);
%     sys_kt = bj(iddata_tension_to_power,  [2 2 2 2 1]);
%     sys_kt = oe(iddata_tension_to_power,  [2 2 1]);

    dt = 1;
%     iddata_tension_to_power     = iddata(furnace_power', tension', dt);
    iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);
    
    for a=2:max_order for b = 3:max_order %for c = 1:4 for d = 1:4
%     sys_kt = armax(iddata_tension_to_power,  [a b c 1]);
%     sys_kt = bj(iddata_tension_to_power,  [a b c d 1]);
    sys_kd = oe(iddata_bfd_to_capstan_speed,  [a b 1]);   
    [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
%     [a b c d fit]
    disp([i a b fit])

    if (30 < abs(fit) && abs(fit) < 200)
        fit_matrix_bool(i, (a-1)*max_order+b) = 1;
        fit_matrix(i, (a-1)*max_order+b) = abs(fit);
        subplot(2,1,1); compare(iddata_bfd_to_capstan_speed, sys_kd);
        title(sprintf('%d, %d%d', i,a,b)); 
        ax = gca; ax.Legend.Location = 'southeast';
        
        subplot(2,1,2);
        plot_fft(capstan_speed, capstan_speed, y_hat.OutputData, y_hat.OutputData)

        cd('_oe');
%         saveas(fig, sprintf('%d%d%d%d',a,b,c,d),'png');
        saveas(fig, sprintf('%d,%d%d',i,a,b),'png');
        cd ..;
    end
    end; end %; end; end

end
disp('Done!')

%% analysis (COEFFS BY SUBBATCH MODALITY)
as = zeros(1,length(XTrainTranspose));
bs = zeros(5,length(XTrainTranspose));
cs = zeros(1,length(XTrainTranspose));
ds = zeros(1,length(XTrainTranspose));
fs = zeros(2,length(XTrainTranspose));

for i = 1:length(XTrainTranspose)
    xs = cell2mat(XTrainTranspose(i)); 
    ys = cell2mat(YTrainTranspose(i)); 
    capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
    bfd = ys(1,:); tension = ys(2,:);
    dt = 1;
    iddata_tension_to_power     = iddata(furnace_power', tension', dt);
    iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);
    
    % replace for different graph
    sys_kd = oe(iddata_bfd_to_capstan_speed,  [4 1 1]); 
    [sys_kd.A sys_kd.B sys_kd.C sys_kd.D sys_kd.F];
    as(i) = sys_kd.A;
    bs(:,i) = sys_kd.B;
    cs(i) = sys_kd.C;
    ds(i) = sys_kd.D;
    fs(:,i) = sys_kd.F;

%     [y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed, sys_kd);
%     figure(1); plot(capstan_speed); hold on; plot(y_hat, 'r'); hold off;
%     stl_plottrainingresults_FFTPSD_function(capstan_speed, capstan_speed, y_hat.OutputData, y_hat.OutputData);

%     figure(1); plot(furnace_power); hold on; plot(y_hat, 'r'); hold off;
%     stl_plottrainingresults_FFTPSD_function(furnace_power, furnace_power, y_hat.OutputData, y_hat.OutputData);
%     pause(1)
end

figure; histogram(fs,100); title('F coefficient distribution'); ylabel('Frequency')
for i = 2:4
    figure; histogram(bs(i,:),100);
    title(sprintf('B(%d) coefficient distribution',i));
    ylabel('Frequency');
end
%% Choose model (multi-experiment) 

% cd 'C:\Users\georg\Desktop\MS Research\fiber-draw'
% experiments = cell(1,length(XTrainTranspose));
% dt = 1;
% 
% for i = 1:length(XTrainTranspose)
%     xs = cell2mat(XTrainTranspose(i));
%     ys = cell2mat(YTrainTranspose(i));
%     capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
%     bfd = ys(1,:); tension = ys(2,:);
% 
%     iddata_tension_to_power     = iddata(furnace_power', tension', dt);
%     iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);
%     experiments{i} = iddata_bfd_to_capstan_speed;
%     if i == 1
%         data = iddata_bfd_to_capstan_speed;
%     else
%         data = merge(data, iddata_bfd_to_capstan_speed);
%     end
% end
% 
% % figure; plot(u_concat)
% % figure; plot(y_concat)
% 
% sys_kd = oe(data, [4 1 1])
% % compare(experiments{16}, sys_kd)
% disp('Done!')

%% Choose model (concat data) PART 1

cd 'C:\Users\georg\Desktop\MS Research\fiber-draw'
experiments = cell(1,length(XTrainTranspose));
dt = 1;

for i = 1:length(XTrainTranspose)
    xs = cell2mat(XTrainTranspose(i));
    ys = cell2mat(YTrainTranspose(i));
    capstan_speed = xs(1,:); furnace_power = xs(2,:); preform_speed = xs(3,:);
    bfd = ys(1,:); tension = ys(2,:);

    iddata_tension_to_power     = iddata(furnace_power', tension', dt);
    iddata_bfd_to_capstan_speed = iddata(capstan_speed', bfd',     dt);
    experiments{i} = iddata_bfd_to_capstan_speed;
    if i == 1
        data = iddata_bfd_to_capstan_speed;
        u_concat = bfd;
        y_concat = capstan_speed;
    else
        data = merge(data, iddata_bfd_to_capstan_speed);
        u_concat = horzcat(u_concat, bfd);
        y_concat = horzcat(y_concat, capstan_speed);
    end
end

figure; plot(u_concat)
figure; plot(y_concat)

iddata_bfd_to_capstan_speed_total = iddata(y_concat', u_concat', dt);

%% part 2
cd 'C:\Users\georg\Desktop\MS Research\fiber-draw'
fig = figure; 

folder_name = '_fit';
if exist(folder_name) == 7
    rmdir("C:\Users\georg\Desktop\MS Research\fiber-draw\"+folder_name, 's');
end
mkdir(folder_name);
for a = 1:4 for b = 1:4 % for c = 1:4
% if (a == 2 && b > 3)
%     continue;
% end

% a = 4;
% b = 1;
% oeopt = oeOptions;
% oeopt.InitialCondition = 'estimate'; % 'backcast';
% oeopt.SearchOptions.MaxIterations = 20;
% oeopt.Display = 'off';
sys_kd = oe(iddata_bfd_to_capstan_speed_total, [a b 1]) % 4 1 1 
% sys_kd = armax(iddata_bfd_to_capstan_speed_total,  [a b c 1]);
[y_hat, fit, x0] = compare(iddata_bfd_to_capstan_speed_total, sys_kd);
figure(1); plot(y_concat); hold on; plot(y_hat, 'r'); hold off;
stl_plottrainingresults_FFTPSD_function(y_concat, y_concat, y_hat.OutputData, y_hat.OutputData);
[a b fit]
% [a b c fit]

figure(2); histogram(y_concat); hold on; histogram(y_hat.OutputData); hold off;
xlim([2000 2600])

cd(folder_name);
saveas(fig, num2str(i), 'png');
cd ..;
end; end %; end
disp('Done!')


%% helper function
function [f] = plot_fft(Y, Y_filt, Ypredict, Ypredict_filt)
    %take the FFT
    FY = fft(Y);
    FY_s = fftshift(FY);
    FY_filt = fft(Y_filt);
    %take the FFT
    FYpredict = fft(Ypredict);
    FYpredict_s = fftshift(FYpredict);
    FYpredict_filt = fft(Ypredict_filt);

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