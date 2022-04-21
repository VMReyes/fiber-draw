clc; clear; close all; 

load("run_results\architecture_experiment.mat", "deep_lstm");
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat", "x_test", "y_test");
load("sys_kd_total_armax.mat", "sys_kd_total")
% load("all_file_data.mat") % TODO: may not run

folder_name = "_closed_loop_sim";
if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

for i = 1:length(x_test)
x_sample = x_test{30}; %44
y_sample = y_test{30}; % 44
net = deep_lstm;

% simulate with known controller and learned system model
bfd_setpoint = 125; % setpoint for controller
time_to_steady_state = 30; % timesteps to replay before simulating controller
curr_bfd_error = 0;
capstan_speed_prev = x_sample(4,1); % set to first value
T = length(y_sample); % usually 8000 with some exceptions;
dt = 0.5;

uarray = zeros(1, T);
earray = zeros(1, T);
nn_output = zeros(1, T);
orders_kd_armax = [4 3 1 3];
orders_kt_armax = [3 2 2 3];

% Preform Velocity Lookup Table
slope = [0 0.01 0.2 0.3 0.5 0.8 1 1.5 2 3 4];
speed100 = [4.2 4 2.5 1.6 0 -3 -5 -8 -10 -10 -10];
speed130 = [2.5 2 1 0.5 -1 -2 -3 -4 -8 -10 -10];

slopefit100 = slope(1:end-2);
speedfit100 = speed100(1:end-2);
slopefit130 = slope(1:end-1);
speedfit130 = speed130(1:end-1);

ft = fittype('a*exp(-b*x)+c','independent','x');
[F100, ~] =  fit(slopefit100', speedfit100', ft);
[F130, ~] =  fit(slopefit130', speedfit130', ft);

lookup_table_100 = @(slope_input) max(F100(slope_input), -10);
lookup_table_130 = @(slope_input) max(F130(slope_input), -10);

% figure; hold on; plot(slope,speed100, 'ro'); plot(slope, speed130, 'k^');
% plot(slope, lookup_table_100(slope), 'r'); plot(slope, lookup_table_130(slope),'k');
% grid minor; title('Exponential Regression of Feed Speed Lookup Table');
% xlabel('Slope of Capstan Speed'); ylabel('Feed Speed @ 2700 mm/min');
% legend({'100mm Preform', '130mm Preform'}); ylim([-12 5]); latexify_plot;

disp("Init Complete!")

%% verify Kd works
% % for i = 1:length(x_test)
%     test_x_input = x_test{i}; % [44 30 17 33 29]
%     test_y_input = y_test{i};
%     test_capstan_speed = test_x_input(1,:);
%     test_bfd = test_y_input;
%     iddata_bfd_to_capstan_speed = iddata(test_capstan_speed', test_bfd'-125,     dt);
%     sys_kd = armax(iddata_bfd_to_capstan_speed, orders_kd_armax);
%     compare(iddata_bfd_to_capstan_speed, sys_kd);
% %     disp(i)
% % end
%% Simulate
for t = 1:time_to_steady_state
    [net, curr_bfd_error] = predictAndUpdateState(net, x_sample(:, t));
    nn_output(t) = curr_bfd_error + 125;
    uarray(t) = x_sample(1, t);
    earray(t) = -curr_bfd_error;
end

for t = time_to_steady_state+1 : T

    %simulate with known controller
    earray(t) = -curr_bfd_error;

%     Kd_output = 1.0 * curr_bfd_error;

%     Kd_output = e * sys_kd_total.B(4) * uarray(end-2) / ...
%         (sys_kd_total.F(1) + sys_kd_total.F(2)*uarray(end)+sys_kd_total.F(3)*uarray(end-1)+sys_kd_total.F(4)*uarray(end-2)+sys_kd_total.F(5)*uarray(end-3));

    Kd_output = lsim(sys_kd_total, earray(1:t), 0:dt:dt*(t-1));
    Kd_output = Kd_output(t);
    
    uarray(t) = Kd_output;
    
    capstan_speed_slope = (Kd_output - capstan_speed_prev)/dt;
    preform_velocity = lookup_table_100(capstan_speed_slope);

    nn_input = [Kd_output; x_sample(2:3, t); preform_velocity];

    [net, curr_bfd_error] = predictAndUpdateState(net, nn_input);
    nn_output(t) = curr_bfd_error + 125;

    capstan_speed_prev = Kd_output;
end
disp('Simulation Done!')

% plots
fig = figure(1);
subplot(3,1,1)
plot(nn_output)
xlabel('t'); ylabel('Simulated BFD');
title('Output of closed loop simulations using learned model')

subplot(3,1,2)
plot(uarray)
xlabel('T'); ylabel('Controller Output');
title('Controller Outputs')

subplot(3,1,3)
plot(earray)
xlabel('T'); ylabel('Controller Input');
title('Controller Input / BFD Error')
latexify_plot;

saveas(fig, sprintf('%s\\%s\\%i', i),'png');  
end