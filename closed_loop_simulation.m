%% Init
clc; clear; close all; 

load("run_results\architecture_experiment.mat", "deep_lstm");
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat", "x_test", "y_test");
load("sys_kd_total_armax.mat", "sys_kd_total")
load("sys_kt_total_armax.mat", "sys_kt_total")

folder_name = "_closed_loop_sim";
if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

% simulate with known controller and learned system model
bfd_setpoint = 125; % setpoint for controller
time_to_steady_state = 9; % timesteps to initialize memory
curr_bfd_error = 0;
dt = 0.5;
net = deep_lstm;

A_kd = sys_kd_total.A;
B_kd = sys_kd_total.B; 
C_kd = sys_kd_total.C;
noise_variance_kd = sys_kd_total.NoiseVariance;

A_kt = sys_kt_total.A;
B_kt = sys_kt_total.B; 
C_kt = sys_kt_total.C;
noise_variance_kt = sys_kt_total.NoiseVariance;

orders_kd_armax = [4 3 1 3];
nk_kd = 3;
orders_kt_armax = [3 2 2 3];
nk_kt = 3;

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

%% verify sys_kd_total works
% for i = 1:length(x_test)
%     test_x_input = x_test{i};
%     test_y_input = y_test{i};
%     test_capstan_speed = test_x_input(1,:);
%     test_bfd_error = test_y_input;
%     iddata_bfd_to_capstan_speed = iddata(test_capstan_speed', test_bfd_error',     dt);
%     [y_hat, fit, ~] = compare(iddata_bfd_to_capstan_speed, sys_kd_total);
%     disp([i fit]); 
% %     figure(1); compare(iddata_bfd_to_capstan_speed, sys_kd_total); pause(1);
% end

%% verify time stepping works
% test_x_input = x_test{16}; 
% test_y_input = y_test{16};
% test_capstan_speed = test_x_input(1,:);
% test_bfd_error = test_y_input;
% 
% 
% for T = 30:length(test_bfd)
%     iddata_bfd_to_capstan_speed = iddata(test_capstan_speed(1:T)', test_bfd_error(1:T)',dt);
%     [y_hat, fit, ~] = compare(iddata_bfd_to_capstan_speed, sys_kd_total);
%     disp([T fit])
%     figure(1); compare(iddata_bfd_to_capstan_speed, sys_kd_total);
% end

%% Simulate

for subbatch = 19 % good subbatches [16 18 19 29] [7 16 18 19 23 29 30 33 44 45]
    x_sample = x_test{subbatch}; % rows: capstan speed, furnace power, He temp, preform velocity
    y_sample = y_test{subbatch}; % just BFD error

    capstan_speed_prev = x_sample(1,1); % set to first value?
    T = length(y_sample); % usually 8000 with some exceptions

    Kd_output_array = zeros(1, T); %uarray
    Kd_input_array = zeros(1, T); %earray
    Kt_input_array = zeros(1, T);
    Kt_output_array = zeros(1, T);
    nn_output = zeros(1, T);
    white_noise = sqrt(noise_variance_kd) .* randn(1,T);

    disp("Simulation Started...")
    for t = 1:time_to_steady_state
        Kd_output_array(t) = x_sample(1, t);
        capstan_speed_slope = (Kd_output_array(t) - capstan_speed_prev)/dt;
        preform_velocity = lookup_table_100(capstan_speed_slope);
        Kt_output_array(t) = x_sample(2, t);
%         Kt_input_array(t) = ???
        [net, curr_bfd_error] = predictAndUpdateState(net, [x_sample(1:3, t); preform_velocity]);
        nn_output(t) = curr_bfd_error + 125;
        Kd_input_array(t) = curr_bfd_error;
        capstan_speed_prev = x_sample(1, t);
    end

    for t = time_to_steady_state+1 : T

        % ----------- Kd: BFD error -> capstan speed -------------
        Kd_input_array(t) = curr_bfd_error;

%         Kd_output = 1.0 * curr_bfd_error; % dummy test case
% % or this:
%         Kd_output = lsim(sys_kd_total, Kd_input_array(1:t), 0:dt:dt*(t-1));
%         Kd_output = Kd_output(end);
% % or this:
        Kd_output = B_kd(4)*Kd_input_array(t-nk_kd-3) + B_kd(5)*Kd_input_array(t-nk_kd-4) + B_kd(6)*Kd_input_array(t-nk_kd-5) ...
                    + white_noise(t) + C_kd(2)*white_noise(t-1) ...
                    - (A_kd(2)*Kd_output_array(t-1) + A_kd(3)*Kd_output_array(t-2) + A_kd(4)*Kd_output_array(t-3) + A_kd(5)*Kd_output_array(t-4));

        Kd_output_array(t) = Kd_output;

        % ----------- Lookup Table  -------------
        capstan_speed_slope = (Kd_output - capstan_speed_prev)/dt;
        preform_velocity = lookup_table_100(capstan_speed_slope);

        % ----------- Kt: tension error -> furnace power  -------------
%         Kt_input_array(t) = ??? % where is tension signal?
        Kt_output = B_kt(4)*Kt_input_array(t-nk_kt-3) + B_kt(5)*Kt_input_array(t-nk_kt-4) ...
                    + white_noise(t) + C_kt(2)*white_noise(t-1) + C_kt(3)*white_noise(t-2) ...
                    - (A_kt(2)*Kt_output_array(t-1) + A_kt(3)*Kt_output_array(t-2) + A_kt(4)*Kt_output_array(t-3));

        Kt_output_array(t) = Kt_output;

        % ----------- Combine  -------------
        nn_input = [Kd_output; x_sample(2:3, t); preform_velocity];
%         nn_input = [x_sample(:,t)]; % dummy test case

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
    ylim([124.8 125.5])

    subplot(3,1,2)
    plot(Kd_output_array)
    xlabel('T'); ylabel('Capstan Speed');
    title('Kd Controller Output')
    ylim([1400 2800])

    subplot(3,1,3)
    plot(Kd_input_array)
    xlabel('T'); ylabel('BFD Error');
    title('Kd Controller Input')
    ylim([-0.2 0.5])
    
%     subplot(4,1,4)
%     plot(Kt_output_array)
%     xlabel('T'); ylabel('Capstan Speed');
%     title('Kt Controller Output')
%     ylim([-0.2 0.5])

%     subplot(4,1,5)
%     plot(Kt_input_array)
%     xlabel('T'); ylabel('Tension Error');
%     title('Kt Controller Input')
%     ylim([-0.2 0.5])
    
    latexify_plot;

    saveas(fig, sprintf('%s\\%i', folder_name, subbatch),'png');
end