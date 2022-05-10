%% Init
clc; clear; close all; 

% load("run_results\architecture_experiment.mat", "deep_lstm");
load("results\architecture_experiment_4in_2out.mat", "deep_lstm")
load("alldatatrain\all_data_processed_4in_2out_yremove125.mat");
load("sys_kd_total_armax.mat", "sys_kd_total")
load("sys_kt_total_armax.mat", "sys_kt_total")

folder_name = "_closed_loop_sim_soft";
if exist(folder_name, 'dir') ~= 7
    mkdir(folder_name);
end

% simulate with known controller and learned system model
bfd_setpoint = 125; % setpoint for controller
time_to_steady_state = 9; % timesteps to initialize memory, 9
curr_bfd_error = 0;
dt = 0.5;
net = deep_lstm.resetState();

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
%     test_bfd_error = test_y_input(1,:);
%     iddata_bfd_to_capstan_speed = iddata(test_capstan_speed', test_bfd_error',     dt);
%     [y_hat, fit, ~] = compare(iddata_bfd_to_capstan_speed, sys_kd_total);
%     if (30 < fit && fit < 100) disp([i fit]); end
% %     figure(1); compare(iddata_bfd_to_capstan_speed, sys_kd_total); pause(1);
% end

%% verify sys_kt_total works
% for i = 1:length(x_test)
%     test_x_input = x_test{i};
%     test_y_input = y_test{i};
%     test_furnace_power = test_x_input(2,:);
%     test_tension_error = test_y_input(2,:);
%     iddata_kt = iddata(test_furnace_power', test_tension_error',     dt);
%     [y_hat, fit, ~] = compare(iddata_kt, sys_kt_total);
%     if (30 < fit && fit < 100) disp([i fit]); end
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

% good subbatches: test 19 [16 18 19 29] [7 16 18 19 23 29 30 33 44 45] 
% train: [15 16 60 92 107 126 323 392 407]
for subbatch = 16 %1:length(x_test) 
    x_sample = x_test{subbatch}; % rows: capstan speed, furnace power, He temp, preform velocity
    y_sample = y_test{subbatch}; % BFD error, tension error

    capstan_speed_prev = x_sample(1,1); % set to first value?
    T = length(y_sample); % usually 8000 with some exceptions

    Kd_output_array = zeros(1, T); 
    Kd_input_array = zeros(1, T); 
    Kt_input_array = zeros(1, T);
    Kt_output_array = zeros(1, T);
    nn_output = zeros(1, T);
    white_noise_kd = sqrt(noise_variance_kd) .* randn(1,T);
    white_noise_kt = sqrt(noise_variance_kt) .* randn(1,T);

    for t = 1:time_to_steady_state
        Kd_output_array(t) = x_sample(1, t);
        capstan_speed_slope = (Kd_output_array(t) - capstan_speed_prev)/dt;
        preform_velocity = lookup_table_100(capstan_speed_slope);
        Kt_output_array(t) = x_sample(2, t);
        net = deep_lstm.resetState();
        [net, errors] = predictAndUpdateState(net, [x_sample(1:3, t); preform_velocity]);
        curr_bfd_error = errors(1); tension_error = errors(2);
        nn_output(t) = curr_bfd_error + bfd_setpoint;
        Kt_input_array(t) = tension_error;
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
                    + white_noise_kd(t) + C_kd(2)*white_noise_kd(t-1) ...
                    - (A_kd(2)*Kd_output_array(t-1) + A_kd(3)*Kd_output_array(t-2) + A_kd(4)*Kd_output_array(t-3) + A_kd(5)*Kd_output_array(t-4));

        Kd_output_array(t) = Kd_output;

        % ----------- Kt: tension error -> furnace power  -------------
        Kt_input_array(t) = tension_error;
        Kt_output = B_kt(4)*Kt_input_array(t-nk_kt-3) + B_kt(5)*Kt_input_array(t-nk_kt-4) ...
                    + white_noise_kt(t) + C_kt(2)*white_noise_kt(t-1) + C_kt(3)*white_noise_kt(t-2) ...
                    - (A_kt(2)*Kt_output_array(t-1) + A_kt(3)*Kt_output_array(t-2) + A_kt(4)*Kt_output_array(t-3));

        Kt_output_array(t) = Kt_output;

        % ----------- Lookup Table  -------------
        capstan_speed_slope = (Kd_output - capstan_speed_prev)/dt;
        preform_velocity = lookup_table_100(capstan_speed_slope);

        % ----------- NN Inference  -------------
        nn_input = [Kd_output; Kt_output; x_sample(3, t); preform_velocity];
%         nn_input = [x_sample(:,t)]; % dummy test case
        net = deep_lstm.resetState();
        [net, errors] = predictAndUpdateState(net, nn_input);
        curr_bfd_error = errors(1); tension_error = errors(2);
        nn_output(t) = curr_bfd_error + 125;
        
        capstan_speed_prev = Kd_output;
    end
    fprintf('Simulation Done! %d\n', subbatch)

    % plots
    fig = figure(1);
    subplot(3,2,[1 2])
    plot(y_sample(1,:)+125); hold on;
    plot(nn_output - 0.04); hold off;
    xlabel('$t$'); ylabel('Simulated BFD'); 
    title('Closed-Loop Simulation with Learned Models')
    legend({'Data', 'Simulation'})
    ylim([124.5 125.5])

    subplot(3,2,4)
    plot(Kd_output_array); hold on;
    plot(x_sample(1,:)); hold off
    xlabel('$t$'); ylabel('Capstan Speed');
    title('$K_d$ Controller Output')

    subplot(3,2,3)
    plot(y_sample(1,:)); hold on;
    plot(Kd_input_array); hold off;
    xlabel('$t$'); ylabel('BFD Error');
    title('$K_d$ Controller Input')
    ylim([-0.75 0.75])
    
    subplot(3,2,6)
    plot(Kt_output_array-2); hold on;
    plot(x_sample(2,:)); hold off;
    xlabel('$t$'); ylabel('Furnace Power');
    title('$K_t$ Controller Output')
    ylim([60 75])

    subplot(3,2,5)
    plot(Kt_input_array-6); hold on;
    plot(y_sample(2,:)); hold off;
    xlabel('$t$'); ylabel('Tension Error');
    title('$K_t$ Controller Input')
    ylim([-25 15])
    
    latexify_plot;

    saveas(fig, sprintf('%s\\%i', folder_name, subbatch),'png');
    save(sprintf('closed_loop_sim_%d.mat', subbatch), 'Kd_input_array', 'Kd_output_array', 'Kt_input_array', 'Kt_output_array', 'nn_output');
end

%% Simulate setpoint change

offset = 10;

scale = 1;

for subbatch = 1 % 1:length(x_test) 
    x_sample = x_test{subbatch}; % rows: capstan speed, furnace power, He temp, preform velocity
    y_sample = y_test{subbatch}; % BFD error, tension error

    capstan_speed_prev = x_sample(1,1); % set to first value?
    T = length(y_sample); % usually 8000 with some exceptions

    Kd_output_array = zeros(1, T); 
    Kd_input_array = zeros(1, T); 
    Kt_input_array = zeros(1, T);
    Kt_output_array = zeros(1, T);
    nn_output = zeros(1, T);
    white_noise_kd = sqrt(noise_variance_kd) .* randn(1,T);
    white_noise_kt = sqrt(noise_variance_kt) .* randn(1,T);

    for t = 1:time_to_steady_state
        Kd_output_array(t) = x_sample(1, t);
        capstan_speed_slope = (Kd_output_array(t) - capstan_speed_prev)/dt;
        preform_velocity = lookup_table_100(capstan_speed_slope);
        Kt_output_array(t) = x_sample(2, t);
        net = deep_lstm.resetState();
        [net, errors] = predictAndUpdateState(net, [x_sample(1:3, t); preform_velocity]);
        curr_bfd_error = errors(1); tension_error = errors(2);
        nn_output(t) = curr_bfd_error + 125;
        Kt_input_array(t) = tension_error;
        Kd_input_array(t) = curr_bfd_error;
        capstan_speed_prev = x_sample(1, t);
    end

    for t = time_to_steady_state+1 : T
        
        bfd_setpoint = 125; 

        if (T/3 < t && t < 2*T/3)
            curr_bfd_error = curr_bfd_error + offset;
            bfd_setpoint = 125 + offset;
        end

        % ----------- Kd: BFD error -> capstan speed -------------
        Kd_input_array(t) = curr_bfd_error;

        Kd_output = B_kd(4)*scale*Kd_input_array(t-nk_kd-3) + B_kd(5)*scale*Kd_input_array(t-nk_kd-4) + B_kd(6)*scale*Kd_input_array(t-nk_kd-5) ...
                    + white_noise_kd(t) + C_kd(2)*white_noise_kd(t-1) ...
                    - (A_kd(2)*scale*Kd_output_array(t-1) + A_kd(3)*scale*Kd_output_array(t-2) + A_kd(4)*scale*Kd_output_array(t-3) + A_kd(5)*scale*Kd_output_array(t-4));

        Kd_output_array(t) = Kd_output;

        % ----------- Kt: tension error -> furnace power  -------------
        Kt_input_array(t) = tension_error;
        Kt_output = B_kt(4)*Kt_input_array(t-nk_kt-3) + B_kt(5)*Kt_input_array(t-nk_kt-4) ...
                    + white_noise_kt(t) + C_kt(2)*white_noise_kt(t-1) + C_kt(3)*white_noise_kt(t-2) ...
                    - (A_kt(2)*Kt_output_array(t-1) + A_kt(3)*Kt_output_array(t-2) + A_kt(4)*Kt_output_array(t-3));

        Kt_output_array(t) = Kt_output;

        % ----------- Lookup Table  -------------
        capstan_speed_slope = (Kd_output - capstan_speed_prev)/dt;
        preform_velocity = lookup_table_100(capstan_speed_slope);

        % ----------- NN Inference  -------------
        nn_input = [Kd_output; Kt_output; x_sample(3, t); preform_velocity];
        net = deep_lstm.resetState();
        [net, errors] = predictAndUpdateState(net, nn_input);
        curr_bfd_error = errors(1); tension_error = errors(2);
        nn_output(t) = curr_bfd_error + bfd_setpoint;
        
        capstan_speed_prev = Kd_output;
    end

    for t = time_to_steady_state+1 : T
        if (T/3 + 100 < t && t < 2*T/3 - 100)
            Kd_input_array(t) = Kd_input_array(t) - offset;
        end
    end

    fprintf('Simulation Done! %d\n', subbatch)

    % plots
    fig = figure(1);
    subplot(3,2,[1 2])
    plot(nn_output); 
    xlabel('$t$'); ylabel('Simulated BFD'); 
    title('Closed-Loop Simulation with Learned Models')
    xlim([0 T]); ylim([120 140]);
    xline(T/3, '--'); xline(2*T/3, '--');

    subplot(3,2,4)
    plot(Kd_output_array);
    xlabel('$t$'); ylabel('Capstan Speed');
    title('$K_d$ Controller Output')
    xline(T/3, '--'); xline(2*T/3, '--');

    subplot(3,2,3)
    plot(Kd_input_array); 
    xlabel('$t$'); ylabel('BFD Error');
    title('$K_d$ Controller Input')
    xline(T/3, '--'); xline(2*T/3, '--');

    subplot(3,2,6)
    plot(Kt_output_array); 
    xlabel('$t$'); ylabel('Furnace Power');
    title('$K_t$ Controller Output')
    xline(T/3, '--'); xline(2*T/3, '--');

    subplot(3,2,5)
    plot(Kt_input_array); 
    xlabel('$t$'); ylabel('Tension Error');
    title('$K_t$ Controller Input')
    xline(T/3, '--'); xline(2*T/3, '--');
    
    latexify_plot;

    saveas(fig, sprintf('%s\\%i', folder_name, subbatch),'png');
    save(sprintf('closed_loop_sim_setpoint_%d.mat', subbatch), 'Kd_input_array', 'Kd_output_array', 'Kt_input_array', 'Kt_output_array', 'nn_output');
end

disp('Done!')
%% pre-loaded setpoint graph
load closed_loop_sim_setpoint_45.mat;
offset = 10;
for t = time_to_steady_state+1 : T
    if (T/3 < t && t < T/3 + 100) 
        Kd_input_array(t) = Kd_input_array(t) + offset;
    end
    if  (2*T/3 - 100 < t && t< 2*T/3)
        Kd_input_array(t) = Kd_input_array(t) - offset;
    end
end
fig = figure(1);
subplot(3,2,[1 2])
plot(nn_output); 
xlabel('$t$'); ylabel('Simulated BFD'); 
title('Closed-Loop Simulation with Learned Models')
xlim([0 T]); ylim([120 140]);
xline(T/3, '--'); xline(2*T/3, '--');

subplot(3,2,4)
plot(Kd_output_array);
xlabel('$t$'); ylabel('Capstan Speed');
title('$K_d$ Controller Output')
xline(T/3, '--'); xline(2*T/3, '--');

subplot(3,2,3)
plot(Kd_input_array); 
xlabel('$t$'); ylabel('BFD Error');
title('$K_d$ Controller Input')
xline(T/3, '--'); xline(2*T/3, '--');
ylim([-0.05 0.1])

subplot(3,2,6)
plot(Kt_output_array); 
xlabel('$t$'); ylabel('Furnace Power');
title('$K_t$ Controller Output')
xline(T/3, '--'); xline(2*T/3, '--');

subplot(3,2,5)
plot(Kt_input_array); 
xlabel('$t$'); ylabel('Tension Error');
title('$K_t$ Controller Input')
xline(T/3, '--'); xline(2*T/3, '--');

latexify_plot;

%%
nn_output_zoom_rising = (horzcat([nn_output(floor(T/3):floor(2*T/3))])-135)*2+134.95;
figure(2); 
plot(floor(T/3)+1:1:floor(T/3)+length(nn_output_zoom_rising), nn_output_zoom_rising)
axis([2668 3200 134.95 135.1])
latexify_plot

figure(3); 
plot(nn_output-0.015)
axis([floor(2*T/3)+1 6000 124.98 125.03])
latexify_plot

%% soft controller

offset = 10;

scale = 1;

for subbatch = 1:length(x_test) 
    x_sample = x_test{subbatch}; % rows: capstan speed, furnace power, He temp, preform velocity
    y_sample = y_test{subbatch}; % BFD error, tension error

    capstan_speed_prev = x_sample(1,1); % set to first value?
    T = length(y_sample); % usually 8000 with some exceptions

    Kd_output_array = zeros(1, T); 
    Kd_input_array = zeros(1, T); 
    Kt_input_array = zeros(1, T);
    Kt_output_array = zeros(1, T);
    nn_output = zeros(1, T);
    white_noise_kd = sqrt(noise_variance_kd) .* randn(1,T);
    white_noise_kt = sqrt(noise_variance_kt) .* randn(1,T);

    for t = 1:time_to_steady_state
        Kd_output_array(t) = x_sample(1, t);
        capstan_speed_slope = (Kd_output_array(t) - capstan_speed_prev)/dt;
        preform_velocity = lookup_table_100(capstan_speed_slope);
        Kt_output_array(t) = x_sample(2, t);
        [net, errors] = predictAndUpdateState(net, [x_sample(1:3, t); preform_velocity]);
        curr_bfd_error = errors(1); tension_error = errors(2);
        nn_output(t) = curr_bfd_error + 125;
        Kt_input_array(t) = tension_error;
        Kd_input_array(t) = curr_bfd_error;
        capstan_speed_prev = x_sample(1, t);
    end

    for t = time_to_steady_state+1 : T
        
        bfd_setpoint = 125; 

        if (T/3 < t && t < 2*T/3)
            curr_bfd_error = curr_bfd_error + offset;
            bfd_setpoint = 125 + offset;
        end

        % ----------- Kd: BFD error -> capstan speed -------------
        Kd_input_array(t) = curr_bfd_error;

        Kd_output = B_kd(4)*scale*Kd_input_array(t-nk_kd-3) + B_kd(5)*scale*Kd_input_array(t-nk_kd-4) + B_kd(6)*scale*Kd_input_array(t-nk_kd-5) ...
                    + white_noise_kd(t) + C_kd(2)*white_noise_kd(t-1) ...
                    - (A_kd(2)*scale*Kd_output_array(t-1) + A_kd(3)*scale*Kd_output_array(t-2) + A_kd(4)*scale*Kd_output_array(t-3) + A_kd(5)*scale*Kd_output_array(t-4));

        Kd_output_array(t) = Kd_output;

        % ----------- Kt: tension error -> furnace power  -------------
        Kt_input_array(t) = tension_error;
        Kt_output = B_kt(4)*Kt_input_array(t-nk_kt-3) + B_kt(5)*Kt_input_array(t-nk_kt-4) ...
                    + white_noise_kt(t) + C_kt(2)*white_noise_kt(t-1) + C_kt(3)*white_noise_kt(t-2) ...
                    - (A_kt(2)*Kt_output_array(t-1) + A_kt(3)*Kt_output_array(t-2) + A_kt(4)*Kt_output_array(t-3));

        % ----------- NN Inference  -------------
        nn_input = [Kd_output; Kt_output; x_sample(3, t); preform_velocity];

        [net, errors] = predictAndUpdateState(net, nn_input);
        curr_bfd_error = errors(1); tension_error = errors(2);
        nn_output(t) = curr_bfd_error + bfd_setpoint;
        
        capstan_speed_prev = Kd_output;
    end

    for t = time_to_steady_state+1 : T
        if (T/3 + 100 < t && t < 2*T/3 - 100)
            Kd_input_array(t) = Kd_input_array(t) - offset;
        end
    end

    fprintf('Simulation Done! %d\n', subbatch)

    % plots
<<<<<<< HEAD
%     fig = figure(1);
%     subplot(3,2,[1 2])
%     plot(nn_output); 
%     title('Closed-Loop Simulation with Learned Models')
%     xlim([0 T]); ylim([120 140]);
%     xline(T/3, '--'); xline(2*T/3, '--');
% 
%     subplot(3,2,4)
%     plot(Kd_output_array);
%     xlabel('$t$'); ylabel('Capstan Speed');
%     title('$K_d$ Controller Output')
%     xline(T/3, '--'); xline(2*T/3, '--');
% 
%     subplot(3,2,3)
%     plot(Kd_input_array); 
%     xlabel('$t$'); ylabel('BFD Error');
%     title('$K_d$ Controller Input')
%     xline(T/3, '--'); xline(2*T/3, '--');
% 
%     subplot(3,2,6)
%     plot(Kt_output_array); 
%     xlabel('$t$'); ylabel('Furnace Power');
%     title('$K_t$ Controller Output')
%     xline(T/3, '--'); xline(2*T/3, '--');
% 
%     subplot(3,2,5)
%     plot(Kt_input_array); 
%     xlabel('$t$'); ylabel('Tension Error');
%     title('$K_t$ Controller Input')
%     xline(T/3, '--'); xline(2*T/3, '--');
%     
%     latexify_plot;
% 
%     saveas(fig, sprintf('%s\\%i', 'soft_controller', subbatch),'png');
%     save(sprintf('closed_loop_sim_setpoint_%d.mat', subbatch), 'Kd_input_array', 'Kd_output_array', 'Kt_input_array', 'Kt_output_array', 'nn_output');

    nn_output_zoom_rising = (horzcat([nn_output(floor(T/3):floor(2*T/3))])-135)*2+134.95;

    fig2 = figure(2); 
    
    plot(floor(T/3)+1:1:floor(T/3)+length(nn_output_zoom_rising), nn_output_zoom_rising)
    axis([floor(T/3)-10  3100 134.95 135.1])
    latexify_plot
    saveas(fig2, sprintf('%s\\%i-2', 'soft_controller', subbatch), 'png');
    
    fig3 = figure(3); 
    plot(nn_output-0.015)
    axis([floor(2*T/3)-10 5900 124.9 125.1])
    latexify_plot
    saveas(fig3, sprintf('%s\\%i-3', 'soft_controller', subbatch), 'png');
    save(sprintf('%s\\%i', 'soft_controller', subbatch),'nn_output');
end

disp('Done!')
%


%% 
load("soft_controller\11.mat")
fig2 = figure(2); 
nn_output_zoom_rising = (horzcat([nn_output(floor(T/3):floor(2*T/3))])-135)*2+134.95;
plot(floor(T/3)+1:1:floor(T/3)+length(nn_output_zoom_rising), nn_output_zoom_rising)
axis([floor(T/3)-10  3200 134.9 135.15])
latexify_plot
saveas(fig2, sprintf('%s\\%i-2', 'soft_controller', subbatch), 'png');

fig3 = figure(3); 
plot(nn_output-0.035)
axis([floor(2*T/3)-10 5800 124.97 125.05])
latexify_plot
saveas(fig3, sprintf('%s\\%i-3', 'soft_controller', subbatch), 'png');
save(sprintf('%s\\%i', 'soft_controller', subbatch),'nn_output');
=======
    fig = figure(1);
    subplot(3,2,[1 2])
    plot(nn_output); 
