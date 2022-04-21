clc; clear; close all; 

load("run_results\architecture_experiment.mat", "deep_lstm");
load("alldatatrain\all_data_processed_4in_1out_yremove125.mat", "x_test", "y_test");
load("sys_kd_total_oe.mat", "m", "sys_kd_total")

x_sample = x_test{1};
y_sample = y_test{1};
net = deep_lstm;

%simulate with known controller and learned system model
bfd_setpoint = 125; % setpoint for controller
time_to_steady_state = 5; % timesteps to replay before simulating controller
curr_bfd = 0;

uarray = [];
earray = [];
YoutSymA_NN = zeros(size(8000));

gainP = 1.0;

% initialize controller state
for ii = 1:time_to_steady_state
    [net, curr_bfd] = predictAndUpdateState(net,x_sample(:,ii));
    curr_bfd = curr_bfd + 125;
    YoutSymA_NN(ii) = curr_bfd;
    uarray{end+1} = x_sample(1:1, ii);
end

for ii = 1:8000-time_to_steady_state
    e = 125 - curr_bfd;
    %e = e*gainP;

    %simulate with known controller

    %Un00 = b0/a0*e;
    
    %Un00 = gainP*e; 

    Un00 = e * sys_kd_total.B(4) * uarray{end-2} / ...
        (sys_kd_total.F(1) + sys_kd_total.F(2)*uarray{end}+sys_kd_total.F(3)*uarray{end-1}+sys_kd_total.F(4)*uarray{end-2}+sys_kd_total.F(5)*uarray{end-3});

    earray{end+1} = e;
    uarray{end+1} = Un00;
    
    %for ref from synthetic: Ynp1 = (Un00*b0 - Yn00*a1 - Ynm1*a2)/a0;
    %Ynp1 = (Un00 - Yn00*coeffY(2) - Ynm1*coeffY(3));
    v = [Un00; x_sample(2:end, ii+1)];
    %v = [Un00 Yn00 Ynm1   ]';
    [net, curr_bfd] = predictAndUpdateState(net,v);
    curr_bfd = curr_bfd + 125;
    %[net,Ynp1] = predictAndUpdateState(net,e);%why not Un00?
    % FIXED ABOVE DO NOT SCALE WITH b0/a0 - that was just in the original synthetic generation / simulation
    %
    YoutSymA_NN(ii+time_to_steady_state) = curr_bfd;
end

figure(1)
subplot(2,1,1)
plot(YoutSymA_NN)

xlabel('T'); ylabel('Simulated BFD');
title('Output of closed loop simulations using learned model')

subplot(2,1,2)
plot(cell2mat(uarray))
xlabel('T'); ylabel('Controller Output');
title('Controller Outputs')
latexify_plot;

