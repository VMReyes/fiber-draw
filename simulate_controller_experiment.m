net = trainedNetworkModel;

%simulate with known controller and learned system model
Ynm1 = inf;
Yn00 = 0;
Ynp1 = 0;
e = 0; 
Unm1 = 0;
uarray = zeros(size(T));uarray = uarray(:);
earray = zeros(size(T));earray = uarray(:);
YoutSymA_NN = zeros(size(T));
for ii = 1:length(T)
    Ynm1 = Yn00;
    Yn00 = Ynp1;

    Unm1 = Un00;
    
    em1 = e;
    e = Ysetpoint_input2closedloopsystem(ii) - Yn00;
    %e = e*gainP;

    %simulate with known controller

    %Un00 = b0/a0*e; 
    Un00 = gainP*e; 

    earray(ii) = e;
    uarray(ii) = Un00;
    
    %for ref from synthetic: Ynp1 = (Un00*b0 - Yn00*a1 - Ynm1*a2)/a0;
    %Ynp1 = (Un00 - Yn00*coeffY(2) - Ynm1*coeffY(3));
    v = [Un00 ]';
    %v = [Un00 Yn00 Ynm1   ]';
    [net,Ynp1] = predictAndUpdateState(net,v);
    %[net,Ynp1] = predictAndUpdateState(net,e);%why not Un00?
    % FIXED ABOVE DO NOT SCALE WITH b0/a0 - that was just in the original synthetic generation / simulation
    %
    YoutSymA_NN(ii) = Yn00;
end



%init state of model
net = trainedNetworkModel;

%simulate with learned controller and learned system model
Ynm1 = inf;
Yn00 = 0;
Ynp1 = 0;
e = 0; 
uarray = zeros(size(T));uarray = uarray(:);
earray = zeros(size(T));earray = uarray(:);
YoutSymB_NN = zeros(size(T));
for ii = 1:length(T)
    Ynm1 = Yn00;
    Yn00 = Ynp1;

    em1 = e;
    e = Ysetpoint_input2closedloopsystem(ii) - Yn00;
    %e = e*gainP;

    %simulate with learned controller
    Un00 = 20 * (e*coeffU(2) + em1*coeffU(3)); %WHY dT gain.. HERE NOW?
%gainP*1/dT 
%gainP*1/dT^2 

    earray(ii) = e;
    uarray(ii) = Un00;

    %for ref from synthetic: Ynp1 = (Un00*b0 - Yn00*a1 - Ynm1*a2)/a0;
    %Ynp1 = (Un00 - Yn00*coeffY(2) - Ynm1*coeffY(3));
    v = [Un00 ]';
    %v = [Un00 Yn00 Ynm1 ]';
    [net,Ynp1] = predictAndUpdateState(net,v);
    %[net,Ynp1] = predictAndUpdateState(net,e);
    %why not Un00?
    %
    YoutSymB_NN(ii) = Yn00;
end


figure(200001)
subplot(2,1,1)
plot(T,Ysetpoint_input2closedloopsystem)
xlabel('T'); ylabel('Input');
title('Input Data for closed loop simulations using learned model')
subplot(2,1,2)
plot(Y)
hold on
plot(YoutSymA_NN,'r')
plot(YoutSymB_NN,'g')
hold off
xlabel('T'); ylabel('Y NN, output');
legend('training Y', ...
    'simulated Y NN with known controller', ...
    'simulated Y NN with learned controller')
title('Output of closed loop simulations using learned model')

