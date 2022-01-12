


% Data Driven approach to Modeling of Real-Time Close Loop Manufacturing Systems 
% --- learning model for System Plant while embedded in a control system
% --- intent is to them use model to develop improved controller.

SystemOption 	= 1; %1 or 2
DataOption      = 1; %1 or 2;

modelcnt2use4sym = 4; % oe (1), arx (2), armax (3), bj (4)
modelcnt2use4symK = 4; % oe (1), arx (2), armax (3), bj (4)

%learn a model for internal Base System but with data from Closed loop operation


% Simulate synthetic data - stepping and lsim to compare
for ff = 1:1
    
%% "Load Data" (simulate synthetic data)
% Here we simulate our system, in reality you would collect data from the
% real physical systen.

% Define a Base model system that is approximately like Fiber Draw - with
% close loop just around angular velocity. (e.g. input is Spindle Ang Vel,
% out of Fiber Diameter)

    num = [1];
    if(SystemOption == 1)
        den = [1 .7 .5];
    end
    if(SystemOption == 2)
        den = [1 2 5];
    end
    sys1 = tf(num,den);
    
% Close loop around diameter:
% Controller is just proportional
    sysK = 1;
% Gain
    gainP = 3;    
% Closed loop fransfer function
    sysFoward = sysK *sys1; 
    sysClose  = gainP*sysFoward / (1 + gainP * sysFoward);

% Define the input sequence(s) for which we will collect data.
% In general these will be the data of normal system operation.
% - in some casese, as we do here, we can do some experiments

% Time vector

    dT = .1;
    T = 0:dT:40;
    
    %input - Option 1
    if(DataOption == 1)
        %INsig = round(sin(2*pi*.2*T))+1;
        INsig = .1*round(sin(2*pi*.2*T));
        INsig = INsig(:);
    end
    %input - Option 2
    if(DataOption == 2)
        INsig = zeros(length(T),1);
        INsig(25:50) = 1;
        INsig(200:250) = -1;
        INsig=idinput(length(T),'PRBS');
        INsig=INsig+1;
    end
    

    
% At this point we can simulate the system using native Maltab tools, but
% that hides some of the data we need to collect from our virtual system.
% Here we simulate for reality check

    sysClose;

    %simulate and add noise
    [Ycheck,ttcheck] = lsim(sysClose,INsig,T);%sys1
    
% Here we simulate the system using difference equation definitions instead
% of the native Maltab functions, as you will want to do this for system you
% learn and in order to collect data internal to the system (alternate
% would be to use native Matlab and state space to make an internal
% paramter an additional output).
    
% grab the coefficients from above
    d0 = den(1);
    d1 = den(2);
    d2 = den(3);

    n0 = num(1);
    
% approximate Base systen on difference equation approximations to the
% differential equation of the such that 
% y[n+1] = u[n]*b0/a0 - y[n]*a1/a0 - y[n-1]*a2/a0;

    a0 = d0/dT^2 + d1/2/dT; 
    a1 = -2*d0/dT^2 + d2;
    a2 = d0/dT^2 - d1/2/dT;

% for the closed loop system then input is the setpoint on Y
    Ysetpoint_input2closedloopsystem = INsig;

% init condtion -- variables storing recent history of output Y
    Ym1 = inf;%read as Y[n-1]
    Y00 = 0; %read as Y[n]
    Yp1 = 0; %read as Y[n+1]

% storage of error signal, command output of controller and output for the
% close loop simulation
    uarray = zeros(size(T));uarray = uarray(:);
    earray = zeros(size(T));earray = earray(:);
    YoutClose = zeros(size(T));
    
% simulate the closed loop system with difference equation definition
%
    for ii = 1:length(T)
        % shift recent history
        Ym1 = Y00;
        Y00 = Yp1;

        %error is diff between setpoint and ouput
        E = Ysetpoint_input2closedloopsystem(ii) - Y00;

        %Here we are modeling a proportional controller, control action is
        %scale of current error
        N_x_U = gainP*E;
        %Un00 = U(ii);

        %log data with 'sensors' on the machine
        earray(ii) = E;
        uarray(ii) = N_x_U;

        %Un00 = Uin(ii);
        %
        Yp1 = (N_x_U*n0 - Y00*a1 - Ym1*a2)/a0;
        %
        YoutClose(ii) = Y00;
    end
    YoutCloseB4noise = YoutClose;
    YoutCloseB4noise = YoutCloseB4noise(:);
    %add some noise
    Ynoise = rand(size(YoutClose))/5 * .1;
    YoutClose = YoutClose+Ynoise;
    YoutClose = YoutClose(:);
    
    earrayClose = earray;
    uarrayClose = uarray;
    %Y = Y - mean(Y);  %for some ML tools, may be helpful

    %plot for reality check confirmation of diff eq sym vs native lsim
    figure(3)
    subplot(2,1,1)
    plot(ttcheck,Ycheck)
    hold on
    plot(T,YoutCloseB4noise,'r')
    plot(T,YoutClose,'g')
    hold off
    legend('Ycheck from lsim','YoutClose from difference eq. loop','YoutClose from difference eq. loop + noise')
    title('Simulation sanity check')
    %
    subplot(2,1,2)
    plot(Ycheck - YoutCloseB4noise,'b')
    hold on
    hold off
    legend('lsim - stepped')
    



end

%% Select the input / output for the block we want to learn - here the PLANT
for ff = 1:1
    %learn a model for internal Base System but with data from Closed loop operation
    INsig4train = uarrayClose;
    Y = YoutClose;


    % 
    %  INsig4train = uarrayBase;
    %  Y = YoutBase;

    figure(100)
    subplot(3,1,1)
    plot(T,INsig)
    xlabel('T'); ylabel('Input');
    title('Input signal to closed loop system ')
    subplot(3,1,2)
    plot(T,INsig4train)
    xlabel('T'); ylabel('Input');
    title('Input Data for training (subsystem) model')
    subplot(3,1,3)
    plot(T,Y)
    xlabel('T'); ylabel('Y');
    title('Output Data for training (subsystem) model')

end


%% Learn 4 models of PLANT with different assumptions on structure
for ff = 1:1 %
    modelcnt = 0;
    for na = [1:4]
        modelcnt = modelcnt + 1;
        % Create Data Structure
        data = iddata(Y,INsig4train,dT);

        if (modelcnt == 1)
            strlabel = 'oe';
            ordrs = [2 2 1]; %[nb nf nk];
            fitsys1 = oe(data, ordrs);
        end
        if (modelcnt == 2)
            strlabel = 'arx'; %[na nb nk]
            ordrs = [ 2 2 1];
            fitsys1 = arx(data, ordrs);
        end
        if (modelcnt == 3)
            strlabel = 'armax';
            ordrs = [ 2 2 2 1]; %[na nb nc nk]
            fitsys1 = armax(data, ordrs);
        end
        if (modelcnt == 4)
            strlabel = 'box-jenkins';
            ordrs = [ 2 2 2 2 1]; %[nb nc nd nf nk]
            fitsys1 = bj(data, ordrs);
        end


        models{modelcnt} =  fitsys1;
        %[YpredictA,ttemp] = lsim(fitsys1,U,T); %works for linear systems
        [Ypredict,ttemp] = sim(fitsys1,INsig4train); %works for linear and nonlinear systems


        figure(10000*na+100)
        subplot(2,3,1)
        plot(T,Y)
        hold on
        plot(T,Ypredict)
        hold off
        xlabel('k'); ylabel('Y, Ypredict');
        title(strlabel)

        % ressidual
        res = Y - Ypredict;
        err = sqrt(sum(res.^2))/length(res);
        errarray(modelcnt) = err;
        %[xcres,lags] = xcorr(res,res, 'coeff');
        [xcres,lags] = xcorr(res,res, floor(length(res)/2), 'unbiased');

        %figure(20)
        subplot(2,3,4)
        plot(res,'.')
        xlabel('k'); ylabel('residual');
        %figure(21)
        subplot(2,3,2)
        hist(res)
        xlabel('bin'); ylabel('Hist res');
        %figure(22)
        subplot(2,3,5)
        normplot(res)
        %figure(23)
        subplot(2,3,6)
        plot(lags,xcres)
        xlabel('lag'); ylabel('auto corr');
        %figure(25)
        subplot(2,3,3)
        plot(Y,Ypredict,'.')
        hold on 
        plot([min(Y) max(Y)],[min(Y) max(Y)],'r')
        hold off
        xlabel('Y'); ylabel('Ypredict');

    end



    %figure(110)
    %bode(sys1)

    figure(200)
    plot(errarray)
    xlabel('model'); ylabel('error');
    title('error vs model type')
    legend('oe (1), arx (2), armax (3), bj (4)')


    figure(300)
    plot(T,Y);
    hold on
    for ii = 1:length(models)
        plot(T,sim(models{ii},INsig4train))
    end
    ylim([min(Y) max(Y)])
    hold off
    legend('training','oe','arx','armax','bj')

end

%% Learn 4 models of CONTROLLER with different assumptions on structure
for ff = 1:1 %
    modelcnt = 0;
    for na = [1:4]
        modelcnt = modelcnt + 1;
        % Create Data Structure
        data = iddata(uarrayClose,earrayClose,dT);

        if (modelcnt == 1)
            strlabel = 'oe';
            ordrs = [2 2 1]; %[nb nf nk];
            fitsys1 = oe(data, ordrs);
        end
        if (modelcnt == 2)
            strlabel = 'arx'; %[na nb nk]
            ordrs = [ 2 2 1];
            fitsys1 = arx(data, ordrs);
        end
        if (modelcnt == 3)
            strlabel = 'armax';
            ordrs = [ 2 2 2 1]; %[na nb nc nk]
            fitsys1 = armax(data, ordrs);
        end
        if (modelcnt == 4)
            strlabel = 'box-jenkins';
            ordrs = [ 2 2 2 2 1]; %[nb nc nd nf nk]
            fitsys1 = bj(data, ordrs);
        end


        modelsK{modelcnt} =  fitsys1;
        %[YpredictA,ttemp] = lsim(fitsys1,U,T); %works for linear systems
        [Upredict,ttemp] = sim(fitsys1,earrayClose); %works for linear and nonlinear systems


        figure(11000*na+100)
        subplot(2,3,1)
        plot(T,uarrayClose)
        hold on
        plot(T,Upredict)
        hold off
        xlabel('k'); ylabel('Y, Ypredict');
        title(strlabel)

        % ressidual
        res = uarrayClose - Upredict;
        err = sqrt(sum(res.^2))/length(res);
        errarray(modelcnt) = err;
        %[xcres,lags] = xcorr(res,res, 'coeff');
        [xcres,lags] = xcorr(res,res, floor(length(res)/2), 'unbiased');

        %figure(20)
        subplot(2,3,4)
        plot(res,'.')
        xlabel('k'); ylabel('residual');
        %figure(21)
        subplot(2,3,2)
        hist(res)
        xlabel('bin'); ylabel('Hist res');
        %figure(22)
        subplot(2,3,5)
        normplot(res)
        %figure(23)
        subplot(2,3,6)
        plot(lags,xcres)
        xlabel('lag'); ylabel('auto corr');
        %figure(25)
        subplot(2,3,3)
        plot(uarrayClose,Upredict,'.')
        hold on 
        plot([min(uarrayClose) max(uarrayClose)],[min(uarrayClose) max(uarrayClose)],'r')
        hold off
        xlabel('U'); ylabel('Upredict');

    end



    %figure(110)
    %bode(sys1)

    figure(210)
    plot(errarray)
    xlabel('model'); ylabel('error');
    title('error vs model type')
    legend('oe (1), arx (2), armax (3), bj (4)')


    figure(310)
    plot(T,uarrayClose);
    hold on
    for ii = 1:length(modelsK)
        plot(T,sim(modelsK{ii},earrayClose))
    end
    ylim([min(uarrayClose) max(uarrayClose)])
    hold off
    legend('training','oe','arx','armax','bj')

end

%% Use the selected learned model for (re)simulation of original known controller and learned plant model
for ff = 1:1
fitsys1 = models{modelcnt2use4sym}
fitsys1K = modelsK{modelcnt2use4symK}

if(modelcnt2use4sym == 1)
    N_coef = fitsys1.B;
    D_coef = fitsys1.F;
end
if(modelcnt2use4sym == 2)
    N_coef = fitsys1.B;
    D_coef = fitsys1.A;
end
if(modelcnt2use4sym == 3)
    N_coef = fitsys1.B;
    D_coef = fitsys1.A;
end
if(modelcnt2use4sym == 4)
    N_coef = fitsys1.B;
    D_coef = fitsys1.F;
end

if(modelcnt2use4symK == 1)
    NK_coef = fitsys1K.B;
    DK_coef = fitsys1K.F;
end
if(modelcnt2use4symK == 2)
    NK_coef = fitsys1K.B;
    DK_coef = fitsys1K.A;
end
if(modelcnt2use4symK == 3)
    NK_coef = fitsys1K.B;
    DK_coef = fitsys1K.A;
end
if(modelcnt2use4symK == 4)
    NK_coef = fitsys1K.B;
    DK_coef = fitsys1K.F;
end

%% Model Use 

    %simulate with KNOWN controller and (partial) learned system model
    Ym1 = inf;
    Y00 = 0;
    Yp1 = 0;
    E = 0; 
    %
    Um2 = 0;
    Um1 = 0;
    U = 0;
    %
    uarray = zeros(size(T));uarray = uarray(:);
    earray = zeros(size(T));earray = earray(:);
    YoutSymA = zeros(size(T));
    for ii = 1:length(T)
        Ym1 = Y00;
        Y00 = Yp1;

        em1 = E;
        E = Ysetpoint_input2closedloopsystem(ii) - Y00;
        %e = e*gainP;

        %Um2 = Um1;
        %Um1 = U;
        %simulate with KNOWN controller
        U = gainP * E ;
        
        %scale based on KNOWN numerator model of plant
        N_x_U =  n0/a0 * U; 
        %numerator times input history - based on LEARNED numerator model of plant
        %N_x_U = (N_coef(2)*U + N_coef(3)*Um1);
        %for ref from synthetic: Ynp1 = (Un00*b0 - Yn00*a1 - Ynm1*a2)/a0;
        Yp1 = (N_x_U - Y00*D_coef(2) - Ym1*D_coef(3));
        %
        earray(ii) = E;
        uarray(ii) = U; %UnScaled;
        YoutSymA(ii) = Y00;
    end
    earraySymA = earray;
    uarraySymA = uarray;
    
    
    %simulate with KNOWN controller and learned system model
    Ym1 = inf;
    Y00 = 0;
    Yp1 = 0;
    E = 0; 
    %
    Um2 = 0;
    Um1 = 0;
    U = 0;
    %
    uarray = zeros(size(T));uarray = uarray(:);
    earray = zeros(size(T));earray = earray(:);
    YoutSymA2 = zeros(size(T));
    for ii = 1:length(T)
        Ym1 = Y00;
        Y00 = Yp1;

        em1 = E;
        E = Ysetpoint_input2closedloopsystem(ii) - Y00;
        %e = e*gainP;

        Um2 = Um1;
        Um1 = U;
        %simulate with KNOWN controller
        U = gainP * E ;
        
        %scale based on KNOWN numerator model of plant
        %N_x_U =  n0/a0 * U; 
        %numerator times input history - based on LEARNED numerator model of plant
        N_x_U = (N_coef(2)*U + N_coef(3)*Um1);
        %for ref from synthetic: Ynp1 = (Un00*b0 - Yn00*a1 - Ynm1*a2)/a0;
        Yp1 = (N_x_U - Y00*D_coef(2) - Ym1*D_coef(3));
        %
        earray(ii) = E;
        uarray(ii) = U; %UnScaled;
        YoutSymA2(ii) = Y00;
    end
    earraySymA2 = earray;
    uarraySymA2 = uarray;


    %simulate with LEARNED controller and learned system model
    Ym1 = inf;  Y00 = 0;    Yp1 = 0;
    Em2 = inf;  Em1 = 0;    E   = 0;
    Um2 = 0;    Um1 = 0;    U = 0;
    %
    uarray = zeros(size(T));uarray = uarray(:);
    earray = zeros(size(T));earray = earray(:);
    YoutSymB = zeros(size(T));
    for ii = 1:length(T)
        Ym1 = Y00;  Y00 = Yp1;
        Em2 = Em1;  Em1 = E;
        Um2 = Um1;  Um1 = U;

        E = Ysetpoint_input2closedloopsystem(ii) - Y00;

        % KNOWN controller
        %U = gainP * E ;
        %
        % LEARNED controller
        %numerator times input history
        NK_x_E  = (NK_coef(2)*E + NK_coef(3)*Em1);
        %and denominator
        U       = NK_x_E - Um1*DK_coef(2) -  Um2*DK_coef(3);

        % PLANT
        %numerator times input history - based on LEARNED numerator model of plant
        N_x_U = (N_coef(2)*U + N_coef(3)*Um1);
        %and denominator.  for ref from synthetic: Ynp1 = (Un00*b0 - Yn00*a1 - Ynm1*a2)/a0;
        Yp1 = N_x_U - Y00*D_coef(2) - Ym1*D_coef(3);
        %

        earray(ii) = E;
        uarray(ii) = U; %N_x_U;
        YoutSymB(ii) = Y00;
    end
    earraySymB = earray;
    uarraySymB = uarray;
    
    figure(100001)
    subplot(2,1,1)
    plot(T,Ysetpoint_input2closedloopsystem)
    xlabel('T'); ylabel('Input');
    title('Input Data for closed loop simulations using learned model')
    subplot(2,1,2)
    plot(Y)
    hold on
    plot(YoutSymA,'m')
    plot(YoutSymA2,'r')
    plot(YoutSymB,'g')
    hold off
    xlabel('T'); ylabel('Y, output');
    legend('training Y','simulated Y with known controller','simulated Y with known controller','simulated Y with learned controller')
    title('Output of closed loop simulations using learned model')

        
    figure(100002)
    subplot(4,1,1)
    plot(T,Ysetpoint_input2closedloopsystem)
    xlabel('T'); ylabel('Input');
    title('Input Data for closed loop simulations using learned model')
    %
    subplot(4,1,2)
    plot(earrayClose)
    hold on
    plot(earraySymA,'m')
    plot(earraySymA2,'r')
    plot(earraySymB,'g')
    hold off
    xlabel('T'); ylabel('E, output');
    legend('training E','partial simulated E with known controller','simulated E with known controller','simulated E with learned controller')
    title('Output of closed loop simulations using learned model')
    %
    subplot(4,1,3)
    plot(uarrayClose)
    hold on
    plot(uarraySymA,'m')
    plot(uarraySymA2,'r')
    plot(uarraySymB,'g')
    hold off
    xlabel('T'); ylabel('U, output');
    legend('training U','partial simulated U with known controller','simulated U with known controller','simulated U with learned controller')
    title('Output of closed loop simulations using learned model')
    %
    subplot(4,1,4)
    plot(Y)
    hold on
    plot(YoutSymA,'m')
    plot(YoutSymA2,'r')
    plot(YoutSymB,'g')
    hold off
    xlabel('T'); ylabel('Y, output');
    legend('training Y','partial simulated Y with known controller','simulated Y with known controller','simulated Y with learned controller')
    title('Output of closed loop simulations using learned model')


end

%return

%%%%%%%%%%%%%%%%%%%%
%% Learn NN model of plant 
%%%%%%%%%%%%%%%%%%%%
for ff = 1:1

%     INsig4train = uarrayClose;
%     Y = YoutClose;
% 
bSubBatchCounter = 1
XTrainTRANSPOSE{bSubBatchCounter}=INsig4train';
%XTrainTRANSPOSE{bSubBatchCounter}=[INsig4train   [0 ; YoutClose(1:end-1) ]   [0 ; 0; YoutClose(2:end-1) ] ]';
YTrainTRANSPOSE{bSubBatchCounter}=YoutClose';

%[rrr,ccc]  = size(INsig4train');
[rrr,ccc]  = size(XTrainTRANSPOSE{1});

    %%
    % Similarly, create a shorter validation signal to use during network training.

    xval = [];%
    yval = [];%

    %% Create and Train Network

    % Create the network architecture.

    numResponses = 1;
    featureDimension = rrr; %1;
    %numHiddenUnits = 10% 100;
    maxEpochs = 200; %500; %1000 %2000
    miniBatchSize = 200;%200;

    Networklayers = [sequenceInputLayer(featureDimension) ...
        lstmLayer(50)...     %50 then  
        lstmLayer(25)... %  
        lstmLayer(10)...
        fullyConnectedLayer(numResponses) ...
        regressionLayer];

    %%
       options = trainingOptions('adam', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'GradientThreshold',10, ...
        'SequenceLength','longest', ...             % padding
        'SequencePaddingDirection', 'right', ...    %
        'SequencePaddingValue',0, ...               %
        'Shuffle','once', ... % 'Plots','training-progress',...     
        'LearnRateSchedule','piecewise',...
        'LearnRateDropPeriod',1900,... %200 %100
        'Verbose',1);%,...
        %'ValidationData',[{xval'} {yval'}]);

    trainedNetworkModel = trainNetwork(XTrainTRANSPOSE,YTrainTRANSPOSE,Networklayers,options);

    %% 
    % Use the trained network to evaluate the system response. Compare the system 
    % and estimated responses in a plot.
    
    Yorig = cell2mat([YTrainTRANSPOSE])';
    Xorig = cell2mat([XTrainTRANSPOSE])';
    
    Ypredict = predict(trainedNetworkModel,Xorig');

    figure(700001)
    title('response estimation')
    plot(Yorig,'k')
    hold on
    plot(Ypredict,'r')
    grid on
    legend('System','Estimated')
    title('Evaluation')

    %%
    %simulate one step at a time preserving state
    % Loop over the time steps in a sequence. 
    % Predict the scores of each time step and update the network state.

    % here we just do a for loop in the same way that predict would using
    % an input which is the original traning inputs and state sequence
    net = trainedNetworkModel;
    numTimeSteps = ccc;
    for ii = 1:numTimeSteps
        v = Xorig(ii,:)';
        [net,score] = predictAndUpdateState(net,v);
        Ypredict_stepped(:,ii) = score;
    end

    % here we just do a for loop in a slightly different way that gets us
    % close to putting a controller in the for loop an input which is the
    % original input and output of  state from previouos predicted state
    net = resetState(trainedNetworkModel);
    numTimeSteps = ccc;
    %
    Ym1=inf;
    Y00 = 0;
    Yp1 = 0;
    for ii = 1:numTimeSteps
        Ym1 = Y00;
        Y00 = Yp1;
        v = [Xorig(ii,1) ]';
        %v = [Xorig(ii,1) Yn00 Ynm1]';
        [net,Yp1] = predictAndUpdateState(net,v);
        Ypredict_stepped2(:,ii) = Yp1;
    end

    figure(700002)
    subplot(2,1,1)
    title('response estimation')
    plot(Yorig,'k')
    hold on
    plot(Ypredict,'r')
    plot(Ypredict_stepped,'b')
    plot(Ypredict_stepped2,'m')
    grid on
    legend('System','Estimated', 'Estimated stepped', 'Estimated stepped closed')
    title('Evaluation')
    subplot(2,1,2)
    plot(Ypredict - Ypredict_stepped,'b+')
    hold on
    plot(Ypredict - Ypredict_stepped2,'m')
    hold off
        legend('Est - Est stepped', 'Est - Est stepped closed')
    title('Evaluation')
    %%

end
%return;

%% Use the selected learned NN model for (re)simulation of original known controller and learned plant model
for ff = 1:1
% Model Use with NN

    %init state of model
    net = resetState(trainedNetworkModel);

    %simulate with LEARNED/KNOWN controller and LEARNED system model
    Ym1 = inf;  Y00 = 0;    Yp1 = 0;
    Em2 = inf;  Em1 = 0;    E   = 0;
    Um2 = 0;    Um1 = 0;    U = 0;
    %
    uarray = zeros(size(T));uarray = uarray(:);
    earray = zeros(size(T));earray = earray(:);
    YoutSymA_NN = zeros(size(T));
    for ii = 1:length(T)
        Ym1 = Y00;  Y00 = Yp1;
        Em2 = Em1;  Em1 = E;
        Um2 = Um1;  Um1 = U;

        E = Ysetpoint_input2closedloopsystem(ii) - Y00;

        % KNOWN controller
        U = gainP * E ;
        %
        % LEARNED controller
        % numerator times input history
        %NK_x_E  = (NK_coef(2)*E + NK_coef(3)*Em1);
        % and denominator
        %U       = NK_x_E - Um1*DK_coef(2) -  Um2*DK_coef(3);

%         % PLANT
%         % numerator times input history - based on LEARNED numerator model of plant
%         N_x_U = (N_coef(2)*U + N_coef(3)*Um1);
%         % and denominator.  for ref from synthetic: Ynp1 = (Un00*b0 - Yn00*a1 - Ynm1*a2)/a0;
%         Yp1 = N_x_U - Y00*D_coef(2) - Ym1*D_coef(3);
%         %

        % PLANT NN
        v = [U ]';
        %v = [Un00 Yn00 Ynm1   ]';
        [net,Yp1] = predictAndUpdateState(net,v);

        earray(ii) = E;
        uarray(ii) = U;
        YoutSymA_NN(ii) = Yp1; %Yp1; %Y00;
    end
    earraySymA_NN = earray;
    uarraySymA_NN = uarray;



    %init state of model
    net = resetState(trainedNetworkModel);

    %simulate with LEARNED/KNOWN controller and LEARNED system model
    Ym1 = inf;  Y00 = 0;    Yp1 = 0;
    Em2 = inf;  Em1 = 0;    E   = 0;
    Um2 = 0;    Um1 = 0;    U = 0;
    %
    uarray = zeros(size(T));uarray = uarray(:);
    earray = zeros(size(T));earray = earray(:);
    YoutSymB_NN = zeros(size(T));
    for ii = 1:length(T)
        Ym1 = Y00;  Y00 = Yp1;
        Em2 = Em1;  Em1 = E;
        Um2 = Um1;  Um1 = U;

        E = Ysetpoint_input2closedloopsystem(ii) - Y00;

        % KNOWN controller
        %U = gainP * E ;
        %
        % LEARNED controller
        % numerator times input history
        NK_x_E  = (NK_coef(2)*E + NK_coef(3)*Em1);
        % and denominator
        U       = NK_x_E - Um1*DK_coef(2) -  Um2*DK_coef(3);

%         % PLANT
%         % numerator times input history - based on LEARNED numerator model of plant
%         N_x_U = (N_coef(2)*U + N_coef(3)*Um1);
%         % and denominator.  for ref from synthetic: Ynp1 = (Un00*b0 - Yn00*a1 - Ynm1*a2)/a0;
%         Yp1 = N_x_U - Y00*D_coef(2) - Ym1*D_coef(3);
%         %

        % PLANT NN
        v = [U ]';
        %v = [Un00 Yn00 Ynm1   ]';
        [net,Yp1] = predictAndUpdateState(net,v);

        earray(ii) = E;
        uarray(ii) = U;
        YoutSymB_NN(ii) = Yp1; %Yp1; %Y00;
    end
    earraySymB_NN = earray;
    uarraySymB_NN = uarray;


    
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

    
end

    %%
save tempsave_dynsymfit
