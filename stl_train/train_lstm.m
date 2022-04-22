function [net, info] = train_lstm(net, x_train, y_train, x_test, y_test, miniBatchSize)

maxEpochs       = 4; %150 ; %150; %500; %1000

options = trainingOptions('adam', ...
                'InitialLearnRate', 0.01, ...
                'OutputNetwork','best-validation-loss', ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'GradientThreshold',10, ...
                'SequenceLength','longest', ...             % padding
                'SequencePaddingDirection', 'right', ...    %
                'SequencePaddingValue',0, ...               %
                'Shuffle','once', ...
                'Plots','training-progress',...
                'LearnRateSchedule','piecewise',...
                'LearnRateDropPeriod',2,... %was 100
                'Verbose',1, ...
                'ValidationData', {x_test, y_test});%,...
[net, info] = trainNetwork(x_train,y_train,net,options);
