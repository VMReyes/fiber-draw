function [Networklayers] = create_simple_fc(xTrain, yTrain)

[featureDimension,~] = size(xTrain{1});
[numResponses, ~] = size(yTrain{1});

Networklayers = [sequenceInputLayer(featureDimension, "Normalization","zscore") ...
    fullyConnectedLayer(numResponses) ...
    regressionLayer];