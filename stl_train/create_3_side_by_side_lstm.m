function [Networklayers] = create_3_side_by_side_lstm(xTrain, yTrain)

[featureDimension,~] = size(xTrain{1});
[numResponses, ~] = size(yTrain{1});

lgraph = layerGraph();

tempLayers = sequenceInputLayer(featureDimension,"Name","sequence");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = lstmLayer(128,"Name","lstm_2");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = lstmLayer(128,"Name","lstm_1");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = lstmLayer(128,"Name","lstm_3");
lgraph = addLayers(lgraph,tempLayers);

tempLayers = [
    concatenationLayer(1,3,"Name","concat")
    fullyConnectedLayer(numResponses,"Name","fc")
    regressionLayer("Name","regressionoutput")];
lgraph = addLayers(lgraph,tempLayers);

% clean up helper variable
clear tempLayers;

lgraph = connectLayers(lgraph,"sequence","lstm_2");
lgraph = connectLayers(lgraph,"sequence","lstm_1");
lgraph = connectLayers(lgraph,"sequence","lstm_3");
lgraph = connectLayers(lgraph,"lstm_2","concat/in2");
lgraph = connectLayers(lgraph,"lstm_1","concat/in1");
lgraph = connectLayers(lgraph,"lstm_3","concat/in3");

Networklayers = lgraph;